//! Parallel batch scheduler for independent Cotton2K worker processes.
//!
//! The binary reads a JSON jobs file, assigns each profile a run directory, and
//! keeps at most the requested number of child workers active. Results are
//! collected into a JSON summary; process creation, filesystem access, and
//! cancellation-file observation are intentional side effects of this CLI.

use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::thread;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

#[derive(Debug, Deserialize)]
struct BatchJobInput {
    profile: PathBuf,
    #[serde(default)]
    run_id: String,
    cancel_file: Option<PathBuf>,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum BatchJobFile {
    Jobs(Vec<BatchJobInput>),
    Wrapped { jobs: Vec<BatchJobInput> },
}

#[derive(Debug)]
struct BatchArgs {
    jobs_path: PathBuf,
    runs_root: PathBuf,
    max_parallel: Option<usize>,
    worker_bin: Option<PathBuf>,
}

#[derive(Debug)]
struct ScheduledJob {
    profile: PathBuf,
    run_id: String,
    run_dir: PathBuf,
    cancel_file: Option<PathBuf>,
}

#[derive(Debug)]
struct ActiveJob {
    child: Child,
    spec: ScheduledJob,
}

#[derive(Debug, Serialize)]
struct BatchJobResult {
    run_id: String,
    profile: PathBuf,
    run_dir: PathBuf,
    status: String,
    exit_code: i32,
}

#[derive(Debug, Serialize)]
struct BatchSummary {
    total: usize,
    succeeded: usize,
    failed: usize,
    cancelled: usize,
    max_parallel: usize,
    results: Vec<BatchJobResult>,
}

fn parse_args() -> Result<BatchArgs, String> {
    let mut args = std::env::args().skip(1);
    match args.next().as_deref() {
        Some("run") => {}
        _ => {
            return Err(
                "usage: cotton2k-batch run --jobs <path> --runs-root <path> [--max-parallel <n>] [--worker-bin <path>]".to_string(),
            )
        }
    }

    let mut jobs_path = None;
    let mut runs_root = None;
    let mut max_parallel = None;
    let mut worker_bin = None;
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--jobs" => {
                let value = args
                    .next()
                    .ok_or_else(|| "missing value for --jobs".to_string())?;
                jobs_path = Some(PathBuf::from(value));
            }
            "--runs-root" => {
                let value = args
                    .next()
                    .ok_or_else(|| "missing value for --runs-root".to_string())?;
                runs_root = Some(PathBuf::from(value));
            }
            "--max-parallel" => {
                let value = args
                    .next()
                    .ok_or_else(|| "missing value for --max-parallel".to_string())?;
                let parsed = value
                    .parse::<usize>()
                    .map_err(|_| format!("invalid --max-parallel value: {value}"))?;
                if parsed == 0 {
                    return Err("--max-parallel must be at least 1".to_string());
                }
                max_parallel = Some(parsed);
            }
            "--worker-bin" => {
                let value = args
                    .next()
                    .ok_or_else(|| "missing value for --worker-bin".to_string())?;
                worker_bin = Some(PathBuf::from(value));
            }
            _ => return Err(format!("unknown argument: {arg}")),
        }
    }

    let jobs_path = jobs_path.ok_or_else(|| "--jobs is required".to_string())?;
    let runs_root = runs_root.ok_or_else(|| "--runs-root is required".to_string())?;
    Ok(BatchArgs {
        jobs_path,
        runs_root,
        max_parallel,
        worker_bin,
    })
}

fn resolve_jobs(path: &Path) -> io::Result<Vec<BatchJobInput>> {
    let raw = fs::read_to_string(path)?;
    let parsed: BatchJobFile = serde_json::from_str(&raw).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("failed to parse jobs file '{}': {e}", path.display()),
        )
    })?;
    let jobs = match parsed {
        BatchJobFile::Jobs(jobs) => jobs,
        BatchJobFile::Wrapped { jobs } => jobs,
    };
    Ok(jobs)
}

fn generated_run_id(index: usize) -> String {
    let now = chrono::Utc::now().format("%Y%m%d-%H%M%S");
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .subsec_nanos();
    format!("{now}-{:04}-{:08x}", index, nanos)
}

fn default_max_parallel() -> usize {
    let physical = num_cpus::get_physical();
    std::cmp::min(4, physical.max(1))
}

fn worker_bin_name() -> &'static str {
    if cfg!(windows) {
        "cotton2k-worker.exe"
    } else {
        "cotton2k-worker"
    }
}

fn resolve_worker_bin(explicit: Option<PathBuf>) -> PathBuf {
    if let Some(path) = explicit {
        return path;
    }
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_cotton2k-worker") {
        return PathBuf::from(path);
    }
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_cotton2k_worker") {
        return PathBuf::from(path);
    }
    if let Ok(current_exe) = std::env::current_exe() {
        let candidate = current_exe.with_file_name(worker_bin_name());
        if candidate.is_file() {
            return candidate;
        }
    }
    PathBuf::from("cotton2k-worker")
}

fn to_scheduled_jobs(jobs: Vec<BatchJobInput>, runs_root: &Path) -> VecDeque<ScheduledJob> {
    jobs.into_iter()
        .enumerate()
        .map(|(index, job)| {
            let run_id = if job.run_id.trim().is_empty() {
                generated_run_id(index)
            } else {
                job.run_id.trim().to_string()
            };
            let run_dir = runs_root.join(&run_id);
            ScheduledJob {
                profile: job.profile,
                run_id,
                run_dir,
                cancel_file: job.cancel_file,
            }
        })
        .collect()
}

fn spawn_worker(worker_bin: &PathBuf, job: &ScheduledJob) -> io::Result<Child> {
    let mut command = Command::new(worker_bin);
    command
        .arg("run")
        .arg("--profile")
        .arg(&job.profile)
        .arg("--run-dir")
        .arg(&job.run_dir)
        .arg("--run-id")
        .arg(&job.run_id);
    if let Some(cancel_file) = &job.cancel_file {
        command.arg("--cancel-file").arg(cancel_file);
    }
    command.stdout(Stdio::null());
    command.stderr(Stdio::inherit());
    command.spawn()
}

fn classify_exit(exit_code: i32) -> &'static str {
    match exit_code {
        0 => "succeeded",
        4 => "cancelled",
        _ => "failed",
    }
}

fn run_batch(args: BatchArgs) -> Result<BatchSummary, String> {
    let jobs = resolve_jobs(&args.jobs_path).map_err(|e| e.to_string())?;
    if jobs.is_empty() {
        return Err("jobs file contains no jobs".to_string());
    }
    fs::create_dir_all(&args.runs_root).map_err(|e| {
        format!(
            "failed to create runs root '{}': {e}",
            args.runs_root.display()
        )
    })?;

    let max_parallel = args
        .max_parallel
        .unwrap_or_else(default_max_parallel)
        .max(1);
    let worker_bin = resolve_worker_bin(args.worker_bin);
    let mut pending = to_scheduled_jobs(jobs, &args.runs_root);
    let total = pending.len();
    let mut active: Vec<ActiveJob> = Vec::new();
    let mut results = Vec::with_capacity(total);

    while !pending.is_empty() || !active.is_empty() {
        while active.len() < max_parallel {
            let Some(job) = pending.pop_front() else {
                break;
            };
            match spawn_worker(&worker_bin, &job) {
                Ok(child) => active.push(ActiveJob { child, spec: job }),
                Err(err) => {
                    results.push(BatchJobResult {
                        run_id: job.run_id,
                        profile: job.profile,
                        run_dir: job.run_dir,
                        status: "failed".to_string(),
                        exit_code: 3,
                    });
                    eprintln!("failed to spawn worker process: {err}");
                }
            }
        }

        let mut completed_any = false;
        let mut index = 0usize;
        while index < active.len() {
            match active[index].child.try_wait() {
                Ok(Some(status)) => {
                    let active_job = active.swap_remove(index);
                    let code = status.code().unwrap_or(3);
                    results.push(BatchJobResult {
                        run_id: active_job.spec.run_id,
                        profile: active_job.spec.profile,
                        run_dir: active_job.spec.run_dir,
                        status: classify_exit(code).to_string(),
                        exit_code: code,
                    });
                    completed_any = true;
                }
                Ok(None) => {
                    index += 1;
                }
                Err(err) => {
                    let active_job = active.swap_remove(index);
                    results.push(BatchJobResult {
                        run_id: active_job.spec.run_id,
                        profile: active_job.spec.profile,
                        run_dir: active_job.spec.run_dir,
                        status: "failed".to_string(),
                        exit_code: 3,
                    });
                    eprintln!("failed to poll worker process: {err}");
                    completed_any = true;
                }
            }
        }

        if !completed_any && !active.is_empty() {
            thread::sleep(Duration::from_millis(50));
        }
    }

    let succeeded = results.iter().filter(|r| r.status == "succeeded").count();
    let cancelled = results.iter().filter(|r| r.status == "cancelled").count();
    let failed = results.len().saturating_sub(succeeded + cancelled);

    Ok(BatchSummary {
        total,
        succeeded,
        failed,
        cancelled,
        max_parallel,
        results,
    })
}

fn main() {
    let args = match parse_args() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            std::process::exit(2);
        }
    };

    let summary = match run_batch(args) {
        Ok(summary) => summary,
        Err(message) => {
            eprintln!("{message}");
            std::process::exit(2);
        }
    };

    if let Err(err) = serde_json::to_writer_pretty(std::io::stdout(), &summary) {
        eprintln!("failed to write batch summary JSON: {err}");
        std::process::exit(3);
    }
    println!();

    if summary.failed > 0 {
        std::process::exit(3);
    }
    if summary.cancelled > 0 {
        std::process::exit(4);
    }
}
