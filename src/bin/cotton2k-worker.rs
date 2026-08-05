//! JSONL worker-process boundary for one cancellable Cotton2K simulation.
//!
//! The CLI parses a profile, run directory, run identifier, and optional cancel
//! file, then forwards engine events as flushed JSONL records on stdout. Exit
//! codes distinguish success, cancellation, failure, and input/engine errors;
//! simulation artifacts are created by the shared `run_job` implementation.

use cotton2k::{run_job, RunErrorCode, RunEvent, RunRequest, RunStatus};
use std::io::{self, Write};
use std::path::PathBuf;

#[derive(Debug)]
struct WorkerArgs {
    profile: PathBuf,
    run_dir: PathBuf,
    run_id: String,
    cancel_file: Option<PathBuf>,
}

fn parse_worker_args() -> Result<WorkerArgs, String> {
    let mut args = std::env::args().skip(1);
    match args.next().as_deref() {
        Some("run") => {}
        _ => {
            return Err(
                "usage: cotton2k-worker run --profile <path> --run-dir <path> --run-id <id> [--cancel-file <path>]".to_string(),
            )
        }
    }

    let mut profile = None;
    let mut run_dir = None;
    let mut run_id = None;
    let mut cancel_file = None;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--profile" => {
                let value = args
                    .next()
                    .ok_or_else(|| "missing value for --profile".to_string())?;
                profile = Some(PathBuf::from(value));
            }
            "--run-dir" => {
                let value = args
                    .next()
                    .ok_or_else(|| "missing value for --run-dir".to_string())?;
                run_dir = Some(PathBuf::from(value));
            }
            "--run-id" => {
                let value = args
                    .next()
                    .ok_or_else(|| "missing value for --run-id".to_string())?;
                run_id = Some(value);
            }
            "--cancel-file" => {
                let value = args
                    .next()
                    .ok_or_else(|| "missing value for --cancel-file".to_string())?;
                cancel_file = Some(PathBuf::from(value));
            }
            _ => {
                return Err(format!("unknown argument: {arg}"));
            }
        }
    }

    let profile = profile.ok_or_else(|| "--profile is required".to_string())?;
    let run_dir = run_dir.ok_or_else(|| "--run-dir is required".to_string())?;
    let run_id = run_id.ok_or_else(|| "--run-id is required".to_string())?;

    Ok(WorkerArgs {
        profile,
        run_dir,
        run_id,
        cancel_file,
    })
}

fn write_event_line(stdout: &mut io::StdoutLock<'_>, event: &RunEvent) -> io::Result<()> {
    serde_json::to_writer(&mut *stdout, event)?;
    stdout.write_all(b"\n")?;
    stdout.flush()
}

fn main() {
    let parsed = match parse_worker_args() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            std::process::exit(2);
        }
    };

    let request = RunRequest::new(
        parsed.profile,
        parsed.run_dir,
        parsed.run_id.clone(),
        parsed.cancel_file,
    );

    let mut stdout = io::stdout().lock();
    let result = run_job(request, |event| {
        if let Err(err) = write_event_line(&mut stdout, &event) {
            eprintln!("failed to emit worker event: {err}");
        }
    });

    match result {
        Ok(summary) => match summary.status {
            RunStatus::Succeeded => std::process::exit(0),
            RunStatus::Cancelled => std::process::exit(4),
            RunStatus::Failed => std::process::exit(3),
        },
        Err(err) => {
            let event = RunEvent::Failed {
                run_id: parsed.run_id,
                code: format!("{:?}", err.code).to_ascii_lowercase(),
                message: err.message.clone(),
            };
            if let Err(write_err) = write_event_line(&mut stdout, &event) {
                eprintln!("failed to emit worker error event: {write_err}");
            }
            match err.code {
                RunErrorCode::Input => std::process::exit(2),
                RunErrorCode::Io | RunErrorCode::Internal => std::process::exit(3),
            }
        }
    }
}
