use cotton2k::{run_job, RunErrorCode, RunEvent, RunRequest, RunStatus};
use serde_json::Value;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Mutex, OnceLock};
use std::thread;
use std::time::{SystemTime, UNIX_EPOCH};

fn global_runner_lock() -> std::sync::MutexGuard<'static, ()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(())).lock().unwrap()
}

fn make_temp_dir(prefix: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before unix epoch")
        .as_nanos();
    let dir =
        std::env::temp_dir().join(format!("cotton2k-{}-{}-{}", prefix, std::process::id(), ts));
    fs::create_dir_all(&dir).expect("failed to create temp directory");
    dir
}

fn copy_fixture_dir(src_dir: &Path, dst_dir: &Path) {
    fs::create_dir_all(dst_dir).expect("failed to create fixture destination");
    for entry in fs::read_dir(src_dir).expect("failed to read fixture directory") {
        let entry = entry.expect("failed to read fixture entry");
        let src = entry.path();
        let dst = dst_dir.join(
            src.file_name()
                .expect("fixture entry should have file name"),
        );
        fs::copy(&src, dst).expect("failed to copy fixture file");
    }
}

fn set_fixture_stop_date(profile_path: &Path, stop_date: &str) {
    let mut profile_toml = fs::read_to_string(profile_path).expect("failed to read profile");
    profile_toml = profile_toml.replace(
        "stop_date = \"2020-03-31\"",
        &format!("stop_date = \"{stop_date}\""),
    );
    fs::write(profile_path, profile_toml).expect("failed to write profile");
}

fn run_job_with_large_stack(
    request: RunRequest,
) -> (
    Result<cotton2k::RunSummary, cotton2k::RunError>,
    Vec<RunEvent>,
) {
    let events = std::sync::Arc::new(Mutex::new(Vec::new()));
    let events_for_thread = events.clone();
    let result = thread::Builder::new()
        .stack_size(32 * 1024 * 1024)
        .spawn(move || {
            run_job(request, |event| {
                events_for_thread.lock().unwrap().push(event);
            })
        })
        .expect("failed to spawn run_job thread")
        .join()
        .expect("run_job thread panicked");

    let events = events.lock().unwrap().clone();
    (result, events)
}

fn binary_path(name: &str) -> PathBuf {
    let env_key = format!("CARGO_BIN_EXE_{}", name);
    if let Ok(path) = std::env::var(&env_key) {
        return PathBuf::from(path);
    }

    let alt_key = env_key.replace('-', "_");
    if let Ok(path) = std::env::var(&alt_key) {
        return PathBuf::from(path);
    }

    let exe_name = if cfg!(windows) {
        format!("{name}.exe")
    } else {
        name.to_string()
    };
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("debug")
        .join(exe_name)
}

#[test]
fn run_job_supports_toml_and_json_profiles() {
    let _guard = global_runner_lock();
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/minimal");
    let temp_dir = make_temp_dir("runner-profile-formats");
    let work_dir = temp_dir.join("work");
    copy_fixture_dir(&fixture_dir, &work_dir);

    let profile_toml_path = work_dir.join("profile.toml");
    set_fixture_stop_date(&profile_toml_path, "2020-04-03");

    let toml_source = fs::read_to_string(&profile_toml_path).expect("failed to read TOML profile");
    let toml_value: toml::Value = toml::from_str(&toml_source).expect("failed to parse TOML value");
    let profile_json_path = work_dir.join("profile.json");
    fs::write(
        &profile_json_path,
        serde_json::to_string_pretty(&toml_value).expect("failed to serialize JSON profile"),
    )
    .expect("failed to write JSON profile");

    let toml_run_dir = temp_dir.join("runs/toml");
    let json_run_dir = temp_dir.join("runs/json");

    let (toml_result, _events) = run_job_with_large_stack(RunRequest::new(
        profile_toml_path,
        toml_run_dir.clone(),
        "toml-run".to_string(),
        None,
    ));
    let (json_result, _events) = run_job_with_large_stack(RunRequest::new(
        profile_json_path,
        json_run_dir.clone(),
        "json-run".to_string(),
        None,
    ));

    let toml_summary = toml_result.expect("toml run should succeed");
    let json_summary = json_result.expect("json run should succeed");
    assert_eq!(toml_summary.status, RunStatus::Succeeded);
    assert_eq!(json_summary.status, RunStatus::Succeeded);

    let toml_output =
        fs::read_to_string(toml_run_dir.join("output.csv")).expect("missing TOML output");
    let json_output =
        fs::read_to_string(json_run_dir.join("output.csv")).expect("missing JSON output");
    assert_eq!(
        toml_output, json_output,
        "TOML and JSON runs should produce identical CSV output"
    );

    fs::remove_dir_all(temp_dir).expect("failed to clean temp directory");
}

#[test]
fn run_job_redirects_output_to_run_directory() {
    let _guard = global_runner_lock();
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/minimal");
    let temp_dir = make_temp_dir("runner-output-dir");
    let work_dir = temp_dir.join("work");
    copy_fixture_dir(&fixture_dir, &work_dir);

    let profile_path = work_dir.join("profile.toml");
    set_fixture_stop_date(&profile_path, "2020-04-03");

    let legacy_output_path = work_dir.join("output.csv");
    if legacy_output_path.exists() {
        fs::remove_file(&legacy_output_path).expect("failed to remove fixture output.csv");
    }

    let run_dir = temp_dir.join("runs/job-one");
    let (result, _events) = run_job_with_large_stack(RunRequest::new(
        profile_path,
        run_dir.clone(),
        "job-one".to_string(),
        None,
    ));

    let summary = result.expect("run should succeed");
    assert_eq!(summary.status, RunStatus::Succeeded);
    assert!(
        run_dir.join("output.csv").exists(),
        "run output should be in run_dir"
    );
    assert!(
        !legacy_output_path.exists(),
        "legacy output path should stay untouched"
    );
    assert!(
        run_dir.join("meta.json").exists(),
        "meta.json should be generated"
    );
    assert!(
        run_dir.join("request.json").exists(),
        "request.json should be generated"
    );
    assert!(
        run_dir.join("input/profile.toml").exists(),
        "input profile snapshot should be generated"
    );

    fs::remove_dir_all(temp_dir).expect("failed to clean temp directory");
}

#[test]
fn run_job_reports_input_errors_with_classification() {
    let _guard = global_runner_lock();
    let temp_dir = make_temp_dir("runner-errors");
    let missing_profile = temp_dir.join("missing-profile.toml");
    let run_dir = temp_dir.join("runs/missing");

    let result = run_job(
        RunRequest::new(
            missing_profile.clone(),
            run_dir,
            "bad-run".to_string(),
            None,
        ),
        |_| {},
    );

    let err = result.expect_err("run_job should reject missing profile path");
    assert_eq!(err.code, RunErrorCode::Input);
    assert!(
        err.message.contains(&missing_profile.display().to_string()),
        "error message should include missing profile path"
    );

    fs::remove_dir_all(temp_dir).expect("failed to clean temp directory");
}

#[test]
fn worker_supports_parallel_jobs_in_separate_run_dirs() {
    let _guard = global_runner_lock();
    let worker_bin = binary_path("cotton2k-worker");
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/minimal");
    let temp_dir = make_temp_dir("worker-parallel");

    let case_a = temp_dir.join("case-a");
    let case_b = temp_dir.join("case-b");
    copy_fixture_dir(&fixture_dir, &case_a);
    copy_fixture_dir(&fixture_dir, &case_b);
    set_fixture_stop_date(&case_a.join("profile.toml"), "2020-04-03");
    set_fixture_stop_date(&case_b.join("profile.toml"), "2020-04-03");

    let run_dir_a = temp_dir.join("runs/a");
    let run_dir_b = temp_dir.join("runs/b");

    let child_a = Command::new(&worker_bin)
        .args([
            "run",
            "--profile",
            case_a.join("profile.toml").to_str().unwrap(),
            "--run-dir",
            run_dir_a.to_str().unwrap(),
            "--run-id",
            "job-a",
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn worker A");

    let child_b = Command::new(&worker_bin)
        .args([
            "run",
            "--profile",
            case_b.join("profile.toml").to_str().unwrap(),
            "--run-dir",
            run_dir_b.to_str().unwrap(),
            "--run-id",
            "job-b",
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn worker B");

    let out_a = child_a
        .wait_with_output()
        .expect("failed waiting for worker A");
    let out_b = child_b
        .wait_with_output()
        .expect("failed waiting for worker B");

    assert!(
        out_a.status.success(),
        "worker A failed: {}",
        String::from_utf8_lossy(&out_a.stderr)
    );
    assert!(
        out_b.status.success(),
        "worker B failed: {}",
        String::from_utf8_lossy(&out_b.stderr)
    );

    assert!(
        run_dir_a.join("output.csv").exists(),
        "worker A should generate output.csv"
    );
    assert!(
        run_dir_b.join("output.csv").exists(),
        "worker B should generate output.csv"
    );

    let events_a: Vec<Value> = String::from_utf8_lossy(&out_a.stdout)
        .lines()
        .map(|line| serde_json::from_str(line).expect("worker A should emit valid JSONL"))
        .collect();
    let events_b: Vec<Value> = String::from_utf8_lossy(&out_b.stdout)
        .lines()
        .map(|line| serde_json::from_str(line).expect("worker B should emit valid JSONL"))
        .collect();
    assert!(events_a.iter().any(|e| e["event"] == "started"));
    assert!(events_a.iter().any(|e| e["event"] == "finished"));
    assert!(events_b.iter().any(|e| e["event"] == "started"));
    assert!(events_b.iter().any(|e| e["event"] == "finished"));

    fs::remove_dir_all(temp_dir).expect("failed to clean temp directory");
}

#[test]
fn worker_marks_run_as_cancelled_when_cancel_file_exists() {
    let _guard = global_runner_lock();
    let worker_bin = binary_path("cotton2k-worker");
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/minimal");
    let temp_dir = make_temp_dir("worker-cancel");
    let work_dir = temp_dir.join("work");
    copy_fixture_dir(&fixture_dir, &work_dir);
    set_fixture_stop_date(&work_dir.join("profile.toml"), "2020-12-31");

    let run_dir = temp_dir.join("runs/cancelled");
    let cancel_file = temp_dir.join("cancel.signal");
    fs::write(&cancel_file, "cancel").expect("failed to write cancel signal file");

    let output = Command::new(&worker_bin)
        .args([
            "run",
            "--profile",
            work_dir.join("profile.toml").to_str().unwrap(),
            "--run-dir",
            run_dir.to_str().unwrap(),
            "--run-id",
            "cancel-job",
            "--cancel-file",
            cancel_file.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run worker cancellation scenario");

    assert_eq!(
        output.status.code(),
        Some(4),
        "worker should exit with code 4 when cancelled"
    );

    let meta: Value = serde_json::from_str(
        &fs::read_to_string(run_dir.join("meta.json")).expect("missing cancellation meta.json"),
    )
    .expect("invalid meta.json");
    assert_eq!(meta["status"], "cancelled");

    fs::remove_dir_all(temp_dir).expect("failed to clean temp directory");
}

#[test]
fn cli_default_output_path_stays_profile_directory() {
    let _guard = global_runner_lock();
    let cli_bin = binary_path("cotton2k");
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/minimal");
    let temp_dir = make_temp_dir("cli-compat");
    let work_dir = temp_dir.join("work");
    copy_fixture_dir(&fixture_dir, &work_dir);

    let profile_path = work_dir.join("profile.toml");
    set_fixture_stop_date(&profile_path, "2020-04-03");

    let output_path = work_dir.join("output.csv");
    if output_path.exists() {
        fs::remove_file(&output_path).expect("failed to remove fixture output.csv before cli run");
    }

    let output = Command::new(&cli_bin)
        .arg(profile_path)
        .output()
        .expect("failed to run cli binary");

    assert!(
        output.status.success(),
        "cli should succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output_path.exists(),
        "cli should keep writing output.csv in profile directory"
    );

    fs::remove_dir_all(temp_dir).expect("failed to clean temp directory");
}
