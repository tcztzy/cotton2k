use cotton2k::{run_job, RunRequest, RunStatus};
use csv::ReaderBuilder;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};
use std::time::{SystemTime, UNIX_EPOCH};

fn global_regression_lock() -> std::sync::MutexGuard<'static, ()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(())).lock().unwrap()
}

fn make_temp_dir() -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before unix epoch")
        .as_nanos();
    let dir =
        std::env::temp_dir().join(format!("cotton2k-regression-{}-{}", std::process::id(), ts));
    fs::create_dir_all(&dir).expect("failed to create temp regression directory");
    dir
}

fn copy_fixture_dir(src_dir: &Path, dst_dir: &Path) {
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

#[test]
fn minimal_fixture_three_day_regression_key_columns() {
    let _guard = global_regression_lock();
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/minimal");
    let temp_dir = make_temp_dir();
    copy_fixture_dir(&fixture_dir, &temp_dir);

    let profile_path = temp_dir.join("profile.toml");
    let mut profile_toml = fs::read_to_string(&profile_path).expect("failed to read profile");
    profile_toml = profile_toml.replace("stop_date = \"2020-03-31\"", "stop_date = \"2020-04-03\"");
    fs::write(&profile_path, profile_toml).expect("failed to write modified profile");

    // Run simulation on a larger-stack thread: this path still uses large legacy state copies.
    let run_profile_path = profile_path.clone();
    std::thread::Builder::new()
        .stack_size(32 * 1024 * 1024)
        .spawn(move || {
            let run_dir = run_profile_path
                .parent()
                .expect("profile should have parent directory")
                .to_path_buf();
            let summary = run_job(
                RunRequest::new(run_profile_path, run_dir, String::new(), None),
                |_| {},
            )
            .expect("simulation should complete");
            assert_eq!(summary.status, RunStatus::Succeeded);
        })
        .expect("failed to spawn simulation thread")
        .join()
        .expect("simulation thread panicked");

    let output_path = temp_dir.join("output.csv");
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .from_path(&output_path)
        .expect("failed to open output csv");
    let headers = reader
        .headers()
        .expect("failed to read csv headers")
        .clone();

    let idx = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .unwrap_or_else(|| panic!("column '{name}' not found in output"))
    };
    let date_idx = idx("date");
    let lai_idx = idx("leaf_area_index");
    let root_weight_idx = idx("root_weight");
    let swc10_idx = idx("swc0-10");
    let swc20_idx = idx("swc0-20");
    let swc30_idx = idx("swc0-30");

    let rows: Vec<csv::StringRecord> = reader
        .records()
        .map(|r| r.expect("failed to read output row"))
        .collect();
    assert_eq!(rows.len(), 3, "expected 3 simulated days");

    let expected = [
        (
            "2020-04-01",
            0.001,
            0.20800000000000002,
            0.37218740888258084,
            0.36759272052047426,
            0.3665730085273955,
        ),
        (
            "2020-04-02",
            0.001,
            0.20800000000000002,
            0.3705127990396385,
            0.36786902980003505,
            0.3662366712173369,
        ),
        (
            "2020-04-03",
            0.001,
            0.20800000000000002,
            0.3536694472153136,
            0.3677432905162993,
            0.3660897425976588,
        ),
    ];
    let tol = 1e-10;
    for (i, row) in rows.iter().enumerate() {
        let (exp_date, exp_lai, exp_root_weight, exp_swc10, exp_swc20, exp_swc30) = expected[i];
        assert_eq!(&row[date_idx], exp_date, "date mismatch at row {}", i + 1);

        let lai = row[lai_idx]
            .parse::<f64>()
            .expect("leaf_area_index should be numeric");
        let root_weight = row[root_weight_idx]
            .parse::<f64>()
            .expect("root_weight should be numeric");
        let swc10 = row[swc10_idx]
            .parse::<f64>()
            .expect("swc0-10 should be numeric");
        let swc20 = row[swc20_idx]
            .parse::<f64>()
            .expect("swc0-20 should be numeric");
        let swc30 = row[swc30_idx]
            .parse::<f64>()
            .expect("swc0-30 should be numeric");

        assert!(
            (lai - exp_lai).abs() <= tol,
            "leaf_area_index mismatch at row {}",
            i + 1
        );
        assert!(
            (root_weight - exp_root_weight).abs() <= tol,
            "root_weight mismatch at row {}",
            i + 1
        );
        assert!(
            (swc10 - exp_swc10).abs() <= tol,
            "swc0-10 mismatch at row {}",
            i + 1
        );
        assert!(
            (swc20 - exp_swc20).abs() <= tol,
            "swc0-20 mismatch at row {}",
            i + 1
        );
        assert!(
            (swc30 - exp_swc30).abs() <= tol,
            "swc0-30 mismatch at row {}",
            i + 1
        );
    }

    fs::remove_dir_all(&temp_dir).expect("failed to clean regression temp directory");
}
