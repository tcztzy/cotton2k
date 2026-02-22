use crate::{Profile, State};
use chrono::{NaiveDate, Utc};
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fmt;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunRequest {
    pub profile_path: PathBuf,
    pub run_dir: PathBuf,
    pub run_id: String,
    pub cancel_file: Option<PathBuf>,
}

impl RunRequest {
    pub fn new(
        profile_path: PathBuf,
        run_dir: PathBuf,
        run_id: String,
        cancel_file: Option<PathBuf>,
    ) -> Self {
        Self {
            profile_path,
            run_dir,
            run_id,
            cancel_file,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunStatus {
    Succeeded,
    Failed,
    Cancelled,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunFailure {
    pub code: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunSummary {
    pub run_id: String,
    pub status: RunStatus,
    pub started_at: String,
    pub finished_at: String,
    pub output_csv_path: PathBuf,
    pub meta_path: PathBuf,
    pub days_simulated: usize,
    pub error: Option<RunFailure>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "event", rename_all = "snake_case")]
pub enum RunEvent {
    Started {
        run_id: String,
        started_at: String,
    },
    Progress {
        run_id: String,
        day_index: usize,
        date: String,
    },
    Finished {
        run_id: String,
        finished_at: String,
    },
    Failed {
        run_id: String,
        code: String,
        message: String,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunErrorCode {
    Input,
    Io,
    Internal,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunError {
    pub code: RunErrorCode,
    pub message: String,
}

impl RunError {
    fn new(code: RunErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.code, self.message)
    }
}

impl Error for RunError {}

#[derive(Debug, Clone)]
struct ResolvedProfilePaths {
    profile_path: PathBuf,
    weather_path: PathBuf,
    soil_impedance_path: PathBuf,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct ExecutionOutcome {
    pub(crate) days_simulated: usize,
    pub(crate) cancelled: bool,
}

#[derive(Debug, Serialize)]
struct RunMetadata {
    run_id: String,
    status: RunStatus,
    started_at: String,
    finished_at: String,
    duration_ms: u64,
    profile_path: PathBuf,
    weather_path: PathBuf,
    soil_impedance_path: PathBuf,
    output_csv: PathBuf,
    error: Option<RunFailure>,
    version: String,
}

fn now_rfc3339() -> String {
    Utc::now().to_rfc3339()
}

fn generated_run_id() -> String {
    let now = Utc::now();
    let prefix = now.format("%Y%m%d-%H%M%S").to_string();
    let suffix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .subsec_nanos();
    format!("{prefix}-{suffix:08x}")
}

fn same_directory(lhs: &Path, rhs: &Path) -> bool {
    if lhs == rhs {
        return true;
    }
    match (fs::canonicalize(lhs), fs::canonicalize(rhs)) {
        (Ok(a), Ok(b)) => a == b,
        _ => false,
    }
}

fn parse_profile_file(profile_path: &Path) -> Result<(Profile, String, String), RunError> {
    let mut file = fs::File::open(profile_path).map_err(|e| {
        RunError::new(
            RunErrorCode::Input,
            format!("failed to open profile '{}': {}", profile_path.display(), e),
        )
    })?;
    let mut contents = String::new();
    file.read_to_string(&mut contents).map_err(|e| {
        RunError::new(
            RunErrorCode::Input,
            format!("failed to read profile '{}': {}", profile_path.display(), e),
        )
    })?;

    let extension = profile_path
        .extension()
        .and_then(|x| x.to_str())
        .unwrap_or("toml")
        .to_ascii_lowercase();

    let profile = if extension == "toml" {
        toml::from_str::<Profile>(&contents).map_err(|e| {
            RunError::new(
                RunErrorCode::Input,
                format!(
                    "failed to parse TOML profile '{}': {}",
                    profile_path.display(),
                    e
                ),
            )
        })?
    } else {
        serde_json::from_str::<Profile>(&contents).map_err(|e| {
            RunError::new(
                RunErrorCode::Input,
                format!(
                    "failed to parse JSON profile '{}': {}",
                    profile_path.display(),
                    e
                ),
            )
        })?
    };

    Ok((profile, contents, extension))
}

fn configure_profile_paths(
    profile: &mut Profile,
    profile_path: &Path,
    run_dir: &Path,
) -> Result<ResolvedProfilePaths, RunError> {
    let profile_dir = profile_path.parent().ok_or_else(|| {
        RunError::new(
            RunErrorCode::Input,
            format!(
                "profile path '{}' has no parent directory",
                profile_path.display()
            ),
        )
    })?;
    let profile_filename = profile_path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("profile.toml");
    profile.path = run_dir.join(profile_filename);

    if profile.weather_path.is_relative() {
        profile.weather_path = profile_dir.join(&profile.weather_path);
    }

    let resolved_soil_impedance = match profile.soil_impedance.clone() {
        Some(path) if path.is_relative() => profile_dir.join(path),
        Some(path) => path,
        None => profile_dir.join("soil_imp.csv"),
    };
    profile.soil_impedance = Some(resolved_soil_impedance.clone());

    Ok(ResolvedProfilePaths {
        profile_path: profile_path.to_path_buf(),
        weather_path: profile.weather_path.clone(),
        soil_impedance_path: resolved_soil_impedance,
    })
}

fn write_json_pretty(path: &Path, value: &impl Serialize) -> Result<(), RunError> {
    let bytes = serde_json::to_vec_pretty(value).map_err(|e| {
        RunError::new(
            RunErrorCode::Internal,
            format!("failed to serialize JSON '{}': {}", path.display(), e),
        )
    })?;
    fs::write(path, bytes).map_err(|e| {
        RunError::new(
            RunErrorCode::Io,
            format!("failed to write file '{}': {}", path.display(), e),
        )
    })
}

fn persist_run_inputs(
    request: &RunRequest,
    profile_contents: &str,
    profile_extension: &str,
) -> Result<(), RunError> {
    let request_path = request.run_dir.join("request.json");
    write_json_pretty(&request_path, request)?;

    let input_dir = request.run_dir.join("input");
    fs::create_dir_all(&input_dir).map_err(|e| {
        RunError::new(
            RunErrorCode::Io,
            format!(
                "failed to create input snapshot directory '{}': {}",
                input_dir.display(),
                e
            ),
        )
    })?;

    let ext = if profile_extension == "toml" {
        "toml"
    } else {
        "json"
    };
    let snapshot_path = input_dir.join(format!("profile.{ext}"));
    fs::write(&snapshot_path, profile_contents).map_err(|e| {
        RunError::new(
            RunErrorCode::Io,
            format!(
                "failed to write profile snapshot '{}': {}",
                snapshot_path.display(),
                e
            ),
        )
    })
}

fn write_metadata(
    request: &RunRequest,
    resolved: &ResolvedProfilePaths,
    summary: &RunSummary,
    duration_ms: u64,
) -> Result<(), RunError> {
    let metadata = RunMetadata {
        run_id: summary.run_id.clone(),
        status: summary.status,
        started_at: summary.started_at.clone(),
        finished_at: summary.finished_at.clone(),
        duration_ms,
        profile_path: resolved.profile_path.clone(),
        weather_path: resolved.weather_path.clone(),
        soil_impedance_path: resolved.soil_impedance_path.clone(),
        output_csv: summary.output_csv_path.clone(),
        error: summary.error.clone(),
        version: env!("CARGO_PKG_VERSION").to_string(),
    };
    let meta_path = request.run_dir.join("meta.json");
    write_json_pretty(&meta_path, &metadata)
}

pub(crate) fn execute_profile(
    profile: &mut Profile,
    mut should_cancel: impl FnMut() -> bool,
    mut on_progress: impl FnMut(usize, NaiveDate),
) -> Result<ExecutionOutcome, Box<dyn Error>> {
    profile.initialize()?;
    profile.output_file_headers()?;
    // Keep legacy globals and model-owned state synchronized while remaining modules
    // still read/write the legacy global storage.
    profile.model_state.legacy.daynum = profile.model_state.legacy.day_start - 1;
    profile.model_state.legacy.b_end = false;
    profile.model_state.legacy.write_to_globals();

    let mut days_simulated = 0usize;
    for _ in profile.model_state.legacy.day_start..(profile.model_state.legacy.day_finish + 1) {
        if should_cancel() {
            return Ok(ExecutionOutcome {
                days_simulated,
                cancelled: true,
            });
        }

        let mut state = if !profile.states.is_empty() {
            let mut new_state = profile.states.last().unwrap().clone();
            new_state.date = new_state.date.succ_opt().unwrap();
            new_state
        } else {
            State::new(
                profile,
                NaiveDate::from_yo_opt(
                    profile.model_state.legacy.iyear,
                    profile.model_state.legacy.day_start as u32,
                )
                .unwrap(),
            )
        };

        // Execute simulation for this day.
        let mut model_state = std::mem::take(&mut profile.model_state);
        let simulation_result = state.simulate_this_day(profile, &mut model_state);
        profile.model_state = model_state;
        profile.model_state.legacy.read_from_globals();

        if let Err(err) = simulation_result {
            if err.level() == 0 {
                println!("{}", err.message());
                break;
            }
            return Err(Box::new(err));
        }

        profile.write_record()?;
        let current_date = state.date;
        profile.states.push(state);
        days_simulated += 1;
        on_progress(days_simulated, current_date);

        if profile.model_state.legacy.b_end {
            break;
        }
    }

    Ok(ExecutionOutcome {
        days_simulated,
        cancelled: false,
    })
}

pub fn run_job(
    mut request: RunRequest,
    mut on_event: impl FnMut(RunEvent),
) -> Result<RunSummary, RunError> {
    if request.run_id.trim().is_empty() {
        request.run_id = generated_run_id();
    }

    let profile_path = request.profile_path.clone();
    if !profile_path.is_file() {
        return Err(RunError::new(
            RunErrorCode::Input,
            format!("profile file '{}' does not exist", profile_path.display()),
        ));
    }

    let profile_dir = profile_path.parent().ok_or_else(|| {
        RunError::new(
            RunErrorCode::Input,
            format!(
                "profile path '{}' has no parent directory",
                profile_path.display()
            ),
        )
    })?;

    if request.run_dir.exists() && !same_directory(&request.run_dir, profile_dir) {
        return Err(RunError::new(
            RunErrorCode::Input,
            format!(
                "run directory '{}' already exists",
                request.run_dir.display()
            ),
        ));
    }

    fs::create_dir_all(&request.run_dir).map_err(|e| {
        RunError::new(
            RunErrorCode::Io,
            format!(
                "failed to create run directory '{}': {}",
                request.run_dir.display(),
                e
            ),
        )
    })?;

    let (mut profile, profile_contents, profile_extension) = parse_profile_file(&profile_path)?;
    let resolved_paths = configure_profile_paths(&mut profile, &profile_path, &request.run_dir)?;
    persist_run_inputs(&request, &profile_contents, &profile_extension)?;

    let started_at = now_rfc3339();
    let timer = Instant::now();
    on_event(RunEvent::Started {
        run_id: request.run_id.clone(),
        started_at: started_at.clone(),
    });

    let output_csv_path = request.run_dir.join("output.csv");
    let meta_path = request.run_dir.join("meta.json");

    let execution = execute_profile(
        &mut profile,
        || {
            request
                .cancel_file
                .as_ref()
                .map(|path| path.exists())
                .unwrap_or(false)
        },
        |day_index, date| {
            on_event(RunEvent::Progress {
                run_id: request.run_id.clone(),
                day_index,
                date: date.format("%F").to_string(),
            });
        },
    );

    let finished_at = now_rfc3339();
    let duration_ms = timer.elapsed().as_millis() as u64;

    let summary = match execution {
        Ok(outcome) if outcome.cancelled => {
            let failure = RunFailure {
                code: "cancelled".to_string(),
                message: "simulation cancelled by cancel-file signal".to_string(),
            };
            on_event(RunEvent::Failed {
                run_id: request.run_id.clone(),
                code: failure.code.clone(),
                message: failure.message.clone(),
            });
            RunSummary {
                run_id: request.run_id.clone(),
                status: RunStatus::Cancelled,
                started_at,
                finished_at,
                output_csv_path,
                meta_path,
                days_simulated: outcome.days_simulated,
                error: Some(failure),
            }
        }
        Ok(outcome) => {
            on_event(RunEvent::Finished {
                run_id: request.run_id.clone(),
                finished_at: finished_at.clone(),
            });
            RunSummary {
                run_id: request.run_id.clone(),
                status: RunStatus::Succeeded,
                started_at,
                finished_at,
                output_csv_path,
                meta_path,
                days_simulated: outcome.days_simulated,
                error: None,
            }
        }
        Err(err) => {
            let failure = RunFailure {
                code: "simulation_error".to_string(),
                message: err.to_string(),
            };
            on_event(RunEvent::Failed {
                run_id: request.run_id.clone(),
                code: failure.code.clone(),
                message: failure.message.clone(),
            });
            RunSummary {
                run_id: request.run_id.clone(),
                status: RunStatus::Failed,
                started_at,
                finished_at,
                output_csv_path,
                meta_path,
                days_simulated: 0,
                error: Some(failure),
            }
        }
    };

    write_metadata(&request, &resolved_paths, &summary, duration_ms)?;
    Ok(summary)
}
