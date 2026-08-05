//! Single-profile command-line entry point for the Cotton2K engine.
//!
//! The binary maps its first positional argument to a [`RunRequest`], runs one
//! profile in that profile's directory, and converts the resulting status into
//! a process exit result. Parsing and simulation errors are reported through
//! the standard error return; model output is written by `run_job`.

use cotton2k::{run_job, RunRequest, RunStatus};
use std::io;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let profile_path = std::env::args().nth(1).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "profile file path should be provided",
        )
    })?;
    let profile_path = std::path::PathBuf::from(profile_path);
    let run_dir = profile_path
        .parent()
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "profile path must have a parent directory",
            )
        })?
        .to_path_buf();

    let summary = run_job(
        RunRequest::new(profile_path, run_dir, String::new(), None),
        |_| {},
    )?;

    match summary.status {
        RunStatus::Succeeded => Ok(()),
        RunStatus::Cancelled => {
            Err(io::Error::new(io::ErrorKind::Interrupted, "simulation cancelled").into())
        }
        RunStatus::Failed => {
            let message = summary
                .error
                .as_ref()
                .map(|e| e.message.clone())
                .unwrap_or_else(|| "simulation failed".to_string());
            Err(io::Error::other(message).into())
        }
    }
}
