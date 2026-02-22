#[pyo3::pymodule]
mod cotton2k {
    use cotton2k::{run_job, RunRequest, RunStatus};
    use pyo3::exceptions::{PyRuntimeError, PyValueError};
    use pyo3::prelude::*;

    #[pyfunction]
    fn run(path: &str) -> PyResult<()> {
        let profile_path = std::path::PathBuf::from(path);
        let run_dir = profile_path
            .parent()
            .ok_or_else(|| PyValueError::new_err("profile path must have a parent directory"))?
            .to_path_buf();

        let summary = run_job(
            RunRequest::new(profile_path, run_dir, String::new(), None),
            |_| {},
        )
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

        match summary.status {
            RunStatus::Succeeded => Ok(()),
            RunStatus::Cancelled => Err(PyRuntimeError::new_err("simulation cancelled")),
            RunStatus::Failed => {
                let message = summary
                    .error
                    .as_ref()
                    .map(|e| e.message.clone())
                    .unwrap_or_else(|| "simulation failed".to_string());
                Err(PyRuntimeError::new_err(message))
            }
        }
    }
}
