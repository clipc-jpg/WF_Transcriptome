
mod utils;



use std::process::Command;


use build_tasks::utils::get_linux_build_dir;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Call redirected correctly to build_tasks main!");

    #[cfg(target_os = "windows")]
    {
        let lx_build_dir = std::env::var("LINUX_TARGET_DIR")
                                    .map_err(|_e| "Mandatory Environment variable not set: 'LINUX_TARGET_DIR'")?;
        //let lx_build_dir = get_linux_build_dir()?;
        // Build and run the Windows runner
        let status = Command::new("wsl")
            .args(&["-d", "Ubuntu"])
            .args(&["cargo", "run", "--package", "build_tasks", "--bin", "linux_build_containerize", "--target-dir", &lx_build_dir])
            .status()
            .map_err(|e| format!("Failed to build or run build_task binary: {e:?}"))?;
        std::process::exit(status.code().unwrap_or(1));
    }

    #[cfg(target_os = "linux")]
    {
        // Build and run the Linux runner
        let status = Command::new("cargo")
            .args(&["run", "--bin", "linux_runner"])
            .status()
            .expect("failed to run linux_runner");
        std::process::exit(status.code().unwrap_or(1));
    }

}
