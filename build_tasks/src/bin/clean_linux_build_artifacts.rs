

use std::process::Command;

use build_tasks::utils::get_linux_build_dir;


fn main() -> Result<(), Box<dyn std::error::Error>> {

    #[cfg(not(any(target_os = "windows", target_os = "linux")))]
    return Err("Cleanup is designed and necessary for Windows WSL Linux builds only and does not support other operating systems.".into());


    #[cfg(target_os = "windows")]
    Command::new("wsl")
        .args(&["-d", "Ubuntu"])
        .args(&["cargo", "cleanlx"])
        .status()
        .map_err(|e| format!("Cleaning linux target directory inside WSL failed: {e:?}"))?;


    // #[cfg(not(target_os = "linux"))]
    // return Err("Cleanup is designed and necessary for WSL Linux builds only.".into());


    //todo!("Test correctness");
    let lx_build_dir = get_linux_build_dir()?;


    let cargo_clean = Command::new("cargo")
        .args(&["clean", "--target-dir", &lx_build_dir])
        .status()
        .map_err(|e| format!("Cleaning target directory {lx_build_dir:?} failed: {e:?}"))?;

    if cargo_clean.success() {
        println!("Successfully cleaned linux target direcoty.");
        return Ok(());
    } else {
        if let Some(err_code) = cargo_clean.code() {
            return Err(format!("Cleaning linux target directory failed with code {err_code}.").into());
        } else {
            return Err("Cleaning linux target directory failed without status code".into());
        }
    }
}



