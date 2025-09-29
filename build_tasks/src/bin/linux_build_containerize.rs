


use std::process::Command;
use std::path::{Path, PathBuf};
use std::env;

use build_tasks::utils::{get_linux_build_dir, effective_target_dir, workspace_root};



fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(not(target_os = "linux"))]
    return Err("Singularity container builds work only on Linux systems.".into());

    let lx_build_dir = get_linux_build_dir()?;

    // let self_config_build = Command::new("cargo")
    //                             .args(&["build", "-p", "self_configurator", "--target-dir", &lx_build_dir, "--release"])
    //                             .status()
    //                             .map_err(|e| format!("Building self configurator failed: {e:?}"))?;

    let self_config_build = Command::new("dx")
                                .args(&["build", "-p", "self_configurator", "--release"])
                                .status()
                                .map_err(|e| format!("Building self configurator failed: {e:?}"))?;

    if self_config_build.success() {
        println!("Successfully built 'self_configurator' binary.");
    } else {
        if let Some(err_code) = self_config_build.code() {
            return Err(format!("Build of 'self_configurator' binary failed with code {err_code}.").into());
        } else {
            return Err("Build of 'self_configurator' binary failed without status code".into());
        }
    }


    let project_dir = workspace_root();

    //todo!("Creating the container");

    let container_build_dir = project_dir.join("container");

    let wf_deps_sif = container_build_dir.join("wf_deps.sif");

    if !wf_deps_sif.is_file() {
        let wf_deps_build = Command::new("sudo")
                                .current_dir(&container_build_dir)
                                .args(&["singularity", "build", "wf_deps.sif", "wf_deps.def"])
                                .status()
                                .map_err(|e| format!("Containerization step failed: {e:?}"))?;

        if wf_deps_build.success() {
            println!("Successfully built container 'wf_deps.sif'.");
        } else {
            if let Some(err_code) = wf_deps_build.code() {
                return Err(format!("Build of container 'wf_deps.sif' failed with code {err_code}.").into());
            } else {
                return Err("Build of container 'wf_deps.sif' failed without status code".into());
            }
        }
    }


    let wf_apps_build = Command::new("sudo")
                            .current_dir(&container_build_dir)
                            .args(&["singularity", "build", "wf_apps_auto.sif", "wf_apps.def"])
                            .status()
                            .map_err(|e| format!("Containerization step failed: {e:?}"))?;

    if wf_apps_build.success() {
        println!("Successfully built container 'wf_apps_auto.sif'.");
    } else {
        if let Some(err_code) = self_config_build.code() {
            return Err(format!("Build of container failed with code {err_code}.").into());
        } else {
            return Err("Build of cotnaienr failed without status code".into());
        }
    }

    return Ok(());
}






