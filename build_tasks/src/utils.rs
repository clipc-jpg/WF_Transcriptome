

use std::path::{Path, PathBuf};
use std::process::Command;

use dirs;
use serde_json;


// rust does not canonicalize '~' in Paths, since '~' may be the actual name of a directory or file
pub fn expand_tilde<P: AsRef<Path>>(path: P) -> PathBuf {
    let p = path.as_ref();

    if let Some(str_path) = p.to_str() {
        if str_path.starts_with("~") {
            if str_path == "~" || str_path.starts_with("~/") {
                if let Some(home_dir) = dirs::home_dir() {
                    return home_dir.join(&str_path[2..]);
                }
            }
        }
    }

    p.to_path_buf()
}

pub fn get_linux_build_dir() -> Result<String, Box<dyn std::error::Error>>{
    let lx_build_dir = std::env::var("LINUX_TARGET_DIR")
        .map_err(|_e| "Mandatory Environment variable not set: 'LINUX_TARGET_DIR'")?;

    let lx_build_dir = expand_tilde(&lx_build_dir)
        .canonicalize()
        .map_err(|e| format!("Canonicalization failed for '{lx_build_dir}': {e:?}"))?;

    let lx_build_dir = lx_build_dir.into_os_string().into_string()
        .map_err(|e| format!("String representation of canonicalized path failed. Path: {e:?}"))?;

    return Ok(lx_build_dir);
}



// dioxus does not expose a --target-dir argument, so always the default directory will be used
pub fn cargo_default_target_dir() -> String {
    let output = Command::new("cargo")
    .args(["metadata", "--format-version", "1", "--no-deps"])
    .output()
    .expect("Failed to run cargo metadata");

    let json: serde_json::Value =
    serde_json::from_slice(&output.stdout).expect("Invalid JSON from cargo metadata");

    json["target_directory"].as_str().unwrap().to_string()
}

pub fn effective_target_dir() -> String {
    // 1. If CARGO_TARGET_DIR is set, that one will be used
    if let Ok(env_dir) = std::env::var("CARGO_TARGET_DIR") {
        return env_dir;
    }

    // 2. Otherwise, cargo metadata will be used
    cargo_default_target_dir()
}


// for finding the singularity container def file
pub fn workspace_root() -> PathBuf {

    let cargo_toml_path = Command::new("cargo")
                            .args(["locate-project", "--workspace", "--message-format=plain"])
                            .output()
                            .expect("Failed to run cargo locate-project");

    let path_str = String::from_utf8(cargo_toml_path.stdout).expect("Cargo locate-project emitted invalid UTF-8");
    let mut path = PathBuf::from(path_str.trim()); path.pop(); // cargo toml always exists and is a file

    return path;
}

