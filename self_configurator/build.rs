

fn main() {
    println!("cargo:rerun-if-changed=build.rs"); // avoid rerunning unnecessarily

    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let host = std::env::var("HOST").unwrap_or_default();
    let target = std::env::var("TARGET").unwrap_or_default();

    println!("Build script: target OS = {}", target_os);
    println!("Build script: host triple = {}", host);
    println!("Build script: target triple = {}", target);
}


