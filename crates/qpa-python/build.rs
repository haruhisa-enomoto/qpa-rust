//! Embed the source checkout's Git identity in the Python extension.

use std::{
    env,
    path::{Path, PathBuf},
    process::Command,
};

const UNKNOWN_COMMIT: &str = "unknown";
const SOURCE_PATHS: &[&str] = &[
    "Cargo.toml",
    "Cargo.lock",
    "crates/qpa-core",
    "crates/qpa-gbnp",
    "crates/qpa-python",
    "python/qpa",
];

fn git_output(root: &Path, arguments: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(root)
        .args(arguments)
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    String::from_utf8(output.stdout)
        .ok()
        .map(|value| value.trim().to_owned())
        .filter(|value| !value.is_empty())
}

fn valid_full_object_id(value: &str) -> bool {
    matches!(value.len(), 40 | 64) && value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

fn watch_engine_sources(root: &Path) {
    for relative in SOURCE_PATHS {
        let path = root.join(relative);
        if path.exists() {
            println!("cargo:rerun-if-changed={}", path.display());
        }
    }
}

fn engine_sources_are_dirty(root: &Path) -> bool {
    let mut command = Command::new("git");
    command
        .arg("-C")
        .arg(root)
        .args(["status", "--porcelain=v1", "--untracked-files=all", "--"]);
    command.args(SOURCE_PATHS);
    command
        .output()
        .ok()
        .filter(|output| output.status.success())
        .is_some_and(|output| !output.stdout.is_empty())
}

fn watch_git_identity(root: &Path) {
    let Some(git_dir) = git_output(root, &["rev-parse", "--absolute-git-dir"]) else {
        return;
    };
    let git_dir = PathBuf::from(git_dir);
    println!("cargo:rerun-if-changed={}", git_dir.join("HEAD").display());

    let common_dir = git_output(
        root,
        &["rev-parse", "--path-format=absolute", "--git-common-dir"],
    )
    .map(PathBuf::from)
    .unwrap_or(git_dir);
    if let Some(reference) = git_output(root, &["symbolic-ref", "--quiet", "HEAD"]) {
        println!(
            "cargo:rerun-if-changed={}",
            common_dir.join(reference).display()
        );
    }
    let packed_refs = common_dir.join("packed-refs");
    if packed_refs.is_file() {
        println!("cargo:rerun-if-changed={}", packed_refs.display());
    }
}

fn main() {
    println!("cargo:rerun-if-env-changed=QPA_PORT_GIT_COMMIT");

    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    watch_git_identity(&root);
    watch_engine_sources(&root);

    let base_commit = match env::var("QPA_PORT_GIT_COMMIT") {
        Ok(value) if valid_full_object_id(&value) => value.to_ascii_lowercase(),
        Ok(value) => panic!(
            "QPA_PORT_GIT_COMMIT must be a full 40- or 64-character Git object ID, got {value:?}"
        ),
        Err(_) => git_output(&root, &["rev-parse", "--verify", "HEAD"])
            .filter(|value| valid_full_object_id(value))
            .unwrap_or_else(|| UNKNOWN_COMMIT.to_owned()),
    };
    let commit = if base_commit != UNKNOWN_COMMIT && engine_sources_are_dirty(&root) {
        format!("{base_commit}-dirty")
    } else {
        base_commit
    };
    if commit == UNKNOWN_COMMIT {
        println!(
            "cargo:warning=QPA port Git commit is unavailable; \
             set QPA_PORT_GIT_COMMIT when building outside a Git checkout"
        );
    }
    println!("cargo:rustc-env=QPA_PORT_GIT_COMMIT={commit}");
}
