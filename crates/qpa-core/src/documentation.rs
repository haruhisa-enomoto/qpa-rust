//! QPA documentation-build metadata.
//!
//! QPA's `MakeQPADocumentation` is a thin `GAPDoc` build helper rather than a
//! representation-theoretic computation. The Rust port exposes the identical
//! inputs as a plan so that a documentation frontend can choose how to invoke
//! `GAPDoc`-compatible tooling without a library call writing process-global
//! files.

use std::path::{Path, PathBuf};

/// Fixed inputs passed by QPA's `MakeQPADocumentation` to `GAPDoc`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QpaDocumentationPlan {
    directory: PathBuf,
    main_file: &'static str,
    additional_files: Vec<PathBuf>,
    book_name: &'static str,
    renderer: &'static str,
    copy_html_style_files: bool,
}

impl QpaDocumentationPlan {
    /// Return QPA's package `doc` directory.
    #[must_use]
    pub fn directory(&self) -> &Path {
        &self.directory
    }

    /// Return QPA's main `GAPDoc` XML filename.
    #[must_use]
    pub const fn main_file(&self) -> &'static str {
        self.main_file
    }

    /// Return QPA's explicitly supplied auxiliary-file list.
    #[must_use]
    pub fn additional_files(&self) -> &[PathBuf] {
        &self.additional_files
    }

    /// Return QPA's `GAPDoc` book name.
    #[must_use]
    pub const fn book_name(&self) -> &'static str {
        self.book_name
    }

    /// Return QPA's selected mathematics renderer.
    #[must_use]
    pub const fn renderer(&self) -> &'static str {
        self.renderer
    }

    /// Return whether QPA copies the GAP HTML style files after the build.
    #[must_use]
    pub const fn copy_html_style_files(&self) -> bool {
        self.copy_html_style_files
    }
}

/// Translate the fixed build inputs of QPA's `MakeQPADocumentation`.
///
/// `package_root` corresponds to GAP's QPA package directory; QPA selects its
/// `doc` child, builds `qpadocumentation.xml` as book `QPA` using `MathJax`,
/// and then copies the HTML style files.
#[must_use]
pub fn make_qpa_documentation(package_root: impl AsRef<Path>) -> QpaDocumentationPlan {
    QpaDocumentationPlan {
        directory: package_root.as_ref().join("doc"),
        main_file: "qpadocumentation.xml",
        additional_files: Vec::new(),
        book_name: "QPA",
        renderer: "MathJax",
        copy_html_style_files: true,
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use serde_json::{Value, json};

    use super::make_qpa_documentation;

    #[test]
    fn pinned_qpa_documentation_plan_matches() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/final_declared_surface.schema-v1.json"
        ))
        .unwrap();
        let plan = make_qpa_documentation("/opt/QPA");
        assert_eq!(plan.directory(), Path::new("/opt/QPA/doc"));
        assert_eq!(
            json!({
                "function_bound": true,
                "directory_count": 1,
                // This is an oracle fact about the pinned QPA package, not a
                // requirement that distributable Rust tests contain its tree.
                "main_exists": true,
                "main": plan.main_file(),
                "files": plan.additional_files(),
                "book_name": plan.book_name(),
                "renderer": plan.renderer(),
                "copy_html_style_files": plan.copy_html_style_files(),
            }),
            fixture["documentation"]
        );
    }
}
