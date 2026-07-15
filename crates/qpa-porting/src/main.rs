//! Generate the authoritative QPA declaration and installed-method inventory.

use std::{
    collections::BTreeMap,
    env,
    fmt::Write as _,
    fs,
    io::{self, Write as _},
    path::{Path, PathBuf},
};

const CONSTRUCTS: &[&str] = &[
    "BindGlobal",
    "DeclareAttribute",
    "DeclareCategory",
    "DeclareCategoryCollections",
    "DeclareCategoryFamily",
    "DeclareFilter",
    "DeclareGlobalFunction",
    "DeclareHandlingByNiceBasis",
    "DeclareInfoClass",
    "DeclareOperation",
    "DeclareOperations",
    "DeclareProperty",
    "DeclareRepresentation",
    "DeclareSynonym",
    "DeclareSynonymAttr",
    "InstallGlobalFunction",
    "InstallHandlingByNiceBasis",
    "InstallImmediateMethod",
    "InstallMethod",
    "InstallOtherMethod",
    "InstallTrueMethod",
];

#[derive(Debug)]
struct Entry {
    construct: Box<str>,
    name: Box<str>,
    source_file: Box<str>,
    source_start: usize,
    source_end: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = env::args_os().skip(1);
    let source_root = PathBuf::from(
        arguments
            .next()
            .unwrap_or_else(|| "upstream/qpa/lib".into()),
    );
    let output = arguments.next().map(PathBuf::from);
    if arguments.next().is_some() {
        return Err("usage: qpa-porting [source-root] [output-file]".into());
    }

    let mut files = Vec::new();
    collect_gap_files(&source_root, &mut files)?;
    files.sort();
    let mut entries = Vec::new();
    for file in files {
        entries.extend(scan_file(&file)?);
    }

    let rendered = render_ledger(&source_root, &entries)?;
    if let Some(output) = output {
        fs::write(output, rendered)?;
    } else {
        io::stdout().write_all(rendered.as_bytes())?;
    }
    Ok(())
}

fn collect_gap_files(directory: &Path, result: &mut Vec<PathBuf>) -> io::Result<()> {
    for entry in fs::read_dir(directory)? {
        let path = entry?.path();
        if path.is_dir() {
            collect_gap_files(&path, result)?;
        } else if matches!(
            path.extension().and_then(|value| value.to_str()),
            Some("gd" | "gi")
        ) {
            result.push(path);
        }
    }
    Ok(())
}

fn scan_file(path: &Path) -> io::Result<Vec<Entry>> {
    let source = fs::read_to_string(path)?;
    let masked = mask_non_code(&source);
    let bytes = masked.as_bytes();
    let newlines = bytes
        .iter()
        .enumerate()
        .filter_map(|(position, byte)| (*byte == b'\n').then_some(position))
        .collect::<Vec<_>>();
    let source_file = path.to_string_lossy().replace('\\', "/").into_boxed_str();
    let mut entries = Vec::new();
    let mut position = 0;

    while position < bytes.len() {
        if !is_identifier_start(bytes[position]) {
            position += 1;
            continue;
        }
        let start = position;
        position += 1;
        while position < bytes.len() && is_identifier_continue(bytes[position]) {
            position += 1;
        }
        let identifier = &masked[start..position];
        if !CONSTRUCTS.contains(&identifier) {
            continue;
        }
        let mut open = position;
        while open < bytes.len() && bytes[open].is_ascii_whitespace() {
            open += 1;
        }
        if bytes.get(open) != Some(&b'(') {
            continue;
        }
        let Some(close) = matching_parenthesis(bytes, open) else {
            continue;
        };
        let first_end = first_argument_end(bytes, open + 1, close);
        let raw_name = source[open + 1..first_end].trim();
        let name = parse_gap_name(raw_name);
        entries.push(Entry {
            construct: identifier.into(),
            name,
            source_file: source_file.clone(),
            source_start: line_at(&newlines, start),
            source_end: line_at(&newlines, close),
        });
    }
    Ok(entries)
}

fn mask_non_code(source: &str) -> String {
    let bytes = source.as_bytes();
    let mut masked = bytes.to_vec();
    let mut position = 0;
    while position < bytes.len() {
        match bytes[position] {
            b'#' => {
                while position < bytes.len() && bytes[position] != b'\n' {
                    masked[position] = b' ';
                    position += 1;
                }
            }
            quote @ (b'"' | b'\'') => {
                masked[position] = b' ';
                position += 1;
                while position < bytes.len() {
                    let byte = bytes[position];
                    masked[position] = if byte == b'\n' { b'\n' } else { b' ' };
                    position += 1;
                    if byte == b'\\' && position < bytes.len() {
                        masked[position] = if bytes[position] == b'\n' {
                            b'\n'
                        } else {
                            b' '
                        };
                        position += 1;
                    } else if byte == quote {
                        break;
                    }
                }
            }
            _ => position += 1,
        }
    }
    String::from_utf8(masked).expect("masking preserves UTF-8 byte boundaries")
}

fn matching_parenthesis(bytes: &[u8], open: usize) -> Option<usize> {
    let mut depth = 0_usize;
    for (offset, byte) in bytes[open..].iter().enumerate() {
        match byte {
            b'(' => depth += 1,
            b')' => {
                depth = depth.checked_sub(1)?;
                if depth == 0 {
                    return Some(open + offset);
                }
            }
            _ => {}
        }
    }
    None
}

fn first_argument_end(bytes: &[u8], start: usize, close: usize) -> usize {
    let mut round = 0_usize;
    let mut square = 0_usize;
    let mut curly = 0_usize;
    for (position, &byte) in bytes.iter().enumerate().take(close).skip(start) {
        match byte {
            b'(' => round += 1,
            b')' => round = round.saturating_sub(1),
            b'[' => square += 1,
            b']' => square = square.saturating_sub(1),
            b'{' => curly += 1,
            b'}' => curly = curly.saturating_sub(1),
            b',' if round == 0 && square == 0 && curly == 0 => return position,
            _ => {}
        }
    }
    close
}

fn parse_gap_name(raw: &str) -> Box<str> {
    if raw.starts_with('"')
        && raw.ends_with('"')
        && let Ok(value) = serde_json::from_str::<String>(raw)
    {
        return value.into_boxed_str();
    }
    raw.split_whitespace().next().unwrap_or("<unknown>").into()
}

fn line_at(newlines: &[usize], position: usize) -> usize {
    newlines.partition_point(|newline| *newline < position) + 1
}

fn is_identifier_start(byte: u8) -> bool {
    byte.is_ascii_alphabetic() || byte == b'_'
}

fn is_identifier_continue(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || byte == b'_'
}

fn render_ledger(
    source_root: &Path,
    entries: &[Entry],
) -> Result<String, Box<dyn std::error::Error>> {
    let mut result = String::new();
    result.push_str("schema_version = 1\n");
    result.push_str("repository = \"qpa\"\n");
    result.push_str("revision = \"8d5f847258d4853700ce3a5d947db2a1441109f1\"\n");
    writeln!(
        result,
        "source_root = {}",
        serde_json::to_string(&source_root.to_string_lossy())?
    )?;
    writeln!(result, "entry_count = {}", entries.len())?;

    let mut occurrences = BTreeMap::<(Box<str>, Box<str>, Box<str>), usize>::new();
    for entry in entries {
        let key = (
            entry.source_file.clone(),
            entry.construct.clone(),
            entry.name.clone(),
        );
        let occurrence = occurrences.entry(key).or_default();
        *occurrence += 1;
        let file_stem = Path::new(entry.source_file.as_ref())
            .file_stem()
            .and_then(|value| value.to_str())
            .unwrap_or("unknown");
        let identifier = format!(
            "{file_stem}:{}:{}:{}",
            entry.construct, entry.name, occurrence
        );

        result.push_str("\n[[entries]]\n");
        push_string(&mut result, "id", &identifier)?;
        push_string(&mut result, "construct", &entry.construct)?;
        push_string(&mut result, "name", &entry.name)?;
        push_string(&mut result, "source_file", &entry.source_file)?;
        writeln!(result, "source_start = {}", entry.source_start)?;
        writeln!(result, "source_end = {}", entry.source_end)?;
        push_string(
            &mut result,
            "classification",
            if entry.construct.starts_with("Declare") {
                "declared_surface"
            } else {
                "installed_implementation"
            },
        )?;
        push_string(&mut result, "status", "inventory")?;
        push_string(&mut result, "rust_symbol", "")?;
        result.push_str("tests = []\n");
        result.push_str("deviations = []\n");
    }
    Ok(result)
}

fn push_string(result: &mut String, key: &str, value: &str) -> Result<(), serde_json::Error> {
    result.push_str(key);
    result.push_str(" = ");
    result.push_str(&serde_json::to_string(value)?);
    result.push('\n');
    Ok(())
}
