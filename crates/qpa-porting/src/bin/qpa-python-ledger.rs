//! Generate the operation-oriented Python exposure ledger from the QPA source ledger.

use std::{
    collections::{BTreeMap, BTreeSet},
    env,
    fmt::Write as _,
    fs,
    io::{self, Write as _},
    path::{Path, PathBuf},
};

#[derive(Debug, Default)]
struct SourceEntry {
    id: Box<str>,
    construct: Box<str>,
    name: Box<str>,
    source_file: Box<str>,
    rust_symbol: Box<str>,
    tests: Vec<Box<str>>,
    deviations: Vec<Box<str>>,
}

#[derive(Debug, Default)]
struct Operation {
    entry_ids: Vec<Box<str>>,
    constructs: BTreeSet<Box<str>>,
    source_files: BTreeSet<Box<str>>,
    rust_symbols: BTreeSet<Box<str>>,
    rust_tests: BTreeSet<Box<str>>,
    deviations: BTreeSet<Box<str>>,
}

#[derive(Debug, Default)]
struct Exposure {
    status: Box<str>,
    target: Box<str>,
    tests: Vec<Box<str>>,
    notes: Box<str>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = env::args_os().skip(1);
    let source = PathBuf::from(
        arguments
            .next()
            .unwrap_or_else(|| "porting/qpa-operations.toml".into()),
    );
    let exposure = PathBuf::from(
        arguments
            .next()
            .unwrap_or_else(|| "python/bindings.toml".into()),
    );
    let output = arguments.next().map(PathBuf::from);
    if arguments.next().is_some() {
        return Err(
            "usage: qpa-python-ledger [qpa-ledger] [exposure-manifest] [output-file]".into(),
        );
    }

    let entries = read_source_entries(&source)?;
    let exposures = read_exposures(&exposure)?;
    let rendered = render(&source, &exposure, &entries, &exposures)?;
    if let Some(output) = output {
        fs::write(output, rendered)?;
    } else {
        io::stdout().write_all(rendered.as_bytes())?;
    }
    Ok(())
}

fn read_source_entries(path: &Path) -> Result<Vec<SourceEntry>, Box<dyn std::error::Error>> {
    let source = fs::read_to_string(path)?;
    let mut result = Vec::new();
    let mut entry: Option<SourceEntry> = None;
    for line in source.lines() {
        if line == "[[entries]]" {
            if let Some(entry) = entry.take() {
                result.push(entry);
            }
            entry = Some(SourceEntry::default());
            continue;
        }
        let Some(entry) = entry.as_mut() else {
            continue;
        };
        if let Some(value) = string_value(line, "id")? {
            entry.id = value;
        } else if let Some(value) = string_value(line, "construct")? {
            entry.construct = value;
        } else if let Some(value) = string_value(line, "name")? {
            entry.name = value;
        } else if let Some(value) = string_value(line, "source_file")? {
            entry.source_file = value;
        } else if let Some(value) = string_value(line, "rust_symbol")? {
            entry.rust_symbol = value;
        } else if let Some(value) = string_array_value(line, "tests")? {
            entry.tests = value;
        } else if let Some(value) = string_array_value(line, "deviations")? {
            entry.deviations = value;
        }
    }
    if let Some(entry) = entry {
        result.push(entry);
    }
    if result.iter().any(|entry| {
        entry.id.is_empty()
            || entry.construct.is_empty()
            || entry.name.is_empty()
            || entry.source_file.is_empty()
    }) {
        return Err("QPA ledger contains an incomplete entry".into());
    }
    Ok(result)
}

fn read_exposures(path: &Path) -> Result<BTreeMap<Box<str>, Exposure>, Box<dyn std::error::Error>> {
    let source = fs::read_to_string(path)?;
    let mut result = BTreeMap::new();
    let mut name: Option<Box<str>> = None;
    let mut exposure: Option<Exposure> = None;
    for line in source.lines() {
        if line == "[[bindings]]" {
            finish_exposure(&mut result, name.take(), exposure.take())?;
            exposure = Some(Exposure::default());
            continue;
        }
        let Some(current) = exposure.as_mut() else {
            continue;
        };
        if let Some(value) = string_value(line, "qpa_name")? {
            name = Some(value);
        } else if let Some(value) = string_value(line, "status")? {
            current.status = value;
        } else if let Some(value) = string_value(line, "target")? {
            current.target = value;
        } else if let Some(value) = string_array_value(line, "tests")? {
            current.tests = value;
        } else if let Some(value) = string_value(line, "notes")? {
            current.notes = value;
        }
    }
    finish_exposure(&mut result, name, exposure)?;
    Ok(result)
}

fn finish_exposure(
    result: &mut BTreeMap<Box<str>, Exposure>,
    name: Option<Box<str>>,
    exposure: Option<Exposure>,
) -> Result<(), Box<dyn std::error::Error>> {
    let Some(exposure) = exposure else {
        return Ok(());
    };
    let name = name.ok_or("Python exposure entry has no qpa_name")?;
    match exposure.status.as_ref() {
        "partial" | "exposed" | "represented" => {
            if exposure.target.is_empty() || exposure.tests.is_empty() {
                return Err(format!("Python exposure for {name} is incomplete").into());
            }
        }
        "unavailable" => {
            if exposure.notes.is_empty() {
                return Err(format!("unavailable Python exposure {name} has no reason").into());
            }
        }
        status => {
            return Err(format!("invalid Python exposure status {status:?} for {name}").into());
        }
    }
    if result.insert(name.clone(), exposure).is_some() {
        return Err(format!("duplicate Python exposure for {name}").into());
    }
    Ok(())
}

fn string_value(line: &str, key: &str) -> Result<Option<Box<str>>, Box<dyn std::error::Error>> {
    let Some(value) = line.strip_prefix(&format!("{key} = ")) else {
        return Ok(None);
    };
    Ok(Some(
        serde_json::from_str::<String>(value)?.into_boxed_str(),
    ))
}

fn string_array_value(
    line: &str,
    key: &str,
) -> Result<Option<Vec<Box<str>>>, Box<dyn std::error::Error>> {
    let Some(value) = line.strip_prefix(&format!("{key} = ")) else {
        return Ok(None);
    };
    Ok(Some(
        serde_json::from_str::<Vec<String>>(value)?
            .into_iter()
            .map(String::into_boxed_str)
            .collect(),
    ))
}

fn render(
    source: &Path,
    exposure_path: &Path,
    entries: &[SourceEntry],
    exposures: &BTreeMap<Box<str>, Exposure>,
) -> Result<String, Box<dyn std::error::Error>> {
    let mut operations = BTreeMap::<Box<str>, Operation>::new();
    let mut seen_ids = BTreeSet::new();
    for entry in entries {
        if !seen_ids.insert(entry.id.clone()) {
            return Err(format!("duplicate QPA source entry {}", entry.id).into());
        }
        let operation = operations.entry(entry.name.clone()).or_default();
        operation.entry_ids.push(entry.id.clone());
        operation.constructs.insert(entry.construct.clone());
        operation.source_files.insert(entry.source_file.clone());
        if !entry.rust_symbol.is_empty() {
            operation.rust_symbols.insert(entry.rust_symbol.clone());
        }
        operation.rust_tests.extend(entry.tests.iter().cloned());
        operation
            .deviations
            .extend(entry.deviations.iter().cloned());
    }
    for name in exposures.keys() {
        if !operations.contains_key(name) {
            return Err(format!("Python exposure names unknown QPA operation {name}").into());
        }
    }

    let exposed = operations
        .keys()
        .filter(|name| exposures.contains_key(*name))
        .count();
    let mut result = String::new();
    result.push_str("schema_version = 1\n");
    push_string(&mut result, "qpa_ledger", &source.to_string_lossy())?;
    push_string(
        &mut result,
        "exposure_manifest",
        &exposure_path.to_string_lossy(),
    )?;
    writeln!(result, "source_entry_count = {}", entries.len())?;
    writeln!(result, "accounted_entry_count = {}", seen_ids.len())?;
    writeln!(result, "operation_count = {}", operations.len())?;
    writeln!(result, "exposed_or_partial_count = {exposed}")?;
    result.push_str(
        "coverage_unit = \"one QPA name with every declaration and installed overload grouped\"\n",
    );

    for (name, operation) in operations {
        let disposition = disposition(&operation.constructs);
        let exposure = exposures.get(&name);
        result.push_str("\n[[operations]]\n");
        push_string(&mut result, "qpa_name", &name)?;
        push_string(&mut result, "python_name", &python_name(&name))?;
        push_string(&mut result, "disposition", disposition)?;
        push_string(
            &mut result,
            "status",
            exposure.map_or("planned", |value| value.status.as_ref()),
        )?;
        push_string(
            &mut result,
            "python_target",
            exposure.map_or("", |value| value.target.as_ref()),
        )?;
        push_string(
            &mut result,
            "notes",
            exposure.map_or("", |value| value.notes.as_ref()),
        )?;
        push_array(&mut result, "entry_ids", operation.entry_ids.iter())?;
        push_array(&mut result, "constructs", operation.constructs.iter())?;
        push_array(&mut result, "source_files", operation.source_files.iter())?;
        push_array(&mut result, "rust_symbols", operation.rust_symbols.iter())?;
        push_array(&mut result, "rust_tests", operation.rust_tests.iter())?;
        push_array(&mut result, "deviations", operation.deviations.iter())?;
        push_array(
            &mut result,
            "python_tests",
            exposure.into_iter().flat_map(|value| value.tests.iter()),
        )?;
    }
    Ok(result)
}

fn disposition(constructs: &BTreeSet<Box<str>>) -> &'static str {
    if constructs.iter().any(|construct| {
        matches!(
            construct.as_ref(),
            "DeclareCategory"
                | "DeclareCategoryCollections"
                | "DeclareCategoryFamily"
                | "DeclareFilter"
                | "DeclareRepresentation"
        )
    }) {
        "runtime_type_or_protocol"
    } else if constructs.contains("DeclareInfoClass") {
        "diagnostic_channel"
    } else if constructs
        .iter()
        .all(|construct| matches!(construct.as_ref(), "DeclareSynonym" | "DeclareSynonymAttr"))
    {
        "alias"
    } else {
        "operation"
    }
}

fn python_name(name: &str) -> String {
    let unquoted = name.trim_matches('"');
    let operator = match unquoted {
        "\\[\\]" => Some("index"),
        "\\^" => Some("power"),
        "\\*" => Some("multiply"),
        "\\+" => Some("add"),
        "\\-" => Some("subtract"),
        "\\." => Some("dot"),
        "\\/" => Some("divide"),
        "\\<" => Some("less_than"),
        "\\=" => Some("equal"),
        _ => None,
    };
    if let Some(operator) = operator {
        return operator.to_owned();
    }
    let name = name
        .strip_prefix("1st")
        .map_or_else(|| name.to_owned(), |rest| format!("First{rest}"));
    let characters = name.chars().collect::<Vec<_>>();
    let mut result = String::new();
    for (index, character) in characters.iter().copied().enumerate() {
        if !character.is_ascii_alphanumeric() {
            if !result.ends_with('_') {
                result.push('_');
            }
            continue;
        }
        let previous = index.checked_sub(1).and_then(|value| characters.get(value));
        let next = characters.get(index + 1);
        if character.is_ascii_uppercase()
            && !result.is_empty()
            && !result.ends_with('_')
            && (previous.is_some_and(|value| value.is_ascii_lowercase() || value.is_ascii_digit())
                || next.is_some_and(char::is_ascii_lowercase))
        {
            result.push('_');
        }
        result.push(character.to_ascii_lowercase());
    }
    result.trim_matches('_').to_owned()
}

fn push_string(result: &mut String, key: &str, value: &str) -> Result<(), serde_json::Error> {
    result.push_str(key);
    result.push_str(" = ");
    result.push_str(&serde_json::to_string(value)?);
    result.push('\n');
    Ok(())
}

fn push_array<'a>(
    result: &mut String,
    key: &str,
    values: impl IntoIterator<Item = &'a Box<str>>,
) -> Result<(), serde_json::Error> {
    let values = values.into_iter().map(Box::as_ref).collect::<Vec<_>>();
    result.push_str(key);
    result.push_str(" = ");
    result.push_str(&serde_json::to_string(&values)?);
    result.push('\n');
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::python_name;

    #[test]
    fn qpa_names_become_stable_python_names() {
        assert_eq!(python_name("PathAlgebra"), "path_algebra");
        assert_eq!(python_name("GBNPGroebnerBasis"), "gbnp_groebner_basis");
        assert_eq!(python_name("1stSyzygy"), "first_syzygy");
        assert_eq!(python_name("IsBound\\[\\]"), "is_bound");
        assert_eq!(python_name("\\[\\]"), "index");
        assert_eq!(python_name("\"\\^\""), "power");
    }
}
