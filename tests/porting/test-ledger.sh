#!/bin/sh

set -eu

script_dir=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH='' cd -- "$script_dir/../.." && pwd)
actual=$(mktemp "${TMPDIR:-/tmp}/qpa-rust-ledger.XXXXXX")
cleanup() {
    rm -f -- "$actual"
}
trap cleanup EXIT HUP INT TERM

cd "$repo_root"
cargo run --quiet -p qpa-porting -- upstream/qpa/lib "$actual"
diff --unified porting/qpa-operations.toml "$actual"
