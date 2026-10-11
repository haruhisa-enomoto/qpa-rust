#!/bin/sh

set -eu

script_dir=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH='' cd -- "$script_dir/../.." && pwd)

actual=$(
    "$repo_root/scripts/gap-reference" \
        "$repo_root/tests/gap-reference/version_guard.g"
)

for expected in \
    '"gap_version":"4.17dev"' \
    '"gap_revision":"4efe877859ab0bc0ac041bf635a996fb897ee0c9"' \
    '"qpa_version":"1.37"' \
    '"qpa_revision":"8d5f847258d4853700ce3a5d947db2a1441109f1"' \
    '"gbnp_version":"1.1.0"' \
    '"gbnp_revision":"e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"'
do
    case "$actual" in
        *"$expected"*) ;;
        *)
            echo "oracle metadata is missing $expected: $actual" >&2
            exit 1
            ;;
    esac
done

printf '%s\n' "$actual"
