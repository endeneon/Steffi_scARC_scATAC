#!/usr/bin/env bash
# Keep the first 4 columns plus all columns whose header ends with "_GT".
set -euo pipefail

usage() {
	cat <<EOF
Usage: $(basename "$0") -i <input.tsv> [-o <output.tsv>]

  -i, --input    Input TSV file (required)
  -o, --output   Output TSV file (default: <input basename>_GT_only.tsv)
  -h, --help     Show this help
EOF
}

input=""
output=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-i | --input)
		input="${2:-}"
		shift 2
		;;
	-o | --output)
		output="${2:-}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "Unknown argument: $1" >&2
		usage >&2
		exit 1
		;;
	esac
done

if [[ -z "$input" ]]; then
	echo "Error: -i/--input is required" >&2
	usage >&2
	exit 1
fi
if [[ ! -f "$input" ]]; then
	echo "Error: input file not found: $input" >&2
	exit 1
fi

if [[ -z "$output" ]]; then
	output="${input%.tsv}_GT_only.tsv"
fi

awk -F'\t' -v OFS='\t' '
NR == 1 {
    n = 0
    for (i = 1; i <= NF; i++) {
        if (i <= 4 || $i ~ /_GT$/) keep[++n] = i
    }
}
{
    line = $(keep[1])
    for (j = 2; j <= n; j++) line = line OFS $(keep[j])
    print line
}
' "$input" >"$output"

echo "Wrote: $output"
