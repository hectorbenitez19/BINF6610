
#!/usr/bin/env bash
set -euo pipefail

SHEET=${1:?Missing samplesheet}
REF=${2:?Missing reference}
REGION=${3:?Missing region}

errors=0
declare -A seen=()

while IFS=, read -r id cond rep lt r1 r2 rest; do
    [[ "$id" == "sample_id" ]] && continue
    [[ -z "$id" ]] && continue

    r1=$(basename "$r1")

    if [[ -n "$r2" ]]; then
        r2=$(basename "$r2")
    fi

    if [[ -z "$id" || ! "$id" =~ ^[A-Za-z0-9._\ -]+$ ]]; then
        echo "Invalid sample ID: $id" >&2
        ((errors+=1))
    fi

    if [[ -n "${seen[$id]:-}" ]]; then
        echo "Duplicate sample ID: $id" >&2
        ((errors+=1))
    fi
    seen[$id]=1

    if [[ "$lt" != "paired" && "$lt" != "single" ]]; then
        echo "$id: invalid library type: $lt" >&2
        ((errors+=1))
    fi

    if [[ ! -s "$r1" ]]; then
        echo "$id: missing R1: $r1" >&2
        ((errors+=1))
    else
        if ! gzip -t "$r1"; then
            echo "$id: invalid R1 gzip" >&2
            ((errors+=1))
        fi

        lines=$(gzip -dc "$r1" | wc -l)
        if (( lines % 4 != 0 )); then
            echo "$id: R1 FASTQ line count not divisible by 4" >&2
            ((errors+=1))
        fi
    fi

    if [[ "$lt" == "paired" ]]; then
        if [[ ! -s "$r2" ]]; then
            echo "$id: missing R2: $r2" >&2
            ((errors+=1))
        else
            if ! gzip -t "$r2"; then
                echo "$id: invalid R2 gzip" >&2
                ((errors+=1))
            fi

            lines=$(gzip -dc "$r2" | wc -l)
            if (( lines % 4 != 0 )); then
                echo "$id: R2 FASTQ line count not divisible by 4" >&2
                ((errors+=1))
            fi
        fi
    fi

done < "$SHEET"

[[ -s "$REF" ]] || {
    echo "Missing reference: $REF" >&2
    ((errors+=1))
}

[[ -s "${REF}.fai" ]] || {
    echo "Missing FASTA index" >&2
    ((errors+=1))
}

dict="${REF%.*}.dict"
[[ -s "$dict" ]] || {
    echo "Missing reference dictionary" >&2
    ((errors+=1))
}

if (( errors > 0 )); then
    echo "Validation failed with $errors error(s)" >&2
    exit 65
fi

echo "Validation passed"
