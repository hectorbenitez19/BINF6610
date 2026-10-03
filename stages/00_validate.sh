# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# 0 · validate — check everything before computing anything.
stage_validate() {
    local id cond rep lt r1 r2 problems=0 n1 n2 r1_ok r2_ok local dict="${REF%.*}.dict"

    while IFS=, read -r id cond rep lt r1 r2; do
        [[ -n "$id" ]] || { log "a row has no sample_id"; problems=$(( problems + 1 )); continue; }

        # the files exist and are not empty
        [[ -s "$r1" ]] || { log "$id: R1 missing or empty: $r1"; problems=$(( problems + 1 )); }
        if [[ "$lt" == paired ]]; then
            [[ -s "$r2" ]] || { log "$id: declared paired but R2 is missing"; problems=$(( problems + 1 )); }
        fi

        # the r2 column and the declared layout agree
        if [[ -z "$r2" && "$lt" == paired ]]; then
            log "$id: r2_fastq is empty but library_type says paired"
            problems=$(( problems + 1 ))
        fi

        # the gzip streams are whole -- both of them. A stream that fails here is
        # not read below: under pipefail, `gzip -dc` on it would end the script
        # instead of letting it report every problem.
        r1_ok=0 r2_ok=0
        if [[ -s "$r1" ]]; then
            if gzip -t "$r1" 2>/dev/null; then r1_ok=1
            else log "$id: R1 is not a valid gzip file"; problems=$(( problems + 1 )); fi
        fi
        if [[ "$lt" == paired && -s "$r2" ]]; then
            if gzip -t "$r2" 2>/dev/null; then r2_ok=1
            else log "$id: R2 is not a valid gzip file"; problems=$(( problems + 1 )); fi
        fi

        # the records are whole, and the mates agree
        if (( r1_ok )); then
            n1=$(gzip -dc "$r1" | wc -l)
            (( n1 % 4 == 0 )) || { log "$id: R1 has $n1 lines, not a whole number of records"
                                   problems=$(( problems + 1 )); }
            if (( r2_ok )); then
                n2=$(gzip -dc "$r2" | wc -l)
                (( n1 == n2 )) || { log "$id: R1 has $(( n1 / 4 )) reads, R2 has $(( n2 / 4 ))"
                                    problems=$(( problems + 1 )); }
            fi
        fi
    done < <(tail -n +2 "$SHEET")

    # duplicate sample ids. sort | uniq -d prints only the repeats.
    local dupes
    dupes=$(awk -F, 'NR>1 { print $1 }' "$SHEET" | sort | uniq -d)
    [[ -z "$dupes" ]] || { log "duplicate sample_id: $dupes"; problems=$(( problems + 1 )); }

    # the reference is where the config says it is
[[ -s "${REF}.fai" ]] || {
    log "reference FASTA index missing: ${REF}.fai"
    problems=$(( problems + 1 ))
}

# check if the reference dictionary exists
[[ -s "$dict" ]] || {
    log "reference dictionary missing: $dict"
    problems=$(( problems + 1 ))
}

# check the bwa index exists
for ext in amb ann bwt pac sa; do
    [[ -s "${REF}.${ext}" ]] || {
        log "BWA index missing: ${REF}.${ext}"
        problems=$(( problems + 1 ))
    }
done

    (( problems == 0 )) || die "validation failed with ${problems} problem(s)"
    log "validation passed"
}

