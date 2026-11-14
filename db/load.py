import csv
import json
import os
import sqlite3
import sys


METADATA_COLUMNS = [
    "sample_id", "phenotype", "sex", "population",
    "family", "father", "mother", "run",
    "read_length", "bases_gb", "study", "native_layout",
    "library_type", "is_synthetic_phenotype"
]


def load_metadata(con, metadata_file):
    """Load cohort sample metadata, skipping comment lines."""

    with open(metadata_file, newline="") as f:
        reader = csv.DictReader(
            (line for line in f if line.strip() and not line.startswith("#")),
            delimiter="\t"
        )

        for row in reader:
            columns = ", ".join(METADATA_COLUMNS)
            placeholders = ", ".join("?" for _ in METADATA_COLUMNS)

            con.execute(
                f"INSERT INTO samples ({columns}) VALUES ({placeholders}) "
                "ON CONFLICT(sample_id) DO NOTHING",
                [row[col] for col in METADATA_COLUMNS]
            )


def load_run(con, folder):
    """Load one pipeline run as an atomic transaction."""

    folder = os.path.abspath(folder)

    with open(os.path.join(folder, "manifest.json")) as f:
        manifest = json.load(f)

    pipeline = manifest["pipeline"]
    run_id = pipeline["run_id"]


    if con.execute(
        "SELECT 1 FROM runs WHERE run_id = ?", (run_id,)
    ).fetchone():
        print(f"SKIPPED: {folder} ({run_id} already loaded)")
        return

    with con:

        con.execute(
            """
            INSERT INTO pipelines (name, version, implementation)
            VALUES (?, ?, ?)
            ON CONFLICT(name, version, implementation) DO NOTHING
            """,
            (
                pipeline["name"],
                pipeline["version"],
                pipeline["implementation"]
            )
        )

        pipeline_id = con.execute(
            """
            SELECT pipeline_id
            FROM pipelines
            WHERE name = ? AND version = ? AND implementation = ?
            """,
            (
                pipeline["name"],
                pipeline["version"],
                pipeline["implementation"]
            )
        ).fetchone()[0]


        params = manifest["params"]
        platform = manifest["platform"]

        con.execute(
            """
            INSERT INTO runs (
                run_id, pipeline_id, started_at, finished_at,
                exit_status, platform, executor, samplesheet,
                reference, region, folder
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                run_id,
                pipeline_id,
                pipeline["started_at"],
                pipeline["finished_at"],
                pipeline["exit_status"],
                platform["kind"],
                platform["executor"],
                params["samplesheet"],
                params["ref"],
                params["region"],
                folder
            )
        )


        for sample in manifest["samples"]:
            sample_id = sample["sample_id"]


            con.execute(
                """
                INSERT INTO samples (sample_id)
                VALUES (?)
                ON CONFLICT(sample_id) DO NOTHING
                """,
                (sample_id,)
            )


            con.execute(
                """
                INSERT INTO run_samples (
                    run_id, sample_id, condition, library_type
                )
                VALUES (?, ?, ?, ?)
                """,
                (
                    run_id,
                    sample_id,
                    sample["condition"],
                    sample["library_type"]
                )
            )


        with open(os.path.join(folder, "variants.tsv"), newline="") as f:
            reader = csv.reader(f, delimiter="\t")

            for row in reader:
                if not row or all(not value.strip() for value in row):
                    continue

                # Skip optional header row
                if row[0].strip().lower() == "chrom" and row[1].strip().lower() == "pos":
                    continue

                chrom, pos, ref, alt, qual, filter_value = row

                con.execute(
                    """
                    INSERT INTO variants (
                        run_id, chrom, pos, ref, alt, qual, filter
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        run_id,
                        chrom,
                        int(pos),
                        ref,
                        alt,
                        float(qual) if qual not in ("", ".") else None,
                        filter_value
                    )
                )



        for output in manifest["outputs"]:
            checksum = output["checksum"]

            if not checksum.startswith("sha256:"):
                raise ValueError(f"Unexpected checksum: {checksum}")

            con.execute(
                """
                INSERT INTO outputs (run_id, path, sha256)
                VALUES (?, ?, ?)
                """,
                (
                    run_id,
                    output["path"],
                    checksum.removeprefix("sha256:")
                )
            )

    print(f"LOADED: {folder} ({run_id})")


def main():
    if len(sys.argv) < 4:
        sys.exit(
            "Usage: python3 db/load.py DATABASE METADATA RUN_FOLDER..."
        )

    database = sys.argv[1]
    metadata = sys.argv[2]
    folders = sys.argv[3:]

    with sqlite3.connect(database) as con:
        con.execute("PRAGMA foreign_keys = ON")
        load_metadata(con, metadata)

        for folder in folders:
            load_run(con, folder)


if __name__ == "__main__":
    main()
