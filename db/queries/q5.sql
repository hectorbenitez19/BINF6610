SELECT
    r.run_id,
    r.folder || '/' || o.path AS file,
    o.sha256
FROM runs AS r
JOIN outputs AS o
    ON r.run_id = o.run_id
WHERE o.path = 'cohort.filtered.vcf.gz';
