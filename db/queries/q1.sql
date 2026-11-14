SELECT
    r.run_id,
    p.implementation,
    r.platform,
    COUNT(rs.sample_id) AS n_samples,
    ROUND(
        (julianday(r.finished_at) - julianday(r.started_at)) * 24 * 60,
        1
    ) AS minutes
FROM runs AS r
JOIN pipelines AS p
    ON r.pipeline_id = p.pipeline_id
LEFT JOIN run_samples AS rs
    ON r.run_id = rs.run_id
GROUP BY r.run_id
ORDER BY julianday(r.started_at) ASC;
