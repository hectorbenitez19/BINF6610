SELECT
    s.sample_id,
    s.sex,
    COUNT(rs.run_id) AS n_runs
FROM samples AS s
LEFT JOIN run_samples AS rs
    ON s.sample_id = rs.sample_id
GROUP BY s.sample_id, s.sex
ORDER BY s.sample_id;
