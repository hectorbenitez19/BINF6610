SELECT
    r.run_id,
    COUNT(v.run_id) AS n_variants,
    SUM(CASE WHEN v.filter = 'PASS' THEN 1 ELSE 0 END) AS n_pass
FROM runs AS r
LEFT JOIN variants AS v
    ON r.run_id = v.run_id
GROUP BY r.run_id;