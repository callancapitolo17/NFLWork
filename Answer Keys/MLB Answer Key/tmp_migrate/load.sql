COPY pbp_all FROM 'tmp_migrate/pbp_all.csv' (FORMAT 'csv', quote '"', delimiter ',', header 0);
COPY mlb_betting_history FROM 'tmp_migrate/mlb_betting_history.csv' (FORMAT 'csv', quote '"', delimiter ',', header 0);
