/*

indices for catalogdb tables, to be run after bulk uploads

psql -f catalogdb.sql -h db.sdss.utah.edu -U sdssdb_admin -p 5432 sdss5db

drop index catalogdb.gaia_dr1_tgas_dec_index;

*/


-- Indices

CREATE INDEX CONCURRENTLY gaia_dr1_source_ra_index ON catalogdb.gaia_dr1_source using BTREE (ra);
CREATE INDEX CONCURRENTLY gaia_dr1_source_dec_index ON catalogdb.gaia_dr1_source using BTREE (dec);
CREATE INDEX CONCURRENTLY gaia_dr1_source_l_index ON catalogdb.gaia_dr1_source using BTREE (l);
CREATE INDEX CONCURRENTLY gaia_dr1_source_b_index ON catalogdb.gaia_dr1_source using BTREE (b);
CREATE INDEX CONCURRENTLY gaia_dr1_source_ecl_lon_index ON catalogdb.gaia_dr1_source using BTREE (ecl_lon);
CREATE INDEX CONCURRENTLY gaia_dr1_source_ecl_lat_index ON catalogdb.gaia_dr1_source using BTREE (ecl_lat);
CREATE INDEX CONCURRENTLY gaia_dr1_source_phot_g_mean_flux_index ON catalogdb.gaia_dr1_source using BTREE (phot_g_mean_flux);
CREATE INDEX CONCURRENTLY gaia_dr1_source_phot_g_mean_mag_index ON catalogdb.gaia_dr1_source using BTREE (phot_g_mean_mag);
CREATE INDEX CONCURRENTLY gaia_dr1_source_solution_id_index ON catalogdb.gaia_dr1_source using BTREE (solution_id);
CREATE INDEX CONCURRENTLY gaia_dr1_source_source_id_index ON catalogdb.gaia_dr1_source using BTREE (source_id);

CREATE INDEX CONCURRENTLY gaia_dr1_tgas_ra_index ON catalogdb.gaia_dr1_tgas using BTREE (ra);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_dec_index ON catalogdb.gaia_dr1_tgas using BTREE (dec);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_l_index ON catalogdb.gaia_dr1_tgas using BTREE (l);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_b_index ON catalogdb.gaia_dr1_tgas using BTREE (b);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_ecl_lon_index ON catalogdb.gaia_dr1_tgas using BTREE (ecl_lon);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_ecl_lat_index ON catalogdb.gaia_dr1_tgas using BTREE (ecl_lat);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_phot_g_mean_flux_index ON catalogdb.gaia_dr1_tgas using BTREE (phot_g_mean_flux);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_phot_g_mean_mag_index ON catalogdb.gaia_dr1_tgas using BTREE (phot_g_mean_mag);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_solution_id_index ON catalogdb.gaia_dr1_tgas using BTREE (solution_id);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_source_id_index ON catalogdb.gaia_dr1_tgas using BTREE (source_id);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_hip_index ON catalogdb.gaia_dr1_tgas using BTREE (hip);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_tycho2_id_index ON catalogdb.gaia_dr1_tgas using BTREE (tycho2_id);
CREATE INDEX CONCURRENTLY gaia_dr1_tgas_duplicated_source_index ON catalogdb.gaia_dr1_tgas using BTREE (duplicated_source);

