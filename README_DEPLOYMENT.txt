Flu-GDB Mutation Explorer - Deployment Package
===============================================

QUICK START:
1. Extract this archive on your Shiny server
2. Run: Rscript install_packages.R
3. Ensure NCBI BLAST+ is installed (blastp, blastx and makeblastdb available)
4. Check flu_gdb_app_data/ is in place, including master_blastdb/ (see below)
5. Configure Shiny Server to serve this directory
6. Restart Shiny Server

REQUIRED DATA:
The BLAST database of reference proteins, flu_gdb_app_data/master_blastdb/
reference_proteins_db, is not optional. Every tab that takes a submitted
sequence identifies it against that one database - Batch Screening reads the
product off an HSP, and the Tree and Adaptation Mutations tabs read the segment,
the product and, for HA, the H3 or H5 numbering scheme off the same hit. The app
refuses to start without it and says so. It is supplied ready-made with the
upstream data, alongside reference_proteins_db.fasta, the FASTA it is built from.

For detailed instructions, see the deployment checklist in the
artifacts directory of the development environment.

SYSTEM REQUIREMENTS:
- R >= 4.0
- NCBI BLAST+ tools
- Shiny Server
- Sufficient RAM (recommend 8GB+)

DATA FOLDERS:
The app reads one generated folder, flu_gdb_app_data/, the single output folder of
the upstream pipeline. Besides the files below it holds trees/ (the midpoint-rooted
cluster trees), blastdb/ (sgt_1 ... sgt_8), master_blastdb/ (the reference protein
database) and IAV_DB_summary.log. Keep it in step with the GLUE server.

flu_gdb_app_data/ holds, among others:
  product_positions.rds, product_columns.rds   all 12 products, keyed seg1 ... seg8_NEP
  reference_numbering.rds   every reference residue with its alignment column; the
                            Tree tab uses it to carry a Position search between
                            references. Optional - without it the Position box keeps
                            its number when the reference changes
  alignment_column_counts.rds   optional; the "All sequences" plot
  cluster_column_residues.rds   optional; the tree pop-up's residue breakdown

The older sequence_positions.rds, alignment_columns.rds and secondary_orf_positions.rds
are no longer read and can be deleted.

VALIDATION:
After deployment, check the logs for:
"[SUCCESS] Verified reference 'XXX' for 'segX'"
for all 8 segments.
