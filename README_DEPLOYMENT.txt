Flu-GDB Mutation Explorer - Deployment Package
===============================================

QUICK START:
1. Extract this archive on your Shiny server
2. Run: Rscript install_packages.R
3. Ensure NCBI BLAST+ is installed (blastp command available)
4. Configure Shiny Server to serve this directory
5. Restart Shiny Server

For detailed instructions, see the deployment checklist in the
artifacts directory of the development environment.

SYSTEM REQUIREMENTS:
- R >= 4.0
- NCBI BLAST+ tools
- Shiny Server
- Sufficient RAM (recommend 8GB+)

VALIDATION:
After deployment, check the logs for:
"[SUCCESS] Verified reference 'XXX' for 'segX'"
for all 8 segments.
