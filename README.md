# calciumIPACL
Code to reproduce analyses in the Calcium IPACL paper: link_to_paper

## Quality Control (QC) using 1_qc.ipynb

QC was done with the following order:
1. Calculate QC metrics.
2. Filter low quality cells based on MADs parameters.
3. Write files for runing SoupX in R.
4. Read output files from SoupX.
5. Filter out low quality genes based on the SoupX-corrected count matrix.
6. Detecting doublets.
7. Write QC-ed count matrix for downstream analyses.


## SoupX (in R) using 1_qc_soupx.ipynb

Read output files generated from the first half of 1_qc.ipynb, run SoupX, write corrected count matrix for further analyses.

## Combine samples
Combine individual samples into the following groups with batch correction using Scanorama for downstream analyses:
1. Each treatment, including: control, condense milk, elevated platform, quinine, social
2. 
