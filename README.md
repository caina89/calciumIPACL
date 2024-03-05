# calciumIPACL
Code to reproduce analyses in the Calcium IPACL paper: link_to_paper

## Quality Control (QC) using 1_qc.ipynb (in Python)

QC was done with the following order:
1. Calculate QC metrics.
2. Filter low quality cells based on MADs parameters.
3. Write files for runing SoupX in R.
4. Read output files from SoupX.
5. Filter out low quality genes based on the SoupX-corrected count matrix.
6. Detecting doublets.
7. Write QC-ed count matrix for downstream analyses.


## SoupX using 1_qc_soupx.ipynb (in R).

Read output files generated from the first half of 1_qc.ipynb, run SoupX, write corrected count matrix for further analyses.

## Combine samples using 2_combine_*.ipynb (in Python)
Combine individual samples into the following groups and perform batch correction using Scanorama for downstream analyses:
1. Each treatment, including:
   (1) control, (2) condense milk, (3) elevated platform, (4) quinine, (5) social.
2. Each treatment with control, including:
   (1) condense milk + control,
   (2) elevated platform + control,
   (3) quinine + control,
   (4) social + control.
4. Two treatments, including:
   (1) elevated platform + condense milk,
   (2) quinine + condense milk,
   (3) quinine + elevated platform,
   (4) social + condense milk,
   (5) social + elevated platform,
   (6) social + quinine.
5. All treatments without control.
6. All treatments with control.
