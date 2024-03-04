# calciumIPACL
Code to reproduce analysee in the Calcium IPACL paper: link_to_paper

## Quality Control (QC)

QC was done using code from 1_QC.ipynb with the following order:
1. Calculate QC metrics.
2. Filter low quality cells based on MADs parameters.
3. Write files for runing SoupX in R.
4. Read output files from SoupX.
5. Filter out low quality genes based on the SoupX-corrected count matrix.
6. Detecting doublets.
7. Write QC-ed count matrix for downstream analyses.


## SoupX (in R)

Run SoupX using 1_qc_soupx.ipynb. Read output files generated from the first half of 1_qc.ipynb, write corrected count matrix for further analyses.

