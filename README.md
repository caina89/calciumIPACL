# calciumIPACL
Code to reproduce analyses in the Calcium IPACL paper: link_to_paper

## Quality Control (QC) using 1_qc.ipynb

Perfrom QC on each individual sample with the following steps:
1. Calculate QC metrics.
2. Filter out low quality cells based on MADs parameters.
3. Write files for runing SoupX in R.
4. Read output files from SoupX.
5. Filter out low quality genes based on the SoupX-corrected count matrix.
6. Detect doublets.
7. Write QC-ed count matrix for further analyses.


## SoupX using 1_qc_soupx.ipynb (in R).

Read output files generated from the first half of 1_qc.ipynb, run SoupX in R, write corrected count matrix for further analyses.


## Combine samples using 2_combine_*.ipynb
Combine individual samples into the following groups and perform batch correction using Scanorama for downstream analyses:
1. Each treatment, including:<br>
   (1) control, (2) condense milk, (3) elevated platform, (4) quinine, (5) social.
2. Each treatment with control, including:<br>
   (1) condense milk + control,<br>
   (2) elevated platform + control,<br>
   (3) quinine + control,<br>
   (4) social + control.
4. Two treatments, including:<br>
   (1) elevated platform + condense milk,<br>
   (2) quinine + condense milk,<br>
   (3) quinine + elevated platform,<br>
   (4) social + condense milk,<br>
   (5) social + elevated platform,<br>
   (6) social + quinine.
5. All treatments without control.
6. All treatments with control.


## Generate plots and tables from MapMyCells annotations using 3_mapmycells_*.ipynb
Generte UMAPs, barplots and tables using annotations from MapMyCells. For each individual sample and combined samples. 


## Differentially expressed genes (DEGs) analyses using 4_DEG_*.ipynb
For each data listed below, generate pseudobulks using decoupleR, run DEG analysis using DESeq2.
1. 1 treatment vs. controls (using the combined data in category 2)
2. 1 treatment vs. 1 treatment (using the combined data in category 4)
3. 1 treatment vs. 3 treatments (using the combined data in category 5)






