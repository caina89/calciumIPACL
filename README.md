# Calcium IPACL

Code to reproduce the analyses presented in the Calcium IPACL paper: link_to_paper

## Quality Control (QC) using 1_qc.ipynb

Perfrom QC on each individual sample using the following steps:
1. Calculate QC metrics.
2. Filter out low-quality cells based on MADs parameters.
3. Generate files for running SoupX in R.
4. Read output files from SoupX.
5. Filter out low-quality genes based on the SoupX-corrected count matrix.
6. Detect doublets using Scrublet.
7. Write QC-ed count matrix for further analyses.


## SoupX using 1_qc_soupx.ipynb (in R).

Read the output files generated from the first half of '1_qc.ipynb', run SoupX in R, and write the corrected count matrix for further analyses.


## Combine samples using 2_combine_*.ipynb

Combine individual samples into the following groups and perform batch correction using Scanorama for downstream analyses:
1. Each treatment, including:<br>
   (1) Control, (2) Condensed milk, (3) Elevated platform, (4) Quinine, (5) Social
2. Each treatment with control, including:<br>
   (1) Condensed milk + Control<br>
   (2) Elevated platform + Control<br>
   (3) Quinine + Control<br>
   (4) Social + Control
3. Two treatments, including:<br>
   (1) Elevated platform + Condensed milk<br>
   (2) Quinine + Condensed milk<br>
   (3) Quinine + Elevated platform<br>
   (4) Social + Condensed milk<br>
   (5) Social + Elevated platform<br>
   (6) Social + Quinine
4. All treatments without control.
5. All treatments with control.

As such, we retain as many overlapping genes in each combination as possible for further analyses.


## Explore MapMyCells annotations using 3_mapmycells_*.ipynb

For each individual sample and combined samples, generate UMAPs, bar plots and tables using annotations from MapMyCells. 


## Differentially expressed genes (DEGs) analyses using 4_DEG_*.ipynb

For each data listed below, generate pseudobulks using decoupleR, and run DEG analysis using DESeq2.
1. 1 treatment vs. controls (using the combined data in category 2)
2. 1 treatment vs. 1 treatment (using the combined data in category 3)
3. 1 treatment vs. 3 treatments (using the combined data in category 4)


## Generate plots using codes from 'plots'

Create visually appealing and consistent plots across the entire project.


