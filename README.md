# Wastewater-Based Genomic Surveillance of Respiratory Syncytial Virus

## Computational Analysis of Wastewater-Derived NGS Data (2022/2023 and 2023/2024 Seasons) of RSV (Subtypes A and B)

---

## Repository Structure  
This repository contains configuration files and reference materials for RSV analysis using **V-pipe**.

### Workflow Configuration Files  
- **Path:** `workflow/configuration/`  
- Configuration files for RSV analysis using **V-pipe**.  
- `samples.tsv` files are available.

### Reference Genomes  
Reference genomes used for RSV analysis:
- **RSV-A:** *hRSV/A/England/397/2017* (GISAID accession no. EPI_ISL_412866)  
- **RSV-B:** *HRSV/B/AUS/VIC-RCH056/2019* (GISAID accession no. EPI_ISL_1653999)  
- **File Location:** `workflow/references/concat_EPI_ISL_412866_EPI_ISL_1653999.fasta`

### Primers & Inserts Genomic Coordinates (BED Files)  
- **RSV-A:** `workflow/references/primers/RSV-A/`  
- **RSV-B:** `workflow/references/primers/RSV-B/`  

### RSV Data Analysis Scripts  
- **General Analysis (Subtype Independent):** `workflow/shared/`  
- **RSV-A Specific Analysis:** `workflow/shared/RSV_A_analysis/`  
- **RSV-B Specific Analysis:** `workflow/shared/RSV_B_analysis/`  

### Figures Published in the Preprint  
- **File Location:** `preprint/plots/`

---

## About This Study  
This study explores the potential of wastewater-based genomic surveillance for **tracking RSV lineages** and **identifying key mutations**. Given the importance of the **RSV F gene** in vaccine and therapeutic strategies, our findings highlight the value of wastewater as a complementary approach for monitoring RSV evolution.

---

## Citation  
If you use this repository in your research, please cite:  
**Wastewater-based sequencing of Respiratory Syncytial Virus enables tracking of lineages and identifying mutations at antigenic sites.**  
Jolinda de Korne-Elenbaas, Auguste Rimaite, Ivan Topolsky, David Dreifuss, Charlyne Bürki, Lara Fuhrmann, Louis du Plessis, William J. Fitzsimmons, Emily E. Bendall, Tanja Stadler, Niko Beerenwinkel, Timothy R Julian; MedRxiv, 2025
