[Supplementary_README.txt](https://github.com/user-attachments/files/31295017/Supplementary_README.txt)
SUPPLEMENTARY MATERIALS README
==============================

Integration of Bioinformatics and Machine Learning to Characterize 
Fusobacterium nucleatum's Pathogenicity

Authors: Zihan Tian, Pietro Lio'

---
CONTENTS
--------------------
1. Supplementary Tables (S1-S15)
2. Supplementary Figures (S1-S3)
3. Supplementary Data (S1-S7)
4. Software Versions
5. File Naming Notes
6. Data Availability

1. SUPPLEMENTARY TABLES
Separated tables can be found on GitHub with the Supplementary Data (S1-S7)
--------------------

Table S1 (a-c): BLAST analysis results for genes in predicted pathogenicity islands PAI1, PAI2, and PAI3. Each table lists query proteins, subject hits, percent identity, alignment length, E-value, and bit score from BLASTP searches against the NCBI non-redundant (nr) and UniProt Swiss-Prot databases.

Table S2: Promoter predictions for F. nucleatum ATCC 25586 genome using 
PePPER. Contains predicted promoters with confidence scores 
(0.0-1.0), -10 box sequences, and promoter sequences.
	-Supplementary_Table_S2_promoters_coordinates.xlsx contains the same information in coordinate format required for running the promoter analysis in Data S4.

Table S3: STRING protein-protein interaction network data for all 3 PAI proteins. Includes interaction edges, combined scores, and evidence channel breakdowns from STRING v12.0.

Table S4: AlphaFold3 Server structure predictions (monomer and multimer) and DeepTMHMM transmembrane topology predictions for all 3 PAI proteins. Includes pTM and ipTM confidence scores for monomer and complex models, respectively, and topology classifications (GLOB, SP, TM, SP+TM, BETA).

Table S5: Co-expression and differential expression analysis (GSE161360). Includes DESeq2 differential expression results across growth phases and genome-wide Spearman co-expression statistics.

Table S6: Comprehensive virulence factor analysis using VFDB and VICTORS databases for all genomic regions of interest. Hits from both databases are shown with BLASTP statistics.

Table S7: Curated iron-related reaction (65) and metabolite (37) sets used for FBA. Includes category assignments (primary iron exchange, alternative iron exchange, siderophore ligand exchange, iron transport, siderophore/heme uptake, iron reduction/dissociation, other) and explanation of 5 reactions and 3 metabolites excluded as keyword-matching false positives.

Table S8: Gene classified as high-CAI (top 20%, CAI ≥ 0.771) 

Table S9: Gene classified as low-CAI (bottom 20%, CAI ≤ 0.670) 

Table S10: Relative Synonymous Codon Usage (RSCU) analysis result, with amino acid, sequence, and RSCU score

Table S11: eggNOG-mapper v2 (Huerta-Cepas et al., 2018) result. Key virulence-relevant findings are listed. Sheet 1 contains raw annotations (seed orthologs, COG categories, KEGG pathway mappings, GO terms). Sheet 2 contains curated virulence-relevant findings with FNxxxx locus tag assignments.

Table S12: Comparative Genomics Analysis, organized from the results of Supplementary_data_S6_PAI_comparative_genomics_pipeline. 
Three components:
  1) gene_content_conservation: All 75 PAI genes with coordinates, product, protein ID, aa length, conservation (n/15 and n/17), and mobility element classification (IS family, status, length, location).
  2) compositional_phylogenetic: GC%, z-scores (GC + dinucleotide), atypicality flags, mobility genes inside/flanking, gene trees tested with nRF values, mean nRF for each PAI.
  3) Comparative Genomics Lists: All 17 genomes used for comparative genomics analysis.

Table S13: a) Comparison between AlphaFold3 monomer models and experimentally determined reference structures from the Protein Data Bank. Includes the query protein, reference structure, TM-score (query-normalized and reference-normalized), RMSD, aligned length, sequence identity across the aligned region, and interpretation of fold agreement.
b) Conserved catalytic or functional motifs mapped between AlphaFold3 models and their corresponding reference structures. For each protein, reports the reference motif, the aligned query residue(s), whether the key residue is conserved, whether the motif aligns at the same structural position, and a brief note on functional interpretation.

Table S14: Lysozyme inhibitors. C-type lysozyme (MliC-type) and LprI-family lysozyme inhibitor. Pairwise sequence comparison and structural alignment with experimentally determined reference structures from the PDB.

Table S15:Prokka annotations validated with RNAseq validation.

---

2. SUPPLEMENTARY FIGURES
---------------------

Figure S1 (a-e): Promoter motif analysis of high and low CAI genes and PAI loci using WebLogo (Crooks et al., 2004; Schneider & Stephens, 1990). 
	(a) Sequence logo of enriched 8-mer motifs in promoters of high CAI genes (CAI > 0.771). 
	(b) Sequence logo of promoter motifs from low CAI genes (CAI < 0.670). 
	(c) Promoter motif analysis of genes in PAI1. 
	(d) Promoter motif analysis of genes in PAI2. 
	(e) Promoter motif analysis of genes in PAI3. 
The motifs were generated from 8-mer sequences extracted from predicted promoter regions, filtered by minimum occurrence, sorted by frequency, and visualized with WebLogo.

Figure S2: Gene-level conservation heatmap for the three candidate PAIs across 15 Fusobacterium genomes and 2 outgroups. Each row represents a PAI-encoded gene; each column represents a genome. Color indicates presence/absence based on TBLASTN homology

Figure S3 (a-c): Synteny analysis of the three candidate PAIs across representative Fusobacterium genomes. (a) PAI1. (b) PAI2. (c)PAI3. Homologous genes are connected by colored blocks; inversions are shown by crossing lines. PAI3 exhibits the most stable synteny; PAI1 and PAI2 show local rearrangements associated with flanking transposase insertions.

---

3. SUPPLEMENTARY DATA
On Github
------------------

Data S1: Genome-scale metabolic model of F. nucleatum ATCC 25586 in SBML format, generated using CarveMe with the BiGG universal model template and Prokka protein annotations as input.

Data S2: F. nucleatum ATCC 25586 genome sequence in FASTA format (NCBI accession: AE009951).

Data S3: Genome annotation in GenBank format (NCBI accession: AE009951).

Data S4: Python analysis pipeline (Jupyter notebook) containing code for:
  - GC skew and cumulative GC skew analysis
  - Codon Adaptation Index (CAI) analysis
  - Promoter analysis (using PePPER v2.0 output)
  - Literature keyword search for protein function annotation
  - BLASTP against the UniProt Swiss-Prot database and NCBI non-redundant
    (nr) protein database (via Biopython v1.79 NCBIWWW.qblast)
  - Targeted BLASTP searches against two curated virulence factor databases:
    the Virulence Factor Database (VFDB) and the Virulence Factors of
    Pathogenic Bacteria database (VICTORS)
  - STRING protein interaction network analysis for 3 PAI protein sets (STRING v12.0 API)
  - Metabolic modeling (FBA under iron-availability conditions using COBRApy)
  - Essential gene and iron homeostasis analysis (NetGenes database +
    GenBank cross-referencing)

Data S5: R script for DESeq2 differential expression analysis across growth phases, genome-wide Spearman co-expression and permutation testing. Python script WIG-to-expression-matrix conversion utility for processing GEO GSE161360 RNA-seq coverage tracks; Prokka annotation validation against RNA-seq data

Data S6: Python analysis pipeline for Comparative Genomic Analysis.
- M1: whole-region BLASTN conservation matrix (query coverage %, best nucleotide identity %) across 15 Fusobacterium genomes and 2 outgroups
- M2: gene-level TBLASTN presence/absence matrix
- M3: compositional analysis (GC content, dinucleotide frequency distances vs. 500 size-matched random genomic windows)
- M4: mobility element annotation (IS elements, transposases)
- M5: targeted HGT tests (composition, mobility context, gene-tree concordance via MAFFT + FastTree, normalized Robinson-Foulds distances)
Requires: Python 3.11, NCBI BLAST+, MAFFT, FastTree, Biopython, pandas, matplotlib, openpyxl.

Data S7: Structure analysis. Convert AlphaFold-downloaded structures from cif to PDB, and TM-align and MAFFT analyses 

---

4. SOFTWARE VERSIONS
-----------------

Python 3.9 (Data S4 notebook)
  - COBRApy 0.26.0
  - Biopython 1.79
  - pandas 1.4.0
  - numpy 1.22.0
  - matplotlib
  - openpyxl

Python 3.11 (Data S6 comparative genomics pipeline)
  - Biopython (SeqIO, AlignIO, Phylo)
  - pandas, matplotlib, openpyxl

R 4.1.0
  - DESeq2 (Bioconductor)
  - bnlearn 4.7
  - tidyverse
  - pheatmap
  - RColorBrewer
  - ggrepel

External tools and web servers:
  - CarveMe 1.5.1 (genome-scale metabolic model reconstruction)
  - IslandViewer 4 (genomic island prediction; web server)
  - Prokka v1.14 (genome annotation)
  - eggNOG-mapper v2.0 (functional annotation; eggNOG 5.0 database)
  - PePPER v2.0 (promoter prediction; web server)
  - STRING v12.0 (protein-protein interaction networks; web API)
  - DeepTMHMM (transmembrane topology prediction; web server)
  - AlphaFold3 Server (structural modeling; alphafoldserver.com)
  - NCBI BLAST+ (BLASTN, TBLASTN for comparative genomics; BLASTP via
    Biopython NCBIWWW.qblast for nr/Swiss-Prot/VFDB/VICTORS searches)
  - MAFFT (multiple sequence alignment for gene trees and structural analysis)
  - FastTree (phylogenetic tree inference for gene-tree concordance tests)
  - TM-align (structural alignment of AlphaFold3 models vs. PDB references)
  - WebLogo (sequence logo visualization for promoter motif analysis)
  - NetGenes database (predicted essentiality; downloaded gene set)

---

5. FILE NAMING NOTES
----------------

Input files referenced in the analysis code use their original names from 
public databases. The supplementary files have been renamed for clarity:

Original Name -> Supplementary Name
EX:
25586_NCBI.gbk -> Supplementary_data_S3_ geneBank.gbk

---

6. DATA AVAILABILITY
-----------------

F. nucleatum ATCC 25586 genome: NCBI accession AE009951.2
Transcriptomic data: GEO accession GSE161360
AlphaFold Server: https://alphafoldserver.com/
STRING database: https://string-db.org/
IslandViewer 4: https://www.pathogenomics.sfu.ca/islandviewer/
NetGenes database: https://tubic.org/netgenes/

All code, supplementary tables, and data files supporting this study are
publicly available at GitHub:
https://github.com/exporterdream/fn-pai-characterization

For questions, contact: [zihantian123@outlook.com]
