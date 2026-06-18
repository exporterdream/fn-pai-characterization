# Computational Characterization of Pathogenicity Islands in Fusobacterium nucleatum ATCC 25586
Code and supplementary materials for: Integration of Bioinformatics and Machine Learning to Characterize Fusobacterium nucleatum's Pathogenicity. Tian &amp; Lio'

The study presents the first systematic computational characterization of pathogenicity islands (PAIs) in F. nucleatum ATCC 25586, integrating genomic island prediction, functional annotation, co-expression network analysis, structural modeling, and genome-scale metabolic modeling to identify three candidate PAIs and propose a mechanistic hypothesis linking hemolysin-mediated iron acquisition to colorectal cancer progression.

Code zenodo DOI: 10.5281/zenodo.20748204
    
SUPPLEMENTARY TABLES
--------------------

Table S1a-c: BLAST analysis results for genes in predicted pathogenicity islands PAI1, PAI2, and PAI3.

Table S2: Promoter predictions for F. nucleatum ATCC 25586 genome using 
PePPER. Contains predicted promoters with confidence scores 
(0.0-1.0), -10 box sequences, and promoter sequences.
	-Supplementary_Table_S2_promoters_coordinates.xlsx is same information, but less informative than Table S2, it is important for running Promoter analysis for Data S4. 

Table S3: STRING protein-protein interaction network data for all 3 PAI proteins.

Table S4: AlphaFold-Multimer structure (monomer and multimer) predictions and TMHMM transmembrane helix predictions for all 3 PAI proteins. Includes pTM and ipTM confidence scores.

Table S5: Co-expression network edges predicted by Bayesian network learning (bnlearn) from F. nucleatum transcriptomic data (n=18 samples) that is significant.

Table S6: Comprehensive virulence factor analysis using databases (VFDB + VICTORS) for all genomic regions of interest. Hits from both database shown.

Table S7: Flux balance analysis (FBA) results showing predicted growth rates under different iron availability conditions.

Table S8: Gene classified as high-CAI (top 20%, CAI ≥ 0.771) 

Table S9: Gene classified as low-CAI (bottom 20%, CAI ≤ 0.670) 

Table S10: Relative Synonymous Codon Usage (RSCU) analysis result, with amino acid, sequence, and RSCU score

Table S11: eggNOG-mapper v2 (Huerta-Cepas et al., 2018) result. Key virulence-relevant findings are listed.

---

SUPPLEMENTARY FIGURES
---------------------

Figure S1 (a-e): Promoter motif analysis of high and low CAI genes 
and PAI loci using WebLogo (Crooks et al., 2004; Schneider & Stephens, 1990). 
	(a) Sequence logo of enriched 8-mer motifs in promoters of high CAI genes 
(CAI > 0.771). 
	(b) Sequence logo of promoter motifs from low CAI genes 
(CAI < 0.670). 
	(c) Promoter motif analysis of genes in PAI1. 
	(d) Promoter motif analysis of genes in PAI2. 
	(e) Promoter motif analysis of genes in PAI3. 
The motifs were generated from 8-mer sequences extracted from predicted promoter regions, filtered by minimum occurrence, sorted by frequency, and visualized with WebLogo.

---

SUPPLEMENTARY DATA
------------------

Data S1: Genome-scale metabolic model of F. nucleatum ATCC 25586 in SBML 
format, generated using CarveMe.

Data S2: F. nucleatum ATCC 25586 genome sequence in FASTA format 
(NCBI accession: AE009951).

Data S3: Genome annotation in GenBank format.

Data S4: Python analysis pipeline (Jupyter notebook) containing code for:
- GC skew
- CAI analysis
-Promotor analysis
- Literature Keyword Search for 1.7 to 1.85Mb protein function
- BLASTP against the UniProt Swiss- Prot database and NCBI non-redundant (nr) protein database
- Targeted BLASTP searches against two curated virulence factor databases: the Virulence Factor Database (VFDB; 4,623 experimentally verified virulence factors) and the Virulence Factors of Pathogenic Bacteria database (VICTORS; 4,964 virulence factors)
- String analysis for 3 PAI proteins
- Metabolic modeling
- Essential gene and iron homeostasis analysis

Data S5: R script for Bayesian network learning (bnlearn) co-expression 
analysis. With readme for code. convert_wig_to_expression is not part of the pipeline, but a tool to achieve an online database for analysis.

---

SOFTWARE VERSIONS
-----------------

Python 3.9
- COBRApy 0.26.0
- Biopython 1.79
- pandas 1.4.0
- numpy 1.22.0


R 4.1.0
- bnlearn 4.7

External tools:
- CarveMe 1.5.1
- IslandViewer 4
- Alien Hunter 1.7
- PePPER (web server)
- STRING v11.5
- AlphaFold-Multimer v2.3
- NCBI BLAST+: brew install blast (Mac)
- networkx: pip install networkx (optional, for network analysis)

---

FILE NAMING NOTE
----------------

Input files referenced in the analysis code use their original names from 
public databases. The supplementary files have been renamed for clarity:

Original Name -> Supplementary Name
EX:
25586_NCBI.gbk -> Supplementary_data_S3_ geneBank.gbk

---

DATA AVAILABILITY
-----------------

F. nucleatum ATCC 25586 genome: NCBI accession AE009951
Transcriptomic data: GEO accession [GSE161360]

For questions, contact: [zihantian123@outlook.com]
