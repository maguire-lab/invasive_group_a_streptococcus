# Invasive Group A Streptococcus Shared Hospital Laboratory

Our goal is to perform genomic analyses of all Streptococcus pyogenes isolates sequenced at the Shared Hospital Laboratory as of 2024-01-22. 

## Bactopia 
To do this we performed a bactopia run on all SHL samples (listed in `20240122_gas_sample_sheet.csv`).

`bactopia --samples 20240122_gas_sample_sheet.csv --outdir 20240122_gas --run_name 20240122_gas -profile docker

We also performed emmtyping and traditional MLST within bactopia.

Then we generated a pangenome of these samples:

`bactopia --wf pangenome --bactopia 20240122_gas --outdir 20240122_gas_pangenome --run_name 20240122_gas_pangenome -profile docker --use_panaroo`

These can be recreated with `run.sh` at the toplevel of the repository.

## Manual Analyses
The following manual analyses were then performed (with a `run.sh` within each folder):

- `assembly_typer/` - m1uk SNP calling

- `poppunk/` - poppunk clustering

- `pyseer/` - genome wide association studies of invasiveness and unitigs

- `virulence_adhesion_factors/` - BLAST based identification of a curated set of viruluence and adhesion factors

## Methods

SHL isolates (n=166) were assembled using Bactopia v3.0.0[^1] workflow with default settings. 
In brief, this involves QC with fastp v0.23.4[^2] followed by assembly using shovill v1.1.0[^3] and annotation using prokka v1.14.6[^4] and amrfinderplus v3.11.18 (with database v2023-08-08.2)[^5].  
Via Bactopia subworkflows/tools the following analyses were also performed: multi-locus sequence typing (MLST) using mlst v2.23.0 [^6] against the 2.23.0-20230907 database and Enrigh et. al. scheme[^7], emm types inferred using emmtyper v0.2.0[^8], and a core genome phylogeny inferred using panaroo v1.3.4[^9], clonalframeml v1.12[^10], and iqtree v2.2.2.7[^11]. 
All SHL isolates were then clustered against a curated reference set of 2,085 assembled Streptococcus pyogenes [^19] genomes with poppunk v2.6.2 and the spyogenes v1.0.0 database[^12].
A minimum spanning tree was then inferred from cluster assignments via poppunk.
Both the core genome phylogeny of SHL samples and the minimum spanning tree of SHL and representative contextual genomes was visualised using microreact[^13].

A customised set of virulence and adhesion factors was (supplemental table X) from Virulence Factor Database (VFDB)[^14] and RefSeq[^15]. These were annotated using BLASTN+ v2.15.0[^16] with a minimum identity and query coverage of 80%.
Alignments were also generated using minimap2[^17] and gofasta[^18] against the MGAS8232 covR and covS reference sequences.

[^1] Petit III, Robert A., and Timothy D. Read. "Bactopia: a flexible pipeline for complete analysis of bacterial genomes." Msystems 5.4 (2020): 10-1128.
[^2] Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor." Bioinformatics 34.17 (2018): i884-i890.
[^3] https://github.com/tseemann/shovill
[^4] Seemann, Torsten. "Prokka: rapid prokaryotic genome annotation." Bioinformatics 30.14 (2014): 2068-2069.
[^5] Feldgarden, Michael, et al. "AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence." Scientific reports 11.1 (2021): 12728.
[^6] https://github.com/tseemann/mlst
[^7] Enright, Mark C., et al. "Multilocus sequence typing of Streptococcus pyogenes and the relationships between emm type and clone." Infection and immunity 69.4 (2001): 2416-2427.
[^8] https://github.com/MDU-PHL/emmtyper
[^9] Tonkin-Hill, Gerry, et al. "Producing polished prokaryotic pangenomes with the Panaroo pipeline." Genome biology 21 (2020): 1-21.
[^10] Didelot, Xavier, and Daniel J. Wilson. "ClonalFrameML: efficient inference of recombination in whole bacterial genomes." PLoS computational biology 11.2 (2015): e1004041.
[^11] Minh, Bui Quang, et al. "IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era." Molecular biology and evolution 37.5 (2020): 1530-1534.
[^12] Lees, John A., et al. "Fast and flexible bacterial genomic epidemiology with PopPUNK." Genome research 29.2 (2019): 304-316.
[^13] Argimón, Silvia, et al. "Microreact: visualizing and sharing data for genomic epidemiology and phylogeography." Microbial genomics 2.11 (2016): e000093.
[^14] Liu, Bo, et al. "VFDB 2022: a general classification scheme for bacterial virulence factors." Nucleic acids research 50.D1 (2022): D912-D917.
[^15] O'Leary, Nuala A., et al. "Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation." Nucleic acids research 44.D1 (2016): D733-D745.
[^16] Camacho, Christiam, et al. "BLAST+: architecture and applications." BMC bioinformatics 10 (2009): 1-9.
[^17] Li, Heng. "Minimap2: pairwise alignment for nucleotide sequences." Bioinformatics 34.18 (2018): 3094-3100.
[^18] Jackson, Ben. "gofasta: command-line utilities for genomic epidemiology research." Bioinformatics 38.16 (2022): 4033-4035.
[^19] Davies, Mark R., et al. "Atlas of group A streptococcal vaccine candidates compiled using large-scale comparative genomics." Nature genetics 51.6 (2019): 1035-1043.
