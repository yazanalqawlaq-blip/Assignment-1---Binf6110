# Assignment-1---Binf6110
Assignment 1 for Binf6110, part 1

# Assignment 1 — Part 1 (Genome assembly plan)

## Introduction
Salmonella enterica is a Gram negative, rod shaped, facultative anaerobic bacterium that belongs to the family Enterobacteriaceae (Jajere, 2019). Salmonella enterica is a major cause of foodborne illness and remains an important public health concern worldwide (Havelaar et al., 2015). According to the WHO global burden estimates, most of the overall foodborne disease burden comes from diarrheal illnesses, and Salmonella enterica accounts for a meaningful share of that burden across the world (Havelaar et al., 2015). In this assignment, I will use Oxford Nanopore R10 long read sequencing data to assemble a consensus genome for Salmonella enterica and then compare it to a reference genome from NCBI through alignment, variant calling, and visualization. The significance of this project is that assembling and comparing the genome provides a clear, practical way to understand how genome assembly quality and reference based analysis affect what we conclude about variation in a clinically important pathogen.
A key challenge in bacterial genome assembly is resolving repeats and mobile genetic elements, which can fragment assemblies or lead to misassemblies if reads do not span repetitive regions (Wick and Holt, 2019). Long reads are useful because they can bridge repeats and support highly contiguous bacterial assemblies, including recovery of plasmids when present (Sereika et al., 2022). The main tradeoff is that nanopore reads can still contain base level errors, and small indels can persist and show up as false differences during reference comparison if polishing and validation are not handled carefully (Luan et al., 2024). Reference based analysis also requires caution because alignment artifacts and low support calls can look like real variants, so quality checks and visualization are needed before interpreting results (Thorvaldsdóttir et al., 2013).
## Proposed Methods

## References
