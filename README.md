# Assignment-1---Binf6110

# Assignment 1 — Part 1 (Genome assembly plan)

## Introduction
Salmonella enterica is a Gram negative, rod shaped, facultative anaerobic bacterium that belongs to the family Enterobacteriaceae (Jajere, 2019). Salmonella enterica is a major cause of foodborne illness and remains an important public health concern worldwide (Havelaar et al., 2015). According to the WHO global burden estimates, most of the overall foodborne disease burden comes from diarrheal illnesses, and Salmonella enterica accounts for a meaningful share of that burden across the world (Havelaar et al., 2015). In this assignment, I will use Oxford Nanopore R10 long read sequencing data to assemble a consensus genome for Salmonella enterica and then compare it to a reference genome from NCBI through alignment, variant calling, and visualization. The significance of this project is that assembling and comparing the genome provides a clear, practical way to understand how genome assembly quality and reference based analysis affect what we conclude about variation in a clinically important pathogen.
A key challenge in bacterial genome assembly is resolving repeats and mobile genetic elements, which can fragment assemblies or lead to misassemblies if reads do not span repetitive regions (Wick and Holt, 2019). Long reads are useful because they can bridge repeats and support highly contiguous bacterial assemblies, including recovery of plasmids when present (Sereika et al., 2022). The main tradeoff is that nanopore reads can still contain base level errors, and small indels can persist and show up as false differences during reference comparison if polishing and validation are not handled carefully (Luan et al., 2024). Reference based analysis also requires caution because alignment artifacts and low support calls can look like real variants, so quality checks and visualization are needed before interpreting results (Thorvaldsdóttir et al., 2013).
For this workflow, I am choosing Flye for de novo assembly because it is designed for long read data and is widely used for prokaryotic genomes (Kolmogorov et al., 2019). A major advantage of Flye is that it often produces highly contiguous bacterial assemblies that are suitable for downstream comparison (Wick and Holt, 2019). A limitation of Flye is that it can require more memory than some alternatives, which matters for compute planning (Wick and Holt, 2019). Using Nanopore R10 long reads is advantageous because long read length supports contiguity, but a limitation is that homopolymer associated indels can still occur and need careful handling (Sereika et al., 2022). To improve base accuracy after assembly, I will polish using medaka because polishing tools have been shown to improve microbial nanopore assemblies and reduce errors that affect downstream interpretation (Lee et al., 2021).
## Software versions for  Proposed Methods
Flye 2.9.6-b1802
minimap2 2.30-r1287
samtools 1.23
bcftools 1.23
Medaka 2.0.1
NanoPlot 1.46.2
QUAST 5.3.0
Bandage 0.9.0
IGV 2.19.7
## Proposed Methods
Starting from the Oxford Nanopore R10 FASTQ reads for Salmonella enterica with expected accuracy Q20+ and an N50 of approximately 5 to 15 kb, I will use the dataset deposited on NCBI SRA under SRR32410565 and document the exact download command and files in the repository. I will first run a basic read quality control using NanoPlot to summarize read length and quality distributions and confirm the data match expectations before assembly (De Coster et al., 2018). If the dataset contains an extreme tail of very short or very low quality reads, I will apply a clearly stated filtering threshold and briefly justify it. I will then assemble the genome using Flye with the Nanopore high quality preset and record all non default parameters, including any genome size estimate and thread settings, as well as the full command used (Kolmogorov et al., 2019). After generating the draft assembly, I will map the raw reads back to the assembly using minimap2 with an ONT mapping preset in order to evaluate support across the contigs and to prepare inputs for polishing (Li, 2018). I will then polish the assembly using medaka to improve base level accuracy and reduce residual errors that could otherwise appear as false differences during reference based analysis, and I will report the exact medaka model used because model choice depends on the sequencing and basecalling configuration (Lee et al., 2021). For reference based comparison, I will download an appropriate Salmonella enterica reference genome from NCBI and report the accession, then align the polished assembly to the reference using minimap2 and map reads to the reference for downstream variant calling (Li, 2018). I will convert, sort, and index alignment files using SAMtools to generate coordinate sorted indexed BAM files suitable for downstream analysis and visualization (Li et al., 2009). I will call SNPs and short indels relative to the reference using BCFtools and apply explicit filtering criteria based on quality and depth before interpretation (Danecek et al., 2021). I will evaluate assembly statistics and reference based metrics using QUAST and use IGV to visualize representative regions and confirm read support for key variants, including screenshots as figures to support interpretation (Gurevich et al., 2013). I will also inspect the assembly graph in Bandage to check for circular contigs and unresolved structures that could affect conclusions, and I will keep all commands, parameters, figures, and software versions in the GitHub repository with multiple commits to document the workflow clearly.
## References
Jajere, S. M. (2019). A review of Salmonella enterica with particular focus on the pathogenicity and virulence factors, host specificity and antimicrobial resistance including multidrug resistance. Veterinary World, 12(4), 504–521. doi:10.14202/vetworld.2019.504-521

Havelaar, A. H., Kirk, M. D., Torgerson, P. R., et al. (2015). World Health Organization global estimates and regional comparisons of the burden of foodborne disease in 2010. PLOS Medicine, 12(12), e1001923. doi:10.1371/journal.pmed.1001923

De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: Visualizing and processing long read sequencing data. Bioinformatics, 34(15), 2666–2669. doi:10.1093/bioinformatics/bty149

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37(5), 540–546. doi:10.1038/s41587-019-0072-8

Wick, R. R., & Holt, K. E. (2019). Benchmarking of long-read assemblers for prokaryote whole genome sequencing. F1000Research, 8, 2138. doi:10.12688/f1000research.21782.4

Sereika, M., Vetcher, A., et al. (2022). Oxford Nanopore  R10.4 long-read sequencing enables the generation of  near-finished bacterial genomes from pure cultures and metagenomes  without short-read or reference polishing. Nature Methods. doi:10.1038/s41592-022-01539-7

Lee, J., Kong, M. , Oh, J., et al. (2021). Comparative evaluation of nanopore polishing tools for microbial genome assembly and polishing strategies for  downstream analysis. Scientific Reports, 11, 20740. doi:10.1038/s41598-021-00178-w

Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. doi:10.1093/bioinformatics/btp352

Danecek, P., Bonfield, J. K., Liddle, J., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi:10.1093/gigascience/giab008

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. doi:10.1093/bioinformatics/bty191

Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072–1075. doi:10.1093/bioinformatics/btt086

Thorvaldsdóttir, H., Robinson, J. T., & Mesirov, J. P. (2013). Integrative Genomics Viewer (IGV): High-performance genomics data visualization and exploration. Briefings in Bioinformatics, 14(2), 178–192. doi:10.1093/bib/bbs017

Wick, R. R., Schultz, M. B., Zobel, J., & Holt, K. E. (2015). Bandage: Interactive visualization of de novo genome assemblies. Bioinformatics, 31(20), 3350–3352. doi:10.1093/bioinformatics/btv383

NCBI Sequence Read Archive (SRA). Oxford Nanopore R10 reads for Salmonella enterica. Run accession: SRR32410565.
