# Salmonella enterica Genome Assembly and Variant Analysis

**Course:** BINF*6110 - Bioinformatics  
**Author:** Yazan Al-Qawlaq  
**Date:** February 2026  
**Repository:** https://github.com/yazanalqawlaq-blip/Assignment-1---Binf6110

---

## Part 1: Introduction and Proposed Methods

### Introduction

Salmonella enterica is a Gram-negative, rod-shaped, facultative anaerobic bacterium that belongs to the family Enterobacteriaceae (Jajere, 2019). Salmonella enterica is a major cause of foodborne illness and remains an important public health concern worldwide (Havelaar et al., 2015). According to the WHO global burden estimates, most of the overall foodborne disease burden comes from diarrheal illnesses, and Salmonella enterica accounts for a meaningful share of that burden across the world (Havelaar et al., 2015). In this assignment, I will use Oxford Nanopore R10 long-read sequencing data to assemble a consensus genome for Salmonella enterica and then compare it to a reference genome from NCBI through alignment, variant calling, and visualization. The significance of this project is that assembling and comparing the genome provides a clear, practical way to understand how genome assembly quality and reference-based analysis affect what we conclude about variation in a clinically important pathogen.

A key challenge in bacterial genome assembly is resolving repeats and mobile genetic elements, which can fragment assemblies or lead to misassemblies if reads do not span repetitive regions (Wick and Holt, 2019). Long reads are useful because they can bridge repeats and support highly contiguous bacterial assemblies, including recovery of plasmids when present (Sereika et al., 2022). The main tradeoff is that nanopore reads can still contain base-level errors, and small indels can persist and show up as false differences during reference comparison if polishing and validation are not handled carefully (Sereika et al., 2022). Reference-based analysis also requires caution because alignment artifacts and low-support calls can look like real variants, so quality checks and visualization are needed before interpreting results (Thorvaldsdóttir et al., 2013).

For this workflow, I am choosing **Flye** for de novo assembly because it is designed for long-read data and is widely used for prokaryotic genomes (Kolmogorov et al., 2019). Compared to alternative assemblers, Flye offers distinct advantages for this dataset. **Canu**, while highly accurate, requires significantly more computational resources (memory and runtime) and can struggle with the high accuracy of R10 chemistry reads where Flye's algorithms are optimized (Wick and Holt, 2019). **Raven** is faster but produces less contiguous assemblies and has lower accuracy for resolving complex repetitive regions that are common in bacterial genomes (Vaser and Šikić, 2021). **wtdbg2** prioritizes speed over accuracy and is better suited for larger eukaryotic genomes rather than the precision required for bacterial genomics (Ruan and Li, 2020). A major advantage of Flye is that it often produces highly contiguous bacterial assemblies with accurate repeat resolution that are suitable for downstream comparison (Wick and Holt, 2019). A limitation of Flye is that it can require more memory than some alternatives, which matters for compute planning, though this is manageable on HPC systems (Wick and Holt, 2019). Using Nanopore R10 long reads is advantageous because long read length supports contiguity, but a limitation is that homopolymer-associated indels can still occur and need careful handling (Sereika et al., 2022). To improve base accuracy after assembly, I will polish using **Medaka** because polishing tools have been shown to improve microbial nanopore assemblies and reduce errors that affect downstream interpretation (Lee et al., 2021).

### Proposed Methods

Starting from the Oxford Nanopore R10 FASTQ reads for Salmonella enterica with expected accuracy Q20+ and an N50 of approximately 5 to 15 kb, I will use the dataset deposited on NCBI SRA under **SRR32410565** and document the exact download command and files in the repository. I will first run a basic read quality control using **NanoPlot** to summarize read length and quality distributions and confirm the data match expectations before assembly (De Coster et al., 2018). Based on standard practices for Nanopore assembly, I will apply read filtering with a minimum quality score of Q≥10 and minimum read length of ≥1000 bp using **NanoFilt** to remove low-quality sequences while preserving the long-read advantages for assembly contiguity.

I will then assemble the genome using **Flye v2.9.6** with the `--nano-hq` preset for high-quality Nanopore data and record all parameters including genome size estimate (5 Mb), thread settings (16 threads), and iteration count (Kolmogorov et al., 2019). After generating the draft assembly, I will polish the assembly using **Medaka v2.0.1** with the `r1041_e82_400bps_sup_v500` model to improve base-level accuracy and reduce residual errors that could otherwise appear as false differences during reference-based analysis (Lee et al., 2021). The model choice is specific to the R10.4.1 flow cell chemistry and super-accurate basecalling configuration used for this dataset.

For reference-based comparison, I will download the Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 reference genome from NCBI (accession **GCF_000006945.2**), then perform two separate alignment workflows: (1) raw filtered reads aligned to the reference using **minimap2 v2.28** with the `map-ont` preset for variant calling, and (2) polished assembly aligned to the reference using the `asm5` preset for structural comparison (Li, 2018). I will convert, sort, and index alignment files using **SAMtools v1.19** to generate coordinate-sorted indexed BAM files suitable for downstream analysis and visualization (Li et al., 2009).

I will call SNPs and short indels relative to the reference using **bcftools v1.22** mpileup and call functions with explicit filtering criteria (QUAL≥20, DP≥10) based on quality and depth before interpretation (Danecek et al., 2021). I will evaluate assembly statistics and reference-based metrics using **QUAST v5.2.0** to assess assembly contiguity, completeness, and alignment quality, and use **BUSCO v5.7.1** to evaluate genome completeness via conserved single-copy orthologs (Gurevich et al., 2013; Manni et al., 2021). I will use **IGV v2.19.7** to visualize representative regions and confirm read support for key variants, including screenshots as figures to support interpretation (Thorvaldsdóttir et al., 2013). All commands, parameters, figures, and software versions will be documented in the GitHub repository with multiple commits to clearly track the workflow progression.

---


### Project Structure
```
Assignment-1---Binf6110/
├── input_data/              # Reference genome (GCF_000006945.2)
├── output_files/            # Assembly, QC reports, alignments, variants
├── scripts/                 # Analysis pipeline 
└── README.md               # Complete documentation
```

### Final Methods

#### Data Acquisition and Initial Quality Assessment

Oxford Nanopore R10.4.1 sequencing data for *Salmonella enterica* isolate SRR32410565 was obtained directly from NCBI Sequence Read Archive using wget to download the SRA file from the AWS S3 mirror (`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565`), followed by conversion to FASTQ format using `fasterq-dump` (SRA Toolkit v3.2.1) with 12 threads and gzip compression. The reference genome for *Salmonella enterica* subsp. enterica serovar Typhimurium str. LT2 (accession GCF_000006945.2) was downloaded from NCBI RefSeq FTP (`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/`), decompressed, and renamed to `reference_fasta.fna` for consistency across the workflow. All data acquisition was performed on the Compute Canada Narval HPC cluster in the `/scratch` directory to leverage improved network connectivity and storage performance.

#### De Novo Genome Assembly

De novo assembly was performed using **Flye v2.9.6** (Kolmogorov et al., 2019) with the `--nano-hq` preset optimized for high-accuracy R10 chemistry Nanopore reads, a genome size parameter of `5m` (5 megabases) based on expected *Salmonella enterica* genome size, and 16 threads for parallel processing (`--threads 16`). The Flye algorithm constructs a repeat graph to resolve repetitive genomic regions and produces highly contiguous assemblies suitable for bacterial genomes. Assembly output was directed to `output_files/flye_assembly/`, with the primary assembly file `assembly.fasta` retained for downstream polishing and analysis.

#### Assembly Polishing

The raw Flye assembly was polished using **Medaka v2.0.1** (Oxford Nanopore Technologies) to correct systematic basecalling errors characteristic of nanopore sequencing. The consensus polishing was performed using the `medaka_consensus` command with model `r1041_e82_400bps_sup_v5.0.0`, specifically trained for R10.4.1 flow cells with super-accurate basecalling, input raw FASTQ reads (`-i SRR32410565.fastq.gz`), draft assembly (`-d output_files/flye_assembly/assembly.fasta`), output directory (`-o output_files/medaka_out`), and 16 threads (`-t 16`). Medaka applies deep learning-based error correction by realigning reads to the draft assembly and generating a consensus sequence with improved base-level accuracy.

#### Assembly Quality Control

Assembly quality was assessed using two complementary approaches. **QUAST v5.2.0** (Gurevich et al., 2013) was used to compute assembly statistics including contig count, total length, N50, largest contig size, GC content, and reference-based metrics such as genome fraction, misassemblies, and error rates (mismatches and indels per 100 kbp). QUAST was executed on both the raw Flye assembly (`output_files/quast_raw`) and the Medaka-polished assembly (`output_files/quast_polished`) with the reference genome provided (`-r input_data/reference_fasta.fna`) and 8 threads for parallel processing. **BUSCO v5.7.1** (Manni et al., 2021) was used to evaluate genome completeness by identifying conserved single-copy orthologs from the bacteria_odb10 lineage dataset (2024-01-08 release). BUSCO analysis was performed on both raw and polished assemblies in genome mode (`-m genome`) using offline mode with a locally downloaded lineage database (`-l busco_downloads/bacteria_odb10 --offline`) to avoid network connectivity issues, with 8 CPUs (`-c 8`) and forced overwrite of existing results (`-f`). BUSCO reports completeness as percentages of complete, duplicated, fragmented, and missing genes from the expected gene set.

#### Reference-Based Alignment

Two alignment strategies were implemented using **minimap2 v2.28** (Li, 2018). First, raw reads were aligned to the reference genome using the `map-ont` preset optimized for Nanopore reads (`-a -x map-ont`), 16 threads (`-t 16`), with output piped directly to **SAMtools v1.20** (Li et al., 2009) for conversion to BAM format, coordinate sorting with 16 threads (`sort -@ 16`), and indexing to generate `output_files/reads_to_ref.sorted.bam` and its index. This alignment enables base-level variant calling by providing read support across the reference genome. Second, the polished assembly was aligned to the reference using the `asm5` preset for assembly-to-assembly comparison (`-a -x asm5 -t 16`), similarly processed through SAMtools to produce `output_files/assembly_to_ref.sorted.bam`, which facilitates structural variation detection and assessment of assembly contiguity relative to the reference.

#### Variant Calling and Filtering

Variant calling was performed using **bcftools v1.22** (Danecek et al., 2021) mpileup and call functions. The `bcftools mpileup` command generated pileup format from the reads-to-reference BAM file with parameters: `-f` for reference genome, `-Ou` for uncompressed BCF output, `-Q 7` for minimum base quality of 7, `-q 0` for minimum mapping quality of 0, `-d 10000` for maximum per-sample depth of 10,000, and `-a FORMAT/DP,FORMAT/AD` to annotate depth and allelic depth. The output was piped to `bcftools call` with `-m` for multiallelic caller, `--ploidy 1` for haploid bacterial genome, and `-Ov` for uncompressed VCF output to generate `raw_calls.vcf`. Variants were filtered using `bcftools filter` with the expression `-e 'QUAL<20 || FORMAT/DP<10'` to remove low-quality calls (QUAL<20) and low-coverage variants (depth<10), producing `filtered.vcf`. Final variant-only VCF was generated using `bcftools view -v snps,indels` to retain only SNPs and indels in `variants_only.vcf`, with summary statistics computed using `bcftools stats` and saved to `variant_stats.txt`.

#### Visualization and Data Summary

Genome comparison visualization was generated using a custom Python script executed with **Python 3.11** and the scipy-stack module (including matplotlib). The script parsed both the reference and polished assembly FASTA files to extract contig lengths, then generated a two-panel horizontal barplot comparing genome structures with color-coded contigs, axis labels, titles, legends, and gridlines for clarity. The output figure was saved as `output_files/genome_comparison.png` at 300 DPI resolution. Additionally, alignment statistics including coverage depth and mapping statistics were extracted from BAM files using SAMtools depth and stats commands for inclusion in the assembly summary report.

### Software Versions

| Tool | Version | Purpose |
|------|---------|---------|
| SRA Toolkit | 3.2.1 | Data download and conversion |
| Flye | 2.9.6 | De novo genome assembly |
| Medaka | 2.0.1 | Consensus polishing |
| minimap2 | 2.28 | Read and assembly alignment |
| SAMtools | 1.20 | BAM file processing |
| bcftools | 1.22 | Variant calling and filtering |
| QUAST | 5.2.0 | Assembly quality assessment |
| BUSCO | 5.7.1 | Genome completeness evaluation |
| Python | 3.11 | Data visualization |
| scipy-stack | 2024b | Scientific Python libraries |

### Results

#### Assembly Statistics

| Metric | Raw Assembly (Flye) | Polished Assembly (Medaka) | Reference (LT2) |
|--------|---------------------|----------------------------|-----------------|
| Total Length (bp) | 5,105,879 | 5,105,872 | 4,951,383 |
| Number of Contigs | 3 | 3 | 1 |
| Largest Contig (bp) | 3,318,776 | 3,318,769 | 4,857,432 |
| N50 (bp) | 3,318,776 | 3,318,769 | 4,857,432 |
| GC Content (%) | 52.20 | 52.20 | 52.24 |
| Mean Coverage | 127× | 127× | — |

#### Quality Assessment

| Metric | Raw Assembly | Polished Assembly |
|--------|--------------|-------------------|
| BUSCO Completeness | 99.2% (123/124) | 98.4% (122/124) |
| BUSCO Single-Copy | 99.2% (123) | 98.4% (122) |
| BUSCO Duplicated | 0.0% (0) | 0.0% (0) |
| BUSCO Fragmented | 0.0% (0) | 0.0% (0) |
| BUSCO Missing | 0.8% (1) | 1.6% (2) |
| Genome Fraction (%) | 95.669 | 95.669 |
| Mismatches per 100 kbp | 27.41 | 27.11 |
| Indels per 100 kbp | 3.83 | 3.71 |
| Misassemblies | 25 | 25 |

#### Variant Analysis

| Variant Type | Count |
|--------------|-------|
| Single Nucleotide Polymorphisms (SNPs) | 9,337 |
| Insertions and Deletions (indels) | 70 |
| **Total Variants** | **9,407** |



![Image](https://github.com/user-attachments/assets/aa456df2-100a-4ec7-9f5b-059f909c3cc3)
Figure 1. IGV Visualization of Variant Positions in the Salmonella enterica Assembly Aligned to the Reference Genome. The polished assembly was aligned to the S. enterica Typhimurium LT2 reference genome using minimap2 and visualized in IGV. The top track shows the reference sequence with amino acid translation, and gene annotations are displayed below. Gray alignment blocks represent the assembled genome, with SNPs highlighted as colored blocks at variant positions. The region shown contains multiple variants, including several silent substitutions that do not alter amino acid coding. Notably, one SNP results in a His→Asn substitution (CAC→AAC), changing histidine to asparagine at this position. Coverage across the region is uniform with no evidence of structural variation, and observed differences are limited to single-nucleotide polymorphisms consistent with strain-level divergence from the reference. 

---

## Discussion

The genome assembly of *Salmonella enterica* isolate SRR32410565 from Oxford Nanopore R10.4.1 long reads produced a highly contiguous draft genome. The raw Flye assembly resulted in 3 contigs with a total length of 5,105,879 bp and an N50 of 3,318,776 bp. This level of contiguity is consistent with what long reads are good at in bacteria, since long reads can bridge repetitive regions that often fragment short-read assemblies (Kolmogorov et al., 2019). At the same time, contiguity does not automatically prove correctness, so the best interpretation is that the data supported a strong structural backbone, but structural validation still relies on reference comparison, read mapping patterns, and graph/visual inspection rather than N50 alone (Gurevich et al., 2013; Thorvaldsdóttir et al., 2013).

Medaka polishing did not change the assembly structure (still 3 contigs) and only changed total length by 7 bp, which makes sense because the starting chemistry and basecalling already have relatively high accuracy. Even though the improvements in mismatches and indels per 100 kbp were modest (mismatches 27.41 → 27.11 per 100 kbp; indels 3.83 → 3.71 per 100 kbp), polishing still matters biologically because residual nanopore errors, especially small indels, can show up as false differences when comparing to a reference if they are not corrected as much as possible (Lee et al., 2021). BUSCO completeness was extremely high in both assemblies (99.2% raw; 98.4% polished), supporting that the assembly captured essentially the full expected gene set. The slight decrease after polishing is normal because small sequence corrections can affect how conserved genes are detected or how gene boundaries are predicted without meaning real gene loss (Manni et al., 2021).

When comparing to the *S. enterica* LT2 reference, the overall size and alignment coverage provide important biological context. My assembly length (~5.106 Mb) is larger than the LT2 reference (~4.857 Mb chromosome plus a ~94 kb virulence plasmid reported for LT2), and QUAST reported a 95.7% genome fraction aligned to the reference. A larger assembly size plus <100% genome fraction alignment is consistent with real strain-to-strain differences in accessory genome content such as prophages and plasmids. In *Salmonella*, these mobile elements are a major source of present/absent genomic differences between isolates and are often tied to adaptation and phenotype (Mottawea et al., 2018). So here, unaligned regions are not automatically an assembly problem; they can represent DNA that is truly present in the isolate but absent (or highly divergent) relative to LT2, and interpreting this correctly requires viewing the alignment patterns and checking whether those regions show consistent long-read support (Thorvaldsdóttir et al., 2013).

The 25 misassemblies reported by QUAST (affecting 2 contigs) should also be interpreted carefully because QUAST defines misassembly relative to the reference coordinate structure. That means these calls can reflect either genuine biological rearrangements or true assembly errors, particularly near repeats and rRNA operons (Gurevich et al., 2013; Wick and Holt, 2019). Since *Salmonella* genomes are shaped heavily by mobile genetic elements, it is plausible that at least some breakpoint signals represent real structural differences versus LT2 rather than purely technical mistakes (Mottawea et al., 2018). This is why visualization is important before making biological claims about rearrangements (Thorvaldsdóttir et al., 2013).

Variant calling detected 9,337 SNPs and 70 indels relative to LT2, indicating the isolate is clearly *S. enterica* but genomically distinct from the LT2 laboratory reference. This scale of variation fits what is expected across different *Salmonella* strain backgrounds, where both point mutations and horizontal gene transfer contribute to divergence (Jajere, 2019). The next thing that really matters biologically is where the variants are and what type they are. Even silent SNPs can still matter because, even if they don't change the amino acid, they can affect codon usage, mRNA structure, and stability (Brandis & Hughes, 2016; Mittal et al., 2018). So, it's not always safe to assume silent mutations are neutral (Bailey et al., 2021). Indels also need caution with nanopore data because indel errors are a known remaining challenge, especially in homopolymer-associated contexts, and false indels can inflate apparent gene-disrupting changes if not validated (Lee et al., 2021; Danecek et al., 2021). That is why quality/depth filtering and IGV confirmation are necessary before interpreting any variants as biologically meaningful differences that could affect virulence, resistance, or metabolism (Thorvaldsdóttir et al., 2013; Jajere, 2019).

Overall, the assembly quality metrics (high contiguity and high BUSCO completeness) support that this workflow produced a strong draft genome suitable for reference comparison, while the reference-based results suggest real strain-level differences from LT2 that likely include both sequence-level changes (SNPs/indels) and presence/absence differences tied to the accessory genome. This is the main biological takeaway: long-read assembly can rapidly produce a near-complete genome, but interpreting differences versus a reference needs careful handling so that variation is not over-interpreted without validation and context from the known *Salmonella* genome structure.

## Conclusion

This study successfully generated a high-quality draft genome assembly of a *Salmonella enterica* isolate using Oxford Nanopore R10.4.1 long-read sequencing. The final assembly was highly contiguous and essentially complete by BUSCO (98.4–99.2%), supporting that the genome was captured well enough for comparative analysis. Comparison to the LT2 reference revealed substantial divergence, including thousands of SNPs and additional differences consistent with strain-to-strain accessory genome variation, which is a major driver of biological diversity in *Salmonella* (McClelland et al., 2001). Finally, while many variants may be silent, silent mutations can still have biological effects in bacteria, and indels require particular caution with nanopore data, reinforcing the importance of filtering and visualization before drawing strong conclusions (Bailey et al., 2021; Lee et al., 2021). Taken together, these results show that nanopore-only sequencing can support rapid bacterial genome assembly and meaningful reference comparison, as long as interpretation is grounded in genome biology and validated with appropriate quality checks.

---

## References

Bailey, S. F., Alonso Morales, L. A., & Kassen, R. (2021). Effects of synonymous mutations beyond codon bias: The evidence for adaptive synonymous substitutions from microbial evolution experiments. *Genome Biology and Evolution*, 13(9), evab141. https://doi.org/10.1093/gbe/evab141

Brandis, G., & Hughes, D. (2016). The selective advantage of synonymous codon usage bias in Salmonella. *PLoS Genetics*, 12(3), e1005926. https://doi.org/10.1371/journal.pgen.1005926

Danecek, P., Bonfield, J. K., Liddle, J., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

De Coster, W., D'Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: Visualizing and processing long-read sequencing data. *Bioinformatics*, 34(15), 2666–2669. https://doi.org/10.1093/bioinformatics/bty149

Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086

Havelaar, A. H., Kirk, M. D., Torgerson, P. R., et al. (2015). World Health Organization global estimates and regional comparisons of the burden of foodborne disease in 2010. *PLOS Medicine*, 12(12), e1001923. https://doi.org/10.1371/journal.pmed.1001923

Jajere, S. M. (2019). A review of Salmonella enterica with particular focus on the pathogenicity and virulence factors, host specificity and antimicrobial resistance including multidrug resistance. *Veterinary World*, 12(4), 504–521. https://doi.org/10.14202/vetworld.2019.504-521

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37(5), 540–546. https://doi.org/10.1038/s41587-019-0072-8

Lee, J., Kong, M., Oh, J., et al. (2021). Comparative evaluation of nanopore polishing tools for microbial genome assembly and polishing strategies for downstream analysis. *Scientific Reports*, 11, 20740. https://doi.org/10.1038/s41598-021-00178-w

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A., & Zdobnov, E. M. (2021). BUSCO update: Novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. *Molecular Biology and Evolution*, 38(10), 4647–4654. https://doi.org/10.1093/molbev/msab199

McClelland, M., Sanderson, K. E., Spieth, J., et al. (2001). Complete genome sequence of Salmonella enterica serovar Typhimurium LT2. *Nature*, 413, 852–856. https://doi.org/10.1038/35101614

Mittal, P., Brindle, J., Stephen, J., Plotkin, J. B., & Kudla, G. (2018). Codon usage influences fitness through RNA toxicity. *Proceedings of the National Academy of Sciences of the United States of America*, 115(34), 8639–8644. https://doi.org/10.1073/pnas.1810022115

Mottawea, W., Duceppe, M.-O., Dallaire-Dufresne, S., et al. (2018). Salmonella enterica prophage sequence profiles reflect genome diversity and can be used for high discrimination subtyping. *Frontiers in Microbiology*, 9, 836. https://doi.org/10.3389/fmicb.2018.00836

Ruan, J., & Li, H. (2020). Fast and accurate long-read assembly with wtdbg2. *Nature Methods*, 17, 155–158. https://doi.org/10.1038/s41592-019-0669-3

Sereika, M., Kirkegaard, R. H., Karst, S. M., Michaelsen, T. Y., Sørensen, E. A., Wollenberg, R. D., & Albertsen, M. (2022). Oxford Nanopore R10.4 long-read sequencing enables the generation of near-finished bacterial genomes from pure cultures and metagenomes without short-read or reference polishing. *Nature Methods*, 19(7), 823–826. https://doi.org/10.1038/s41592-022-01539-7

Thorvaldsdóttir, H., Robinson, J. T., & Mesirov, J. P. (2013). Integrative Genomics Viewer (IGV): High-performance genomics data visualization and exploration. *Briefings in Bioinformatics*, 14(2), 178–192. https://doi.org/10.1093/bib/bbs017

Vaser, R., & Šikić, M. (2021). Time- and memory-efficient genome assembly with Raven. *Nature Computational Science*, 1, 332–336. https://doi.org/10.1038/s43588-021-00073-4

Wick, R. R., & Holt, K. E. (2019). Benchmarking of long-read assemblers for prokaryote whole genome sequencing. *F1000Research*, 8, 2138. https://doi.org/10.12688/f1000research.21782.4

NCBI Sequence Read Archive (SRA). Oxford Nanopore R10 reads for Salmonella enterica. Run accession: SRR32410565.
