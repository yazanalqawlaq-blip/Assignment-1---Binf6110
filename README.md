# Salmonella enterica Genome Assembly and Variant Analysis

**Course:** BINF*6110 - Bioinformatics  
**Author:** Yazan Al-Qawlaq  
**Date:** February 2026  
**Repository:** https://github.com/yazanalqawlaq-blip/Assignment-1---Binf6110

---


### Introduction

Salmonella enterica is a Gram-negative, rod-shaped, facultative anaerobic bacterium that belongs to the family Enterobacteriaceae (Jajere, 2019). Salmonella enterica is a major cause of foodborne illness and remains an important public health concern worldwide (Havelaar et al., 2015). According to the WHO global burden estimates, most of the overall foodborne disease burden comes from diarrheal illnesses, and Salmonella enterica accounts for a meaningful share of that burden across the world (Havelaar et al., 2015). In this assignment, I will use Oxford Nanopore R10 long-read sequencing data to assemble a consensus genome for Salmonella enterica and then compare it to a reference genome from NCBI through alignment, variant calling, and visualization. The significance of this project is that assembling and comparing the genome provides a clear, practical way to understand how genome assembly quality and reference-based analysis affect what we conclude about variation in a clinically important pathogen.

A key challenge in bacterial genome assembly is resolving repeats and mobile genetic elements, which can fragment assemblies or lead to misassemblies if reads do not span repetitive regions (Wick and Holt, 2019). Long reads are useful because they can bridge repeats and support highly contiguous bacterial assemblies, including recovery of plasmids when present (Sereika et al., 2022). The main tradeoff is that nanopore reads can still contain base-level errors, and small indels can persist and show up as false differences during reference comparison if polishing and validation are not handled carefully (Sereika et al., 2022). Reference-based analysis also requires caution because alignment artifacts and low-support calls can look like real variants, so quality checks and visualization are needed before interpreting results (Thorvaldsdóttir et al., 2013).

For this workflow, I am choosing **Flye** for de novo assembly because it is designed for long-read data and is widely used for prokaryotic genomes (Kolmogorov et al., 2019). Compared to alternative assemblers, Flye offers distinct advantages for this dataset. **Canu**, while highly accurate, requires significantly more computational resources (memory and runtime) and can struggle with the high accuracy of R10 chemistry reads where Flye's algorithms are optimized (Wick and Holt, 2019). **Raven** is faster but produces less contiguous assemblies and has lower accuracy for resolving complex repetitive regions that are common in bacterial genomes (Vaser and Šikić, 2021). **wtdbg2** prioritizes speed over accuracy and is better suited for larger eukaryotic genomes rather than the precision required for bacterial genomics (Ruan and Li, 2020). A major advantage of Flye is that it often produces highly contiguous bacterial assemblies with accurate repeat resolution that are suitable for downstream comparison (Wick and Holt, 2019). A limitation of Flye is that it can require more memory than some alternatives, which matters for compute planning, though this is manageable on HPC systems (Wick and Holt, 2019). Using Nanopore R10 long reads is advantageous because long read length supports contiguity, but a limitation is that homopolymer-associated indels can still occur and need careful handling (Sereika et al., 2022). To improve base accuracy after assembly, I will polish using **Medaka** because polishing tools have been shown to improve microbial nanopore assemblies and reduce errors that affect downstream interpretation (Lee et al., 2021).


---


### Project Structure
```
Assignment-1---Binf6110/
├── input_data/              # Reference genome (GCF_000006945.2)
├── output_files/            # Assembly, QC reports, alignments, variants
├── scripts/                 # Analysis pipeline 
└── README.md               # Complete documentation
```
---

## Final Methods

### Data Acquisition

Oxford Nanopore R10.4.1 sequencing data for a *Salmonella enterica* isolate (SRA accession SRR32410565) was downloaded directly from the NCBI Sequence Read Archive AWS S3 mirror using `wget`, then converted to FASTQ format using `fasterq-dump` (SRA Toolkit v3.2.1) with 4 threads and compressed with `gzip`. The reference genome for *Salmonella enterica* subsp. enterica serovar Typhimurium str. LT2 (accession GCF_000006945.2) was downloaded from NCBI RefSeq, decompressed, and stored as `reference_fasta.fna`. All analyses were performed on the Compute Canada Narval HPC cluster using SLURM job scheduling and Apptainer containers for software reproducibility.

No pre-assembly read filtering or quality control was performed, as the R10.4.1 chemistry with super-accurate basecalling produces reads at Q20+ accuracy, and Flye is designed to handle uncorrected long reads directly (Kolmogorov et al., 2019).

### De Novo Genome Assembly

De novo assembly was performed using Flye v2.9.6 (Kolmogorov et al., 2019) with the `--nano-hq` preset, which is optimized for high-accuracy Nanopore reads with error rates below 5%. Assembly was run with 16 threads (`--threads 16`). Flye was selected over alternative long-read assemblers such as Canu and miniasm because it constructs repeat graphs that effectively resolve repetitive regions in bacterial genomes while maintaining computational efficiency, and it produces polished consensus sequences as part of its core pipeline (Kolmogorov et al., 2019; Wick & Holt, 2021).

### Assembly Polishing

The raw Flye assembly was polished using Medaka v2.2.0 (Oxford Nanopore Technologies) to correct residual basecalling errors. Polishing was performed with `medaka_consensus` using the `r1041_e82_400bps_sup_v5.0.0` model, which is specifically trained for R10.4.1 flow cell chemistry with super-accurate basecalling. The original FASTQ reads were provided as input alongside the draft assembly, and polishing was run with 16 threads. Medaka applies neural network-based error correction by aligning reads to the draft assembly and generating an improved consensus sequence.

### Assembly Quality Assessment

Assembly quality was evaluated at two stages—before and after polishing—using two complementary tools. QUAST v5.2.0 (Gurevich et al., 2013) was used to compute contiguity and reference-based metrics including contig count, total length, N50, genome fraction, misassemblies, and error rates (mismatches and indels per 100 kbp), with the reference genome provided via the `-r` flag and 8 threads for parallel processing. BUSCO v5.7.1 (Manni et al., 2021) was used to assess genome completeness by identifying conserved single-copy orthologs from the `bacteria_odb10` lineage dataset (2024-01-08 release) in genome mode (`-m genome`), run in offline mode with a locally downloaded database to avoid network connectivity issues on the HPC cluster.

### Reference-Based Alignment

Two alignment strategies were implemented using minimap2 v2.28 (Li, 2018) and processed with SAMtools v1.20 (Li et al., 2009). First, raw reads were aligned to the reference genome using the `map-ont` preset (`-a -x map-ont -t 16`), with output piped directly through SAMtools for BAM conversion, coordinate sorting (`sort -@ 16`), and indexing. This read-level alignment provides the per-base coverage required for variant calling. Second, the polished assembly was aligned to the reference using the `asm5` preset (`-a -x asm5 -t 16`) for whole-genome structural comparison, similarly processed into a sorted, indexed BAM file.

### Variant Calling and Filtering

Variant calling was performed on the reads-to-reference alignment using bcftools v1.22 (Danecek et al., 2021). The `bcftools mpileup` command was run with a minimum base quality of 10 (`-Q 10`), minimum mapping quality of 5 (`-q 5`), maximum per-sample depth of 8,000 (`-d 8000`), and annotation of per-sample depth and allelic depth (`--annotate FORMAT/DP,FORMAT/AD,INFO/AD`). Output was piped to `bcftools call` with the multiallelic caller (`-m`) and haploid ploidy (`--ploidy 1`). Variants were then normalized using `bcftools norm` to left-align indels and split multi-allelic sites before filtering. Quality filtering was applied using `bcftools filter` with the expression `QUAL<30 || FORMAT/DP<15` to exclude low-confidence and low-coverage calls. The final variant set was extracted using `bcftools view -v snps,indels`, and summary statistics were generated with `bcftools stats`.

### Visualization

Three figures were generated to characterize the assembly and its relationship to the reference genome. An IGV screenshot (Thorvaldsdóttir et al., 2013) was produced by loading the reference genome, gene annotations (GFF), polished assembly alignment, and variant calls into IGV v2.19.7 on a local machine to visualize variant positions with gene context. A Circos-style circular genome plot was generated using the Circa web application (https://circa.omgenomics.com/) to display the genome-wide distribution of 9,407 variants across the reference chromosome and plasmid. A grouped bar chart comparing key quality metrics (N50, total length, BUSCO completeness, mismatches per 100 kbp, and indels per 100 kbp) between the raw and polished assemblies was generated using a custom Python script with matplotlib, executed on Narval via the `scipy-stack/2024b` module.

### Computational Environment

All analyses were executed on the Compute Canada Narval HPC cluster. Bioinformatics tools (SRA Toolkit, Flye, Medaka, QUAST, BUSCO) were run via Apptainer v1.2.4 containers to ensure software reproducibility. Alignment and variant calling tools (minimap2 v2.28, SAMtools v1.20, bcftools v1.22) were loaded as environment modules (`StdEnv/2023`). All scripts, parameters, and outputs are documented in the project GitHub repository.

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

---

![Image](https://github.com/user-attachments/assets/aa456df2-100a-4ec7-9f5b-059f909c3cc3)

**Figure 1. IGV Visualization of Variant Positions in the *Salmonella enterica* Assembly Aligned to the Reference Genome.** The polished assembly was aligned to the *S. enterica* Typhimurium LT2 reference genome using minimap2 and visualized in IGV. The top track shows the reference sequence with amino acid translation, and gene annotations are displayed below. Gray alignment blocks represent the assembled genome, with SNPs highlighted as colored blocks at variant positions. The region shown contains multiple variants, including several silent substitutions that do not alter amino acid coding. Notably, one SNP results in a His→Asn substitution (CAC→AAC), changing histidine to asparagine at this position. Coverage across the region is uniform with no evidence of structural variation, and observed differences are limited to single-nucleotide polymorphisms consistent with strain-level divergence from the reference.


---

![Image](https://github.com/user-attachments/assets/4ea114f9-44df-4554-9dea-84d21d437988)

**Figure 2. Circos Plot Showing Genome-Wide Variant Distribution in *Salmonella enterica* Assembly.** The circular plot displays the reference genome structure of *S. enterica* Typhimurium LT2 with the main chromosome (NC_003197.2) and virulence plasmid (NC_003277.2) shown as outer segments. The inner track shows the distribution of 9,407 variants (9,337 SNPs and 70 indels) identified by aligning the polished assembly to the reference genome. Variant positions are displayed as radial lines, with density reflecting the frequency of genetic differences across genomic regions. The widespread distribution of variants around the genome indicates strain-level divergence rather than localized assembly errors or hotspots, supporting that the observed differences represent genuine biological variation between the sequenced isolate and the LT2 reference strain. Circos plot was generated using Circa (https://omgenomics.com/circa).


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
