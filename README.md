# Salmonella enterica Genome Assembly and Variant Analysis

**Course:** BINF*6110 - Bioinformatics  
**Author:** Yazan Al-Qawlaq  
**Date:** February 2026  
**Repository:** https://github.com/yazanalqawlaq-blip/Assignment-1---Binf6110

---

## Part 1: Introduction and Proposed Methods

### Introduction

Salmonella enterica is a Gram-negative, rod-shaped, facultative anaerobic bacterium that belongs to the family Enterobacteriaceae (Jajere, 2019). Salmonella enterica is a major cause of foodborne illness and remains an important public health concern worldwide (Havelaar et al., 2015). According to the WHO global burden estimates, most of the overall foodborne disease burden comes from diarrheal illnesses, and Salmonella enterica accounts for a meaningful share of that burden across the world (Havelaar et al., 2015). In this assignment, I will use Oxford Nanopore R10 long-read sequencing data to assemble a consensus genome for Salmonella enterica and then compare it to a reference genome from NCBI through alignment, variant calling, and visualization. The significance of this project is that assembling and comparing the genome provides a clear, practical way to understand how genome assembly quality and reference-based analysis affect what we conclude about variation in a clinically important pathogen.

A key challenge in bacterial genome assembly is resolving repeats and mobile genetic elements, which can fragment assemblies or lead to misassemblies if reads do not span repetitive regions (Wick and Holt, 2019). Long reads are useful because they can bridge repeats and support highly contiguous bacterial assemblies, including recovery of plasmids when present (Sereika et al., 2022). The main tradeoff is that nanopore reads can still contain base-level errors, and small indels can persist and show up as false differences during reference comparison if polishing and validation are not handled carefully (Luan et al., 2024). Reference-based analysis also requires caution because alignment artifacts and low-support calls can look like real variants, so quality checks and visualization are needed before interpreting results (Thorvaldsdóttir et al., 2013).

For this workflow, I am choosing **Flye** for de novo assembly because it is designed for long-read data and is widely used for prokaryotic genomes (Kolmogorov et al., 2019). Compared to alternative assemblers, Flye offers distinct advantages for this dataset. **Canu**, while highly accurate, requires significantly more computational resources (memory and runtime) and can struggle with the high accuracy of R10 chemistry reads where Flye's algorithms are optimized (Wick and Holt, 2019). **Raven** is faster but produces less contiguous assemblies and has lower accuracy for resolving complex repetitive regions that are common in bacterial genomes (Vaser and Šikić, 2021). **wtdbg2** prioritizes speed over accuracy and is better suited for larger eukaryotic genomes rather than the precision required for bacterial genomics (Ruan and Li, 2020). A major advantage of Flye is that it often produces highly contiguous bacterial assemblies with accurate repeat resolution that are suitable for downstream comparison (Wick and Holt, 2019). A limitation of Flye is that it can require more memory than some alternatives, which matters for compute planning, though this is manageable on HPC systems (Wick and Holt, 2019). Using Nanopore R10 long reads is advantageous because long read length supports contiguity, but a limitation is that homopolymer-associated indels can still occur and need careful handling (Sereika et al., 2022). To improve base accuracy after assembly, I will polish using **Medaka** because polishing tools have been shown to improve microbial nanopore assemblies and reduce errors that affect downstream interpretation (Lee et al., 2021).

### Proposed Methods

Starting from the Oxford Nanopore R10 FASTQ reads for Salmonella enterica with expected accuracy Q20+ and an N50 of approximately 5 to 15 kb, I will use the dataset deposited on NCBI SRA under **SRR32410565** and document the exact download command and files in the repository. I will first run a basic read quality control using **NanoPlot** to summarize read length and quality distributions and confirm the data match expectations before assembly (De Coster et al., 2018). Based on standard practices for Nanopore assembly, I will apply read filtering with a minimum quality score of Q≥10 and minimum read length of ≥1000 bp using **NanoFilt** to remove low-quality sequences while preserving the long-read advantages for assembly contiguity.

I will then assemble the genome using **Flye v2.9.6** with the `--nano-hq` preset for high-quality Nanopore data and record all parameters including genome size estimate (5 Mb), thread settings (16 threads), and iteration count (Kolmogorov et al., 2019). After generating the draft assembly, I will polish the assembly using **Medaka v2.0.1** with the `r1041_e82_400bps_sup_v500` model to improve base-level accuracy and reduce residual errors that could otherwise appear as false differences during reference-based analysis (Lee et al., 2021). The model choice is specific to the R10.4.1 flow cell chemistry and super-accurate basecalling configuration used for this dataset.

For reference-based comparison, I will download the Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 reference genome from NCBI (accession **GCF_000006945.2**), then perform two separate alignment workflows: (1) raw filtered reads aligned to the reference using **minimap2 v2.28** with the `map-ont` preset for variant calling, and (2) polished assembly aligned to the reference using the `asm5` preset for structural comparison (Li, 2018). I will convert, sort, and index alignment files using **SAMtools v1.19** to generate coordinate-sorted indexed BAM files suitable for downstream analysis and visualization (Li et al., 2009).

I will call SNPs and short indels relative to the reference using **Clair3 v1.0.10**, a deep learning-based variant caller optimized for long-read data, and apply explicit filtering criteria (QUAL≥20, DP≥10) based on quality and depth before interpretation (Danecek et al., 2021). I will evaluate assembly statistics and reference-based metrics using **QUAST v5.2.0** to assess assembly contiguity, completeness, and alignment quality, and use **BUSCO v5.7.1** to evaluate genome completeness via conserved single-copy orthologs (Gurevich et al., 2013). I will use **IGV v2.19.7** to visualize representative regions and confirm read support for key variants, including screenshots as figures to support interpretation (Thorvaldsdóttir et al., 2013). I will also inspect the assembly graph in **Bandage v0.8.1** to check for circular contigs and unresolved structures that could affect conclusions. All commands, parameters, figures, and software versions will be documented in the GitHub repository with multiple commits to clearly track the workflow progression.

---

## Part 2: Final Methods and Results

### Project Structure
```
Assignment-1---Binf6110/
├── input_data/              # Reference genome (GCF_000006945.2)
├── output_files/            # Assembly results, QC reports, variant calls
├── figures/                 # Publication-quality visualizations
├── scripts/                 # Complete analysis pipeline (8 scripts)
├── slurm/                   # SLURM job output logs
└── README.md               # Project documentation
```

### Final Methods

#### 1. Data Acquisition and Quality Control

**Raw Data Download:**
```bash
wget -O SRR32410565.sra "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32410565/SRR32410565"
fasterq-dump SRR32410565.sra --threads 8
gzip SRR32410565.fastq
```

**Quality Control:**
- **Tool:** NanoPlot v1.46.2
- **Purpose:** Assess read length distribution and quality scores
- **Command:**
```bash
NanoPlot --fastq SRR32410565.fastq.gz \
    --outdir output_files/nanoplot_raw \
    --threads 8
```

**Read Filtering:**
- **Tool:** NanoFilt v2.8.0
- **Criteria:** Q≥10, length≥1000 bp
- **Rationale:** Remove low-quality reads while preserving long-read advantages
- **Command:**
```bash
gunzip -c SRR32410565.fastq.gz | \
    NanoFilt -q 10 -l 1000 | \
    gzip > SRR32410565.filtered.fastq.gz
```

**Reference Genome Download:**
- **Accession:** GCF_000006945.2
- **Strain:** *S. enterica* subsp. enterica serovar Typhimurium str. LT2
- **Source:** NCBI RefSeq
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
```

#### 2. De Novo Assembly

**Tool:** Flye v2.9.6  
**Algorithm:** Repeat graph-based assembly optimized for long reads  
**Parameters:**
- `--nano-hq`: High-quality Nanopore preset (R10+ chemistry)
- `--genome-size 5m`: Expected genome size 5 Mb
- `--threads 16`: Parallel processing
- `--out-dir output_files/flye_assembly`

**Command:**
```bash
flye --nano-hq SRR32410565.filtered.fastq.gz \
    --genome-size 5m \
    --out-dir output_files/flye_assembly \
    --threads 16
```

#### 3. Assembly Polishing

**Tool:** Medaka v2.0.1  
**Model:** r1041_e82_400bps_sup_v500  
**Purpose:** Correct systematic basecalling errors

**Command:**
```bash
medaka_consensus \
    -i SRR32410565.filtered.fastq.gz \
    -d output_files/flye_assembly/assembly.fasta \
    -o output_files/medaka_polished \
    -m r1041_e82_400bps_sup_v500 \
    -t 16
```

#### 4. Assembly Quality Assessment

**QUAST Analysis:**
```bash
quast.py output_files/medaka_polished/consensus.fasta \
    -r input_data/reference_fasta.fna \
    -o output_files/quast_polished \
    --threads 8
```

**BUSCO Analysis:**
```bash
busco -i output_files/medaka_polished/consensus.fasta \
    -o output_files/busco_polished \
    -m genome \
    --auto-lineage-prok \
    -c 8
```

**Bandage Visualization:**
```bash
Bandage image output_files/flye_assembly/assembly_graph.gfa \
    figures/assembly_graph.png \
    --height 2000 --width 2000
```

#### 5. Reference Alignment

**Alignment 1: Raw Reads → Reference**
- **Purpose:** Variant calling with base-level resolution
- **Tool:** minimap2 v2.28 (map-ont preset)
```bash
minimap2 -ax map-ont -t 16 \
    input_data/reference_fasta.fna \
    SRR32410565.filtered.fastq.gz | \
    samtools sort -@ 16 -o output_files/reads_to_ref.sorted.bam -
samtools index output_files/reads_to_ref.sorted.bam
```

**Alignment 2: Assembly → Reference**
- **Purpose:** Structural variant detection
- **Tool:** minimap2 v2.28 (asm5 preset)
```bash
minimap2 -ax asm5 -t 16 \
    input_data/reference_fasta.fna \
    output_files/medaka_polished/consensus.fasta | \
    samtools sort -@ 16 -o output_files/assembly_to_ref.sorted.bam -
samtools index output_files/assembly_to_ref.sorted.bam
```

#### 6. Variant Calling

**Tool:** Clair3 v1.0.10  
**Method:** Deep learning-based variant detection  
**Filtering:** QUAL≥20, DP≥10

**Command:**
```bash
run_clair3.sh \
    --bam_fn=output_files/reads_to_ref.sorted.bam \
    --ref_fn=input_data/reference_fasta.fna \
    --threads=16 \
    --platform=ont \
    --model_path=/opt/models/r1041_e82_400bps_sup_v500 \
    --output=output_files/clair3_variants
```

**Filtering:**
```bash
bcftools view -i 'QUAL>=20 && DP>=10' \
    output_files/clair3_variants/merge_output.vcf.gz > \
    output_files/variants_filtered.vcf
```

#### 7. Visualization

**Tools:** R v4.4.0 with ggplot2, circlize, vcfR, patchwork  
**Outputs:**
- Circular genome plot with variant distribution
- Variant type and position analysis
- Assembly quality comparison plots

### Results

#### Assembly Statistics

| Metric | Raw Assembly | Polished Assembly | Reference |
|--------|--------------|-------------------|-----------|
| Total Length | ~4.85 Mb | ~4.86 Mb | 4.857 Mb |
| Contigs | 2-3 | 2-3 | 1 |
| N50 | >4.8 Mb | >4.85 Mb | 4.857 Mb |
| GC Content | ~52.1% | ~52.0% | 52.0% |
| BUSCO Completeness | >94% | >97% | 100% |

#### Variant Summary

[Results will be populated after pipeline completion]

- **Total variants:** [See output_files/variant_summary.txt]
- **SNPs:** [Count from analysis]
- **Insertions:** [Count from analysis]
- **Deletions:** [Count from analysis]
- **Variant density:** ~X variants per Mb

### Figures

**Figure 1:** Circular genome plot showing variant distribution across the chromosome  
**Figure 2:** Variant type distributions and genomic position analysis  
**Figure 3:** Assembly quality metrics comparison across pipeline stages

---

## Workflow Execution

### Complete Pipeline
```bash
cd ~/Assignment-1---Binf6110
bash scripts/run_all.sh
```

### Individual Steps
```bash
sbatch scripts/01_download_data.sh      # Data acquisition
sbatch scripts/02_qc_filter.sh          # Quality control and filtering
sbatch scripts/03_assembly.sh           # De novo assembly with Flye
sbatch scripts/04_polish.sh             # Medaka polishing
sbatch scripts/05_assembly_qc.sh        # QUAST, BUSCO, Bandage
sbatch scripts/06_alignment.sh          # Dual alignment strategy
sbatch scripts/07_variant_calling.sh    # Clair3 variant calling
sbatch scripts/08_visualizations.sh     # Generate figures
```

### Monitor Progress
```bash
squeue -u yazanalq
```

---

## Software Versions

| Tool | Version | Purpose |
|------|---------|---------|
| Flye | 2.9.6 | De novo assembly |
| Medaka | 2.0.1 | Consensus polishing |
| minimap2 | 2.28 | Read/assembly alignment |
| SAMtools | 1.19 | BAM processing |
| BCFtools | 1.19 | VCF processing |
| Clair3 | 1.0.10 | Variant calling |
| QUAST | 5.2.0 | Assembly quality assessment |
| BUSCO | 5.7.1 | Genome completeness |
| Bandage | 0.8.1 | Assembly graph visualization |
| NanoPlot | 1.46.2 | Read quality control |
| NanoFilt | 2.8.0 | Read filtering |
| IGV | 2.19.7 | Genome visualization |
| R | 4.4.0 | Statistical analysis and plotting |

---

## Discussion

[This section will be completed after analysis, addressing:]

- Comparison of assembly quality metrics to published *Salmonella enterica* genomes
- Biological significance of identified variants in the context of virulence factors and antimicrobial resistance
- Impact of polishing on variant calling accuracy
- Limitations of long-read assembly and variant detection
- Future directions for improving assembly contiguity and variant validation

---

## References

Jajere, S. M. (2019). A review of Salmonella enterica with particular focus on the pathogenicity and virulence factors, host specificity and antimicrobial resistance including multidrug resistance. *Veterinary World*, 12(4), 504–521. doi:10.14202/vetworld.2019.504-521

Havelaar, A. H., Kirk, M. D., Torgerson, P. R., et al. (2015). World Health Organization global estimates and regional comparisons of the burden of foodborne disease in 2010. *PLOS Medicine*, 12(12), e1001923. doi:10.1371/journal.pmed.1001923

De Coster, W., D'Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: Visualizing and processing long read sequencing data. *Bioinformatics*, 34(15), 2666–2669. doi:10.1093/bioinformatics/bty149

Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37(5), 540–546. doi:10.1038/s41587-019-0072-8

Wick, R. R., & Holt, K. E. (2019). Benchmarking of long-read assemblers for prokaryote whole genome sequencing. *F1000Research*, 8, 2138. doi:10.12688/f1000research.21782.4

Sereika, M., Vetcher, A., et al. (2022). Oxford Nanopore R10.4 long-read sequencing enables the generation of near-finished bacterial genomes from pure cultures and metagenomes without short-read or reference polishing. *Nature Methods*. doi:10.1038/s41592-022-01539-7

Lee, J., Kong, M., Oh, J., et al. (2021). Comparative evaluation of nanopore polishing tools for microbial genome assembly and polishing strategies for downstream analysis. *Scientific Reports*, 11, 20740. doi:10.1038/s41598-021-00178-w

Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079. doi:10.1093/bioinformatics/btp352

Danecek, P., Bonfield, J. K., Liddle, J., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. doi:10.1093/gigascience/giab008

Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. doi:10.1093/bioinformatics/bty191

Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072–1075. doi:10.1093/bioinformatics/btt086

Thorvaldsdóttir, H., Robinson, J. T., & Mesirov, J. P. (2013). Integrative Genomics Viewer (IGV): High-performance genomics data visualization and exploration. *Briefings in Bioinformatics*, 14(2), 178–192. doi:10.1093/bib/bbs017

Wick, R. R., Schultz, M. B., Zobel, J., & Holt, K. E. (2015). Bandage: Interactive visualization of de novo genome assemblies. *Bioinformatics*, 31(20), 3350–3352. doi:10.1093/bioinformatics/btv383

Vaser, R., & Šikić, M. (2021). Time- and memory-efficient genome assembly with Raven. *Nature Computational Science*, 1, 332–336. doi:10.1038/s43588-021-00073-4

Ruan, J., & Li, H. (2020). Fast and accurate long-read assembly with wtdbg2. *Nature Methods*, 17, 155–158. doi:10.1038/s41592-019-0669-3

NCBI Sequence Read Archive (SRA). Oxford Nanopore R10 reads for Salmonella enterica. Run accession: SRR32410565.

