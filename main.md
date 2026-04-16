```bash
conda activate fastqc
fastqc *.fastq.gz 
```

following https://github.com/deepalivasoya/MHCtyping from https://link.springer.com/article/10.1007/s00251-016-0945-7#Sec2
In folder
```bash
conda create --name antibodies
source activate antibodies
conda install -c conda-forge perl
conda install -c bioconda blast
conda install -c bioconda flash
conda install -c bioconda sickle-trim
conda install -c bioconda fastqc
conda install -c bioconda snakemake
```

In folder custom
```bash
conda install -c bioconda fastqc trimmomatic flash
trimmomatic PE -phred33 \
  Gradient-C7_S272_L001_R1_001.fastq.gz Gradient-C7_S272_L001_R2_001.fastq.gz \
  r1_paired.fq.gz r1_unpaired.fq.gz \
  r2_paired.fq.gz r2_unpaired.fq.gz \
  SLIDINGWINDOW:4:28 MINLEN:50
flash r1_paired.fq.gz r2_paired.fq.gz -m 10 -M 150 -o bovine_merged
```

Spades after trimming:
```bash
conda install bioconda::spades
mkdir spades
spades.py -1 r1_paired.fq.gz -2 r2_paired.fq.gz -o spades
```

======= SPAdes pipeline finished WITH WARNINGS!

=== Error correction and assembling warnings:
 * 0:00:01.339     2M / 256M  WARN    General                 (kmer_coverage_model.cpp   : 219)   Too many erroneous kmers, the estimates might be unreliable
 * 0:00:10.281     3M / 323M  WARN    General                 (kmer_coverage_model.cpp   : 328)   Valley value was estimated improperly, reset to 1
 * 0:00:10.282     3M / 323M  WARN    General                 (kmer_coverage_model.cpp   : 367)   Failed to determine erroneous kmer threshold. Threshold set to: 1
 * 0:00:03.899     4M / 255M  WARN    General                 (kmer_coverage_model.cpp   : 367)   Failed to determine erroneous kmer threshold. Threshold set to: 17
======= Warnings saved to /mnt/harddisk/biostar/NCB/antibodies/custom/spades/warnings.log

SPAdes log can be found here: /mnt/harddisk/biostar/NCB/antibodies/custom/spades/spades.log

Thank you for using SPAdes! If you use it in your research, please cite:

  Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A. and Korobeynikov, A., 2020. Using SPAdes de novo assembler. Current protocols in bioinformatics, 70(1), p.e102.
  doi.org/10.1002/cpbi.102
