```bash
conda activate fastqc
fastqc *.fastq.gz 
```

following https://github.com/deepalivasoya/MHCtyping from https://link.springer.com/article/10.1007/s00251-016-0945-7#Sec2
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

```bash
conda install -c bioconda fastqc trimmomatic flash
trimmomatic PE -phred33 \
  Gradient-C7_S272_L001_R1_001.fastq.gz Gradient-C7_S272_L001_R2_001.fastq.gz \
  r1_paired.fq.gz r1_unpaired.fq.gz \
  r2_paired.fq.gz r2_unpaired.fq.gz \
  SLIDINGWINDOW:4:28 MINLEN:50
flash r1_paired.fq.gz r2_paired.fq.gz -m 10 -M 150 -o bovine_merged
```
