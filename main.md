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
