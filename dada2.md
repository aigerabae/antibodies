The primer sequences need to be without adapter/overhang

Original primers:
Primers: 
1 fragment 
Inf_alpha-73F-3 TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGCTGCAATGCCATCTGCTCT 
Inf_alpha_rev415 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCTCCTCCTGCCTCAGACAG  

2 fragment   
Inf_alpha_f186-2 TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGCCTGCAGGACAGAAATGACTT  
Inf_alpha_R601-3 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCTTTCTCCTGAAACTCTCCTGC  


fwd1
zcat Gradient-C7_S272_L001_R1_001.fastq.gz | grep "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGCTGCAATGCCATCTGCTCT" | wc -l   # 0
zcat Gradient-C7_S272_L001_R1_001.fastq.gz | grep                                 "GGCTGCAATGCCATCTGCTCT" | wc -l   #17
zcat Gradient-C7_S272_L001_R1_001.fastq.gz | grep                                  "GCTGCAATGCCATCTGCTCT" | wc -l   #79395
zcat Gradient-C7_S272_L001_R1_001.fastq.gz | grep                                   "CTGCAATGCCATCTGCTCT" | wc -l   #79951

rev1
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCTCCTCCTGCCTCAGACAG" | wc -l   # 0
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep                                  "GCCTCCTCCTGCCTCAGACAG" | wc -l   #41
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep                                   "CCTCCTCCTGCCTCAGACAG" | wc -l   #82721
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep                                    "CTCCTCCTGCCTCAGACAG" | wc -l   #83574

fwd2
zcat Gradient-C7_S272_L001_R1_001.fastq.gz | grep "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGCCTGCAGGACAGAAATGACTT" | wc -l   # 0
zcat Gradient-C7_S272_L001_R1_001.fastq.gz | grep                                 "GGCCTGCAGGACAGAAATGACTT" | wc -l   #9
zcat Gradient-C7_S272_L001_R1_001.fastq.gz | grep                                  "GCCTGCAGGACAGAAATGACTT" | wc -l   #121873
zcat Gradient-C7_S272_L001_R1_001.fastq.gz | grep                                   "CCTGCAGGACAGAAATGACTT" | wc -l   #122397

rev2
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCTTTCTCCTGAAACTCTCCTGC" | wc -l   # 0
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep                                  "GCCTTTCTCCTGAAACTCTCCTGC" | wc -l   # 3
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep                                   "CCTTTCTCCTGAAACTCTCCTGC" | wc -l   # 44233
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep                                    "CTTTCTCCTGAAACTCTCCTGC" | wc -l   # 44761
zcat Gradient-C7_S272_L001_R2_001.fastq.gz | grep                                           "CTGAAACTCTCCTGC" | wc -l   # 46512

To remove:
fwd1:  TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
check: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
fwd2:  TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG

rev1:  GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
check: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
rev2:  GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

Removing TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG and GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

following tutorials: https://benjjneb.github.io/dada2/tutorial.html
this seems newer: https://benjjneb.github.io/dada2/ITS_workflow.html
might want to do it separately for 2 fragments and then combine

Pulling
```bash
docker pull blekhmanlab/dada2:1.26.0a
```

Running
```bash
docker run -it -v ~/biostar/NCB/antibodies/dada2:/home blekhmanlab/dada2:1.26.0a bash
```

Going to directory and opening R
```bash
cd home
R
```

Code in R
```R
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

path <- "/home"  ## CHANGE ME to the directory containing the fastq files.
fnFs <- sort(list.files(path, pattern = "1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "2_001.fastq.gz", full.names = TRUE))

FWD1 <- "GCTGCAATGCCATCTGCTCT"
FWD2 <- "GCCTGCAGGACAGAAATGACTT"
REV1 <- "CCTCCTCCTGCCTCAGACAG"
REV2 <- "CCTTTCTCCTGAAACTCTCCTGC"

# checking orientation of primers
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD1.orients <- allOrients(FWD1)
REV1.orients <- allOrients(REV1)
FWD2.orients <- allOrients(FWD2)
REV2.orients <- allOrients(REV2)
FWD1.orients

# removing N nucleotides from sequences
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# We are now ready to count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations. Identifying and counting the primers on one set of paired end FASTQ files is sufficient, assuming all the files were created using the same library preparation, so we’ll just process the first sample.
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.filtN[[1]]), FWD1.ReverseReads = sapply(FWD1.orients,
    primerHits, fn = fnRs.filtN[[1]]), REV1.ForwardReads = sapply(REV1.orients, primerHits,
    fn = fnFs.filtN[[1]]), REV1.ReverseReads = sapply(REV1.orients, primerHits, fn = fnRs.filtN[[1]]))

rbind(FWD2.ForwardReads = sapply(FWD2.orients, primerHits, fn = fnFs.filtN[[1]]), FWD2.ReverseReads = sapply(FWD2.orients,
    primerHits, fn = fnRs.filtN[[1]]), REV2.ForwardReads = sapply(REV2.orients, primerHits,
    fn = fnFs.filtN[[1]]), REV2.ReverseReads = sapply(REV2.orients, primerHits, fn = fnRs.filtN[[1]]))
```

This table is a Primers:  
                  Forward Complement Reverse RevComp
FWD1.ForwardReads   77538          0       0       0
FWD1.ReverseReads     389          0       0     423
REV1.ForwardReads     588          0       0   28012
REV1.ReverseReads   82311          0       0      43
                  Forward Complement Reverse RevComp
FWD2.ForwardReads  118973          0       0     305
FWD2.ReverseReads     525          0       0   21436
REV2.ForwardReads     307          0       0      16
REV2.ReverseReads   44223          0       0       0


Installing cutadapt via conda:
```bash
#apt install curl
#curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# i copied file from another repo 
bash Miniconda3-latest-Linux-x86_64.sh
conda
```

in another window:
```bash
docker commit 20277fe72a02 dada2:v2.0
docker run -it -v ~/biostar/NCB/antibodies/dada2:/home dada2:v2.0 bash
cd home/
conda create --name cutadapt python=3.7
conda activate cutadapt
conda install bioconda::cutadapt
which cutadapt
```

in R:
```R
cutadapt <- "/home/conda/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R
```

Cutting - separate for 2 fragments
```R
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# primer 1
FWD1.RC <- dada2:::rc(FWD1)
REV1.RC <- dada2:::rc(REV1)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD1, "-a", REV1.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV1, "-A", FWD1.RC) 
# Run Cutadapt 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                              "--discard-untrimmed",         # i do it to separate 2 amplicons
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
```

I manually moved these to primer1 folder inside cutadapt

```bash
# primer 2
FWD2.RC <- dada2:::rc(FWD2)
REV2.RC <- dada2:::rc(REV2)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD2, "-a", REV2.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV2, "-A", FWD2.RC) 
# Run Cutadapt 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                              "--discard-untrimmed",         # i do it to separate 2 amplicons
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
```

Bases preceding removed adapters:
  A: 87.7%
  C: 7.4%
  G: 4.0%
  T: 0.8%
  none/other: 0.0%
WARNING:
    The adapter is preceded by "A" extremely often.
    The provided adapter sequence could be incomplete at its 3' end.

The primers are definitely correct so must be the biological reason for it

Sanity check primer1 files:
```bash
path.cut <- file.path(path, "cutadapt/primer1")
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
# primer 1 in primer1 files
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.cut[[1]]), FWD1.ReverseReads = sapply(FWD1.orients,
    primerHits, fn = fnRs.cut[[1]]), REV1.ForwardReads = sapply(REV1.orients, primerHits,
    fn = fnFs.cut[[1]]), REV1.ReverseReads = sapply(REV1.orients, primerHits, fn = fnRs.cut[[1]]))
# primer 2 in primer1 files
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.cut[[1]]), FWD1.ReverseReads = sapply(FWD1.orients,
    primerHits, fn = fnRs.cut[[1]]), REV1.ForwardReads = sapply(REV1.orients, primerHits,
    fn = fnFs.cut[[1]]), REV1.ReverseReads = sapply(REV1.orients, primerHits, fn = fnRs.cut[[1]]))
```


Sanity check primer1 files:
```bash
path.cut <- file.path(path, "cutadapt/primer2")
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
# primer 2 in primer1 files
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.cut[[1]]), FWD1.ReverseReads = sapply(FWD1.orients,
    primerHits, fn = fnRs.cut[[1]]), REV1.ForwardReads = sapply(REV1.orients, primerHits,
    fn = fnFs.cut[[1]]), REV1.ReverseReads = sapply(REV1.orients, primerHits, fn = fnRs.cut[[1]]))
# primer 2 in primer2 files
rbind(FWD1.ForwardReads = sapply(FWD1.orients, primerHits, fn = fnFs.cut[[1]]), FWD1.ReverseReads = sapply(FWD1.orients,
    primerHits, fn = fnRs.cut[[1]]), REV1.ForwardReads = sapply(REV1.orients, primerHits,
    fn = fnFs.cut[[1]]), REV1.ReverseReads = sapply(REV1.orients, primerHits, fn = fnRs.cut[[1]]))
```
