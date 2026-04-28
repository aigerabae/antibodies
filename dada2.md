The primer sequences need to be without adapter/overhang

Removing TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG and GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

following tutorials: https://benjjneb.github.io/dada2/tutorial.html
this seems newer: https://benjjneb.github.io/dada2/ITS_workflow.html
might want to do it separately for 2 fragments and then combine

```bash
docker pull blekhmanlab/dada2:1.26.0a
docker run -it -v ~/biostar/NCB/antibodies/dada2:/home blekhmanlab/dada2:1.26.0a bash
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

Cutting - didn't do yet bc im not sure about both fragments:
```R
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
```
