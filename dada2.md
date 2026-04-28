The primer sequences need to be without adapter/overhang

Removing TCGTCGGCAGCGTCAGATGTGTATAAGAGACA and GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

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

FWD1 <- "GGCTGCAATGCCATCTGCTCT"
FWD2 <- "GGCCTGCAGGACAGAAATGACTT"
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

This table is a Primer Orientation Matrix:  
                  Forward Complement Reverse RevComp  
FWD1.ForwardReads      17          0       0       0  
FWD1.ReverseReads       1          0       0     108  
REV1.ForwardReads     588          0       0   28012  
REV1.ReverseReads   82311          0       0      43  

                  Forward Complement Reverse RevComp  
FWD2.ForwardReads       9          0       0       8  
FWD2.ReverseReads       1          0       0    1935  
REV2.ForwardReads     307          0       0      16  
REV2.ReverseReads   44223          0       0       0  

Key Takeaways:
- The "Reverse" Bias: For some reason, your library preparation or the way the MiSeq was set up heavily favored the Reverse Primer side. This is why you couldn't find your FWD primers earlier with grep—they simply aren't in the reads in significant numbers.
- Read-Through: Since you see REV1 in both R1 and R2 (one as Forward, one as RevComp), it means your DNA fragment might be shorter than the 301 cycles you programmed. The sequencer read all the way through the fragment and started reading the primer on the other side.
- Primer Dimer vs. Real Data: If these fragments are very short, you might be sequencing "Primer Dimers" (primers just sticking to each other without any actual gene in between).

What to check next:
- Look at the length of your reads in the FASTQ file. If your fragments were supposed to be 500bp but your reads are only 300bp, you will only see one primer. If the fragments are 150bp, you will see both primers, but they will be "flipped" in one of the reads.
- Does the expected size of your DNA fragments (the "legit" part) match the 301bp length of the MiSeq run?
