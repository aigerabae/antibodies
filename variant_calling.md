I will try variant calling with bos taurus referent genome and cruk-ci ampliseq pipeline:
Didn't run anything yet
Alignment:
```bash
while read -r line; do
    echo "Processing file: $line"
    bwa mem -M -t 20 \
        ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        fastp/${line}R1_trimmed.fastq.gz \
        fastp/${line}R2_trimmed.fastq.gz \
    | samtools view -b -o bams/${line}.bam
done < pairs.txt
```

sorting
```bash
conda activate samtools
for i in *.bam; do samtools index "$i"; done
```

Config file - needs fixing
```
params {
    samples               = "input/sample_sheet.tsv"
    amplicons             = "input/amplicons.tsv"
    referenceGenomeFasta  = "../ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    vepAnnotation         = true
    vepCacheDir           = "vep_cache/"
    vepSpecies            = "homo_sapiens"
    vepAssembly           = "GRCh38"
    outputDir             = "results"
    variantCaller         = "HaplotypeCaller"
    minimumAlleleFraction = 0.01
}

profiles {
    myprofile {
        process.executor = 'local'
        executor {
            cpus = 20
            memory = 32.GB
        }
        docker.enabled = true
    }
}
```
