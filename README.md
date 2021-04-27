For rare-disease, the [best practices and expected number of candidate variants for each inheritance mode are known](https://www.biorxiv.org/content/10.1101/2020.08.13.249532v3). The
actual filtering is easily done with a tool like [slivar](https://github.com/brentp/slivar/wiki/rare-disease). 
This is a necessary first step with the following limitations:

 1. it leaves an analyst or clinician with choices on how to prioritize the 10-15 candidates variants or ~100 for autosomal (non *de novo*) dominant.
    - This is quite a small number, but the prioritization after this is highly variable across tools and analysts.
 2. it is limited text/spreadsheet output
 3. it assumes a high-quality, jointly-called VCF is already available

The aim of this project is a more complete workflow that can be run something like:

```
nextflow run rare-disease.nf --ped $ped \
	--cohort cakut-cohort --fasta $fasta \
	--alignment-list cram-paths.txt
```

This will:

 1. Run [DeepVariant](https://github.com/google/deepvariant) and [GLNexus](https://github.com/dnanexus-rnd/GLnexus) (we have shown these tools to give higher quality results for trios) in an efficient nextflow workflow that can be easily run in the cloud or on a cluster.
 1. Decompose and normalize variants.
 1. Annotate with [bcftools csq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870570/) and [snpEff](https://pcingola.github.io/SnpEff/adds/SnpEff_paper.pdf) 
 1. Annotate with allele frequency and inheritance modes using [slivar](https://github.com/brentp/slivar)
 1. Annotate with gene-based annotations:
    - clinvar-gene-phenotype
    - loss-of-function intolerance
 1. Output QC with [somalier](https://github.com/brentp/somalier) and other tools to be shown in [multiQC](https://multiQC.info)


Development and research is underway so that it will:

 1. Auto-generate per-family IGV images or interactive igv.js pages of candidate variants.
 1. Add a high-quality set of SV/CNVs
    - [Manta](https://github.com/Illumina/manta) + SVchannels and [duphold](https://github.com/brentp/duphold) filtering
 1. Add some **prioritization** of variants
    - For example, lower priority to variants filtered in gnomAD
 1. Integrate SV/CNV calls with the snp/indels to find, for example compound heterozygotes with a snp:SV pair.
 1. Evaluate use of [octopus](https://github.com/luntergroup/octopus) to find large indels (and/or SNPs and indels).
 1. Use GTex + phenotypes to further prioritize variants in a family and phenotype-specific way, such that, for example
    variants in genes that are not expressed in relevant tissues are down-weighted.
 1. Provide a graphical-user-interface so that sorting, filtering, note-taking, sharing is simplified
    - this might be adding custom scripts and colors (:barf:) to a google-sheets document.
