For rare-disease, the [best practices and expected number of candidate variants for each inheritance mode are known](https://www.biorxiv.org/content/10.1101/2020.08.13.249532v3). The
actual filtering is easily done with a tool like [slivar](https://github.com/brentp/slivar/wiki/rare-disease). 
This is a necessary first step with the following limitations:

 1. it leaves an analyst or clinician with choices on how to prioritize the 10-15 candidates variants or ~100 for autosomal (non *de novo*) dominant.
    - This is quite a small number, but the prioritization after this is highly variable across tools and analysts.
 2. it is limited text/spreadsheet output
 3. it assumes a high-quality, jointly-called VCF is already available
 4. it leaves the analyst with the chore of getting IGV set up, and browsing each candidate for each family.


## Quickstart 

Note, it is early days for the project. It will produce high-quality SNP/indel candidates
but you may need experience with [nextflow](https://nextflow.io) to run it easily.

This project **currently** has workflow that can be run as:

```
# NOTE that you need to remove everything after \ on each line for the command to work
# the comments here are just for documentation purposes.
nextflow run -resume -profile slurm rare-disease.nf \
    -config nextflow.config \    # a starting config is included in this repo. adjust from there.
    --xams "/path/to/*/*.cram" \ # NOTE that this is a string glob
    --ped $pedigree_file \       # see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format
    --fasta $reference_fasta \
    --gff $gff \                   # e.g. from: ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
    --slivarzip gnomad.hg38.zip  \  # from: https://github.com/brentp/slivar#gnotation-files
    --cohort_name my_rare_disease
```

This does:

 1. Run [DeepVariant](https://github.com/google/deepvariant) and [GLNexus](https://github.com/dnanexus-rnd/GLnexus) (we have shown these tools to give higher quality results for trios) in an efficient nextflow workflow that can be easily run in the cloud or on a cluster.
 1. Decompose and normalize variants.
 1. Annotate with [bcftools csq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870570/) and [snpEff](https://pcingola.github.io/SnpEff/adds/SnpEff_paper.pdf) 
 1. Annotate with allele frequency and inheritance modes using [slivar](https://github.com/brentp/slivar)
 1. Annotate with gene-based annotations:
    - clinvar-gene-phenotype
    - loss-of-function intolerance
 1. Output high-quality calls from [slivar](https://github.com/brentp/slivar) for recessive, dominant, x-linked, compound-het and other
    inheritance modes.

And the key output will be in: `results-rare-disease/slivar.candidates.tsv` which is something one can easily view in excel or other spreadsheet
software.

In coming releases, this will:

 1. Provide pre-made IGV outputs for each candidate.
 1. Output QC with [somalier](https://github.com/brentp/somalier) and other tools to be shown in [multiQC](https://multiQC.info)
 1. Output high-quality SVs (using manta-> graphtyper)


## Future Development

Development and research is underway so that it will:

 1. Auto-generate per-family IGV images or interactive igv.js pages of candidate variants.
 1. Automatically detect trios in the given pedigree file (or via [somalier](https://github.com/brentp/somalier)) and run [DeepTrio](https://www.biorxiv.org/content/10.1101/2021.04.05.438434v1.full)
 1. Add a high-quality set of SV/CNVs
    - [Manta](https://github.com/Illumina/manta) + SVchannels and [duphold](https://github.com/brentp/duphold) filtering
 1. Add some **prioritization** of variants
    - For example, lower priority to variants filtered in gnomAD
 1. Integrate SV/CNV calls with the snp/indels to find, for example compound heterozygotes with a snp:SV pair.
 1. Evaluate use of [octopus](https://github.com/luntergroup/octopus) to find large indels (and/or SNPs and indels).
 1. Use GTex + phenotypes to further prioritize variants in a family and phenotype-specific way, such that, for example
    variants in genes that are not expressed in relevant tissues are down-weighted.
 1. Provide a graphical-user-interface so that sorting, filtering, note-taking, sharing is simplified
