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

## Output

See [this wiki page](https://github.com/brentp/rare-disease-wf/wiki/Workflow-Output) for more information about how to use the output.

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
 1. Generates and links pre-made, standalone [igv.js](https://github.com/igvteam/igv.js)/[jigv](https://github.com/brentp/jigv) outputs for each candidate.

And the key output will be in: `results-rare-disease/${cohort_name}.slivar.candidates.tsv` which is something one can easily view in excel or other spreadsheet software.
In addition, it will create: `results-rare-disease/${cohort_name}.jigv.html` and `results-rare-disease/jigv_plots/*` which together provide an HTML table and interactive [igv.js](https://github.com/brentp/igvteam/igv.js) views of each variant and associated alignments that do not rely on the original alignment files.

In coming releases, this will:

 1. Output QC with [somalier](https://github.com/brentp/somalier) and other tools to be shown in [multiQC](https://multiQC.info)
 1. Output high-quality SVs (using manta-> graphtyper)

## Octopus

currently, [octopus](https://doi.org/10.1038/s41587-021-00861-3) is included as
a separate workflow. This octopus.nf pipeline will detect trios and families
and run them together and then iteratively merge across families using the
`n+1` schema [described in the octopus
docs](https://luntergroup.github.io/octopus/docs/guides/models/population)
Finally, the workflow will do the forest filtering as recommended by the
octopus documentation.
We plan to integrate the octopus and deepvariant calls in the future.


## Future Development

Development and research is underway so that it will:

 1. Add a high-quality set of SV/CNVs
    - [Manta](https://github.com/Illumina/manta) + SVchannels and [duphold](https://github.com/brentp/duphold) filtering
 1. Add some **prioritization** of variants
    - For example, lower priority to variants filtered in gnomAD
 1. Integrate SV/CNV calls with the snp/indels to find, for example compound heterozygotes with a snp:SV pair.
 1. Evaluate use of [octopus](https://github.com/luntergroup/octopus) to find large indels (and/or SNPs and indels).
 1. Use GTex + phenotypes to further prioritize variants in a family and phenotype-specific way, such that, for example
    variants in genes that are not expressed in relevant tissues are down-weighted.
 1. Provide a graphical-user-interface so that sorting, filtering, note-taking, sharing is simplified


## Software Used

+ [DeepVariant](https://github.com/google/deepvariant) Variant Calling with Deep Learning. https://doi.org/10.1038/nbt.4235
+ [GLNexus](https://github.com/dnanexus-rnd/GLnexus) Joint variant calling. http://dx.doi.org/10.1101/343970
+ [octopus](https://github.com/luntergroup/octopus) haplotype-based mutation caller. https://doi.org/10.1038/s41587-021-00861-3
+ [bcftools](https://github.com/samtools/bcftools) BCF/VCF manipulation. https://doi.org/10.1093/gigascience/giab008
+ [bcftools csq](https://github.com/samtools/bcftools) variant consequence annotation. https://doi.org/10.1093/bioinformatics/btx100
+ [htslib](https://github.com/samtools/htslib) C libary for genomics data. https://doi.org/10.1093/gigascience/giab007
+ [slivar](https://github.com/brentp/slivar) variant filtering and annotation. https://doi.org/10.1101/2020.08.13.249532
+ [igv.js](https://github.com/igvteam/igv.js/). javascript genomics viewer. https://doi.org/10.1101/2020.05.03.075499
+ [nextflow](https://nextflow.io/) scientific workflows. https://doi.org/10.1038/nbt.3820
+ [manta](https://github.com/Illumina/manta) structural variant caller. https://doi.org/10.1093/bioinformatics/btv710
+ [dysgu](https://github.com/kcleal/dysgu) structural variant caller. https://doi.org/10.1101/2021.05.28.446147 
+ [paragraph](https://github.com/Illumina/paragraph) structural variant genotyper. https://doi.org/10.1186/s13059-019-1909-7
+ [jasmine](https://github.com/mkirsche/Jasmine) structural variant merging. https://doi.org/10.1101/2021.05.27.445886
+ [duphold](https://github.com/brentp/duphold) structural variant depth annotation. https://doi.org/10.1093/gigascience/giz040
+ [snpEff](http://pcingola.github.io/SnpEff/) variant consequence annotation. https://doi.org/10.4161/fly.19695
+ [svpack](https://github.com/amwenger/svpack/) structural variant annotation.
