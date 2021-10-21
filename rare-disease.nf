nextflow.enable.dsl=2

include  { find_index } from './nf/common'

process DeepVariant {
    label "DeepVariant"
    container = 'docker://gcr.io/deepvariant-docker/deepvariant:1.1.0'
    publishDir "${params.output_dir}/gvcfs/", mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(sample_id), path(aln_file), path(aln_index))
        path(fasta)
        path(fai)

    output:
      //tuple(file("${sample_id}.gvcf.gz"), file("${sample_id}.gvcf.gz.csi"))
      tuple(file("${sample_id}.gvcf.gz"), file("${sample_id}.gvcf.gz.tbi"))


    script:
    is_cram = aln_file.toString().endsWith(".cram")
    """
# faster to convert cram to bam, but easier to use original DV image.
#if [[ "$is_cram" == "true" ]]; then
#    samtools view --write-index --threads ${task.cpus} -bT $fasta -o ${sample_id}.bam $aln_file
#else
#    ln -sf $aln_file ${sample_id}.bam
#    ln -sf $aln_index ${sample_id}.bam.bai
#fi
#--reads=${sample_id}.bam \

echo "TMPDIR:\$TMPDIR"

/opt/deepvariant/bin/run_deepvariant \
    --reads=${aln_file} \
    --intermediate_results_dir=\$TMPDIR \
    --model_type=${params.model_type} \
    --output_gvcf=${sample_id}.gvcf.gz \
    --output_vcf=${sample_id}.vcf.gz \
    --num_shards=${task.cpus} \
    --ref=$fasta
 
#rm -f ${sample_id}.bam
#bcftools index --threads 6 ${sample_id}.gvcf.gz
    """
}

include { split; split_by_size } from "./split"

process glnexus_anno_slivar {
    container = 'docker://brentp/rare-disease:v0.2.2'
    shell = ['/bin/bash', '-euo', 'pipefail']
    publishDir "${params.output_dir}/joint-by-chrom/", mode: 'copy'

    input:
        tuple(path(gvcfs), val(chrom))
        path(fasta)
        path(fai)
        path(gff)
        path(slivar_zip)
        val(cohort_name)

    output: tuple(path(output_file), path(output_csi))

    script:
        output_file = "${cohort_name}.${chrom}.glnexus.anno.bcf"
        output_csi = "${output_file}.csi"
        file("$workDir/file.list.${cohort_name}.${chrom}").withWriter { fh ->
            gvcfs.each { gvcf ->
                fh.write(gvcf.toString()); fh.write("\n")
            }
        }
        """
# GRCh38.99
# GRCh37.75

glnexus_cli \
    -t ${task.cpus} \
    --mem-gbytes 128 \
    --config DeepVariant${params.model_type} \
    --list $workDir/file.list.${cohort_name}.${chrom} \
| bcftools norm --threads 3 -m - -w 10000 -f $fasta -O u \
| bcftools csq --threads 3 -s - --ncsq 50 -g $gff -l -f $fasta - -o - -O v \
| snpEff eff -noStats -dataDir $projectDir GRCh38.99 \
| slivar expr -g $slivar_zip -o $output_file --vcf -


bcftools index --threads 6 $output_file
    """

}

process slivar_rare_disease {
  container = 'docker://brentp/rare-disease:v0.2.2'
  publishDir "${params.output_dir}/joint-by-chrom-slivar/", mode: 'copy'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input: tuple(path(bcf), path(csi))
         path(ped)
         path(gnomad_zip)

  output: tuple(path(slivar_bcf), path(slivar_ch_bcf), path(slivar_bcf_csi), path(slivar_ch_bcf_csi), path(slivar_tsv), path(slivar_counts))

  script:
  slivar_bcf = bcf.getBaseName() + ".slivar.bcf"
  slivar_ch_bcf = bcf.getBaseName() + ".slivar.ch.bcf"
  slivar_bcf_csi = slivar_bcf + ".csi"
  slivar_ch_bcf_csi = slivar_ch_bcf + ".csi"
  slivar_tsv = bcf.getBaseName() + ".slivar.tsv"
  slivar_counts = bcf.getBaseName() + ".slivar.counts.txt"

  """
# NOTE: we do *NOT* limit to impactful so that must be reported and used by slivar tsv

export SLIVAR_SUMMARY_FILE=${slivar_counts}.main
slivar expr --vcf $bcf \
    --ped $ped \
    -o $slivar_bcf \
    --pass-only \
    -g $gnomad_zip \
    --info 'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"' \
    --js /opt/slivar/slivar-functions.js \
    --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \
    --family-expr 'recessive:fam.every(segregating_recessive)' \
    --family-expr 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001' \
    --family-expr 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)' \
    --family-expr 'dominant:INFO.gnomad_popmax_af < 0.005 && fam.every(segregating_dominant)' \
    --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10'

bcftools index --threads 3 $slivar_bcf &

export SLIVAR_SUMMARY_FILE=${slivar_counts}.ch
slivar compound-hets -v $slivar_bcf \
    --sample-field comphet_side --sample-field denovo -p $ped \
  | bcftools view -O b -o $slivar_ch_bcf

bcftools index --threads 3 $slivar_ch_bcf &

tiwih combine_slivar_counts ${slivar_counts}.main ${slivar_counts}.ch > $slivar_counts

slivar tsv \
  -s denovo \
  -s x_denovo \
  -s recessive \
  -s x_recessive \
  -s dominant \
  -i gnomad_popmax_af -i gnomad_popmax_af_filter -i gnomad_nhomalt \
  -i impactful -i genic \
  -c ANN \
  -g /opt/slivar/pli.lookup \
  -g /opt/slivar/clinvar_gene_desc.txt \
  -p $ped \
  $slivar_bcf > $slivar_tsv

slivar tsv \
  -s slivar_comphet \
  -i gnomad_popmax_af -i gnomad_popmax_af_filter -i gnomad_nhomalt \
  -i impactful -i genic \
  -c ANN \
  -g /opt/slivar/pli.lookup \
  -g /opt/slivar/clinvar_gene_desc.txt \
  -p $ped \
  $slivar_ch_bcf \
  | { grep -v ^# || true; } >> $slivar_tsv # || true avoids error if there are no compound hets.


wait
  """
}

process slivar_merge_tsvs {
  // TSVS are by chromosome. just concatentate them here to get a single file.
  container = 'docker://brentp/rare-disease:v0.2.2'
  publishDir "${params.output_dir}/", mode: 'copy'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input: path(tsvs)
         val(cohort_name)

  output: tuple(path("${output_file}"), path("${html_output}"))

  script:
    output_file = "${cohort_name}.slivar.candidates.tsv"
    html_output = "${cohort_name}.jigv.html"
    """
# get header from first file and drop it from other files
# and make sure slivar_comphet id is unique
    awk 'NR == FNR || FNR > 1 { sub(/^slivar_comphet/, "slivar_comphet_"NR, \$0); print; }' $tsvs > ${output_file}

    tiwih slivar_jigv_tsv \
        --html-template /opt/rare-disease/tmpl.html \
        --prefix 'jigv_plots/\${family_id}/\${mode}/\${family_id}.\${mode}.\${site}.js' \
        ${output_file} > ${html_output}

    """
}

process slivar_sum_counts {
  container = 'docker://brentp/rare-disease:v0.2.2'
  publishDir "${params.output_dir}/", mode: 'copy'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input: path(counts)

  output: path("${output_file}")

  script:
    output_file = "slivar.counts.txt"
    """
    tiwih sum_slivar_counts -z $counts > ${output_file}
    """
}

process slivar_split_by_fam {
  container = 'docker://brentp/rare-disease:v0.2.2'
  publishDir "${params.output_dir}/slivar_split_by_fam_mode", mode: 'copy'
  shell = ['/bin/bash', '-euo', 'pipefail']
  input: 
    path(ped)
    path(bcfs)

  output: 
    path("*.fam.*.bcf{,.csi}")

  script:
    """
# NOTE: if --template is changed, must change the mode= part of generate_jigv_pages

tiwih slivar_split_fam --ped $ped \
  --fields denovo,recessive,x_denovo,x_recessive,dominant,slivar_comphet \
  --template 'slivar.rd-byfam.\${field}.fam.\${fam}.bcf' \
  $bcfs

for f in *.fam.*.bcf; do
    # out of order because comphet + other vars.
    bcftools sort -O b -o tmp.bcf \$f
    mv tmp.bcf \$f
    bcftools index \$f
done

    """

}

process generate_jigv_pages {
  container = 'docker://brentp/rare-disease:v0.2.2'
  publishDir "${params.output_dir}/jigv_plots/", mode: 'copy'
  shell = ['/bin/bash', '-euo', 'pipefail']
  cache false


  input:
    tuple(val(family_id), path(vcfs), path(csis))
    path(xams)
    path(indexes)
    path(ped)
    path(fasta)
    path(fai)

  output:
    path("${outdir}")

  script:
  outdir = "./$family_id"
    """

awk '\$1 == "$family_id"' $ped > fam.ped

for vcf in ${vcfs}; do
    mode=\$(basename \$vcf | cut -d. -f 3)

    jigv \
      --prefix "$outdir/\${mode}/${family_id}.\${mode}." \
      --ped fam.ped \
      --sites \$vcf \
      --flank 100 \
      --fasta $fasta \
      $xams &
done
wait

    """

}

// largely cribbed from Joe: https://github.com/brwnj/smoove-nf/blob/master/main.nf
params.help = false
if (params.help) {
    log.info """
----------------------------------------------------------------
rare-disease-nf

Required Arguments:
-------------------

   --xams            A glob or string with **comma separatd paths to
                     bams or crams. Sample names should match those in the
                     --ped argument
   
   --ped             A pedigree (or .fam) file with columns of:
                         $family_id\t$sample_id\t$paternal_id\t$maternal_id\t$sex\t$phenotype
                     see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format


   --fasta           Path to reference fasta 

   --modeltype       "WGS" or "WES" (default is WGS)

   --gff             Path to gff3 file for annotation with bcftools csq.
                     can be downloaded from Ensembl. e.g. for human:
                     ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
                     ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/
  
   --slivarzip       Path to gnotate zip file for slivar see: https://github.com/brentp/slivar#gnotation-files
                     likely:
                      https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip
                     or:
                      https://s3.amazonaws.com/slivar/gnomad.hg37.zip

Optional Arguments:
-------------------

   --cohort_name     optional name for the cohort (default: "rare-disease")
   --output_dir      optional name for where to place results (default: "results-rare-disease")

    """
}
params.xams = false
if(!params.xams) { exit 1, "--xams is required" }
params.ped = false
if(!params.ped) { exit 1, "--ped is required" }
params.fasta = false
if(!params.fasta) { exit 1, "--fasta is required" }
params.gff = false
if(!params.gff) { exit 1, "--gff is required" }
params.slivarzip = false
if(!params.slivarzip) { exit 1, "--slivarzip is required" }
params.model_type = "WGS"
if(!params.model_type) { exit 1, "--model_type ('WGS' or 'WES') is required" }
params.cohort_name = "rare-disease"
params.output_dir = "results-rare-disease"



workflow {

    ped_file = "${file(params.ped).toAbsolutePath()}"
    slivarzip = "${file(params.slivarzip).toAbsolutePath()}"


    input = channel.fromPath(params.xams, checkIfExists: true)
                  | map { it -> tuple(it.baseName, it, find_index(it)) }

    gvcfs_tbis = DeepVariant(input, params.fasta, params.fasta + ".fai") 
    //  something.$chrom.split.gvcf.gz
    sp = split(gvcfs_tbis, params.fasta + ".fai")
    gr_by_chrom = sp.flatMap { it }
         | map { it -> [it, file(file(file(file(it).baseName).baseName).baseName).getExtension() ] } 
         | groupTuple(by: 1) 

    joint_by_chrom = glnexus_anno_slivar(gr_by_chrom, params.fasta, params.fasta + ".fai", params.gff, slivarzip, params.cohort_name)
    // temporary hack since slivar 0.2.1 errors on no usable comphet-side sites.
    jbf = joint_by_chrom | filter { !(it[0].toString() ==~ /.*(MT|Y).glnexus.*/) }

    slivar_results = slivar_rare_disease(jbf, ped_file, slivarzip)

    slivar_tsvs = slivar_results | map { it -> it[4] } | collect
    slivar_counts = slivar_results | map { it -> it[5] } | collect
    // get the regular and comphet bcf
    slivar_vcfs = slivar_results | map { it -> [it[0], it[1]] } | flatten | collect


    slivar_merge_tsvs(slivar_tsvs, params.cohort_name)
    slivar_sum_counts(slivar_counts)

    byFamBcfIdx = slivar_split_by_fam(ped_file, slivar_vcfs) | flatMap { it }

    byFamBcf = byFamBcfIdx | filter { "${it}".endsWith(".bcf") } 
        | map { it -> [ file(file(it).baseName).getExtension(), file(it), file(it + ".csi") ] }
        | groupTuple(by: 0)

    //byFamBcf | view

    xams = input | map { it -> it[1] } | collect
    indexes = input | map { it -> it[2] } | collect

    generate_jigv_pages(byFamBcf, xams, indexes, ped_file, params.fasta, params.fasta + ".fai")
    // TODO: get bams grouped by family. dont want to stage all crams for each family.

}
