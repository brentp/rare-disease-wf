nextflow.enable.dsl=2

process DeepVariant {
    label "DeepVariant"
    container = 'docker://gcr.io/deepvariant-docker/deepvariant:1.1.0'
    publishDir "results-rare-disease/gvcfs/", mode: 'copy'

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
#    samtools view --write-index --threads ${params.cpus} -bT $fasta -o ${sample_id}.bam $aln_file
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
    container = 'docker://brentp/rare-disease:v0.0.4'
    shell = ['/bin/bash', '-euo', 'pipefail']
    publishDir "results-rare-disease/joint-by-chrom/", mode: 'copy'

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
# TODO: can't yet get snpEff working 
# | snpEff ann -noStats -dataDir {params.snpeff_data_dir} GRCh37.75

glnexus_cli \
    -t ${task.cpus} \
    --mem-gbytes 128 \
    --config DeepVariant${params.model_type} \
    --list $workDir/file.list.${cohort_name}.${chrom} \
| bcftools norm --threads 3 -m - -w 10000 -f $fasta -O u \
| bcftools csq --threads 3 -s - --ncsq 50 -g $gff -l -f $fasta - -o - -O u \
| slivar expr -g $slivar_zip -o $output_file --vcf -

bcftools index --threads 6 $output_file
    """

}

process slivar_rare_disease {
  container = 'docker://brentp/rare-disease:v0.0.5'
  publishDir "results-rare-disease/joint-by-chrom-slivar/", mode: 'copy'
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
  -c BCSQ \
  -g /opt/slivar/pli.lookup \
  -g /opt/slivar/clinvar_gene_desc.txt \
  -p $ped \
  $slivar_bcf > $slivar_tsv

slivar tsv \
  -s slivar_comphet \
  -i gnomad_popmax_af -i gnomad_popmax_af_filter -i gnomad_nhomalt \
  -i impactful -i genic \
  -c BCSQ \
  -g /opt/slivar/pli.lookup \
  -g /opt/slivar/clinvar_gene_desc.txt \
  -p $ped \
  $slivar_ch_bcf \
  | { grep -v ^# || true; } >> $slivar_tsv # || true avoids error if there are no compound hets.


wait
  """
}

process slivar_merge_tsvs {
  container = 'docker://brentp/rare-disease:v0.0.6'
  publishDir "results-rare-disease/", mode: 'copy'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input: path(tsvs)

  output: path("${output_file}")

  script:
    output_file = "slivar.candidates.tsv"
    """
# get header from first file and drop it from other files
    awk 'NR == FNR || FNR > 1' $tsvs > ${output_file}
    """
}

process slivar_sum_counts {
  container = 'docker://brentp/rare-disease:v0.0.6'
  publishDir "results-rare-disease/", mode: 'copy'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input: path(counts)

  output: path("${output_file}")

  script:
    output_file = "slivar.counts.txt"
    """
    tiwih sum_slivar_counts -z $counts > ${output_file}
    """
}


workflow {

  //  split(["/data/human/HG002_SVs_Tier1_v0.6.DEL.vcf.gz", "/data/human/HG002_SVs_Tier1_v0.6.DEL.vcf.gz.tbi", "/data/human/g1k_v37_decoy.fa.fai"]) | view
   // DeepVariant(["HG002", "/data/human/hg002.cram", "/data/human/hg002.cram.crai", "/data/human/g1k_v37_decoy.fa", "/data/human/g1k_v37_decoy.fa.fai"]) | view
    fasta = "/hpc/cog_bioinf/GENOMES.old/NF-IAP-resources//GRCh37/Sequence/genome.fa"
    gff = "/home/cog/bpedersen/src/rare-disease-wf/Homo_sapiens.GRCh37.87.chr.gff3.gz"
    slivar_zip = "/hpc/compgen/users/bpedersen/gnomad.hg37.zip"
    samples = [
        ["150424",
        "/hpc/cog_bioinf/ubec/useq/processed_data/external/REN5302/REN5302_5/BAMS/150424_dedup.bam",
        "/hpc/cog_bioinf/ubec/useq/processed_data/external/REN5302/REN5302_5/BAMS/150424_dedup.bai"],
        ["150423",
        "/hpc/cog_bioinf/ubec/useq/processed_data/external/REN5302/REN5302_5/BAMS/150423_dedup.bam",
        "/hpc/cog_bioinf/ubec/useq/processed_data/external/REN5302/REN5302_5/BAMS/150423_dedup.bai"],
        ["150426",
        "/hpc/cog_bioinf/ubec/useq/processed_data/external/REN5302/REN5302_5/BAMS/150426_dedup.bam",
        "/hpc/cog_bioinf/ubec/useq/processed_data/external/REN5302/REN5302_5/BAMS/150426_dedup.bai"]
    ]
    //samples = [
    // [ "hg002",
    //   "/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/NA24385/HG002.hs37d5.2x250.bam"
    //   "/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/NA24385/HG002.hs37d5.2x250.bam.bai"
    // ]]


    input = channel.fromList(samples)

    gvcfs_tbis = DeepVariant(input, fasta, fasta + ".fai") 
    //  something.$chrom.split.gvcf.gz
    sp = split(gvcfs_tbis, fasta + ".fai")
    gr_by_chrom = sp.flatMap { it }
         | map { it -> [it, file(file(file(file(it).baseName).baseName).baseName).getExtension() ] } 
         | groupTuple(by: 1) 

    joint_by_chrom = glnexus_anno_slivar(gr_by_chrom, fasta, fasta + ".fai", gff, slivar_zip, params.cohort)
    // temporary hack since slivar 0.2.1 errors on no usable comphet-side sites.
    jbf = joint_by_chrom | filter { !(it[0].toString() ==~ /.*(MT|Y).glnexus.*/) }

    slivar_results = slivar_rare_disease(jbf, params.ped, slivar_zip)

    slivar_tsvs = slivar_results | map { it -> it[4] } | collect
    slivar_counts = slivar_results | map { it -> it[5] } | collect

    slivar_tsvs | view

    slivar_merge_tsvs(slivar_tsvs)
    slivar_sum_counts(slivar_counts)

}
