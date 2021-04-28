nextflow.enable.dsl=2

process DeepVariant {
    label "DeepVariant"
    container = 'docker://gcr.io/deepvariant-docker/deepvariant:1.1.0'
    publishDir "results-rare-disease", mode: 'copy'
    // container = '/hpc/compgen/users/bpedersen/deepvariant_1_1_0.sif'

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
    --num_shards=${params.cpus} \
    --ref=$fasta
 
#rm -f ${sample_id}.bam
#bcftools index --threads 6 ${sample_id}.gvcf.gz
    """
}

process split {
    container = 'docker://brentp/rare-disease:v0.0.1'

    input: tuple(path(gvcf), path(tbi))
           path(fai)
    output: file("*.split.gvcf.gz")

    script:
    """
    # get large chroms and chrM in one file each
    for chrom in \$(awk '\$2 > 40000000 || \$1 ~/(M|MT)\$/' $fai | cut -f 1); do
        bcftools view $gvcf --threads 3 -O b -o \$(basename $gvcf .gvcf.gz).\${chrom}.split.gvcf.gz \$chrom
    done
    # small HLA and gl chroms all to go single file
    awk '!(\$2 > 40000000 || \$1 ~/(M|MT)\$/) { print \$1"\\t0\\t"\$2+1 }' $fai > other_chroms
    # if it's non-empty then create the extras / other split file
    if [ -s other_chroms ]; then
        bcftools view $gvcf --threads 3 -R other_chroms -O b -o \$(basename $gvcf .gvcf.gz).OTHER.split.gvcf.gz
    fi
    """

}

process glnexus_anno_slivar {
    tag {"GLNexus BCFtools Slivar ${cohort_name}"}
    label {"GLNexus"}
    container = 'docker://brentp/rare-disease:v0.0.1'

    input: tuple(val(gvcfs), val(cohort_name), val(chrom))
    output: tuple(path(output_file), path(output_csi))


    script:
    file("$workDir/file.list.${cohort_name}").withWriter { fh ->
        gvcfs.each { gvcf ->
            fh.write(gvcf.toString()); fh.write("\n")
        }
    }
    output_file = "${cohort_name}.glnexus.anno.bcf"
    output_csi = "${output_file}.csi"
    """
glnexus_cli \
    -t ${params.cpus} \
    --mem-gbytes 128 \
    --config DeepVariant${params.model_type} \
    --list $workDir/file.list.${cohort_name} \
| bcftools norm --threads 3 -m - -w 10000 -f $fasta -O v - \
| snpEff -Xmx4G eff -noStats GRCh38.99 \
| bcftools csq --threads 3 -s - --ncsq 40 -g /data/Homo_sapiens.GRCh38.95.chr.prefix.gff3.gz -l -f $fasta - -o - -O u \
| slivar expr -g /data/gnomad.hg38.genomes.v3.fix.zip -o $output_file

bcftools index --threads 6 $output_file
    """

}

process slivar_rare_disease {
  input: path(bcf)
         path(ped)
         path(gnomad_zip)
  output: tuple(path(slivar_bcf), path(slivar_ch_bcf) path(slivar_bcf_csi), path(slivar_ch_bcf_csi), path(slivar_tsv))

  script:
  slivar_bcf = getBaseName(bcf) + ".slivar.bcf"
  slivar_ch_bcf = getBaseName(bcf) + ".slivar.ch.bcf"
  slivar_bcf_csi = slivar_bcf + ".csi"
  slivar_ch_bcf_csi = slivar_ch_bcf + ".csi"
  slivar_tsv = getBaseName(bcf) + ".slivar.tsv"

  script:
  """
# NOTE: we do *NOT* limit to impactful so that must be reported and used by slivar tsv
slivar expr --vcf $slivar_bcf \
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

slivar compound-hets -v $slivar_bcf \
    --sample-field comphet_side --sample-field denovo -p $ped \
  | bcftools view -O b -o $slivar_ch_bcf

bcftools index --threads 3 $slivar_ch_bcf &

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
  | grep -v ^# >> $slivar_tsv

wait
  """
}

workflow {

  //  split(["/data/human/HG002_SVs_Tier1_v0.6.DEL.vcf.gz", "/data/human/HG002_SVs_Tier1_v0.6.DEL.vcf.gz.tbi", "/data/human/g1k_v37_decoy.fa.fai"]) | view
   // DeepVariant(["HG002", "/data/human/hg002.cram", "/data/human/hg002.cram.crai", "/data/human/g1k_v37_decoy.fa", "/data/human/g1k_v37_decoy.fa.fai"]) | view
    fasta = "/hpc/cog_bioinf/GENOMES/NF-IAP-resources//GRCh37/Sequence/genome.fa"
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
    cinput = channel.fromList(samples)

    gvcfs_tbis = DeepVariant(cinput, fasta, fasta + ".fai") 
    //split(gvcfs_tbis, fasta + ".fai") | view

}
