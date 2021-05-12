
process manta {
    errorStrategy 'terminate' // TODO: change after debugging is done

    container = 'docker://brentp/manta-graphtyper:v0.0.9'
    publishDir "results-rare-disease/manta-vcfs/", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']
    input:
        tuple(val(sample_name), path(bam), path(index))
        path(fasta)
        path(fai)

    output:
        path("*.vcf.gz")

    script:
    """
# limit to larger chroms ( > 10MB)
awk '\$2 > 10000000 || \$1 ~/(M|MT)\$/ { print \$1"\t0\t"\$2 }' $fai | bgzip -c > cr.bed.gz
D=\$TMPDIR
# manta hits the ref a lot so copy to local tmp
cp $fasta \$D/ref.fa
cp $fai \$D/ref.fa.fai
# NOTE: we are copying the bam to local TMP to save on network
# TODO: make this depend on e.g param.copy_bam_local
cp $bam \$D/${sample_name}.bam
cp $index \$D/${sample_name}.bam.bai


tabix cr.bed.gz
configManta.py --bam \$D/${sample_name}.bam --referenceFasta \$D/ref.fa --runDir \$D --callRegions cr.bed.gz
python2 \$D/runWorkflow.py -j ${task.cpus}
mv \$D/results/variants/diploidSV.vcf.gz ${sample_name}.diploidSV.vcf.gz
mv \$D/results/variants/candidateSV.vcf.gz ${sample_name}.candidateSV.vcf.gz
mv \$D/results/variants/candidateSmallIndels.vcf.gz ${sample_name}.candidateSmallIndels.vcf.gz
rm -rf \$D/results/
    """
}

process svimmer {
    container = 'docker://brentp/manta-graphtyper:v0.0.9'
    publishDir "results-rare-disease/manta-merged/", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        val(sample_vcfs)
        path(fai)
    output: tuple(file("${output_file}"), file("${output_file}.tbi"))

    script:
    output_file = "svimmer.merged.vcf.gz"
    file("$workDir/vcfs.list").withWriter { fh ->
            sample_vcfs.each { vcf ->
                fh.write(vcf.toString()); fh.write("\n")
            }
    }
    """
    cat $workDir/vcfs.list | xargs -I{} -P ${task.cpus} bcftools index -f --threads 2 {}
    chroms=\$(awk '\$2 > 10000000 { printf("%s ", \$1) }' $fai) 
    svimmer --max_distance 100 --max_size_difference 60 $workDir/vcfs.list \$chroms \
        | bgzip --threads 3 > $output_file
    tabix $output_file
    """
}

process graphtyper_sv {
    shell = ['/bin/bash', '-euo', 'pipefail']
    container = 'docker://brentp/manta-graphtyper:v0.0.9'

    publishDir "results-rare-disease/svs-joint/", mode: 'copy'
    // TODO: make task.cpus depend on n-samples as well

    input:
        val(region)
        val(samples_bams_indexes)
        tuple(path(merged_sv_vcf), path(vcf_index))
        path(fasta)
        path(fai)
    output:
        tuple(file("genotyped-svs.${reg}.bcf"), file("genotyped-svs.${reg}.bcf.csi"))


    script:
    chrom_pos = "${region}".split(":")
    pos = chrom_pos[1].split("-")
    reg = "${chrom_pos[0]}_${pos[0].padLeft(12, '0')}_${pos[1].padLeft(12, '0')}"
    file("$workDir/bams.list.${reg}").withWriter { fh ->
            samples_bams_indexes.each { bi ->
                fh.write(bi[1].toString()); fh.write("\n")
            }
    }
    """
    bcftools index --threads 4 $merged_sv_vcf # TODO: pass in index
    cat $workDir/bams.list.$reg | parallel -j ${task.cpus} -k "tiwih meandepth --scale-by-read-length {1}" > $workDir/avg.cov.$reg

    graphtyper genotype_sv $fasta $merged_sv_vcf \
        --sams=$workDir/bams.list.$reg \
        --threads=${task.cpus} \
        --avg_cov_by_readlen=$workDir/avg.cov.$reg \
        --force_use_input_ref_for_cram_reading \
        --region $region \
        --output=graphtyper_sv_results/

    ls graphtyper_sv_results/*/*.vcf.gz > file.list
    bcftools concat --threads 3 -O u -o - --file-list file.list \
       | bcftools sort -m 2G -T \$TMPDIR -o genotyped-svs.${reg}.bcf -O b -
    bcftools index --threads 4 genotyped-svs.${reg}.bcf
    """
}

workflow {

    fasta = "/hpc/cog_bioinf/GENOMES.old/NF-IAP-resources//GRCh37/Sequence/genome.fa"
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
    input = channel.fromList(samples)


    manta_results = manta(input, fasta, fasta + ".fai")

    mr = manta_results 
        | map { it -> it.find { it =~ /diploidSV.vcf.gz/ } } | collect


    sv_merged = svimmer(mr, fasta + ".fai")

    regions = split_by_size(fasta + ".fai", 2000000).splitText()
       | map { it -> it.trim() }

    graphtyper_sv(regions, input.toList(), sv_merged, fasta, fasta + ".fai")

}
