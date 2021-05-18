nextflow.enable.dsl=2

include { split_by_size } from "./split"


process manta {
    errorStrategy 'terminate' // TODO: change after debugging is done
    stageInMode "copy"

    container = 'docker://brentp/manta-graphtyper:v0.1.8'
    publishDir "results-rare-disease/manta-vcfs/", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(sample_name), path(bam), path(index))
        path(fasta)
        path(fai)

    output:
        tuple(file("${sample_name}.diploidSV.vcf.gz"), file("${sample_name}.meandepth.txt"))

    script:
    """
# limit to larger chroms ( > 10MB)
awk '\$2 > 10000000 || \$1 ~/(M|MT)\$/ { print \$1"\t0\t"\$2 }' $fai | bgzip -c > cr.bed.gz
tiwih meandepth --scale-by-read-length $bam > ${sample_name}.meandepth.txt

tabix cr.bed.gz
configManta.py --bam ${bam} --referenceFasta $fasta --runDir . --callRegions cr.bed.gz
python2 ./runWorkflow.py -j ${task.cpus}
# TODO: run convert inversions.
mv results/variants/diploidSV.vcf.gz ${sample_name}.diploidSV.vcf.gz
mv results/variants/candidateSV.vcf.gz ${sample_name}.candidateSV.vcf.gz
mv results/variants/candidateSmallIndels.vcf.gz ${sample_name}.candidateSmallIndels.vcf.gz
rm -rf results/
    """
}

process svimmer {
    container = 'docker://brentp/manta-graphtyper:v0.1.8'
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
    # NOTE: removing BNDs
    svimmer --max_distance 100 --max_size_difference 60 $workDir/vcfs.list \$chroms \
        | tiwih setsvalt --drop-bnds /dev/stdin  -o $output_file
    tabix $output_file
    """
}

process graphtyper_sv {
    errorStrategy 'terminate' // TODO: change after debugging is done
    shell = ['/bin/bash', '-euo', 'pipefail']
    container = 'docker://brentp/manta-graphtyper:v0.1.8'

    publishDir "results-rare-disease/svs-joint/", mode: 'copy'
    // TODO: make task.cpus depend on n-samples as well

    input:
        val(region)
        val(samples_bams_indexes)
        tuple(path(merged_sv_vcf), path(vcf_index))
        path(fasta)
        path(fai)
        val(mdps)
    output:
        tuple(file("genotyped-svs.${reg}.bcf"), file("genotyped-svs.${reg}.bcf.csi"))


    script:
    chrom_pos = "${region}".split(":")
    pos = chrom_pos[1].split("-")
    reg = "${chrom_pos[0]}_${pos[0].padLeft(12, '0')}_${pos[1].padLeft(12, '0')}"
    // TODO: how to move these outside of here so we just send in the files that we've written once.
    file("$workDir/bams.list.${reg}").withWriter { fh ->
            samples_bams_indexes.each { bi ->
                fh.write(bi[1].toString()); fh.write("\n")
            }
    }

    file("$workDir/covs.list.${reg}").withWriter { fh ->
        mdps.each { f ->
            fh.write(file(f).text)
        }
    }

    """
    bcftools index $merged_sv_vcf # TODO: pass in index

    graphtyper genotype_sv $fasta $merged_sv_vcf \
        --sams=$workDir/bams.list.$reg \
        --threads=${task.cpus} \
        --avg_cov_by_readlen=$workDir/covs.list.${reg} \
        --force_use_input_ref_for_cram_reading \
        --region $region \
        --output=graphtyper_sv_results/

    ls graphtyper_sv_results/*/*.vcf.gz > file.list
    bcftools concat --threads 3 -O u -o - --file-list file.list \
       | bcftools sort -m 2G -T \$TMPDIR -o genotyped-svs.${reg}.bcf -O b -
    bcftools index --threads 4 genotyped-svs.${reg}.bcf
    """
}
process paragraph {
  errorStrategy 'terminate' // TODO: change after debugging is done
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'docker://brentp/manta-graphtyper:v0.1.8'

  publishDir "results-rare-disease/svs-paragraph/", mode: 'copy'

  input:
     tuple(path(site_vcf), path(site_vcf_index))
     tuple(val(sample), path(bam), path(index))
     path(fasta)
     path(fai)

  output: path("*.paragraph.vcf.gz")

  script:
  """
samplename=$sample
dp=\$(tiwih meandepth $bam)
echo "id\tpath\tdepth\tread length" > sample.manifest
echo "$sample\t$bam\t\$dp\t150" >> sample.manifest
M=\$((dp * 7))
cat sample.manifest

multigrmpy.py -i $site_vcf \
    -m sample.manifest \
    -r $fasta \
    -o t \
    -t ${task.cpus} \
    -M \$M

mv t/genotypes.vcf.gz \${samplename}.paragraph.vcf.gz

  """

}

process merge_svs {
    errorStrategy 'terminate' 
    shell = ['/bin/bash', '-euo', 'pipefail']
    container = 'docker://brentp/manta-graphtyper:v0.1.8'

    publishDir "results-rare-disease/", mode: 'copy'
    // TODO: make task.cpus depend on n-samples as well

    input: path(vcfs)
           path(fasta)
           path(fai)
    output: tuple(file("${output_file}"), file("${output_file}.csi"))

    script:
    output_file = "rare-disease.manta-graphtyper.bcf"
    """
      ls *.bcf | sort -V > file.list
      bcftools concat --threads 2 --naive -O b -f file.list \
          | bcftools view --threads 2 -i "SVMODEL='AGGREGATED'" -O b \
          | tiwih setref -u -o $output_file /dev/stdin $fasta
      bcftools index --threads 2 $output_file
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

    mr = manta_results | map { it -> it[0] }  | collect
    mds = manta_results | map { it -> it[1] }  | collect

    sv_merged = svimmer(mr, fasta + ".fai")

    paragraph(sv_merged, input, fasta, fasta + ".fai")



    /*
    regions = split_by_size(fasta + ".fai", 2000000).splitText()
       | map { it -> it.trim() }

    gt_vcfs = graphtyper_sv(regions, input.toList(), sv_merged, fasta, fasta + ".fai", mds)
       | map { it -> file(it[0]) }  | collect
    merge_svs(gt_vcfs, fasta, fasta + ".fai")
    */

}
