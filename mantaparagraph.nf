nextflow.enable.dsl=2

include { split_by_size } from "./split"


process manta {
    errorStrategy 'terminate' // TODO: change after debugging is done

    container = 'docker://brentp/manta-paragraph:v0.2.4'
    publishDir "results-rare-disease/manta-sample-vcfs/", mode: 'copy'
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

convertInversion.py \$(which samtools) $fasta results/variants/diploidSV.vcf.gz \
    | bcftools view -O z -o ${sample_name}.diploidSV.vcf.gz

mv results/variants/candidateSV.vcf.gz ${sample_name}.candidateSV.vcf.gz
mv results/variants/candidateSmallIndels.vcf.gz ${sample_name}.candidateSmallIndels.vcf.gz
rm -rf results/
    """
}

process jasmine {
    container = 'docker://brentp/manta-paragraph:v0.2.4'
    publishDir "results-rare-disease/jasmine-merged-sites/", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        val(sample_vcfs)
        path(fai)
    output: tuple(file("${output_file}"), file("${output_file}.tbi"))

    script:
    output_file = "jasmine.merged.vcf.gz"
    file("$workDir/vcfs.list").withWriter { fh ->
            sample_vcfs.each { vcf ->
                fh.write(vcf.toString()); fh.write("\n")
            }
    }
    // we don't merge if we only have a single sample.
    if(sample_vcfs.size() > 1) {
        """
	# jasmine can't do gzip.
        cat $workDir/vcfs.list | xargs -I{} -P ${task.cpus} sh -c 'bcftools view -O v {} > {}.un &&  mv {}.un {}'

        chroms=\$(awk '\$2 > 10000000 { printf("%s ", \$1) }' $fai) 
        # NOTE: removing BNDs and setting min start to > 150 as paragraph fails if start < readlength
        jasmine --file_list=${workDir}/vcfs.list out_file=${output_file}.tmp.vcf.gz
        # jasmine currently doesn't output PRECISE which causes problems.
	bcftools view -h ${output_file}.tmp.vcf.gz > header.txt
	head -n -1 header.txt > h.txt
	echo '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">' >> h.txt
	tail -n 1 header.txt >> h.txt
	bcftools reheader -h h.txt -o o.vcf.gz ${output_file}.tmp.vcf.gz
	bcftools sort -m 2G -O z -o ${output_file}.tmp.vcf.gz o.vcf.gz
# NOTE: changed the path to local for debugging
        ~/bin/tiwih setsvalt --drop-bnds --inv-2-ins ${output_file}.tmp.vcf.gz  -o $output_file
        tabix $output_file
        """
    } else {
        """
        tiwih setsvalt --drop-bnds ${sample_vcfs[0]} -o $output_file
        tabix $output_file
        """
    }

}

process paragraph_duphold {
  errorStrategy 'terminate' // TODO: change after debugging is done
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'docker://brentp/manta-paragraph:v0.2.5'

  publishDir "results-rare-disease/paragraph-genotyped-sample-vcfs/", mode: 'copy'

  input:
     tuple(path(site_vcf), path(site_vcf_index))
     tuple(val(sample), path(bam), path(index))
     path(fasta)
     path(fai)

  output: tuple(path("${output_file}"), path("${output_file}.csi"))

  script:
  output_file = "${sample}.paragraph.vcf.gz"

  """
dp=\$(tiwih meandepth $bam)
echo "id\tpath\tdepth\tread length" > sample.manifest
echo "$sample\t$bam\t\$dp\t150" >> sample.manifest
M=\$((dp * 5))
cat sample.manifest

multigrmpy.py -i $site_vcf \
    -m sample.manifest \
    -r $fasta \
    -o t \
    -t ${task.cpus} \
    -M \$M

duphold -d -v t/genotypes.vcf.gz -b $bam -f $fasta -t 4 -o $output_file
bcftools index --threads 3 $output_file
  """

}

process square {
  errorStrategy 'terminate' // TODO: change after debugging is done
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'docker://brentp/manta-paragraph:v0.2.4'
  publishDir "results-rare-disease/", mode: 'copy'

  input: val(sample_vcfs)
  output: tuple(file("${output_file}"), file("${output_file}.csi"))

  script:
  file("$workDir/joint.vcfs.list").withWriter { fh ->
        sample_vcfs.each { vcf ->
            fh.write(vcf[0].toString()); fh.write("\n")
        }
  }
  output_file = "mantaparagraph.vcf.gz"
  """
  bcftools merge -m none --threads 3 -O u --file-list $workDir/joint.vcfs.list \
    | bcftools annotate --threads 3 -x "INFO/GRMPY_ID" -O z -o $output_file
  bcftools index --threads 3 $output_file
  """
}

workflow {

    //fasta = "/media/brentp/transcend/data/human/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    //samples = [
    //    ["HG01051",
    //   "/media/brentp/transcend/data/1kg-pur-trio/HG01051.final.cram",
    //    "/media/brentp/transcend/data/1kg-pur-trio/HG01051.final.cram.crai"],
    //    ["HG01052",
    //    "/media/brentp/transcend/data/1kg-pur-trio/HG01052.final.cram",
    //    "/media/brentp/transcend/data/1kg-pur-trio/HG01052.final.cram.crai"],
    //    ["HG01053",
    //    "/media/brentp/transcend/data/1kg-pur-trio/HG01051.final.cram",
    //    "/media/brentp/transcend/data/1kg-pur-trio/HG01051.final.cram.crai"]
    //]

    samples = [
     ["HG01051", "/hpc/compgen/users/bpedersen/1kg-hc-trio/HG01051.final.cram", "/hpc/compgen/users/bpedersen/1kg-hc-trio/HG01051.final.cram.crai"],
     ["HG01052", "/hpc/compgen/users/bpedersen/1kg-hc-trio/HG01052.final.cram", "/hpc/compgen/users/bpedersen/1kg-hc-trio/HG01052.final.cram.crai"],
     ["HG01053", "/hpc/compgen/users/bpedersen/1kg-hc-trio/HG01053.final.cram", "/hpc/compgen/users/bpedersen/1kg-hc-trio/HG01053.final.cram.crai"],
     ["NA19256", "/hpc/compgen/users/bpedersen/1kg-hc-trio/NA19256.final.cram", "/hpc/compgen/users/bpedersen/1kg-hc-trio/NA19256.final.cram.crai"],
     ["NA19257", "/hpc/compgen/users/bpedersen/1kg-hc-trio/NA19257.final.cram", "/hpc/compgen/users/bpedersen/1kg-hc-trio/NA19257.final.cram.crai"],
     ["NA19258", "/hpc/compgen/users/bpedersen/1kg-hc-trio/NA19258.final.cram", "/hpc/compgen/users/bpedersen/1kg-hc-trio/NA19258.final.cram.crai"]
    ]

    fasta = "/hpc/cog_bioinf/GENOMES.old/1KP_GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa"


    input = channel.fromList(samples)

    manta_results = manta(input, fasta, fasta + ".fai")

    mr = manta_results | map { it -> it[0] }  | collect
    mds = manta_results | map { it -> it[1] }  | collect

    sv_merged = jasmine(mr, fasta + ".fai")

    genotyped = paragraph_duphold(sv_merged, input, fasta, fasta + ".fai")

    square(genotyped.toList()) | view

}
