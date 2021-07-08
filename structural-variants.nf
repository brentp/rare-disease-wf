nextflow.enable.dsl=2

include  { find_index } from './nf/common'

process manta {
    errorStrategy 'terminate' // TODO: change after debugging is done

    container = 'docker://brentp/rare-disease-sv:v0.0.5'
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

process dysgu {
    errorStrategy 'terminate' // TODO: change after debugging is done

    container = 'docker://brentp/rare-disease-sv:v0.0.5'
    publishDir "results-rare-disease/dysgu-sample-vcfs/", mode: 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(sample_name), path(bam), path(index))
        path(fasta)
        path(fai)

    output:
        tuple(file("${output_file}"), file("${output_file}.tbi"))

script:
    output_file = "${sample_name}.dysgu.vcf.gz"
    """
M=\$((dp * 5))
dysgu run --clean \
    --pl pe --mode pe \
    -o dysgu.${sample_name}.vcf \
    -p ${task.cpus} \
    --max-cov \$M \
    --thresholds 0.25,0.25,0.25,0.25,0.25 \
    $fasta \${TMPDIR}/dysgu.${sample_name} ${bam}

#$fasta \${TMPDIR}/dysgu.${sample_name} x.${sample_name}.bam
#cp dysgu.vcf ~/dysgu.${sample_name}.vcf

bcftools sort --temp-dir \$TMPDIR -m 2G -O z -o ${output_file} dysgu.${sample_name}.vcf
bcftools index --tbi ${output_file}
    """

}

process jasmine {
    container = 'docker://brentp/rare-disease-sv:v0.0.5'
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
                // write .vcf even though it's really .vcf.gz since jasmine doesn't accept gz
                // and we change the file below.
                fh.write(vcf.toString().take(vcf.toString().lastIndexOf('.'))); fh.write("\n")
            }
    }
    // we don't merge if we only have a single sample.
    if(sample_vcfs.size() > 1) {
        """
        # jasmine can't do gzip.
        cat $workDir/vcfs.list | xargs -I{} -P ${task.cpus} sh -c 'bcftools view -O v {}.gz -o {}'

        # NOTE: removing BNDs and setting min start to > 150 as paragraph fails if start < readlength
        jasmine --file_list=${workDir}/vcfs.list out_file=${output_file}
        tiwih setsvalt --drop-bnds --inv-2-ins -o ${output_file}.tmp.vcf.gz $output_file
        bcftools sort --temp-dir \$TMPDIR -m 2G -O z -o ${output_file} ${output_file}.tmp.vcf.gz
        bcftools index --tbi $output_file
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
  container = 'docker://brentp/rare-disease-sv:v0.0.5'

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
tsample=\$(tiwih samplename $bam)
echo "id\tpath\tdepth\tread length" > sample.manifest
echo "\$tsample\t$bam\t\$dp\t150" >> sample.manifest
M=\$((dp * 5))
cat sample.manifest

# this is the main paragraph entrypoint
multigrmpy.py -i $site_vcf \
    -m sample.manifest \
    -r $fasta \
    -o t \
    -t ${task.cpus} \
    -M \$M


# duphold adds depth annotations looking at coverage fold-change around Svs
duphold -d -v t/genotypes.vcf.gz -b $bam -f $fasta -t 4 -o $output_file
bcftools index --threads 3 $output_file
  """

}

process square_svcsq {
  errorStrategy 'terminate' // TODO: change after debugging is done
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'docker://brentp/rare-disease-sv:v0.0.5'
  publishDir "results-rare-disease/", mode: 'copy'

  input: val(sample_vcfs)
         file(gff)
         val(cohort_name)
  output: tuple(file("${output_file}"), file("${output_file}.csi"))

  script:
  file("$workDir/joint.vcfs.list").withWriter { fh ->
        sample_vcfs.each { vcf ->
            fh.write(vcf[0].toString()); fh.write("\n")
        }
  }
  output_file = "${cohort_name}.svs.vcf.gz"
  """
  bcftools merge -m none --threads 3 -O u --file-list $workDir/joint.vcfs.list \
    | bcftools annotate --threads 3 -x "INFO/GRMPY_ID" -O z -o ${output_file}.tmp.vcf.gz
  svpack consequence ${output_file}.tmp.vcf.gz $gff | bgzip -c > $output_file
  bcftools index --threads 3 $output_file
  """
}

params.help = false
if (params.help) {
    log.info("""
-------------------------------------------------------------------------
rare-disease-wf sv

Required Arguments:

   --xams            A glob or string with **comma separatd paths to
                     bams or crams. Sample names should match those in the
                     --ped argument

   --ped             A pedigree (or .fam) file with columns of:
                         $family_id\t$sample_id\t$paternal_id\t$maternal_id\t$sex\t$phenotype
                     see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format


   --fasta           Path to reference fasta

   --gff             Path to gff3 file for annotation with bcftools csq.
                     can be downloaded from Ensembl. e.g. for human:
                     ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
                     ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/

Optional Arguments:
-------------------

   --cohort_name     optional name for the cohort (default: "rare-disease")
   --output_dir      optional name for where to place results (default: "results-rare-disease")
   """)

}

params.xams = false
if(!params.xams) { exit 1, "--xams is required" }
params.ped = false
if(!params.ped) { exit 1, "--ped is required" }
params.fasta = false
if(!params.fasta) { exit 1, "--fasta is required" }
params.gff = false
if(!params.gff) { exit 1, "--gff is required" }
params.cohort_name = "rare-disease"
params.output_dir = "results-rare-disease"


workflow {

    input = channel.fromPath(params.xams, checkIfExists: true) 
                  | map { it -> tuple(it.baseName, it, find_index(it)) }

    fasta = file(params.fasta, checkIfEsits: true)
    fai = file(params.fasta + ".fai", checkIfEsits: true)

    manta_results = manta(input, fasta, fai)
    dysgu_results = dysgu(input, fasta, fai)

    mr = manta_results | map { it -> it[0] } 
    dr = dysgu_results | map { it -> it[0] } 

    svs = mr.concat(dr) | collect

    //mds = manta_results | map { it -> it[1] }  | collect

    sv_merged = jasmine(svs, fasta + ".fai")

    genotyped = paragraph_duphold(sv_merged, input, fasta, fasta + ".fai")

    square_svcsq(genotyped.toList(), params.gff, params.cohort_name) | view

}
