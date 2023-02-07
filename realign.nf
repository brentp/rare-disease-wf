nextflow.enable.dsl=2

process realign {
    label "realign"
    container = 'docker://brentp/realign:v0.1.0'
    publishDir "${params.output_dir}/crams/", mode: 'copy'

    shell = ['/bin/bash', '-euo', 'pipefail']
    cpus = 16

    input:
        tuple(val(sample_id), path(aln_file), path(aln_index))
        path(old_fasta)
        path(old_fai)
        path(fasta)
        path(fai)
        path(alt)
        path(amb)
        path(ann)
        path(bwt)
        path(pac)
        path(sa)

    output:
      tuple(file("${sample_id}.cram"), file("${sample_id}.cram.crai"))


    script:
    """

echo "TMPDIR:\$TMPDIR"

# http://lh3.github.io/2021/07/06/remapping-an-aligned-bam
samtools collate --reference $old_fasta -Oun128 ${aln_file} | samtools fastq -OT RG,BC - \
  | bwa mem -pt 16 -CH <(samtools view -H ${aln_file} |grep ^@RG) $fasta - \
  | samblaster --addMateTags \
  | samtools sort --write-index -@4 -m4g -o ${sample_id}.cram -

    """
}


params.help = false
if (params.help) {
    log.info """
----------------------------------------------------------------
realing-nf: realign from one genome (bams/crams) to another with few intermediate files

Required Arguments:
-------------------

   --xams            A glob or string with **comma separatd paths to
                     bams or crams. Sample names should match those in the
   
   --old_fasta       path to fasta to which input bams/crams were aligned.    
   --fasta           Path to reference fasta for realignment (must be indexed with bwa-mem)


Optional

    --output_dir      default: results-realigned

    """
}
params.xams = false
if(!params.xams) { exit 1, "--xams is required" }
params.old_fasta = false
if(!params.old_fasta) { exit 1, "--old_fasta is required" }
params.fasta = false
if(!params.fasta) { exit 1, "--fasta is required" }
params.output_dir = "results-realigned"

workflow {

    old_fasta = "${file(params.old_fasta).toAbsolutePath()}"
    fasta = "${file(params.fasta).toAbsolutePath()}"


    input = channel.fromPath(params.xams, checkIfExists: true)
                  | map { it -> tuple(it.baseName, it, find_index(it)) }

    bams = realign(input,
                   params.old_fasta, params.old_fasta + ".fai",
                   params.fasta,
                   params.fasta + ".fai",
                   params.fasta + ".amb",
                   params.fasta + ".ann",
                   params.fasta + ".bwt",
                   params.fasta + ".pac",
                   params.fasta + ".sa")

    bams | view
}
