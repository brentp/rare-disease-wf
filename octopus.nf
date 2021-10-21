nextflow.enable.dsl=2

include  { find_index } from './nf/common'
params.help = false
if(params.help) {

    log.info("""
------------------------------------------
octopus workflow

Required Arguments
------------------

    --ped       pedigree file with the 6 required columns and a 7th column
                indicating the path to a bam or cram file for each sample.

    --fasta           Path to reference fasta
""")

}

params.ped = false
if(!params.ped) { exit 1, "--ped is required" }
params.fasta = false
if(!params.fasta) { exit 1, "--fasta reference is required" }
file(params.fasta, checkIfExists:true)
file(params.fasta + ".fai", checkIfExists:true)
file(params.ped, checkIfExists:true)
params.chunk_size = 100000000

process forest_filter {
    input: tuple(val(region), path(vcfs), path(crams), path(indexes))
           path(ref)
           path(fai)
    output: path("${output_path}")

    script:
       reg=region.replaceAll(":", "_")
       output_path="octopus.filtered.${reg}.vcf.gz"
       file("$workDir/${reg}.crams.list").withWriter { fh ->
            crams.each { cram ->
                fh.write(cram.toString()); fh.write("\n")
            }
        }
       """
octopus -R $ref -i $workDir/${reg}.crams.list --filter-vcf ${vcfs[0]} \
    --forest-model /opt/germline.v0.7.4.forest.gz \
   -o ${output_path}
      """
}

process octopus_trio {
    input: each(region)
           tuple(val(family_id), val(sample_id), val(dad_id), val(mom_id), path(kid_bam), path(dad_bam), path(mom_bam), path(kid_index), path(dad_index), path(mom_index))
           path(ref)
           path(fai)
    output: tuple(val("${family_id}"), val("${region}"), path("${output_path}"))
    script:
       output_path="${sample_id}.${region.replaceAll(':','_')}.trio.vcf.gz"
       """
hostname
which octopus
octopus -R $ref -I ${kid_bam} ${dad_bam} ${mom_bam} -M  ${mom_id} -F ${dad_id} \
    -p Y=2 chrY=2 chrM=1 MT=1 chrMT=1 --threads ${task.cpus} --one-based-indexing -T ${region} \
    -o ${output_path}
#--bamout "${sample_id}.realigned.bams/" \
       """
}

process octopus_fam_or_single {
    input: each(region)
           tuple(val(family_id), path(bams), path(indexes))
           path(ref)
           path(fai)
    output: tuple(val("${family_id}"), val("${region}"), path("${output_path}"))
    script:
       output_path="${family_id}.${region.replaceAll(':', '_')}.notrio.vcf.gz"
       bamout="${family_id}.realigned.bams.fam/"
       if (bams.size() == 1 ) {
         bamout += "octopus.${family_id}.bam"
       }
       """
hostname
which octopus
octopus -R $ref -I $bams \
    -p MT=1, chrM=1 chrMT=1 Y=2 chrY=2 --threads ${task.cpus} --one-based-indexing -T ${region} \
    -o ${output_path}
#--bamout ${bamout} \
       """
}

include { split_by_size } from "./split"
include { octopus_population as octopus_population1 } from "./octopus_pop" // 20
include { octopus_population as octopus_population2 } from "./octopus_pop" // 400
include { octopus_population as octopus_population3 } from "./octopus_pop" // 8000
include { octopus_population as octopus_population4 } from "./octopus_pop" // 160000

@groovy.transform.EqualsAndHashCode(includeFields=true, allProperties=true, allNames=true, excludes=["dad", "mom"])
@groovy.transform.ToString(includeNames=true, ignoreNulls=true, excludes=["dad", "mom"])
public class Sample {
    String id
    String family_id
    String maternal_id
    String paternal_id
    java.nio.file.Path path
    Sample mom
    Sample dad
}

workflow {

    def trios = []
    def in_trio = [:]
    def sample_by_id = [:]
    def samples = []
    file(params.ped, checkIfExists:true).eachLine { line ->
        if(line[0] == '#' ) { return }
        def toks = line.split("\t")
        if(toks.size() != 7 ){
            println("ERROR: expecting 7 columns in pedigree file; found ${toks.size()} in line\n${line}")
        }
        s = new Sample(id: toks[1], family_id:toks[0], paternal_id: toks[2], maternal_id: toks[3], 
                       path: file(toks[6], checkIfExists: true))
        sample_by_id[s.id] = s
        samples.add(s)
    }
    // get trios
    samples.each { s ->
      s.mom = sample_by_id[s.maternal_id]
      s.dad = sample_by_id[s.paternal_id]
      if(s.dad && s.mom) {
         trios.add(s)
         in_trio[s.id] = true
         in_trio[s.maternal_id] = true
         in_trio[s.paternal_id] = true
      }
    }
    // now collect other samples that are not in a trio and group them by
    // family for calling
    def non_trio_fam = [:]
    def by_fam = [:]
    samples.each { s ->
        if(!by_fam.containsKey(s.family_id)) {
            by_fam[(s.family_id)] = []
        }
        by_fam[(s.family_id)].add(s)

        if(in_trio[s.id]) { return }
        if(!non_trio_fam.containsKey(s.family_id)) {
            non_trio_fam[(s.family_id)] = []
        }
        non_trio_fam[(s.family_id)].add(s)
    }

    regions = split_by_size(params.fasta + ".fai", params.chunk_size).splitText() | map { s -> s.replaceAll("\\s", "") }

    trs = channel.fromList(trios) | map { it -> [it.family_id, it.id, it.dad.id, it.mom.id, it.path, it.dad.path, it.mom.path, find_index(it.path), find_index(it.dad.path), find_index(it.mom.path)] } 
    trio_ch = octopus_trio(regions, trs, params.fasta, params.fasta + ".fai")

    // now add families to list 
    fams = []
    non_trio_fam.each { li -> {
          fam = []
          idx = []
           
          li.value.each { it ->
            fam.add(it.path) 
            idx.add(find_index(it.path))
          }
          fams.add(tuple(li.value[0].family_id, fam, idx))
       }
    }

    group_size = 20

    // now do joint-calling with size 20 as per: https://luntergroup.github.io/octopus/docs/guides/models/population
    by_region = octopus_fam_or_single(regions, channel.fromList(fams), params.fasta, params.fasta + ".fai").concat(trio_ch) 
               | groupTuple(by: 1, size: group_size, remainder:true) \
               | map { it -> [it[0].unique(), it[1], it[2]] } // drop duplicate family ids.

    // by_region is: [[fam1_id, fam2_id, ...], $region, [vcf1, vcf2, ...]] .. vcfs might be longer than fams.

    pop_input = by_region.map { it -> [it[1], it[2], 
                                     it[0].collect(f -> by_fam[f].collect(s -> s.path)).flatten(),
                                     it[0].collect(f -> by_fam[f].collect(s -> find_index(s.path))).flatten() ] }

    // this code below iteratively calls octopus matching the instructions here:
    // https://luntergroup.github.io/octopus/docs/guides/models/population/
    // most often only the first 2 calls will be used, but octopus_population3 and
    // 4 for larger cohorts

    op1 = octopus_population1(pop_input, params.fasta, params.fasta + ".fai", "1") 
           | groupTuple(by: 0, size: group_size, remainder: true) \
           | map { it -> [it[0], it[1], it[2].flatten(), it[3].flatten() ] }

    final1 = op1 | filter { it[1].size() == 1 }
    op1 = op1 | filter { it[1].size() > 1 }


    op2 = octopus_population2(op1, params.fasta, params.fasta + ".fai", "2") 
           | groupTuple(by: 0, size: group_size, remainder: true) \
           | map { it -> [it[0], it[1], it[2].flatten(), it[3].flatten() ] }

    final2 = op2 | filter { it[1].size() == 1 }
    op2 = op2 | filter { it[1].size() > 1 }

    op3 = octopus_population3(op2, params.fasta, params.fasta + ".fai", "3") 
           | groupTuple(by: 0, size: group_size, remainder: true) \
           | map { it -> [it[0], it[1], it[2].flatten(), it[3].flatten() ] }

    final3 = op3 | filter { it[1].size() == 1 }
    op3 = op3 | filter { it[1].size() > 1 }

    op4 = octopus_population4(op3, params.fasta, params.fasta + ".fai", "4") 
           | groupTuple(by: 0, size: group_size, remainder: true) \
           | map { it -> [it[0], it[1], it[2].flatten(), it[3].flatten() ] }

    final4 = op4 | filter { it[1].size() == 1 }

    final_called = final1 | concat(final2, final3, final4) 
    final_called | view
    final_result = forest_filter(final_called, params.fasta, params.fasta + ".fai")
    final_result | view

}
