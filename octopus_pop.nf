process octopus_population {
    input: tuple(val(region), path(vcfs), path(crams), path(indexes))
           path(ref)
           path(fai)
           val(index) // 1, 2, 3, or 4 just for file naming since we run this
           // process iteratively
    output: tuple(val("${region}"), path("${output_path}"), path(crams), path(indexes))
    // TODO: use --bamout bams from previous? instead of original crams?
    script:
       reg=region.replaceAll(":", "_")
       output_path="${reg}.${index}.population.vcf.gz"
       file("$workDir/${reg}.${index}.vcfs.list").withWriter { fh ->
            vcfs.each { vcf ->
                fh.write(vcf.toString()); fh.write("\n")
            }
        }
       file("$workDir/${reg}.${index}.crams.list").withWriter { fh ->
            crams.each { cram ->
                fh.write(cram.toString()); fh.write("\n")
            }
        }
       if(vcfs.size() > 1) {

       """
which bcftools
bcftools --version
octopus --version
echo ${workflow.projectDir}
while read path; do
    bcftools index --threads 4 \$path
done < $workDir/${reg}.${index}.vcfs.list
octopus -R $ref \
    -p Y=2 chrY=2 chrM=1 chrMT=1 MT=1 --threads ${task.cpus} --one-based-indexing \
    -T ${region} \
    --disable-denovo-variant-discovery \
    -i $workDir/${reg}.${index}.crams.list \
    --source-candidates-file $workDir/${reg}.${index}.vcfs.list \
    -o ${output_path}
       """
     } else {
      """
ln -s ${vcfs[0]} ${output_path}
      """
     }
}
