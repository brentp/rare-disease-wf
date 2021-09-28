process octopus_population {
    input: tuple(val(region), path(vcfs), path(crams))
           path(ref)
           path(fai)
           val(index) // 1, 2, 3, or 4 just for file naming since we run this
           // process iteratively
    output: tuple(val("${region}"), path("${output_path}"), path(crams))
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
echo octopus -R $ref \
    -p Y=2 chrY=2 -w \$TMPDIR --threads ${task.cpus} --one-based-indexing \
    --disable-denovo-variant-discovery \
    -i $workDir/${reg}.${index}.crams.list \
    --source-candidates-file $workDir/${reg}.${index}.vcfs.list \
    -o ${output_path} > ${output_path}
       """
     } else {
      """
ln -s ${vcfs[0]} ${output_path}
      """
     }
}
