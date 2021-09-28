
process split {
    container = 'docker://brentp/rare-disease:v0.0.3'

    input: tuple(path(gvcf), path(tbi))
           path(fai)
    output: file("*.split.gvcf.gz")

    script:
    """
    # get large chroms and chrM in one file each
    for chrom in \$(awk '\$2 > 40000000 || \$1 ~/(M|MT)\$/' $fai | cut -f 1); do
        bcftools view $gvcf --threads 3 -O z -o \$(basename $gvcf .gvcf.gz).\${chrom}.split.gvcf.gz \$chrom
    done
    # small HLA and gl chroms all to go single file
    awk '!(\$2 > 40000000 || \$1 ~/(M|MT)\$/) { print \$1"\\t0\\t"\$2+1 }' $fai > other_chroms
    # if it's non-empty then create the extras / other split file
    if [ -s other_chroms ]; then
        bcftools view $gvcf --threads 3 -R other_chroms -O z -o \$(basename $gvcf .gvcf.gz).OTHER.split.gvcf.gz
    fi
    """

}

process split_by_size {
    container = 'docker://brentp/rare-disease:v0.2.1'

    input: path(fai)
           val(chunk_size)
    output: stdout

    script:
    """
tiwih fairegions $fai --region_size $chunk_size
    """
}
