
def find_index(xam_path) {
    base = "${xam_path}".take("${xam_path}".lastIndexOf('.'))
    isbam = xam_path.name.endsWith(".bam")
    if(isbam) {
        if(file(base + ".bai").exists()){
            return file(base + ".bai")
        }
        // .bam.bai
        base = "${xam_path}"
        if(file(base + ".bai").exists()){
            return file(base + ".bai")
        }
        return "INDEX NOT FOUND"
    }
    if(file(base + ".crai").exists()){
        return file(base + ".crai")
    }
    // .cram.crai
    base = "${xam_path}"
    if(file(base + ".crai").exists()){
        return file(base + ".crai")
    }
    return "CRAM INDEX NOT FOUND"
}


