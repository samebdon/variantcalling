process BEDTOOLS_SUBTRACT {
        cpus 1
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(a_bed)
        tuple val(species), path(b_bed)

        output:
        tuple val(species), path("${species}.callable.bed")

        script:
        """
        bedtools subtract -a ${a_bed} -b ${b_bed} > ${species}.callable.bed
        """
}
