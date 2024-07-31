process MERGE_SAMBAMBA {
        memory '8G'

        input:
        path(bams)
        val(species)

        output:
        tuple val(species), path("${species}.bam"), path("${species}.bam.bai")

        script:
        """
        sambamba merge -t ${task.cpus} ${species}.bam ${bams.join(" ")}
        sambamba index -t ${task.cpus} ${species}.bam         
        """
}
