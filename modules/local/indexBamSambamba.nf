process INDEX_BAM_SAMBAMBA{

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f}.bai")

        script:
        """
        sambamba index -t ${task.cpus} ${bam_f}             
        """
}
