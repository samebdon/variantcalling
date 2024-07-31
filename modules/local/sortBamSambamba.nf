process SORT_BAM_SAMBAMBA {
        publishDir params.outdir, mode:'copy'
        memory '4G'

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.coord_sorted.bam")

        script:
        avail_mem = (task.memory.mega*1).intValue()
        """
        sambamba sort -t ${task.cpus} -m ${avail_mem}MB -o ${bam_f.baseName}.coord_sorted.bam ${bam_f}
        """
}
