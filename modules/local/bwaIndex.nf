process BWA_INDEX {

        input:
        path(genome_f)

        output:
        path("${genome_f}.*")

        script:
        """
        bwa index ${genome_f}
        """
}
