process BCFTOOLS_INDEX {
        cpus 1
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(vcf_f)

        output:
        tuple val(meta), path("${vcf_f}.csi")

        script:
        """
        bcftools index -c ${vcf_f} -o ${vcf_f}.csi
        """
}

