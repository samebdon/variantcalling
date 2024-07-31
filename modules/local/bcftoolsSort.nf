process BCFTOOLS_SORT {
        cpus 1
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${species}.hard_filtered.sorted.vcf.gz")
        
        script:
        """
        bcftools sort -Oz ${vcf_f} > ${species}.hard_filtered.sorted.vcf.gz
        """
}     
