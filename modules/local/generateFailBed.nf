process GENERATE_FAIL_BED {
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${species}.vcf_filter_fails.bed")

        script:
        """
        bcftools view --threads ${task.cpus} -H -i "FILTER!='PASS'" ${vcf_f} | \
        perl -lane '\$pad=0; print(\$F[0]."\\t".(\$F[1]-1)."\\t".((\$F[1]-1)+length(\$F[3]))."\\t".\$F[6])' | \
        bedtools sort | \
        bedtools merge > ${species}.vcf_filter_fails.bed
        """
}
