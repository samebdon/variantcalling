process INTERSECT_BEDS{
        publishDir params.outdir, mode:'copy'

        input:
        path(beds, stageAs: "inputs/*")
        val(species)

        output:
        tuple val(species), path("${species}.callable.all.bed"), emit: all
        tuple val(species), path("${species}.callable.freebayes.bed"), emit: overlap

        script:
        """
        N_FILES="\$(ls inputs/*.bed | wc -l)"
        bedtools multiinter -i $beds | cut -f1-5 > ${species}.callable.all.bed
        cat ${species}.callable.all.bed | awk -v var=\$N_FILES '\$4==var'  | cut -f1-3 > ${species}.callable.freebayes.bed
        """
}
