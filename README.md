# WIP

## Introduction

**obscuromics/variantcalling** is a bioinformatics best-practice analysis pipeline for variant calling using [freebayes](https://github.com/freebayes/freebayes) ...

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

## Pipeline summary

The pipeline takes unaligned paired sample reads (FASTQ files) from a CSV file and the reference file in FASTA format, and then uses freebayes to call variants.

Steps involved:

- Align the reads using bwa mem.
- Locate callable sites using mosdepth.
- Merge bam files.
- Run freebayes.
- Apply filters.

