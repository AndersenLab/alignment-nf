

/* 
    Coverage
*/
process coverage {

    tag {"${grouping} -> ${name}" }

//    publishDir "${params.output}/coverage/${grouping}", mode: 'copy'

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 4 * task.attempt }
    memory = { 16.GB * task.attempt }

    input:
        tuple val(grouping), val(name), path("in.bam"), file("in.bam.bai")

    output:
        tuple path("${name}.mosdepth.summary.txt"), \
              path("${name}.mosdepth.global.dist.txt"), \
              path("${name}.per-base.bed.gz"), \
              path("${name}.per-base.bed.gz.csi")

    """
        export MOSDEPTH_PRECISION=5
        mosdepth --threads ${task.cpus} \\
                 ${name} \\
                 in.bam
    """
}

/* 
    samtools idx stats
*/

process idxstats {
    
    tag {"${grouping} -> ${name}" }

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
        tuple val(grouping), val(name), path("in.bam"), file("in.bam.bai")
    output:
        tuple val(grouping), path("${name}.idxstats")

    """
        samtools idxstats in.bam > ${name}.idxstats
    """
}

process stats {
    
    tag {"${grouping} -> ${name}" }

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
        tuple val(grouping), val(name), path("in.bam"), file("in.bam.bai")
    output:
        tuple val(grouping), path("${name}.stats")

    """
        samtools stats in.bam > ${name}.stats
    """
}


process flagstat {
    
    tag {"${grouping} -> ${name}" }

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
        tuple val(grouping), val(name), path("in.bam"), file("in.bam.bai")
    output:
        tuple val(grouping), path("${name}.flagstat")

    """
        samtools flagstat in.bam > ${name}.flagstat
    """

}

/*
    Kmers
*/

process kmer_counting {

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    tag { "${data.strain}" }
    when params.kmers.toString() == "true"

    input:
        tuple val(data), path(genome_path), val(genome_basename), path(fq1), path(fq2)

    output:
        tuple val(data.strain), val(data), file("${data.id}.kmer.tsv")

    """
        # fqs will have same number of lines
        export OFS="\t"
        fq_wc=`zcat fq1 | awk 'NR % 4 == 0' | wc -l`
        zcat fq1 fq2 | \\
        fastq-kmers -k 6 | \\
        awk -v OFS="\t" -v ID=${data.id} -v SM=${data.strain} -v fq_wc="\${fq_wc}" 'NR > 1 { print \$0, SM, ID, fq_wc }' - > ${data.id}.kmer.tsv
    """

}


process aggregate_kmer {

    tag "aggregate-kmer"

    executor 'local'
    container null

    publishDir "${params.output}/_aggregate", mode: 'copy'
    when params.kmers.toString() == "true"

    input:
        file("kmer*.tsv")
    output:
        file("kmers.tsv")

    """
        cat <(echo "kmer\tfrequency\tstrain\tid\treads ") *.tsv > kmers.tsv
    """

}


process validatebam {

    tag {"${grouping} -> ${name}" }

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
        tuple val(grouping), val(name), path("in.bam"), file("in.bam.bai")
    output:
        tuple val(grouping), path("${name}.validatesamfile.txt")

    """
        picard ValidateSamFile I=in.bam MODE=SUMMARY > ${name}.validatesamfile.txt
    """
}


/* MULTI-QC */
process multiqc {

    label 'multiqc'

    publishDir "${params.output}/_aggregate/multiqc", mode: 'copy'

    errorStrategy 'ignore'


    input:
        tuple val(grouping), file("*")

    output:
        tuple path("${grouping}_multiqc_report.html"), path("${grouping}_data/*")

    """
        multiqc . --data-format tsv \\
                  --config multiqc_config.yaml \\
                  --title ${grouping} \\
                  --flat \\
                  --profile-runtime
        # mv data folder to reduce size
        mv ${grouping}_multiqc_report_data/ ${grouping}_data/
    """

}
