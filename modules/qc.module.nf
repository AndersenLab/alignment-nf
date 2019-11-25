

/* 
    Coverage
*/
process coverage {

    tag {"${grouping} -> ${name}" }
    publishDir "${params.bamdir}/WI/strain/coverage/${grouping}", mode: 'copy'
    conda { System.properties['os.name'] != "Mac OS X" ? 'bioconda::mosdepth=0.2.6' : "" }
    label 'mid'

    input:
        tuple grouping, name, path("in.bam"), file("in.bam.bai")

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

    input:
        tuple grouping, name, path("in.bam"), file("in.bam.bai")
    output:
        tuple grouping, path("${name}.idxstats")

    """
        samtools idxstats in.bam > ${name}.idxstats
    """
}

process stats {
    
    tag {"${grouping} -> ${name}" }

    input:
        tuple grouping, name, path("in.bam"), file("in.bam.bai")
    output:
        tuple grouping, path("${name}.stats")

    """
        samtools stats in.bam > ${name}.stats
    """
}


process flagstat {
    
    tag {"${grouping} -> ${name}" }

    input:
        tuple grouping, name, path("in.bam"), file("in.bam.bai")
    output:
        tuple grouping, path("${name}.flagstat")

    """
        samtools flagstat in.bam > ${name}.flagstat
    """
}

/*
    Kmers
*/

process kmer_counting {

    conda 'fastq-tools=0.8'
    tag { "${row.strain}" }

    input:
        tuple row, file("fq1.fq.gz"), file("fq2.fq.gz")
    output:
        tuple val(row.strain), row, file("${row.id}.kmer.tsv")

    """
        # fqs will have same number of lines
        export OFS="\t"
        fq_wc=`zcat fq1.fq.gz | awk 'NR % 4 == 0' | wc -l`
        zcat fq1.fq.gz fq2.fq.gz | \\
        fastq-kmers -k 6 | \\
        awk -v OFS="\t" -v ID=${row.id} -v SM=${row.strain} -v fq_wc="\${fq_wc}" 'NR > 1 { print \$0, SM, ID, fq_wc }' - > ${row.id}.kmer.tsv
    """
}


process aggregate_kmer {

    publishDir "${params.bamdir}/WI/strain/_aggregate", mode: 'copy'

    input:
        file("kmer*.tsv")
    output:
        file("kmers.tsv")

    """
        cat <(echo "kmer\tfrequency\tstrain\tid\treads ") *.tsv > kmers.tsv
    """

}


/* Aggregation */
process multiqc {

    tag { "multiqc" }

    publishDir "${params.bamdir}/WI/strain/_aggregate/multiqc", mode: 'copy'
    conda 'multiqc=1.8'

    input:
        file("*")

    output:
        path("${params.grouping}_multiqc_report.html")
        path("${params.grouping}_multiqc_report_data/*")

    // --config ${params.multiqc_config}
    """
        multiqc . --data-format tsv \\
                  --title ${params.grouping} \\
                  --flat
    """

}