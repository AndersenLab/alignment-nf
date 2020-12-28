#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Alignment Pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
    - Ye Wang <yewangfaith@gmail.com>
*/
nextflow.preview.dsl=2

/* 
    Params
*/

date = new Date().format( 'yyyyMMdd' )
parse_conda_software = file("${workflow.projectDir}/scripts/parse_conda_software.awk")


// Debug
if (params.debug) {
    println """

        *** Using debug mode ***

    """
    params.output = "alignment-${date}-debug"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.fq_prefix = "${workflow.projectDir}/test_data"
    params.species = "ce"

} else {
    // The strain sheet that used for 'production' is located in the root of the git repo
    params.output = "alignment-${date}"
    params.sample_sheet = "${workflow.launchDir}/sample_sheet.tsv"
    params.fq_prefix = ""
}


// Define which genome to map to
// Genome location see config
if (params.species == "ce") {
    params.reference = "${params.reference_ce}"
} else if (params.species == "cb") {
    params.reference = "${params.reference_cb}"
} else if (params.species == "ct") {
    params.reference = "${params.reference_ct}"
} else if (params.species == null) {
    if (params.reference == null) {
        if (params.help) {
        } else { 
        println """

        Please specify a species: ce cb ct with option --species, or a ref genome with --reference"

        """
        exit 1
        }
    }
}



// For now, this pipeline requires NXF_VER 20.01.0
// Prefix this version when running
// e.g.
// NXF_VER=20.01.0 nextflow run ...
//assert System.getenv("NXF_VER") == "20.01.0"

def log_summary() {
/*
    Generates a log
*/

out = '''
                                                            
             ▗▖ ▝▜   ▝                       ▗      ▗▖ ▖▗▄▄▖
             ▐▌  ▐  ▗▄   ▄▄ ▗▗▖ ▗▄▄  ▄▖ ▗▗▖ ▗▟▄     ▐▚ ▌▐
             ▌▐  ▐   ▐  ▐▘▜ ▐▘▐ ▐▐▐ ▐▘▐ ▐▘▐  ▐      ▐▐▖▌▐▄▄▖
             ▙▟  ▐   ▐  ▐ ▐ ▐ ▐ ▐▐▐ ▐▀▀ ▐ ▐  ▐   ▀▘ ▐ ▌▌▐
            ▐  ▌ ▝▄ ▗▟▄ ▝▙▜ ▐ ▐ ▐▐▐ ▝▙▞ ▐ ▐  ▝▄     ▐ ▐▌▐
                         ▖▐
                         ▝▘
''' + """
To run the pipeline:

nextflow main.nf --debug
nextflow main.nf --sample_sheet=name_of_sample_sheet.tsv --species=ce

    parameters              description                                 Set/Default
    ==========              ===========                                 ========================
    --debug                 Set to 'true' to test                       ${params.debug}
    --sample_sheet          sample_sheet (see help)                     ${params.sample_sheet}
    --species               species to map: 'ce', 'cb' or 'ct'          ${params.species}
    --fq_prefix             fastq prefix if not in sample_sheet         ${params.fq_prefix}
    --kmers                 count kmers                                 ${params.kmers}
    --reference             to use in place of default ce cb ct .fa.gz  ${params.reference}
    --output                Location for output                         ${params.output}

    username                                                            ${"whoami".execute().in.text}

"""
out
}


log.info(log_summary())


if (params.help) {
    exit 1
}

// Includes
include coverage as coverage_id from './modules/qc.module.nf' params(params)
include coverage as coverage_strain from './modules/qc.module.nf' params(params)

include idxstats as idxstats_id from './modules/qc.module.nf' params(params)
include idxstats as idxstats_strain from './modules/qc.module.nf' params(params)

include stats as stats_id from './modules/qc.module.nf' params(params)
include stats as stats_strain from './modules/qc.module.nf' params(params)

include flagstat as flagstat_id from './modules/qc.module.nf' params(params)
include flagstat as flagstat_strain from './modules/qc.module.nf' params(params)

include kmer_counting from './modules/qc.module.nf' params(params)
include aggregate_kmer from './modules/qc.module.nf' params(params)

include validatebam as validatebam_strain from './modules/qc.module.nf' params(params)

include multiqc as multiqc_id from './modules/qc.module.nf' params(output: params.output, grouping: "id")
include multiqc as multiqc_strain from './modules/qc.module.nf' params(output: params.output, grouping: "strain")


// Read sample sheet
sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                      .ifEmpty { exit 1, "sample sheet not found" }
                      .splitCsv(header:true, sep: "\t")

workflow {
    
    // check software
    summary(Channel.from("run"))

    aln_in = sample_sheet.map { row -> row.fq1 = params.fq_prefix ? row.fq1 = params.fq_prefix + "/" + row.fq1 : row.fq1; row }
                .map { row -> row.fq2 = params.fq_prefix ? row.fq2 = params.fq_prefix + "/" + row.fq2 : row.fq2; row }
                .map { row -> [row, file(row.fq1), file(row.fq2)]}
    aln_in | (alignment & kmer_counting)
    kmer_counting.out | aggregate_kmer

    /* Merge Bams */
    merge_in = alignment.out.map { row -> [row[0].strain] + row }
                        .groupTuple()
                        .map { strain, row, bam, bai -> [strain, row.sort()[0], bam, bai, row.size()] }
    merge_in | merge_bam | mark_dups

    /* ID Level Stats and multiqc */
    alignment.out.map { row, bam, bai -> ["id", row.id, bam, bai ]} | \
        (coverage_id & idxstats_id & flagstat_id & stats_id)

    coverage_id.out.concat(idxstats_id.out,
                           flagstat_id.out,
                           stats_id.out).collect() | multiqc_id

    /* Strain Level Stats and multiqc */
    mark_dups.out.bams.map { row, bam, bai -> ["strain", row.strain, bam, bai] } | \
        (coverage_strain & idxstats_strain & flagstat_strain & stats_strain & validatebam_strain)

    mark_dups.out.markdups.concat(coverage_strain.out,
                                    idxstats_strain.out,
                                    flagstat_strain.out,
                                    stats_strain.out,
                                    validatebam_strain.out).collect() | multiqc_strain                 

    /* Generate a bam file summary for the next step */
    mark_dups.out.strain_sheet.map { row, bam, bai -> [row.strain, "${params.output}/bam/${row.strain}.bam","${params.output}/bam/${row.strain}.bam.bai"].join("\t") } \
                 .collectFile(name: 'strain_summary.tsv',
                              newLine: true,
                              storeDir: "${params.output}")

}


process summary {
    
    executor 'local'

    // conda 'fd-find'
    publishDir "${params.output}", mode: 'copy'
    
    input:
        val(run)

    output:
        path("sample_sheet.tsv")
        path("summary.txt")
        path("software_versions.txt")

    """
        echo '''${log_summary()}''' > summary.txt
        fd "\\.nf\$" ${workflow.projectDir} --exclude 'work' --exec awk -f ${parse_conda_software} > software_versions.txt
        cat ${params.sample_sheet} > sample_sheet.tsv
    """

}

/* 
    Alignment
*/

process alignment {

    tag { row.id }
    
    label 'md'
    // container "andersenlab/alignment"

    input:
        tuple row, path(fq1), path(fq2)
        
    output:
        tuple row, file("${row.id}.bam"), file("${row.id}.bam.bai")

	script:
		// Construct read group
		RG = ["@RG",
			  "ID:${row.id}",
			  "SM:${row.strain}",
			  "LB:${row.lb}",
			  "PL:illumina"].join("\\t")

    """
        bwa mem -t ${task.cpus} -R '${RG}' ${params.reference} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${task.cpus} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${task.cpus} --show-progress --tmpdir=. --out=${row.id}.bam /dev/stdin
        sambamba index --nthreads=${task.cpus} ${row.id}.bam
        if [[ ! \$(samtools view ${row.id}.bam | head -n 10) ]]; then
            exit 1;
        fi
    """
}


/* 
    Merge ID Bams by Strain
*/
process merge_bam {

    tag { row.strain }

    label 'lg'
    // container "andersenlab/alignment"

    input:
        tuple strain, row, path(bam), path(bai), val(n_count)

    output:
        tuple strain, row, file("${row.strain}.bam"), file("${row.strain}.bam.bai")

    script:
        if (n_count == 1)
            """
                mv ${bam} ${strain}.bam
                mv ${bai} ${strain}.bam.bai
            """
        else
            """
                sambamba merge --nthreads=${task.cpus} --show-progress ${strain}.bam ${bam}
                sambamba index --nthreads=${task.cpus} ${strain}.bam
            """
}

process mark_dups {

    tag { "${strain}" }

    label 'lg'
    publishDir "${params.output}/bam", mode: 'copy', pattern: '*.bam*'
    // container "andersenlab/alignment"

    input:
        tuple val(strain), row, path("${strain}.in.bam"), path("${strain}.in.bam.bai")
    output:
        tuple row, path("${strain}.bam"), path("${strain}.bam.bai"), emit: "bams"
        tuple row, path("${strain}.bam"), path("${strain}.bam.bai"), emit: "strain_sheet"
        path "${strain}.duplicates.txt", emit: "markdups"

    """
        picard -Xmx${task.memory.toGiga()}g -Xms1g MarkDuplicates I=${strain}.in.bam \\
                              O=${strain}.bam \\
                              M=${strain}.duplicates.txt \\
                              VALIDATION_STRINGENCY=SILENT \\
                              REMOVE_DUPLICATES=false \\
                              TAGGING_POLICY=All \\
                              REMOVE_SEQUENCING_DUPLICATES=TRUE \\
                              SORTING_COLLECTION_SIZE_RATIO=0.1

        sambamba index --nthreads=${task.cpus} ${strain}.bam
    """
}
