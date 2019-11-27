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
params.debug = false
params.email = ""
params.reference = "${workflow.projectDir}/WS245/WS245.fa.gz"
parse_conda_software = file("${workflow.projectDir}/scripts/parse_conda_software.awk")

// Debug
if (params.debug.toString() == "true") {
    params.output = "alignment-${date}-debug"
    params.fqs = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.fq_file_prefix = "${workflow.projectDir}/test_data"

} else {
    // The strain sheet that used for 'production' is located in the root of the git repo
    params.output = "alignment-${date}"
    params.fqs = "sample_sheet.tsv"
    params.fq_file_prefix = "";
}

// For now, this pipeline requires NXF_VER 19.09.0
// Prefix this version when running
// e.g.
// NXF_VER=19.09.0-edge nextflow run ...
assert System.getenv("NXF_VER") == "19.09.0-edge"

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
'''

if (params.debug.toString() == "true") {
    out += "\n\n============================= DEBUG =============================\n\n"
}

out += """
    parameters              description                    Set/Default
    ==========              ===========                    ========================
    --debug                 Set to 'true' to test          ${params.debug}
    --fqs                   fastq file (see help)          ${params.fqs}
    --fq_file_prefix        fastq prefix                   ${params.fq_file_prefix}
    --reference             Reference Genome (w/ .gz)      ${params.reference}
    --output                Location for output            ${params.output}
    --email                 Email to be sent results       ${params.email}

    HELP: http://andersenlab.org/dry-guide/pipeline-alignment/
"""
log.info(out)

out
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

include multiqc as multiqc_id from './modules/qc.module.nf' params(output: params.output, grouping: "id")
include multiqc as multiqc_strain from './modules/qc.module.nf' params(output: params.output, grouping: "strain")

process summary {
    
    executor 'local'

    //conda 'fd'
    publishDir "${params.output}", mode: 'copy'
    
    input:
        val(run)

    output:
        path("summary.txt")
        path("software_versions.txt")

    """
        echo '''${log_summary()}''' > summary.txt
        fd "\\.nf\$" ${workflow.projectDir} --exclude 'work' --exec awk -f ${parse_conda_software} > software_versions.txt
    """

}

// process config {

//     """
//         fd "\\.nf$" . --exec awk -f parse_versions.awk > software_versions.txt
//     """
// }

/* 
    Alignment
*/

process alignment {

    tag { row.id }
    
    label 'md'
    conda "bwa=0.7.17 sambamba=0.7.0"

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

    label 'md'
    conda "bwa=0.7.17 sambamba=0.7.0"

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
    publishDir "${params.output}/all", mode: 'copy', pattern: '*.bam*'
    conda 'picard=2.21.3'

    input:
        tuple val(strain), row, path("in.bam"), path("in.bam.bai")
    output:
        tuple row, path("${strain}.bam"), path("${strain}.bam.bai"), emit: "bams"
        path "${strain}.duplicates.txt", emit: "markdups"

    """
        picard MarkDuplicates I=in.bam O=${strain}.bam M=${strain}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
        sambamba index --nthreads=${task.cpus} ${strain}.bam
    """
}

// process symlink_ref_strains {
//     /*
//         This process creates a link for any
//         strain labeled as a 'reference strain'
//         in a new reference folder.

//         So output will be:
//             all/ → All strains
//             reference_strain/ → reference strains

//         The use of symlinks saves on disk space.
//     */
//     tag { "${row.strain}" }

//     executor 'local'
//     when:
//         row.reference_strain == "TRUE"
//     input:
//         tuple row, path("${row.strain}.bam"), path("${row.strain}.bam.bai")

//     script:
//         // check if outdir is abs. path
//         outdir_path = file(params.output).exists() ? "${workflow.projectDir}/${params.outdir}" : params.outdir
//     """
//         mkdir -p ${outdir_path}/reference_strain/
//         if [ ! -L ${outdir_path}/reference_strain/${row.strain}.bam ]; then
//             ln -s ${outdir_path}/all/${row.strain}.bam ${outdir_path}/reference_strain/${row.strain}.bam
//         fi;
//         if [ ! -L ${outdir_path}/reference_strain/${row.strain}.bam.bai ]; then
//             ln -s ${outdir_path}/all/${row.strain}.bam.bai ${outdir_path}/reference_strain/${row.strain}.bam.bai
//         fi;
//     """
// }

// Read sample sheet
sample_sheet = Channel.fromPath(params.fqs, checkIfExists: true)
                      .ifEmpty { exit 1, "sample sheet not found" }
                      .splitCsv(header:true, sep: "\t")

workflow {
    
    // check software
    summary(Channel.from("run"))

    aln_in = sample_sheet.map { row -> row.fq1 = params.fq_file_prefix ? row.fq1 = params.fq_file_prefix + "/" + row.fq1 : row.fq1; row }
                .map { row -> row.fq2 = params.fq_file_prefix ? row.fq2 = params.fq_file_prefix + "/" + row.fq2 : row.fq2; row }
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
        (coverage_strain & idxstats_strain & flagstat_strain & stats_strain)

    mark_dups.out.markdups.concat(coverage_strain.out,
                                    idxstats_strain.out,
                                    flagstat_strain.out,
                                    stats_strain.out).collect() | multiqc_strain                 
}
