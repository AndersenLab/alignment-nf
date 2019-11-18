nextflow.preview.dsl=2 
/*
 * Authors:
 * - Daniel Cook <danielecook@gmail.com>
 * - Ye Wang <yewangfaith@gmail.com>
 *
 */

/* 
    ======
    Params
    ======
*/

date = new Date().format( 'yyyyMMdd' )
params.out = "Alignment-${date}"
params.debug = false
params.tmpdir = "tmp/"
params.email = ""
params.reference = "(required)"
params.env_container = ""

// Debug
if (params.debug == true) {
    log.info '''
        *** Using debug mode ***
    '''
    params.bamdir = "${params.out}/bam"
    params.fqs = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.fq_file_prefix = "${workflow.projectDir}/test_data"

} else {
    // The SM sheet that is used is located in the root of the git repo
    params.bamdir = "(required)"
    params.fqs = "sample_sheet.tsv"
    params.fq_file_prefix = "";
}

log.info '''

             ▗▖ ▝▜   ▝                       ▗      ▗▖ ▖▗▄▄▖
             ▐▌  ▐  ▗▄   ▄▄ ▗▗▖ ▗▄▄  ▄▖ ▗▗▖ ▗▟▄     ▐▚ ▌▐
             ▌▐  ▐   ▐  ▐▘▜ ▐▘▐ ▐▐▐ ▐▘▐ ▐▘▐  ▐      ▐▐▖▌▐▄▄▖
             ▙▟  ▐   ▐  ▐ ▐ ▐ ▐ ▐▐▐ ▐▀▀ ▐ ▐  ▐   ▀▘ ▐ ▌▌▐
            ▐  ▌ ▝▄ ▗▟▄ ▝▙▜ ▐ ▐ ▐▐▐ ▝▙▞ ▐ ▐  ▝▄     ▐ ▐▌▐
                         ▖▐
                         ▝▘
'''

log.info """
    parameters              description                    Set/Default
    ==========              ===========                    =======
    --debug                 Set to 'true' to test          ${params.debug}
    --fqs                   fastq file (see help)          ${params.fqs}
    --fq_file_prefix        fastq prefix                   ${params.fq_file_prefix}
    --reference             Reference Genome (w/ .gz)      ${params.reference}
    --bamdir                Location for bams              ${params.bamdir}
    --tmpdir                A temporary directory          ${params.tmpdir}
    --email                 Email to be sent results       ${params.email}

    HELP: http://andersenlab.org/dry-guide/pipeline-alignment/
"""


// fqs.into {
//     fqs_kmer
//     fqs_align
// }

/* 
    =========
    Alignment
    =========
*/

process alignment {

    label 'high_cpu'
    tag { row.id }
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
  Merge - Merge ID Bams by Strain
*/

process merge_bam {

    tag { row.strain }
    conda "bwa=0.7.17 sambamba=0.7.0"
    label 'high_cpu'

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

    tag { row.strain }
    label 'high_cpu'
    publishDir "${params.bamdir}/WI/strain/all", mode: 'copy', pattern: '*.bam*'
    conda 'picard=2.21.3'

    input:
        tuple strain, row, file("in.bam"), file("in.bam.bai")
    output:
        tuple strain, row, file("${row.strain}.bam"), file("${row.strain}.bam.bai"), emit: "bams"
        path "${strain}.duplicates.txt", emit: "markdups"

    """
    picard MarkDuplicates I=in.bam O=${strain}.bam M=${strain}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${task.cpus} ${strain}.bam
    """
}

process symlink_ref_strains {
    /*
        This process creates a link for any
        strain labeled as a 'reference strain'
        in a new reference folder.

        So output will be:
            all/ → All strains
            reference_strain/ → reference strains

        The use of symlinks saves on disk space.
    */
    tag { row.strain }
    executor 'local'
    when:
        row.reference_strain == "TRUE"
    input:
        tuple strain, row, file("${row.strain}.bam"), file("${row.strain}.bam.bai")

    script:
        // check if bamdir is abs. path
        bamdir_path = file(params.bamdir).exists() ? "${workflow.projectDir}/${params.bamdir}" : params.bamdir
    """
        mkdir -p ${bamdir_path}/WI/strain/reference_strain/
        if [ ! -L ${bamdir_path}/WI/strain/reference_strain/${strain}.bam ]; then
            ln -s ${bamdir_path}/WI/strain/all/${strain}.bam ${bamdir_path}/WI/strain/reference_strain/${strain}.bam
        fi;
        if [ ! -L ${bamdir_path}/WI/strain/reference_strain/${strain}.bam.bai ]; then
            ln -s ${bamdir_path}/WI/strain/all/${strain}.bam.bai ${bamdir_path}/WI/strain/reference_strain/${strain}.bam.bai
        fi;
    """
}

// process kmer_counting {

//     cpus params.cores

//     tag { ID }

//     input:
//         set SM, ID, LB, fq1, fq2, seq_folder from fqs_kmer
//     output:
//         file("${ID}.kmer.tsv") into kmer_set
//     when:
//         params.goal == "isotype"

//     """
//         # fqs will have same number of lines
//         export OFS="\t"
//         fq_wc="`zcat ${fq1} | awk 'NR % 4 == 0' | wc -l`"
        
//         zcat ${fq1} ${fq2} | \\
//         fastq-kmers -k 6 | \\
//         awk -v OFS="\t" -v ID=${ID} -v SM=${SM} -v fq_wc="\${fq_wc}" 'NR > 1 { print \$0, SM, ID, fq_wc }' - > ${ID}.kmer.tsv
//     """
// }


// process merge_kmer {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/phenotype", mode: "copy"

//     input:
//         file("kmer*.tsv") from kmer_set.collect()
//     output:
//         file("kmers.tsv")
//     when:
//         params.goal == "isotype"

//     """
//         cat <(echo "kmer\tfrequency\tSM\tID\twc") *.tsv > kmers.tsv
//     """

// }

// fq_bam_set.into { fq_cov_bam; fq_stats_bam; fq_idx_stats_bam; bam_strain_stats}

// bam_set.into { bam_set1; bam_set2 }

// bam_set1.into { 
//                merged_bams_for_coverage;
//                merged_bams_individual;
//                merged_bams_union;
//                bams_idxstats;
//                bams_stats;
//                fq_concordance_bam
//              }

// bam_set2.into {
//                   bam_idxstats;
//                   bam_stats;
//                   bam_coverage;
//                   bam_telseq;
//                   bam_isotype_stats;
//     }

/*
=======================================================
               Strain bams information
=======================================================
*/

/* 
    ========
    Coverage
    ========
*/
// process coverage_fq {

//     tag { ID }

//     label 'bam_coverage'

//     input:
//         set val(ID), file("${ID}.bam"), file("${ID}.bam.bai") from fq_cov_bam
//     output:
//         file("${ID}.coverage.tsv") into fq_coverage
//     when:
//         params.goal == "strain"


//     """
//         bam coverage ${ID}.bam > ${ID}.coverage.tsv
//     """
// }


// process coverage_fq_merge {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/fq", mode: 'copy', overwrite: true

//     input:
//         file fq_set from fq_coverage.toSortedList()

//     output:
//         file("fq_coverage.full.tsv")
//         file("fq_coverage.tsv")
//     when:
//         params.goal == "strain"

//     """
//         echo -e 'fq\\tcontig\\tstart\\tend\\tproperty\\tvalue' > fq_coverage.full.tsv
//         cat ${fq_set} >> fq_coverage.full.tsv
//         cat <(echo -e 'fq\\tcoverage') <( cat fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > fq_coverage.tsv
//     """
// }

/* 
    ==============
    fq index stats
    ==============
*/

// process fq_idx_stats {
    
//     tag { ID }

//     input:
//         set val(ID), file("${ID}.bam"), file("${ID}.bam.bai") from fq_idx_stats_bam
//     output:
//         file fq_idxstats into fq_idxstats_set
//         file fq_idxstats into strain_bam_idxstats_multiqc
//     when:
//         params.goal == "strain"

//     """
//         samtools idxstats ${ID}.bam | awk '{ print "${ID}\\t" \$0 }' > fq_idxstats
//     """
// }

// process fq_combine_idx_stats {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/fq", mode: 'copy', overwrite: true

//     input:
//         file("?.stat.txt") from fq_idxstats_set.toSortedList()

//     output:
//         file("fq_bam_idxstats.tsv")
        
//     when:
//         params.goal == "strain"

//     """
//         echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > fq_bam_idxstats.tsv
//         cat *.stat.txt >> fq_bam_idxstats.tsv
//     """

// }

/* 
    ============
    fq bam stats
    ============
*/
// process fq_bam_stats {

//     tag { ID }

//     input:
//         set val(ID), file("${ID}.bam"), file("${ID}.bam.bai") from fq_stats_bam

//     output:
//         file 'bam_stat' into fq_bam_stat_files
//     when:
//         params.goal == "strain"

//     """
//         cat <(samtools stats ${ID}.bam | grep ^SN | cut -f 2- | awk '{ print "${ID}\t" \$0 }' | sed 's/://g') > bam_stat
//     """
// }

// process combine_fq_bam_stats {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/fq", mode: 'copy', overwrite: true

//     input:
//         file("*.stat.txt") from fq_bam_stat_files.toSortedList()

//     output:
//         file("fq_bam_stats.tsv")
//     when:
//         params.goal == "strain"

//     """
//         echo -e "ID\\tvariable\\tvalue\\tcomment" > fq_bam_stats.tsv
//         cat *.stat.txt >> fq_bam_stats.tsv
//     """
// }

/*
    SM_idx_stats
*/

// process SM_idx_stats {
    
//     tag { SM }

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bams_idxstats
//     output:
//         file('bam_idxstats.txt') into bam_idxstats_set
//     when:
//         params.goal == "strain"

//     """
//         samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > bam_idxstats.txt
//     """
// }

// process isotype_SM_combine_idx_stats {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/SM", mode: 'copy', overwrite: true

//     input:
//         file("*.stat.txt") from bam_idxstats_set.toSortedList()

//     output:
//         file("SM_bam_idxstats.tsv")
//     when:
//         params.goal == "strain"

//     """
//         echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
//         cat *.stat.txt | sort >> SM_bam_idxstats.tsv
//     """

// }

/*
    SM bam stats
*/

// process SM_bam_stats {

//     tag { SM }

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bams_stats

//     output:
//         file('bam_stat.txt') into SM_bam_stat_files
//     when:
//         params.goal == "strain"

//     """
//         cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat.txt
//     """
// }

// process combine_SM_bam_stats {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/SM", mode: 'copy', overwrite: true

//     input:
//         file("?.stat.txt") from SM_bam_stat_files.toSortedList()

//     output:
//         file("SM_bam_stats.tsv")
//     when:
//         params.goal == "strain"

//     """
//         echo -e "ID\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
//         cat *.stat.txt | sort >> SM_bam_stats.tsv
//     """
// }

/*
     Publish duplicates
*/

// process format_duplicates {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/SM", mode: 'copy', overwrite: true

//     input:
//         val duplicates_set from duplicates_file.toSortedList()

//     output:
//         file("bam_duplicates.tsv")
//     when:
//         params.goal == "strain"


//     """
//         echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > bam_duplicates.tsv
//         for i in ${duplicates_set.join(" ")}; do
//             f=\$(basename \${i})
//             cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> bam_duplicates.tsv
//         done;
//     """
// }

/*
    Coverage Bam
*/

// process coverage_SM {

//     tag { SM }

//     label 'bam_coverage'

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_for_coverage

//     output:
//         file("${SM}.coverage.tsv") into SM_coverage
//         file("${SM}.1mb.coverage.tsv") into SM_1mb_coverage
//         file("${SM}.100kb.coverage.tsv") into SM_100kb_coverage
//     when:
//         params.goal == "strain"


//     """
//         bam coverage ${SM}.bam > ${SM}.coverage.tsv
//         bam coverage --window=1000000 ${SM}.bam > ${SM}.1mb.coverage.tsv
//         bam coverage --window=100000 ${SM}.bam > ${SM}.100kb.coverage.tsv
//     """
// }


// process coverage_SM_merge {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/SM", mode: 'copy', overwrite: true

//     input:
//          file(sm_set) from SM_coverage.toSortedList()

//     output:
//         file("SM_coverage.full.tsv")
//         file("SM_coverage.tsv") into SM_coverage_merged
//     when:
//         params.goal == "strain"

//     """
//         echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
//         cat ${sm_set.join(" ")} >> SM_coverage.full.tsv
//         # Generate condensed version
//         cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6 | sort) > SM_coverage.tsv
//     """
// }

// process coverage_bins_merge {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/SM", mode: 'copy', overwrite: true

//     input:
//         file(mb) from SM_1mb_coverage.toSortedList()
//         val kb_100 from SM_100kb_coverage.toSortedList()

//     output:
//         file("SM_coverage.mb.tsv.gz")
//     when:
//         params.goal == "strain"

//     """
//         echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.mb.tsv
//         cat ${mb.join(" ")} >> SM_coverage.mb.tsv
//         gzip SM_coverage.mb.tsv
//     """
// }

// process bam_strain_stats {

//     cpus params.cores

//     tag { SM }

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_strain_stats
//         file("${params.genome}/*") from reference2.collect()

//     output:
//         file("${SM}.samtools.txt") into strain_SM_samtools_stats_set
//         file("${SM}.bamtools.txt") into strain_SM_bamtools_stats_set
//         file("${SM}_fastqc.zip") into strain_SM_fastqc_stats_set
//         file("${SM}.picard.*") into strain_SM_picard_stats_set
//     when:
//          params.goal == "strain"

//     """
//         samtools stats --threads=${task.cpus} ${SM}.bam > ${SM}.samtools.txt
//         bamtools -in ${SM}.bam > ${SM}.bamtools.txt
//         fastqc --threads ${task.cpus} ${SM}.bam
//         picard CollectAlignmentSummaryMetrics R=${params.genome}/${params.genome}.fa.gz I=${SM}.bam O=${SM}.picard.alignment_metrics.txt
//         picard CollectInsertSizeMetrics I=${SM}.bam O=${SM}.picard.insert_metrics.txt H=${SM}.picard.insert_histogram.txt
//     """

// }



/*
=======================================================
               Isotype bams information
=======================================================
*/

// process bam_isotype_stats {

//     cpus params.cores

//     tag { SM }

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_isotype_stats
//         file("${params.genome}/*") from reference3.collect()
//     output:
//         file("${SM}.samtools.txt") into SM_samtools_stats_set
//         file("${SM}.bamtools.txt") into SM_bamtools_stats_set
//         file("${SM}_fastqc.zip") into SM_fastqc_stats_set
//         file("${SM}.picard.*") into SM_picard_stats_set
//     when:
//          params.goal == "isotype"

//     """
//         samtools stats --threads=${task.cpus} ${SM}.bam > ${SM}.samtools.txt
//         bamtools -in ${SM}.bam > ${SM}.bamtools.txt
//         fastqc --threads ${task.cpus} ${SM}.bam
//         picard CollectAlignmentSummaryMetrics R=${params.genome}/${params.genome}.fa.gz I=${SM}.bam O=${SM}.picard.alignment_metrics.txt
//         picard CollectInsertSizeMetrics I=${SM}.bam O=${SM}.picard.insert_metrics.txt H=${SM}.picard.insert_histogram.txt
//     """

// }

/*
    =================
        inx_stats
    =================
*/

// process isotype_SM_idx_stats {

//     tag { SM }

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_idxstats
//     output:
//         file("${SM}.bam_idxstats") into isotype_bam_idxstats_set
//         file("${SM}.bam_idxstats") into bam_idxstats_multiqc
//     when:
//         params.goal == "isotype"

//     """
//         samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > ${SM}.bam_idxstats
//     """
// }

// process SM_combine_idx_stats {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/SM", mode: "copy"

//     input:
//         file(bam_idxstats) from isotype_bam_idxstats_set.toSortedList()

//     output:
//         file("isotype_bam_idxstats.tsv")
//     when:
//         params.goal == "isotype"

//     """
//         echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > isotype_bam_idxstats.tsv
//         cat ${bam_idxstats.join(" ")} >> isotype_bam_idxstats.tsv
//     """
// }

/*
    =================
    Isotype BAM stats
    =================
*/

// process isotype_bam_stats {

//     tag { SM }

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_stats

//     output:
//         file 'bam_stat' into isotype_SM_bam_stat_files
//     when:
//         params.goal == "isotype"

//     """
//         samtools stats ${SM}.bam | \\
//         grep ^SN | \\
//         cut -f 2- | \\
//         awk '{ print "${SM}\t" \$0 }' | \\
//         sed 's/://g' > bam_stat
//     """
// }

// process combine_isotype_bam_stats {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/SM", mode: "copy"

//     input:
//         file(stat_files) from isotype_SM_bam_stat_files.collect()

//     output:
//         file("isotype_bam_stats.tsv")
//     when:
//         params.goal == "isotype"

//     """
//         echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > isotype_bam_stats.tsv
//         cat ${stat_files.join(" ")} >> SM_stats.tsv
//     """
// }

/*
    ============
    Coverage BAM
    ============
*/
// process isotype_coverage_SM {

//     tag { SM }

//     label 'bam_coverage'

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_coverage

//     output:
//         val SM into isotype_coverage_sample
//         file("${SM}.coverage.tsv") into isotype_coverage
//     when:
//         params.goal == "isotype"


//     """
//         bam coverage ${SM}.bam > ${SM}.coverage.tsv
//     """
// }

// process isotype_coverage_SM_merge {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/SM", mode: 'copy'

//     input:
//         file(sm_set) from isotype_coverage.toSortedList()

//     output:
//         file("isotype_coverage.full.tsv") into mt_content
//         file("isotype_coverage.tsv") into isotype_coverage_merged
//     when:
//         params.goal == "isotype"

//     """
//         echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > isotype_coverage.full.tsv
//         cat ${sm_set.join(" ")} >> isotype_coverage.full.tsv
//         # Generate condensed version
//         cat <(echo -e 'strain\\tcoverage') <(cat isotype_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > isotype_coverage.tsv
//     """
// }

/*
    ==========
    MT content
    ==========
*/

// process output_mt_content {


//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/phenotype", mode: "copy"

//     input:
//         file("isotype_coverage.full.tsv") from mt_content

//     output:
//         file("MT_content.tsv")
//     when:
//         params.goal == "isotype"

//     """
//         cat <(echo -e 'isotype\\tmt_content') <(cat isotype_coverage.full.tsv | awk '/mt_nuclear_ratio/' | cut -f 1,6) > MT_content.tsv
//     """
// }

/*
    ======
    telseq
    ======
*/

// process call_telseq {

//     tag { SM }

//     input:
//         set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_telseq
//     output:
//         file("telseq_out.txt") into telseq_results
//     when:
//         params.goal == "isotype"

//     """
//         telseq -z TTAGGC -H ${SM}.bam > telseq_out.txt
//     """
// }

// process combine_telseq {


//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/phenotype", mode: 'copy'

//     input:
//         file("ind_telseq?.txt") from telseq_results.toSortedList()

//     output:
//         file("telseq.tsv")
//     when:
//         params.goal == "isotype"

//     '''
//         telseq -h > telseq.tsv
//         cat ind_telseq*.txt | egrep -v '\\[|BAMs' >> telseq.tsv
//     '''
// }

// process isotype_multiqc_report {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/report", mode: 'copy'

//     input:
//         file(samtools_stats) from SM_samtools_stats_set.toSortedList()
//         file(fastqc) from SM_fastqc_stats_set.toSortedList()
//         file("picard*.stats.txt") from SM_picard_stats_set.collect()
//         file(bamtools_stats) from SM_bamtools_stats_set.toSortedList()
//         file("bam*.idxstats") from bam_idxstats_multiqc.toSortedList()

//     output:
//         file("multiqc_data/*.json") into multiqc_json_files_isotype
//         file("multiqc.html")

//     when:
//         params.goal == "isotype"

//     """
//         multiqc -k json --filename multiqc.html .
//     """
// }


// process strain_multiqc_report {

//     publishDir "${params.bamdir}/WI/${params.goal}/${params.out}/report", mode: 'copy'

//     input:
//         file(samtools_stats) from strain_SM_samtools_stats_set.toSortedList()
//         file(fastqc) from strain_SM_fastqc_stats_set.toSortedList()
//         file("picard*.stats.txt") from strain_SM_picard_stats_set.collect()
//         file(bamtools_stats) from strain_SM_bamtools_stats_set.toSortedList()
//         file("bam*.idxstats") from strain_bam_idxstats_multiqc.toSortedList()

//     output:
//         file("multiqc_data/*.json") into multiqc_json_files_strain optional true
//         file("multiqc.html")

//     when:
//         params.goal == "strain"

//     """
//         multiqc -k json --filename multiqc.html .
//     """
// }


// // Pipeline work summary, there would be a email notification if email has been given.
// workflow.onComplete {

//     summary = """   
//     Pipeline execution summary
//     ---------------------------
//     Completed at: ${workflow.complete}
//     Duration    : ${workflow.duration}
//     Success     : ${workflow.success}
//     workDir     : ${workflow.workDir}
//     exit status : ${workflow.exitStatus}
//     Error report: ${workflow.errorReport ?: '-'}
//     Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
//     """

//     println summary

//     // mail summary
//     if (params.email) {
//         ['mail', '-s', 'Alignment-nf', params.email].execute() << summary
//     }


// }

// Read sample sheet
sample_sheet = Channel.fromPath(params.fqs, checkIfExists: true)
                      .ifEmpty { exit 1, "sample sheet not found" }
                      .splitCsv(header:true, sep: "\t")



workflow {
    
    aln_in = sample_sheet.map { row -> row.fq1 = params.fq_file_prefix ? row.fq1 = params.fq_file_prefix + "/" + row.fq1 : row.fq1; row }
                .map { row -> row.fq2 = params.fq_file_prefix ? row.fq2 = params.fq_file_prefix + "/" + row.fq2 : row.fq2; row }
                .map { row -> [row, file(row.fq1), file(row.fq2)]}
    aln_in | alignment
    
    merge_in = alignment.out.map { row -> [row[0].strain] + row }
                        .groupTuple()
                        .map { strain, row, bam, bai -> [strain, row.sort()[0], bam, bai, row.size()] }
    merge_in | merge_bam | mark_dups
    mark_dups.out.bams | symlink_ref_strains
                 
}

//tuple row, file("${row.id}.bam"), file("${row.id}.bam.bai")