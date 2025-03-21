#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Alignment Pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
    - Ye Wang <yewangfaith@gmail.com>
*/

// Includes

include { coverage as coverage_id } from './modules/qc.module'
include { coverage as coverage_strain } from './modules/qc.module'
include { idxstats as idxstats_id } from './modules/qc.module'
include { idxstats as idxstats_strain } from './modules/qc.module'
include { stats as stats_id } from './modules/qc.module'
include { stats as stats_strain } from './modules/qc.module'
include { flagstat as flagstat_id } from './modules/qc.module'
include { flagstat as flagstat_strain } from './modules/qc.module'
include { kmer_counting } from './modules/qc.module'
include { aggregate_kmer } from './modules/qc.module'
include { validatebam as validatebam_strain } from './modules/qc.module'
include { multiqc as multiqc_id } from './modules/qc.module'
include { multiqc as multiqc_strain } from './modules/qc.module'


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
} else {
    // The strain sheet that used for 'production' is located in the root of the git repo
    params.output = "alignment-${date}"
    params.sample_sheet = "${workflow.launchDir}/sample_sheet.tsv"
    params.fq_prefix = "${params.data_path}/${params.species}/WI/fastq/dna/"
}
params.ncbi = "${params.data_path}/other/ncbi_blast_db/"

// set project and build defaults for CE, CB, and CT, can always change with argument.
if(params.species == "c_elegans") {
    params.project="PRJNA13758"
    params.ws_build="WS283"
} else if(params.species == "c_briggsae") {
    params.project="QX1410_nanopore"
    params.ws_build="Feb2020"
} else if(params.species == "c_tropicalis") {
    params.project="NIC58_nanopore"
    params.ws_build="June2021"
}

// Define the genome
if(params.species == "c_elegans" | params.species == "c_briggsae" | params.species == "c_tropicalis") {
    params.reference = "${params.data_path}/${params.species}/genomes/${params.project}/${params.ws_build}/${params.species}.${params.project}.${params.ws_build}.genome.fa.gz"
    if (params.bam_dir != null) {
        bam_dir = params.bam_dir
    } else  {
        bam_dir = "${params.data_path}/${params.species}/WI/alignments"
    }
} else if (params.species == null) {
    if (params.reference == null) {
        if (params.help) {
        } else { 
        println """

        Please specify a species: c_elegans c_brigssae c_tropicalis with option --species, or a ref genome with --reference"

        """
        exit 1
        }
    }
    if (params.skip_external_bams == false & params.bam_dir == null & params.help == false) {
        println """

        Please specify a species: c_elegans c_brigssae c_tropicalis with option --species, a pre-existing bam directory with --bam_dir, or skip merging with pre-existing bams with --skip_external_bams"

        """
        exit 1
    }
    bam_dir = params.bam_dir
}



// For now, this pipeline requires NXF_VER 23.0
// Prefix this version when running
// e.g.
// NXF_VER=23.0 nextflow run ...
// assert System.getenv("NXF_VER") >= "23.0"

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
nextflow main.nf --sample_sheet=name_of_sample_sheet.tsv --species=c_elegans

    parameters              description                                 Set/Default
    ==========              ===========                                 ========================
    --debug                 Use --debug to indicate debug mode          ${params.debug}
    --sample_sheet          See test_data/sample_sheet for example      ${params.sample_sheet}
    --species               Species to map: 'c_elegans', 'c_briggsae'   ${params.species}
                              or 'c_tropicalis'                         
    --fq_prefix             Path to fastq if not in sample_sheet        ${params.fq_prefix}
    --bam_dir               Path to pre-existing strain bams            ${bam_dir}
    --skip_kmers            Whether to skip counting kmers              ${params.skip_kmers}
    --skip_external_bams    Whether to skip merging with pre-existing
                            strain bam files                            ${params.skip_external_bams}
    --reference             genome.fasta.gz to use in place of default  ${params.reference}
    --output                Output folder name.                         ${params.output}

    username                                                            ${"whoami".execute().in.text}
    ----------------------------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

"""
out
}


log.info(log_summary())


if (params.help) {
    exit 1
}

workflow {
    
    // Read sample sheet
    sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                        .ifEmpty { exit 1, "sample sheet not found" }
                        .splitCsv(header:true, sep: "\t")

    // check software
    // summary(Channel.from("run"))

    aln_in = sample_sheet.map { row -> row.fq1 = params.fq_prefix ? row.fq1 = params.fq_prefix + "/" + row.fq1 : row.fq1; row }
                .map { row -> row.fq2 = params.fq_prefix ? row.fq2 = params.fq_prefix + "/" + row.fq2 : row.fq2; row }
                .map { row -> [row,
                               "$params.reference".substring(0, "$params.reference".lastIndexOf("/")),
                               "$params.reference".substring("$params.reference".lastIndexOf("/") + 1),
                               file(row.fq1),
                               file(row.fq2)]}
    aln_in | alignment
    aln_in | kmer_counting
    kmer_counting.out | aggregate_kmer

    /* Merge Bams */
    if (params.skip_external_bams == false) {
        inventory_existing_bams( bam_dir, Channel.fromPath(bam_dir) )
        existing_bams = inventory_existing_bams.out.splitCsv(sep: "\t")
            .join( alignment.out
                       .map { row -> [row[0].strain, row[1], row[2]] }
                       .groupTuple() )
            .map{ row -> [row[0], row[1], row[2]] }
        
        /* only stage strains with multiple bams for merging */
        alignment.out.map { row -> [row[0].strain] + [row[1], row[2]] }
            .concat( existing_bams)
            .groupTuple()
            .map { strain, bam, bai -> [strain, bam, bai, bam.size()] }
            .branch { strain, bam, bai, count ->
                need_merge: count > 1
                no_merge: count < 2 }
            .set{ merging }
    } else {
        /* only stage strains with multiple bams for merging */
        alignment.out.map { row -> [row[0].strain] + [row[1], row[2]] }
            .groupTuple()
            .map { strain, bam, bai -> [strain, bam, bai, bam.size()] }
            .branch { strain, bam, bai, count ->
                need_merge: count > 1
                no_merge: count < 2 }
            .set{ merging }
    }

    /* only merge strains with multiple bams */
    merge_bam( merging.need_merge )

    /* recombine merged bams with bams not needing merging */
    merged_bams = merge_bam.out
        .mix( merging.no_merge.map { strain, bam, bai, count -> [strain, bam, bai] } )

    mark_dups( merged_bams )

    /* ID Level Stats and multiqc */
    alignment.out.map { row, bam, bai -> ["id", row.id, bam, bai ]} | \
        (coverage_id & idxstats_id & flagstat_id & stats_id)

    coverage_id.out.concat(idxstats_id.out,
                           flagstat_id.out,
                           stats_id.out).collect()
                           .combine(Channel.fromPath("${workflow.projectDir}/scripts/multiqc_config.yaml"))
                           .map{ row -> ["id", row] }
                            // | multiqc_id

    /* Strain Level Stats and multiqc */
    mark_dups.out.bams.map { strain, bam, bai -> ["strain", strain, bam, bai] } | \
        (coverage_strain & idxstats_strain & flagstat_strain & stats_strain & validatebam_strain)

    mark_dups.out.markdups.concat(coverage_strain.out,
                                  idxstats_strain.out,
                                  flagstat_strain.out,
                                  stats_strain.out,
                                  validatebam_strain.out).collect()
                                 .combine(Channel.fromPath("${workflow.projectDir}/scripts/multiqc_config.yaml"))
                                 .map{ row -> ["strain", row] }
                                //   | multiqc_strain

    /* Generate a bam file summary for the next step */
    strain_summary = mark_dups.out.strain_sheet.map { strain, bam, bai -> [strain, "${strain}.bam","${strain}.bam.bai"].join("\t") } \
                 .collectFile(name: 'strain_summary.tsv',
                              newLine: true,
                              storeDir: "${params.output}")

    // Find final coverage of each strain bam
    find_final_coverage( mark_dups.out.bams )

    // summarize coverage
    // multiqc_strain.out
    //     .combine(strain_summary) 
    //     .combine(Channel.fromPath(params.sample_sheet))
    //     .combine(Channel.fromPath("${workflow.projectDir}/scripts/low_map_cov_for_seq_sheet.Rmd")) | coverage_report
        //.combine(summary.out) 

    // check for npr-1 allele
    if(params.species == "c_elegans") {
        merge_bam.out.combine(Channel.fromPath("${params.reference}")) | npr1_allele_check

        npr1_allele_check.out.collect() | npr1_allele_count
    }
    
    // blobtools
    if(params.blob) {
        coverage_report.out.low_strains
            .splitCsv(sep: '\n', strip: true) | blob_align | blob_assemble | blob_unmapped
    
        blob_unmapped.out
            .combine(Channel.fromPath("${params.ncbi}")) | blob_blast | blob_plot
    }
}


process summary {
    
    executor 'local'
    container null

    // conda 'fd-find'
    publishDir "${params.output}", mode: 'copy'
    
    input:
        val(run)

    output:
        tuple path("sample_sheet.tsv"), path("summary.txt"), path("software_versions.txt")

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

    array 100
    tag { data.id }
    
    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 4 * task.attempt }
    memory = { 16.GB * task.attempt }


    input:
        tuple val(data), path(genome_path), val(genome_basename), path(fq1), path(fq2)
        
    output:
        tuple val(data), file("${data.id}.bam"), file("${data.id}.bam.bai")

	script:
		// Construct read group
		RG = ["@RG",
			  "ID:${data.id}",
			  "SM:${data.strain}",
			  "LB:${data.lb}",
			  "PL:illumina"].join("\\t")

    """
    INDEX=`find -L ${genome_path} -name "${genome_basename}.amb" | sed 's/\\.amb\$//'`
    bwa mem -t ${task.cpus} -R '${RG}' \${INDEX} ${fq1} ${fq2} | \\
    sambamba view --nthreads=${task.cpus} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
    sambamba sort --nthreads=${task.cpus} --show-progress --tmpdir=. --out=${data.id}.bam /dev/stdin
    sambamba index --nthreads=${task.cpus} ${data.id}.bam
    # if [[ ! \$(samtools view ${data.id}.bam | head -n 10) ]]; then
    #     exit 1;
    # fi
    """
}

process inventory_existing_bams {
    
    executor 'local'
    container null
    
    input:
        val(bam_path)
        path(bam_dir)

    output:
        path("bam_list.txt")

    script:
    """
    for I in ${bam_dir}/*.bam; do
        NAME=\$(basename \${I})
        echo -e "\${NAME/.bam/}\\t${bam_path}/\${NAME}\t${bam_path}/\${NAME}.bai" >> bam_list.txt
    done
    """
}


/* 
    Merge ID Bams by Strain
*/
process merge_bam {

    array 20
    tag { "${strain}" }

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
        tuple val(strain), path(bam), path(bai), val(n_count)

    output:
        tuple val(strain), file("new_${strain}.bam"), file("new_${strain}.bam.bai")

    script:
    """
    sambamba merge --nthreads=${task.cpus} --show-progress new_${strain}.bam ${bam}
    sambamba index --nthreads=${task.cpus} new_${strain}.bam
    """
}

process mark_dups {

    array 20
    tag { "${strain}" }

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    publishDir "${params.output}/bam", mode: 'copy', pattern: '*.bam*'

    input:
        tuple val(strain), path("${strain}.in.bam"), path("${strain}.in.bam.bai")
    output:
        tuple val(strain), path("${strain}.bam"), path("${strain}.bam.bai"), emit: "bams"
        tuple val(strain), path("${strain}.bam"), path("${strain}.bam.bai"), emit: "strain_sheet"
        path "${strain}.duplicates.txt", emit: "markdups"
        tuple path("${strain}.bam"), path("${strain}.bam.bai"), emit: "npr"

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

/* 
    Report for low coverage strains
*/

process coverage_report {

    label 'xs'
    label 'R'

    publishDir "${params.output}", mode: 'copy'

    //errorStrategy 'ignore'

    input:
        tuple path("report.html"), path("strain_data/*"), path("strain_summary"), path("sample_sheet"), path("original_low_map_cov_for_seq_sheet.Rmd")

    output:
        path("low_map_cov_for_seq_sheet.html")
        path("low_map_cov_for_seq_sheet.Rmd")
        // path "strains_with_low_values.tsv", emit: "low_strains"
        path "fastq_for_blobtools.tsv", emit: "low_strains"
        path("*.tsv")


    """
    cat original_low_map_cov_for_seq_sheet.Rmd | \\
        sed -e 's/read.delim("sample_sheet.tsv"/read.delim("${sample_sheet}"/g' | \\
        sed -e 's/strain_summary.tsv/${strain_summary}/g' | \\
        sed -e 's+_aggregate/multiqc/strain_data+strain_data+g' > low_map_cov_for_seq_sheet.Rmd
    Rscript -e "rmarkdown::render('low_map_cov_for_seq_sheet.Rmd')"

    """

}

/* 
    Quick check for npr-1 allele (used to be part of concordance)
*/

process npr1_allele_check {
    
    array 20
    label 'sm'
    label 'postgatk'

    input:
        tuple val(strain), path("${strain}.in.bam"), path("${strain}.in.bam.bai"), path("reference")

    output:
        path("${strain}.npr1.bcf")

    """
    bcftools mpileup -f ${reference} \\
    -r X:4768788 \\
    ${strain}.in.bam |  \\
    bcftools call -mv -Ob -o ${strain}.npr1.bcf
    """

}

// Probably not the prettiest way to do this, but gets the job done

process npr1_allele_count {
    
    array 20
    label 'sm'
    label 'postgatk'

    publishDir "${params.output}/", mode: 'copy'

    input:
        path("npr_bcf")

    output:
        path("npr1_allele_strain.tsv")

    """
    # index bcf
    for v in `ls | grep npr_bcf`;
    do
    bcftools index \$v;
    done

    echo -e 'problematic_strain\\tgt' > npr1_allele_strain.tsv
    bcftools merge ${npr_bcf} --missing-to-ref | \\
    bcftools query -f '[%SAMPLE\\t%GT\\n]' | awk '\$2 != "1/1"' >> npr1_allele_strain.tsv

    """

}

process find_final_coverage {

    array 20
    
    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }

    publishDir "${params.output}/bam", mode: 'copy', pattern: '*.coverage'

    input:
        tuple val(strain), path(bam), path(bai)

    output:
        tuple val(strain), path("${strain}.coverage")

    script:
    """
    mosdepth --threads ${task.cpus} ${strain} ${bam}
    awk -v STRAIN="${strain}" 'BEGIN{N=0} \
                               { if (NR > 1){N=N+1; CHROMS[N]=\$1; COVERAGE[\$1]=\$4;}} \
                               END{printf "strain\\tcoverage"; for (I=1;I<N;I++) printf "\\t%s", CHROMS[I]; printf "\\n"; \
                                   printf "%s\\t%f", STRAIN, COVERAGE[CHROMS[N]]; for (I=1;I<N;I++) printf "\\t%s", COVERAGE[CHROMS[I]]; \
                                   printf "\\n";}' ${strain}.mosdepth.summary.txt > ${strain}.coverage
    rm ${strain}.mosdepth*
    rm ${strain}.per-base*
    """
}

/* 
    Blobtools on low coverage strains
*/

process blob_align {

    label 'md'
    label 'blob'

    input:
        val(STRAIN)

    output:
        tuple val(STRAIN), path("Unmapped.out.mate1.step1.fq"), path("Unmapped.out.mate2.step1.fq")


    script:
    def avail_mem = (task.memory.byte).intValue()
    """
    # get fastq pair
    # st=`echo ${STRAIN} | sed 's/\\[//' | sed 's/\\]//'`
    # p1=`cat ${params.sample_sheet} | awk -v st="\$st" '\$0 ~ st { print "${params.fq_prefix}"\$4 }'`
    # p2=`cat ${params.sample_sheet} | awk -v st="\$st" '\$0 ~ st { print "${params.fq_prefix}"\$5 }'`

    # in case of multiple fastq... combine with comma
    # p1=`echo \$p1 | sed 's/ /,/g'`
    # p2=`echo \$p2 | sed 's/ /,/g'`

    p1=`echo "${params.fq_prefix}/"${STRAIN} | sed 's/\\[//' | sed 's/\\]//'`
    p2=`echo \$p1 | sed 's/1P/2P/'`

    # change ref to unzip
    ref=`echo ${params.reference} | sed 's/.gz//'`

    STAR \\
    --runThreadN ${task.cpus} \\
    --runMode genomeGenerate \\
    --limitGenomeGenerateRAM ${avail_mem} \\
    --genomeDir . \\
    --genomeFastaFiles \$ref \\
    --genomeSAindexNbases 12 
    STAR \\
    --runThreadN ${task.cpus} \\
    --genomeDir . \\
    --outSAMtype BAM Unsorted SortedByCoordinate \\
    --outReadsUnmapped Fastx \\
    --twopassMode Basic \\
    --readFilesCommand zcat \\
    --readFilesIn \$p2 \$p1 

    mv Unmapped.out.mate1 Unmapped.out.mate1.step1.fq
    mv Unmapped.out.mate2 Unmapped.out.mate2.step1.fq


    """

}

process blob_assemble {

    label "blob"
    label "md"

    input:
        tuple val(STRAIN), path("Unmapped_mate1_step1.fq"), path("Unmapped_mate2_step1.fq")

    output:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam")


    script:
    def avail_mem = (task.memory.byte).intValue()
    """
    
    spades.py --pe1-1 Unmapped_mate1_step1.fq  --pe1-2 Unmapped_mate2_step1.fq -t 24 -m 160 --only-assembler -o UM_assembly

    STAR \\
    --runThreadN ${task.cpus} \\
    --runMode genomeGenerate \\
    --limitGenomeGenerateRAM ${avail_mem} \\
    --genomeDir . \\
    --genomeFastaFiles UM_assembly/scaffolds.fasta \\
    --genomeSAindexNbases 5 \\

    STAR \\
    --runThreadN ${task.cpus} \\
    --genomeDir . \\
    --outSAMtype BAM Unsorted SortedByCoordinate \\
    --outReadsUnmapped Fastx \\
    --twopassMode Basic \\
    --readFilesIn Unmapped_mate1_step1.fq Unmapped_mate2_step1.fq

    """ 
}

process blob_unmapped {

    label 'alignment'

    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }

    input:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam")

    output:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam"), path("Aligned.sortedByCoord.out.bam.bai")


    """
    samtools index Aligned.sortedByCoord.out.bam
    """ 
}

process blob_blast {

    label "blob"
    label "md"

    input:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam"), path("Aligned.sortedByCoord.out.bam.bai")

    output:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam"), path("Aligned.sortedByCoord.out.bam.bai"), path("assembly.1e25.megablast.out"), \
        path("ncbi_nt")


    """
    blastn \\
    -task megablast \\
    -query UM_assembly/scaffolds.fasta \\
    -db ${params.ncbi}/nt \\
    -outfmt '6 qseqid staxids bitscore std' \\
    -max_target_seqs 1 \\
    -max_hsps 1 \\
    -num_threads 24 \\
    -evalue 1e-25 \\
    -out assembly.1e25.megablast.out

    """ 
}

process blob_plot {

    label "blob"
    label "sm"

    publishDir "${workflow.launchDir}/${params.output}/blobtools/", mode: 'copy'

    input:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam"), path("Aligned.sortedByCoord.out.bam.bai"), path("assembly.1e25.megablast.out"), \
        path("ncbi_nt")

    output:
        tuple file("*.png"), file("*blobplot.stats.txt")

    """
    st=`echo ${STRAIN} | sed 's/\\[//' | sed 's/\\]//'`
    blobtools create -i UM_assembly/scaffolds.fasta -b Aligned.sortedByCoord.out.bam -t assembly.1e25.megablast.out -o \$st --names ${params.ncbi}/names.dmp --nodes ${params.ncbi}/nodes.dmp --db nodesDB.txt
    blobtools plot -i \$st.blobDB.json -o \$st.plot

    """ 
}

