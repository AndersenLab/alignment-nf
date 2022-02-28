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
params.R_libpath = "/projects/b1059/software/R_lib_3.6.0"
params.species = "c_elegans"
params.ncbi = "/projects/b1059/data/other/ncbi_blast_db/"
params.blob = true

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
    params.fq_prefix = "/projects/b1059/data/${params.species}/WI/fastq/dna/"
}

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
    params.reference = "/projects/b1059/data/${params.species}/genomes/${params.project}/${params.ws_build}/${params.species}.${params.project}.${params.ws_build}.genome.fa.gz"
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
    --debug                 Use --debug to indicate debug mode          ${params.debug}
    --sample_sheet          See test_data/sample_sheet for example      ${params.sample_sheet}
    --species               Species to map: 'c_elegans', 'c_briggsae'   ${params.species}
                              or 'c_tropicalis'                         
    --fq_prefix             Path to fastq if not in sample_sheet        ${params.fq_prefix}
    --kmers                 Whether to count kmers                      ${params.kmers}
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
    // summary(Channel.from("run"))

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
    strain_summary = mark_dups.out.strain_sheet.map { row, bam, bai -> [row.strain, "${row.strain}.bam","${row.strain}.bam.bai"].join("\t") } \
                 .collectFile(name: 'strain_summary.tsv',
                              newLine: true,
                              storeDir: "${params.output}")

    // summarize coverage
    multiqc_strain.out
        .combine(strain_summary)
        .combine(summary.out) | coverage_report

    // check for npr-1 allele
    if(params.species == "c_elegans") {
        merge_bam.out.combine(Channel.fromPath("${params.reference}")) | npr1_allele_check

        npr1_allele_check.out.collect() | npr1_allele_count
    }
    
    // blobtools
    if(params.blob) {
        coverage_report.out.low_strains
            .splitCsv(sep: '\n', strip: true)
            .combine(Channel.fromPath(params.sample_sheet))
            .combine(Channel.fromPath("${params.reference}")) | blob_align | blob_assemble | blob_unmapped
    
        blob_unmapped.out
            .combine(Channel.fromPath("${params.ncbi}")) | blob_blast | blob_plot
    }
}


process summary {
    
    executor 'local'

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

    conda "/projects/b1059/software/conda_envs/cegwas2-nf_env"

    publishDir "${workflow.launchDir}/${params.output}/", mode: 'copy'

    //errorStrategy 'ignore'

    input:
        tuple path("report.html"), path("strain_data/*"), path("strain_summary"), path("sample_sheet"), path("summary"), path("software_versions")

    output:
        path("low_map_cov_for_seq_sheet.html")
        path("low_map_cov_for_seq_sheet.Rmd")
        // path "strains_with_low_values.tsv", emit: "low_strains"
        path "fastq_for_blobtools.tsv", emit: "low_strains"
        path("*.tsv")


    """
    cat "${workflow.projectDir}/scripts/low_map_cov_for_seq_sheet.Rmd" | \\
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
    conda "/projects/b1059/software/conda_envs/cegwas2-nf_env"

    input:
        tuple val(strain), row, path("${strain}.in.bam"), path("${strain}.in.bam.bai"), path("reference")

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
    conda "/projects/b1059/software/conda_envs/cegwas2-nf_env"

    publishDir "${workflow.launchDir}/${params.output}/", mode: 'copy'

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

/* 
    Blobtools on low coverage strains
*/

process blob_align {

    container 'andersenlab/blobtools:v2.1'

    cpus 12

    conda "/projects/b1059/software/conda_envs/blobtools"

    input:
        val(STRAIN)

    output:
        tuple val(STRAIN), path("Unmapped.out.mate1.step1.fq"), path("Unmapped.out.mate2.step1.fq")


    """
    # get fastq pair
    # st=`echo ${STRAIN} | sed 's/\\[//' | sed 's/\\]//'`
    # p1=`cat ${params.sample_sheet} | awk -v st="\$st" '\$0 ~ st { print "${params.fq_prefix}"\$4 }'`
    # p2=`cat ${params.sample_sheet} | awk -v st="\$st" '\$0 ~ st { print "${params.fq_prefix}"\$5 }'`

    # in case of multiple fastq... combine with comma
    # p1=`echo \$p1 | sed 's/ /,/g'`
    # p2=`echo \$p2 | sed 's/ /,/g'`

    p1=`echo "${params.fq_prefix}"${STRAIN} | sed 's/\\[//' | sed 's/\\]//'`
    p2=`echo \$p1 | sed 's/1P/2P/'`

    # change ref to unzip
    ref=`echo ${params.reference} | sed 's/.gz//'`


    STAR \\
    --runThreadN 12 \\
    --runMode genomeGenerate \\
    --limitGenomeGenerateRAM 600000000000 \\
    --genomeDir . \\
    --genomeFastaFiles \$ref \\
    --genomeSAindexNbases 12 
    STAR \\
    --runThreadN 12 \\
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

    memory '160 GB'
    cpus 24
    container 'andersenlab/blobtools:v2.1'

    conda "/projects/b1059/software/conda_envs/blobtools"

    input:
        tuple val(STRAIN), path("Unmapped_mate1_step1.fq"), path("Unmapped_mate2_step1.fq")

    output:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam")


    """
    spades.py --pe1-1 Unmapped_mate1_step1.fq  --pe1-2 Unmapped_mate2_step1.fq -t 24 -m 160 --only-assembler -o UM_assembly

    STAR \\
    --runThreadN 24 \\
    --runMode genomeGenerate \\
    --limitGenomeGenerateRAM 600000000000 \\
    --genomeDir . \\
    --genomeFastaFiles UM_assembly/scaffolds.fasta \\
    --genomeSAindexNbases 5 \\

    STAR \\
    --runThreadN 24 \\
    --genomeDir . \\
    --outSAMtype BAM Unsorted SortedByCoordinate \\
    --outReadsUnmapped Fastx \\
    --twopassMode Basic \\
    --readFilesIn Unmapped_mate1_step1.fq Unmapped_mate2_step1.fq

    """ 
}

process blob_unmapped {

    cpus 4

    conda "/projects/b1059/software/conda_envs/samtools"

    input:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam")

    output:
        tuple val(STRAIN), path("UM_assembly/scaffolds.fasta"), path("Aligned.sortedByCoord.out.bam"), path("Aligned.sortedByCoord.out.bam.bai")


    """
    samtools index Aligned.sortedByCoord.out.bam

    """ 
}

process blob_blast {

    cpus 4
    container 'andersenlab/blobtools:v2.1'

    conda "/projects/b1059/software/conda_envs/blast"

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

    executor 'local'

    publishDir "${workflow.launchDir}/${params.output}/blobtools/", mode: 'copy'
    container 'andersenlab/blobtools:v2.1'
    // container 'genomehubs/blobtoolkit:1.1'

    conda "/projects/b1059/software/conda_envs/blobtools"

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

