![Build Docker (env/align.Dockerfile)](https://github.com/AndersenLab/wi-gatk/workflows/Build%20Docker%20(env/align.Dockerfile)/badge.svg)


# alignment-nf

[TOC]

The [alignment-nf](https://github.com/AndersenLab/alignment-nf) pipeline performs alignment for wild isolate sequence data __at the strain level__, and outputs BAMs and related information. Those BAMs can be used for downstream analysis including variant calling, [concordance analysis](http://andersenlab.org/dry-guide/pipeline-concordance/), [wi-gatk-nf (variant calling)](http://andersenlab.org/dry-guide/pipeline-wi/) and other analyses.

This page details how to run the pipeline and how to add new wild isolate sequencing data. You can also find more information on the Andersen Lab [dry guide](http://andersenlab.org/dry-guide/latest/pipeline-alignment/).


# Pipeline overview

```

             ▗▖ ▝▜   ▝                       ▗      ▗▖ ▖▗▄▄▖
             ▐▌  ▐  ▗▄   ▄▄ ▗▗▖ ▗▄▄  ▄▖ ▗▗▖ ▗▟▄     ▐▚ ▌▐
             ▌▐  ▐   ▐  ▐▘▜ ▐▘▐ ▐▐▐ ▐▘▐ ▐▘▐  ▐      ▐▐▖▌▐▄▄▖
             ▙▟  ▐   ▐  ▐ ▐ ▐ ▐ ▐▐▐ ▐▀▀ ▐ ▐  ▐   ▀▘ ▐ ▌▌▐
            ▐  ▌ ▝▄ ▗▟▄ ▝▙▜ ▐ ▐ ▐▐▐ ▝▙▞ ▐ ▐  ▝▄     ▐ ▐▌▐
                         ▖▐
                         ▝▘
    parameters              description                                 Set/Default
    ==========              ===========                                 ========================
    --debug                 Use --debug to indicate debug mode          null
    --sample_sheet          See test_data/sample_sheet for example      null
    --species               Species to map: 'ce', 'cb' or 'ct'          null
    --fq_prefix             Path to fastq if not in sample_sheet        /projects/b1059/data/{species}/WI/fastq/dna/
    --kmers                 Whether to count kmers                      false
    --reference             genome.fasta.gz to use in place of default  defaults for c.e, c.b, and c.t
    --output                Output folder name.                         alignment-{date}

    HELP: http://andersenlab.org/dry-guide/pipeline-alignment/
```

![Pipeline-overview](alignment-nf.drawio.svg)

## Software requirements

* Nextflow v20.01+ (see the dry guide on Nextflow [here](quest-nextflow.md) or the Nextflow documentation [here](https://www.nextflow.io/docs/latest/getstarted.html)). On QUEST, you can access this version by loading the `nf20` conda environment prior to running the pipeline command:

```
module load python/anaconda3.6
source activate /projects/b1059/software/conda_envs/nf20_env
```

* Currently only runs on Quest with conda environments installed at `/projects/b1059/software/conda_envs/`


# Usage

## Testing on Quest

*This command uses a test dataset*

```
nextflow run main.nf --debug -profile quest
```

## Running on Quest

You should run this in a screen session.

```
nextflow run main.nf --sample_sheet <path_to_sample_sheet> --species c_elegans -profile quest -resume
```

# Parameters

## -profile

There are three configuration profiles for this pipeline.

* `local` - Used for local development.
* `quest` - Used for running on Quest.
* `gcp` - For running on Google Cloud (not currently active?).

## --sample_sheet

The `sample sheet` for alignment is the output from the [trim-fq-nf](https://github.com/AndersenLab/trim-fq-nf) pipeline. the `sample sheet` has the following columns:

* __strain__ - the name of the strain. Multiple sequencing runs of the same strain are merged together.
* __id__ - A unique ID for each sequencing run. This must be unique for every single pair of FASTQs.
* __lb__ - A library ID. This should uniquely identify a DNA sequencing library.
* __fq1__ - The path to FASTQ1
* __fq2__ - The path to FASTQ2

![Sample_sheet](img/alignment_sample_sheet.png)


The `library` column is a useful tool for identifying errors by variant callers. For example, if the same library is sequenced twice, and a variant is only observed in one sequencing run then that variant may be excluded as a technical / PCR artifact depending on the variant caller being used.

The alignment pipeline will merge multiple sequencing runs of the same strain into a single bam. However, summary output is provided at both the `strain` and `id` level. In this way, if there is a poor sequencing run it can be identified and removed from a collection of sequencing runs belonging to a strain. **For this reason, it is important that each id be unique and not just the strain name**

## --debug (optional)

You should use `--debug true` for testing/debugging purposes. This will run the debug test set (located in the `test_data` folder) using your specified configuration profile (e.g. local / quest / gcp).

For example:

```
nextflow run main.nf -profile quest --debug -resume
```

Using `--debug` will automatically set the sample sheet to `test_data/sample_sheet.tsv`

### --species (optional)

Defaults to "c_elegans", change to "c_briggsae" or "c_tropicalis" to select correct reference file.

### --fq_prefix (optional)

Within a sample sheet you may specify the locations of FASTQs using an absolute directory or a relative directory. If you want to use a relative directory, you should use the `--fq_prefix` to set the path that should be prefixed to each FASTQ.

### --kmers (optional)

__default__ = false

Toggles kmer-analysis

### --reference (optional)

A fasta reference indexed with BWA. WS245 is packaged with the pipeline for convenience when testing or running locally.

On Quest, the default references are here:

```
c_elegans: /projects/b1059/data/c_elegans/genomes/PRJNA13758/WS276/c_elegans.PRJNA13758.WS276.genome.fa.gz
c_briggsae: /projects/b1059/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome.fa.gz
c_tropicalis: /projects/b1059/data/c_tropicalis/genomes/NIC58_nanopore/June2021/c_tropicalis.NIC58_nanopore.June2021.genome.fa.gz
```

A different `--project` and `--wsbuild` can be used with the `--species` parameter to generate the path to other reference genomes such as:

```
nextflow run main.nf --species c_elegans --project PRJNA13758 --wsbuild WS280
```

### --ncbi (optional)

__Default__ - `/projects/b1059/data/other/ncbi_blast_db/`

Path to the NCBI blast database used for blobtool analysis. Should not need to change.

### --output (optional)

__Default__ - `alignment-YYYYMMDD`

A directory in which to output results. If you have set `--debug true`, the default output directory will be `alignment-YYYYMMDD-debug`.


# Output

```
├── _aggregate
│   ├── kmers.tsv
│   └── multiqc
│       ├── strain_data/
│       │   ├── mqc_mosdepth-coverage-dist-id_1.txt
│       │   ├── mqc_mosdepth-coverage-per-contig_1.txt
│       │   ├── mqc_mosdepth-coverage-plot-id_1.txt
│       │   ├── mqc_picard_deduplication_1.txt
│       │   ├── mqc_samtools-idxstats-mapped-reads-plot_Counts.txt
│       │   ├── mqc_samtools-idxstats-mapped-reads-plot_Normalised_Counts.txt
│       │   ├── mqc_samtools_alignment_plot_1.txt
│       │   ├── multiqc.log
│       │   ├── multiqc_data.json
│       │   ├── multiqc_general_stats.txt
│       │   ├── multiqc_picard_dups.txt
│       │   ├── multiqc_qualimap_bamqc_genome_results.txt
│       │   ├── multiqc_samtools_flagstat.txt
│       │   ├── multiqc_samtools_idxstats.txt
│       │   ├── multiqc_samtools_stats.txt
│       │   └── multiqc_sources.txt
│       ├── strain_multiqc_report.html
│       ├── id_data/
│       │   └──... (same as strain_data/)
│       └── id_multiqc_report.html
├── bam
│   ├── [strain].bam
│   └── [strain].bam.bai
├── blobtools
│   ├── {strain}.*.blobplot.bam0.png
│   ├── {strain}.*.blobplot.read_cov.bam0.png
│   └── {strain}.*.blobplot.stats.txt
├── software_versions.txt
├── sample_sheet.tsv
├── strain_summary.tsv
├── stats_strain_all.tsv
├── stats_strains_with_low_values.tsv
├── sample_sheet_for_seq_sheet.tsv
├── sample_sheet_for_seq_sheet_ALL.tsv
├── low_map_cov_for_seq_sheet.Rmd
├── low_map_cov_for_seq_sheet.html
└── summary.txt
```

Most files should be obvious. A few are detailed below.

* __software_versions.txt__ - Outputs the software versions used for every process (step) of the pipeline.
* __summary.txt__ - Outputs a summary of the parameters used.
* __sample_sheet.tsv__ - The sample sheet (input file) that was used to produce the alignment directory.
* __strain_summary.tsv__ - A summary of all strains and bams in the alignment directory.
* __aggregate__ - Stores data that has been aggregated across all strains or sequencing IDs. 
* __coverage__ - Contains coverage data at the strain or id level, presented in a variety of ways.
* __low_map_cov_for_seq_sheet.(Rmd/html)__ - Report showing low coverage or problematic strains to remove.
* __stats_strain_all.tsv__ - contains stats for all strains, with all replicates combined
* __stats_strains_with_low_values.tsv__ - contains stats for strains with either (1) low number of reads, (2) low mapping rate, and/or (3) low coverage
* __sample_sheet_for_seq_sheet.tsv__ - sample sheet to be added to google sheet, filtered to remove low coverage strains
* __sample_sheet_for_seq_sheet_ALL.tsv__ - sample sheet to be added to google sheet, contains all strains (use this one)
* __blobplot/__ - contains plots for low coverage strains to see if they show contamination issues and if they should be resequenced.
