# alignment-nf

### Typical use for debugging:

```
nextflow main.nf --debug
```

### Typical use for new fastq:

```
nextflow main.nf --sample_sheet=name_of_sample_sheet.tsv --species=c_elegans
```

### Parameters
    parameters              description                                 Set/Default
    ==========              ===========                                 ========================
    --debug                 Use --debug to indicate debug mode          (optional)
    --sample_sheet          See test_data/sample_sheet for example      (required)
    --species               Species to map: 'ce', 'cb' or 'ct'          (required)
    --fq_prefix             Path to fastq if not in sample_sheet        (optional)
    --kmers                 Whether to count kmers                      false
    --reference             genome.fasta.gz to use in place of default  (optional)
    --output                Output folder name.                         alignment-date in current folder


### Overview

![Overview of alignment-nf](http://andersenlab.org/dry-guide/img/alignment.png)


