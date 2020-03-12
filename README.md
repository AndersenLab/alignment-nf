# alignment-nf

[Documentation](http://andersenlab.org/dry-guide/pipeline-alignment/)

## Overview

![Overview of alignment-nf](http://andersenlab.org/dry-guide/img/alignment.png)

## Debugging

This pipeline runs very fast locally with the test set. You can run the following to make changes and evaluate output:

__debug locally__

```bash
NXF_VER=19.09.0-edge nextflow run main.nf -profile local --debug true -resume
```

!! Important

If debugging locally, it is important that you increase the amount of memory docker has available or tools such as picard MarkDups will crash.