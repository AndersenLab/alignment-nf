# alignment-nf

[Documentation](http://andersenlab.org/dry-guide/pipeline-alignment/)

## Quick Reference

This pipeline runs very fast locally with the test set. You can run the following to make changes and evaluate output:

__debug locally__

```bash
NXF_VER=19.09.0-edge nextflow run main.nf -profile local --debug true -resume
```