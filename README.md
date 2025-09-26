---
description: >-
  Pipeline for analyzing variants from research and clinical ACCESS and IMPACT
  samples, and creating per-patient summaries of variants across the different
  samples.
---

# ACCESS Data Analysis Nextflow Pipeline

## Dependencies

* Nextflow ≥ 24.04.2
* Python ≥ 3.9
* pandas
* numpy

## Usage

```
nextflow run main.nf -c nextflow.config -profile conda,juno,accessv1
```

## Flowchart (proposed pipeline)

<figure><img src=".gitbook/assets/Flowchart (1).png" alt=""><figcaption></figcaption></figure>





