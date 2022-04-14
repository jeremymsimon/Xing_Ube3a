---
Title: Xing_Ube3a
Author: Jeremy M. Simon
Date: 4/14/22
output: github_document
---

## Xing et al. 2022 (in preparation)

This repository contains `R` code to analyze scRNA-seq data consistent with the figures presented in Xing _et al._ 2022 (_in preparation_). 

We analyzed scRNA-seq (drop-seq) data from the cerebral cortex of embryonic (E14.5) and early postnatal (P0) mice to explore the consequences of Ube3a-T485A homozygous mutation. 

### Abstract
The E3 ubiquitin ligase Ube3a is biallelically expressed in mitotic cells, including neural progenitors and glial cells, raising the possibility that UBE3A gain-of-function mutations might cause neurodevelopmental disorders irrespective of parent-of-origin. To test this possibility, we engineered a mouse line that harbors an autism-linked UBE3AT485A (T508A in mouse) gain-of-function mutation and evaluated phenotypes in animals that inherited the mutant allele paternally, maternally, or from both parents. We found that both paternally and maternally expressed UBE3AT485A resulted in elevated UBE3A activity in neural progenitors and glial cells where Ube3a is biallelically expressed. Expression of UBE3AT485A from the maternal allele, but not the paternal one, led to a persistent elevation of UBE3A activity in postmitotic neurons. Maternal, paternal, or biparental inheritance of the mutant allele promoted embryonic expansion of Zcchc12 lineage interneurons which mature into Sst and Calb2 expressing interneurons, and caused a spectrum of behavioral phenotypes that differed by parent-of-origin. Phenotypes were distinct from those observed in Angelman syndrome model mice that harbor a Ube3a maternal loss-of-function allele. Our study shows that the UBE3AT485A gain-of-function mutation causes distinct neurodevelopmental phenotypes when inherited maternally or paternally. These findings have clinical implications for a growing number of disease-linked UBE3A gain-of-function mutations. 

### What you'll find here

* `Seurat.R`: Import counts matrices as quantified by [Alevin](https://salmon.readthedocs.io/en/latest/alevin.html), perform normalization using SCTransform, data integration and exploratory analyses using [Seurat](https://satijalab.org/seurat/index.html)
	+ Analysis for the publication was performed using `R` v3.5.2

* `PropTesting.R`: Execute statistical tests comparing cell proportions within each identified cluster using [speckle](https://github.com/Oshlack/speckle)

* `Slingshot.R`: Pseudotiming analysis using [Slingshot](https://github.com/kstreet13/slingshot)
	+ Analysis for the publication was performed using `R` v4.1.0
	+ We updated the Seurat object to be compatible with `Seurat` v4.1 using `UpdateSeuratObject()`

### What you'll find elsewhere

* All data have been deposited to GEO under accession GSEXXXXXXX

* The raw data matrices from Alevin are also available for download from Zenodo [here](https://zenodo.org/record/6459852)

* The Ube3a WT data presented here [has been published previously](https://www.nature.com/articles/s41467-018-08079-9) and the original code for analysis was deposited [here](https://github.com/jeremymsimon/MouseCortex). We also created an [interactive webtool](https://zylkalab.org/datamousecortex) for exploration of the WT data

