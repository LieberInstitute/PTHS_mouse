# A myelin-related transcriptomic profile is shared by Pitt–Hopkins syndrome models and human autism spectrum disorder.

## [Citation](https://doi.org/10.1038/s41593-019-0578-x)
```Phan, B.N., Bohlen, J.F., Davis, B.A. et al. A myelin-related transcriptomic profile is shared by Pitt–Hopkins syndrome models and human autism spectrum disorder. Nat Neurosci 23, 375–385 (2020). https://doi.org/10.1038/s41593-019-0578-x```

# Description of data & code

## Github
The Rscripts and affiliated tables or Rdata objects in this github repositories contain roughly half of what files are required to analyze mouse transcriptomic signatures of PTHS and the relation to other mouse models of syndromic ASD and human idiopathic ASD. The github file size limits prevent uploading of larger file objects. For those larger files, please see the Globus part of this README. 

1) [Sample sheet](tcf4_mouse/tables/stable1_RNAseq_sample_info.xlsx) describing files, genotype, sex, batch, etc between the 2 sets of RNA-seq from the Maher lab and Philpot lab. An aggregated phenotype table of these samples are at [here](tcf4_mouse/tables/maher_philpot_PTHS_mouse_rnaseq_pheno_relativePaths.csv)

2) DEG analyses
- [Maher lab Tcf4<sup>+/tr</sup> mouse DEG](tcf4_mouse/analyze_tcf4_mouse_DESeq2.R)
- [Maher lab Tcf4<sup>+/tr</sup> mouse DEG by age group](tcf4_mouse/analyze_tcf4_mouse_by_ages.R)
- [Mega DEG analyses across multiple Tcf4<sup>+/mut</sup> mouse models by age group](fig1c_analyze_mega_tcf4_ages.R)

3) Cell type-specific Expression Analysis (CSEA)
- [Across mega DEGs by age](fig2a_cell_type_enrichment.R)

4) [CIBERSORT analyses](fig2c_analyze_mega_tcf4_CIBERSORT.R)
Note: CIBERSORT cell type quantifications were estimated w/ the CIBERSORT web portal. These Rscripts only assess the group statistical differences in cell type abundance between animals and cohorts.

5) [Comparison to MeCP2 and Pten mutant mouse models](fig3_compare_mecp2_and_pten.R)
JHPCE path: `/users/bphan/tcf4/PTHS_mouse`

## Globus
Intermediates and and larger genomic data files generated along the RNA-seq data processing pipelines for this paper are deposited on Globus. These Globus would contain the files that would be referenced in the GitHub Rscripts pointing to files within , `/dcl01/lieber/ajaffe/Brady/philpot`, and . To access Globus, a user would have to [make a Globus ID](https://www.globusid.org/create) and [follow instructions](https://docs.globus.org/how-to/get-started/) to grab these data. Below are a brief description of the file types for various aggregates of data. The Globus links are provided below for different datasets and main file objects. The contents of each Globus roughly contains roughly the following:
- raw sequencing files within the `FASTQ` folder
- FastQC reports for each corresponding set of sequencing file with the `FastQC` folder
- alignment `.bam` files within `HISAT2_out` folder
- coverage `.bw` within `Coverage` folder
- featureCounts objects at the gene, exon, or junction levels in the `Counts` folder
- aggregated rawCount objects (# of times a unique mRNA fragment align to a gene)
- aggregated rpkmCount object (rawCounts normalized by gene size and sequencing depth).

1) Maher lab Tcf4<sup>+/tr</sup> mouse processed files:
- counts: `rawCounts_mouse_tcf4_n36_rerun_paired_stranded.rda`
- rpkm: `rpkmCounts_mouse_tcf4_n36_rerun_paired_stranded.rda`
- DESeq2 object from Maher lab Tcf4<sup>+/tr</sup> mouse DEG analyses: `mouse_tcf4_ages_DESeq2_svaAdj.rda`
- DESeq2 object from Mega DEG analyses: `mega_tcf4_ages_DESeq2_svaAdj.rda`
- JHPCE path: `/dcl01/lieber/ajaffe/Brady/mouseRNAseq`
- Globus path: https://app.globus.org/file-manager?origin_id=eee07044-9e5c-11ed-b579-33287ee02ec7&origin_path=%2F
Note: the Mega DEG DESeq2 objects will contain gene counts from Maher, Philpot, and Sweatt lab mouse models of PTHS.

[Philpot lab processed files](to add)
- Philpot lab Tcf4<sup>+/mut</sup> mouse counts: `rawCounts_philpot_OCT20_n58.rda`
- Philpot lab Tcf4<sup>+/mut</sup> mouse rpkm: `rpkmCounts_philpot_OCT20_n58.rda`
- DESeq object of Tcf4 mutation across Philpot mouse samples: `philpot_by_age_DESeq2_svaAdj.rda`
- JHPCE path: `/dcl01/lieber/ajaffe/Brady/philpot`
- Globus path: 
 

[Sweatt lab processed files](to add)
- Philpot lab Tcf4<sup>+/mut</sup> mouse counts: `rawCounts_philpot_OCT20_n58.rda`
- Philpot lab Tcf4<sup>+/mut</sup> mouse rpkm: `rpkmCounts_philpot_OCT20_n58.rda`
- JHPCE path: `/dcl01/lieber/ajaffe/Brady/sweatt`
