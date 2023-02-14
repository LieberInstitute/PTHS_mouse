# A myelin-related transcriptomic profile is shared by Pitt–Hopkins syndrome models and human autism spectrum disorder.

## [Citation](https://doi.org/10.1038/s41593-019-0578-x)
```Phan, B.N., Bohlen, J.F., Davis, B.A. et al. A myelin-related transcriptomic profile is shared by Pitt–Hopkins syndrome models and human autism spectrum disorder. Nat Neurosci 23, 375–385 (2020). https://doi.org/10.1038/s41593-019-0578-x```

# Description of data & code


## Github
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

## [Globus](https://app.globus.org/file-manager?origin_id=eee07044-9e5c-11ed-b579-33287ee02ec7&origin_path=%2F)

Intermediates and final summary files describing the transcriptomics analyses are deposited on Globus. These are all files that would be found at `/dcl01/lieber/ajaffe/Brady/mouseRNAseq` and subdirectories within the Github Rscripts. To access Globus, a user would have to [make a Globus ID](https://www.globusid.org/create) and [follow instructions](https://docs.globus.org/how-to/get-started/) to grab these data.





Internal JHPCE location: 
Raw data files: `/dcl01/lieber/ajaffe/Brady/mouseRNAseq`
Code files: `/users/bphan/tcf4/PTHS_mouse`
