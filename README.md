# scRNA-Seq fGSEA  

A good introduction on gene set enrichment analysis can be found there  

https://www.sc-best-practices.org/conditions/gsea_pathway.html  

The demo consists in integratating (merging) datasets from two conditions of cells (control and stimulated by interferon).    
After integration, a fast implementation of the Wilcoxon test is performed by Presto and a fast GSEA is done by fGSEA.

For the integration step, see  https://satijalab.org/seurat/archive/v3.2/immune_alignment.html  
For Presto, see https://github.com/immunogenomics/presto/blob/master/vignettes/getting-started.Rmd
For fGSEA, see  https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html  
