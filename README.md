# Transcriptomic

Here, I analyzed  expression profiling by high throughput sequencing in accession number GSE121884. 

Above is presented transcriptional analysis of the pancreatic islets of adult Paupar KO mice versus Paupar WT mice.
Aligned reads were assigned to genes using annotations from Ensembl (Mus_musculus.GRCm38.73.gtf) and HTseq-count v0.10.0 (Anders et al., 2015). Differential expression across cohorts was assessed using DESeq2 v1.18 (Love, Huber, and Anders, 2014). Among the downregulated genes with a p-value < 0.05 genes important for alpha cell identity and function.



Firstly, we downloaded data and run FASTQC and MULTIQC. Even though quality of reads was not high, I proceeded to work with row reads. 

Then I run kallisto for quantifying abundances of transcripts. 

In R, I сreated a phenotypic table of the experiment for splitting the data. It was divided in 4 groups. Then reads were filtered to ignore low-expr essed genes. I took wild sample and compared with mutated one (plotCount plot)

Volcano plot shows genes with large fold changes. Ut is probably the most biologically significant genes. Towards the right, there are the most upregulated genes, and towards the left there are the most downregulated genes (volcano plot)  Then the consistency tests were performed. I got a heatmap of 100 most expressed genes (heatmap plot). Also I build a plot of distance between samples (PCA plot), which reduces the dimensionality and allows to cluster data into groups.  Then I investigated enrichments of genes in our samples. It turns out that "0 enriched terms found"... It sounds strange... Later I got the plot of "Starch and sucrose metabolism - Mus musculus" (metabolism_mus plot). Then I analyzed sets according to the content of three sub-ontologies: biological process, molecular function, and cellular component. Obviosly, mutated genes have affected to helicase activity, chromatin DNA binding, lamin binding and ATP-dependent actuvity. 


Lastly, I run the set enrichment analysis ----------------------------------


