# CoCoCoNet User Manual

CoCoCoNet is a simple to use webserver that allows the user to build, view and analyze co-expression networks without having to input experimental data. Because gene expression data is inherently noisy, meta-analysis provides a method to significantly improve the quality of data when assessed using the Guilt by Association (GBA) principle. This data is curated by aggregating the expression reads of many RNA-seq experiments obtained through the NCBI's SRA database [1] using guidelines outlined in [2]. 

To use CoCoCoNet, simply input a list of genes (or a single gene) as gene symbols, EntrezID or EnsemblID and the corresponding species to be used in the construction of the network. Given these, you can select optional parameters, described below, generate the results and view the distribution of co-expression values as well as the network. CoCoCoNet also allows you to compare the network with another species using a 1 to 1 ortholog mapping [3] as well as perform a GBA analysis using EGAD [4].

CoCoCoNet is free to use and available at https://milton.cshl.edu/CoCoCoNet/

## Step 1: Initialization
![initialize](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/initialize.png)
Preloaded genes can be used to test the server’s functionality with the “Example settings”. Users can choose either the top 231 co-expressed yeast genes studied by Eisen et al (5) or the top 102 co-expressed genes associated with Autism Spectrum Disorder (ASD) from Satterstrom et al (6) . More on these examples in the use cases section of the manuscript and this document. In the steps below, we have selected the yeast set.

“__Select Input method:__” allows users to choose how to load in a gene list. Currently, users can select genes from a drop down list, paste a comma separated list, or upload a file with genes listed in new lines. Toggling between these will change the adjacent “Upload gene list” box. 

Users can extend their analyses beyond their gene set with the “__Using genes__” options. 

•	“__That I provide only__" will construct the network using only selected genes. This is the default setting.

•	"__That I provide plus more highly co-expressed genes__" will select the most closely related genes not in the provided set. Here we define the "relation" as having the largest weighted degree of edges connected to the provided genes.



“__Compare my genes to__” options allow users to limit the genes to a high confidence set or all genes based on expression levels. 

•	“__A high confidence gene set__” will match input genes to only a subset of genes filtered on minimum expression level across experiments. This is the default setting.

•	“__Almost all genes__” will match your genes to a very lightly filtered gene set.

  
## Visualizing Results

Once your genes have been properly loaded, the "Generate Results" and "Clear Results" options will appear. Selecting "Generate Results" will display the network, the distribution of co-expression values, and a sliding threshold bar that allows the user to filter the network to include only connections greater than this threshold. The user also has the option to highlight genes along with their direct connections or highlight genes with specified enriched GO terms. Using the example setting of "Highly co-expressed yeast genes" and selecting genes with GO term "ribosome" from the drop down box gives the following network – highlighting genes involved in the ribosome. 

![part1](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/part1.png)

## Ortholog Mapped genes

The next section allows the user to input a second species to compare to the first. Doing so will select genes of the second species with a 1 to 1 ortholog of genes in the provided gene set. After selection, "Generate" and "Clear Section" buttons will appear where selecting "Generate" will again display the network, the co-expression value distribution, and a sliding threshold bar. Again the user has the option to highlight genes along with its nearest neighbors or highlight genes by GO term. Here, we select mouse and selecting by gene, which also highlights the gene in orange and the order 1 connections in green. Red Stars indicate hub genes which are defined to be genes with the largest node degree.

![compare](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/compare.png)


## Extending Guilt by Association by Degree (EGAD)

The final section allows the user to perform Guilt by Association (GBA) analysis on the input genes and the corresponding 1 to 1 ortholog using EGAD [4]. EGAD analyzes enriched GO terms of each species using either neighbor voting or by node degree and reports the corresponding area under the receiver operating curve (AUROC) or the precision recall curve (AUPRC) across 3 cross validation folds. Four figures are shown, the top two figures give the distribution of AUC scores for each species, the bottom figures compares the two speciess. On the left, we compare the AUC scores of common GO terms and on the right, we compare the enrichment p-values for similar GO terms. Optionally, the user can elect to obtain the gene set score which is a measure of connectivity of the input genes to the "High Confidence" gene set or the "Almost all" gene set chosen in step 1. 

![assess](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/assess.png)

## Downloadables

Networks can be downloaded by right clicking and saving the image, histograms and scatter plots can be saves using the "download plot as a png" utility and relevant data can be downloaded at each stage by toggling the "Download" switch.
![exporting](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/exporting.png)


Raw data used to generate results can be downloaded at ![ftp://milton.cshl.edu/data](ftp://milton.cshl.edu/data). Here, you can find: rank standardized co-expression networks, gene to GO relationships, 1-to-1 ortholog maps, gene conversion tables, the null distributions of AUROC scores (measured by EGAD) and node degree, sample input and R scripts used to generate results as well as a CoCoCoNetlite version. 

![exporting](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/download.png)

## Browser Support

| OS | Version | Chrome | FireFox | Microsoft  Edge | Safari | 
| ----- | ----- | ----- | ----- | -----  | ----- | 
| Linux | Ubuntu | 79.0.3945 | 71.0 | N/A | N/A | 
| MacOS | Catalina | 79.0.3945 | 71.0 | N/A | 10.1.2 | 
| Windows | 10 | 79.0.3945 | 71.0 | 44.18362.449.0 | N/A | 


## References

[1] Leinonen, R., Sugawara, H., Shumway, M. on behalf of the International Nucleotide Sequence Database Collaboration (2011) The Sequence Read Archive, _Nucleic Acids Research_, __(39)__, suppl_1, D19–D21, https://doi.org/10.1093/nar/gkq1019

[2] Ballouz, S., Verleyen, W., Gillis, J. (2015) Guidance for RNA-seq co-expression network construction and analysis: safety in numbers, _Bioinformatics_, __(31)__, 13, 2123–2130, https://doi.org/10.1093/bioinformatics/btv118

[3] Kriventseva, E. V., Kuznetsov, D., Tegenfeldt, F., Manni, M., Dias, R., Simão, F. A., Zdobnov, E. M. (2019) OrthoDB v10: sampling the diversity of animal, plant, fungal, protist, bacterial and viral genomes for evolutionary and functional annotations of orthologs, _Nucleic Acids Research_, __(47)__, D1, D807–D811, https://doi.org/10.1093/nar/gky1053

[4] Ballouz, S., Weber, M., Pavlidis, P., Gillis, J. (2017) EGAD: ultra-fast functional analysis of gene networks, _Bioinformatics_, __(33)__, 4, 612–614, https://doi.org/10.1093/bioinformatics/btw695

[5] Eisen, M. B., Spellman, P. T., Brown, P. O., Botstein, D. (1998) Cluster analysis and display of genome-wide expression patterns, _Proceedings of the National Academy of Sciences_, __(95)__, 25, 14863-14868; https://www.pnas.org/content/95/25/14863

[6] Satterstrom F., et. al.(2020) Large-Scale Exome Sequencing Study Implicates Both Developmental and Functional Changes in the Neurobiology of Autism, _Cell_, __180__: 568-584; https://www.cell.com/cell/fulltext/S0092-8674(19)31398-4#%20
