# CoCoCoNet User Manual (UNDER CONSTRUCTION)

CoCoCoNet is a simple to use webserver that allows the user to build, view and analyze co-expression networks without having to input experimental data. Because gene expression data is inherently noisy, meta-analysis provides a method to significantly improve the quality of data when assessed using the Guilt by Association principle. This data is curated by aggregating the expression reads of many RNA-seq experiments obtained through the NCBI's SRA database using guidelines outlined in 

To use CoCoCoNet, simply input a list of genes (or a single gene) and the corresponding species to be used in the construction of the network. Given these, you can select optional parameters, described below, generate the results and view the distribution of co-expression values as well as the network. 



## The Parameters
![param](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/main.png)

**The Species**
  
  Select the species that you wish to match your input gene list to. 

**Select Input method:**
  
  Allows you to choose how to load in a gene list. Currently, you can select genes from a drop down list, paste a comma separated list, or upload a file with genes listed in new lines. 


**Priority vs Meta**
  
  Selecting Priority will only use a subset of the gene set filtered on expresionality across experiments.
  
  Meta on the otherhand will build the network using all available genes in the gene set.
  

**Using genes ...**
  
  "That I provide only" will construct a network using only selected genes
  
  "That I provide plus __ more" will select the most closely related genes not in the provided set. Here we define the "relation" as having the largest weighted degree of edges connected to the provided genes.
  
  
## Generating the results

Once your genes have been properly loaded, the "Generate Results" and "Clear Results" options will appear. Selecting "Generate Results" will display the network, the distribution of Co-expression values, and a sliding threshold bar that allows the user to filter the network to include only connections greater than this threshold. The user also has the option to highlight genes along with their direct connections or highlight genes with specified enriched GO terms. Using the example setting of "Top Co-expressed Yeast Genes" and selecting genes with GO term "translation initiation" gives the folowwing network.
![yeast](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/yeast.png)


## Ortholog Mapped genes

The next section allows the user to input a second species to compare to the first. Doing so will select genes of the second species with a 1 to 1 ortholog of genes in the provided gene set. After selection, "Generate" and "Clear Section" buttons will appear where selecting "Generate" will again display the network, the co-expression value distribution, and a sliding threshold bar. Again the user has the option to highlight genes along with its nearest neighbors or highlight genes by GO term.  
![C_elegans](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/roundworm.png)


## Extending Guilt by Association by Degree (EGAD)
The final section allows the user to perform Guilt by Association (GBA) analysis on each of the above networks using EGAD and visualize the distribution of AUROC  or AUPRC using either neighbor voting or node degreee.


## Downloadables
![download](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/download.png)

You can download the data used to generate the 


## Case Study
