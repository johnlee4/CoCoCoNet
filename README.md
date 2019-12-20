# CoCoCoNet User Manual (UNDER CONSTRUCTION)

CoCoCoNet is a simple to use web server that allows the user to build and view and analyze co-expression networks without having to input experimental data. The data used is curated by aggregating the expression reads of RNA-seq experiments obtained through the NCBI database and is thus a better representation of the average gene co-expression than any individual experiment. 

To use CoCoCoNet, simply input a list of genes to be used in the construction of the network. Given a list of genes, you can then select specific parameters which are described below, generate the results and view the distribution of co-expression values as well as the network. 



## The Parameters
![param](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/parameters.png)

**The Species**
  
  Select the species that you wish to match your input gene list to. 
  
**Priority vs Meta**
  
  Priority will only use a subset of the gene set filtered on expresionality across experiments.
  
  Meta on the otherhand will build the network using all genes in the gene set.
  

**Using genes ...**
  
  "That I provide only" will construct a network using only selected genes
  "That I provide plus more" will search for additional genes that have a co-expression value greater than the threshold for at least one gene in the input list.
  
**Select Input method**
  
  Allows you to choose how to load in a gene list. Currently, you can select genes from a drop down list, paste a comma separated list, or upload a file with genes listed in new lines. 
  
## Generating the results
![yeast](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/yeast.png)


Once your genes have been properly loaded, the "Generate Results" and "Clear Results" options will appear. Selecting "Generate Results" will display the distribution of Co-expression values, the network and a sliding threshold bar that allows the user to display co-expression values greater than this threshold. The user also has the option to highlight genes along with their nearest neighbor or highlight genes with specified GO terms. Note that the GO terms listed are deemed overrepresented in the network by a gene set enrichment analysis. Users can also zoom in and out of the network to view gene symobls and connection weights. 

## Ortholog Mapped genes
![C_elegans](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/C_elegans.png)

Expanding the next section using the (+) symbol allows the user to input a second species to compare to the first. After selection, "Generate" and "Clear" buttons will appear where selecting "Generate" will again display the co-expression value distribution, the network and a sliding threshold bar. Again the user has the option to highlight genes along with its nearest neighbors or highlight genes by GO term.  

## Extending Guilt by Association by Degree (EGAD)
Expanding the final sectionusing the (+) symbbol gives the user the option to perform EGAD analysis on each of the above networks and visualize the distribution of AUROC  or AUPRC using either neighbor voting or node degreee. Also displayed is a summary of the top performing gene - GO term pair predicted using EGAD. 


## Downloadables
![download](https://github.com/johnlee4/CoCoCoNet/blob/master/figures/download.png)

You can download the data used to generate the 
