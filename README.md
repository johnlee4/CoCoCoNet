# OrthoCoNet User Manual

OrthoCoNet is a simple to use web server that allows the user to build and view RNA-seq coexpression networks without having to input experimental data. The data used is currated by aggregating the expression reads of RNA-seq experiments obtained through NCBI and GEMMA (?) and is thus a better representation of the average gene co-expression than any individual experiment. 

To use OrthoCoNet, simply input a list of genes to be used in the construction of the network. Given a list of genes, you can then select specific parameters which are described below, generate the results and view the distribution of co-expression values as well as the network. 

## Parameters
**The Species**
  
  Select the species that you wish to match your input gene list to. 
  
**Priority vs Meta**
  
  Priority will search for genes that appear in more than 50% of the experiments used to generate the aggregate. This ensures that expression levels are similarly powered.
  
  Meta on the otherhand will build the network for all genes in the genome. 
  

**Using genes ...**
  
  "That I provide only" will construct a network using only selected genes
  "That I provide plus more" will search for additional genes that have a co-expression value greater than the threshold for at least one gene in the input list.
  
**Select Input method**
  
  Allows you to choose how to load in a gene list. Currently, you can select genes from a drop down list, paste a comma separated list, or upload a file with genes listed in new lines. 
  
## Generating the results
Once your genes have been properly loaded, 
