# load minimal dependencies
library(EGAD)
library(rhdf5)
library(Matrix):


getEnsemblID <- function(species, gene_list){
# Input gene_list is allowed to be anything. EnsemblID , entrezID, genesymbol, or any synonyms
# input should be an array/vector
# GeneInfo file is read from the ftp - change to read from where you store it.
# outputs
# [[1]] ensemble genes
# [[2]] input genes not found
# [[3]] original gene set

    # get gene info 
# adjust file location as needed
####################################################################################################################################
    geneInfo = read.delim(paste0("/ftp/data/geneInfo/",species,"_geneInfo.tab"), stringsAsFactors=F, sep="\t")
####################################################################################################################################

    gene_list = as.character(gene_list)


    # search for genes in EnsemblID
    m               = match(toupper(gene_list) , toupper(geneInfo$EnsemblID) ) 
    ind_found       = m[!is.na(m)]
        ind_found   = ind_found[! duplicated(ind_found )]   # remove duplicates. 

    ensemblGenes    = geneInfo$EnsemblID[ind_found]
    genes_found     = gene_list [ !is.na(m)  ] 
    genes_not_found = gene_list [  is.na(m)  ] 

    # search for genes in symbols 
    m               = match(toupper(genes_not_found), toupper(geneInfo$GeneSymbol) )
    ind_found       = m[!is.na(m)]
        ind_found   = ind_found[! duplicated(ind_found )]   # remove duplicates. 

    ensemblGenes    = c(ensemblGenes , geneInfo$EnsemblID[ind_found])

    genes_found     = c(genes_found , genes_not_found [ !is.na(m)  ]  )
    genes_not_found = genes_not_found [  is.na(m)  ] 


    # search for genes in EntrezID
    m               = match(genes_not_found, geneInfo$EntrezID)
    ind_found       = m[!is.na(m)]
        ind_found   = ind_found[! duplicated(ind_found )]   # remove duplicates. 
    
    ensemblGenes    = c(ensemblGenes , geneInfo$EnsemblID[ind_found])
    
    genes_found     = c(genes_found , genes_not_found [ !is.na(m) ] )
    genes_not_found = genes_not_found [  is.na(m)  ] 

    # search in synonyms
    temp       = paste0("\\|", genes_not_found ,"\\|")
    ind_found       = lapply(temp , grep , geneInfo$Synonyms,ignore.case=TRUE )
    
    lens            = lengths(ind_found)

    ind_found       = unlist(ind_found)
        ind_found   = ind_found[! duplicated(ind_found )]   # remove duplicates. 

    ensemblGenes    = c(ensemblGenes , geneInfo$EnsemblID[ind_found])
    
    genes_found     = c(genes_found , rep(genes_not_found, lens ) )
    genes_not_found = genes_not_found[ lens==0]  

    # tidy up
    genes_found     = genes_found[!is.na(genes_found)]
    noENS           = genes_found[ is.na(ensemblGenes)  ]
    genes_not_found = c(genes_not_found,noENS)
    ensemblGenes    = ensemblGenes[! duplicated(ensemblGenes)]

    
    N = length(ensemblGenes)
    if (N == 0 ){
        print("None of your input genes were found or cannot be converted to the correct format!")
    }
    
    count = length(genes_not_found)
    if (count > 1 ) {
        
        if (count > 5)  {
            temp = paste(genes_not_found[1:5],collapse=", ")
            temp = paste(temp ,"\n and",  count -5 , "others." )
            print(paste("The following input genes were not found! \n ",temp) )
        } else {
            temp = paste(genes_not_found,collapse=", ")
            print(paste("The following input genes were not found! \n ",temp) )
        }
    }


    if (N > 1000 ){
        
        ind = sample(1:length(ensemblGenes) , 1000 )
        ensemblGenes = ensemblGenes[ind]
        genes_found  = genes_found[ind]
        print(paste("We found",N,"genes. Requesting too many genes will slow our server down for other users - so we took a random sample of 1000 genes instead." ) )
    }

    # ensemblGenes should be used to match to the data, 
    # genes_not_found should be reported to the user
    # genes_found should be used as display names
    return(list(ensemblGenes, genes_not_found, make.unique(genes_found) ) )  
}

getGOA <- function(species){
# reads in the hdf5 file related genes to GO terms and stores as a sparse matrix

# adjust file location as needed
########################################################################################    
    filename = paste0("/ftp/data/gene2go/",species,"_gene2go.hdf5" )
########################################################################################
    # read in gene2go matrices (hdf5) in COO format
    go      = h5read(filename, name ="GO")
    genes   = h5read(filename, name ="genes")
    ind     = h5read(filename, name ="ind")
    
    h5closeAll()

    # store as a sparsematrix
    GO_annots       = sparseMatrix(i=ind[,1] ,j= ind[,2] ,dimnames=list(toupper(genes),toupper(go)))
    csums = colSums( GO_annots)
    GO_annots       = GO_annots[, (csums <=1000)&(csums>=10) ]
    return(GO_annots)
}



geneSetScore <- function( species, gene_list ,flag1= "prio",op = "AUROC"){
    # get score of geneset using EGAD
    # apply egad to all genes 
    # takes in the ensemblIDs
    filename = paste0("/ftp/data/networks/" ,species,"_",flag1, "AggNet.hdf5")

    
    row = h5read(filename, name= "row"  )
    col = h5read(filename, name= "col"  )
 
    # just in case, make sure it's symmetric
    genes = intersect(toupper(row) , toupper(col))

    ind = match(gene_list, genes)
    ind = ind[!is.na(ind)]
    network = h5read(filename, name= "agg" )
    rownames(network) = toupper(row)
    colnames(network) = toupper(col)
    h5closeAll()

    network=network[genes,genes]

    labels = matrix(0, length(genes) , 1)
    rownames(labels ) = genes
    colnames(labels ) = "Gene Set"

    labels[ rownames(labels) %in%  toupper(gene_list),1] =1



    roc <- neighbor_voting(labels, network, 3, output =op)
    #   preds <- predictions(labels, network.sub)

    return(roc)
  
}

#
networkMaker <- function(species,gene_list,flag1="prio"){
    # genes from file are listed as ENSEMBL type - need to make sure to convert input genes to this format first.


    # get indices of genes in hdf5 file
    filename        = paste0("/ftp/data/networks/",species,"_",flag1,"AggNet.hdf5")
    allGenes        = h5read(filename, name= "row")

    genesInDB       = intersect(gene_list , allGenes )
    genesNotInDB    = setdiff(gene_list, allGenes)
    ind             = match(genesInDB , allGenes)

    net             = h5read(filename, name = "agg" ,index=list(ind,ind))
    rownames(net)   = genesInDB
    colnames(net)   = genesInDB

    # symmetric part removed - along with diagonal which is always 1.0

    net[ lower.tri(net,diag=TRUE)] = 0 
    ind             = which(net > 0, arr.ind = TRUE )

    df              = data.frame(
                        from = rownames(net)[ind[,1]],
                        to = colnames(net)[ind[,2]] ,
                        value = net[ind],
                        stringsAsFactors=F
                        )

    h5closeAll()


    
    return ( list(df,net ))
}


gene_set_enrichment <- function(gene_list, species){
## Perform gene set enrichment to look for over represented GO terms in the sample (gene_list) 
## Sample should be preprocessed in function <filterGeneList>
#   Input:  gene_list       - User defined gene list to consider
#           GOA             - ngCMatrices of species gene x GO 
#           voc             - Data frame containing <GO term, description, cc/bp/mf> global variable
#   Output: results         - Data frame containing 


    GOA = getGOA(species)

    genes.names     = rownames(GOA)            
    labels.names    = colnames(GOA)
    genes.counts    = rowSums(GOA)
    labels.counts   = colSums(GOA)      # p - gives number of genes in a GO group
    # look only at genes with known annotations
    m               = match ( toupper(gene_list), toupper(genes.names) )
    filt.genes      = !is.na(m)
    filt.labels     = m[filt.genes]


    labels.counts.set   = rep( sum(filt.genes), length(labels.counts) )  # q - for broadcasting later.  # of genes involved 
    m               = match (labels.names,  voc[,1] )
    v.f             = !is.na(m)
    v.g             = m[v.f]

    if(  length(filt.labels) == 1 )  { genes.counts.set = GOA[filt.labels,] }
    else { genes.counts.set = colSums(GOA[filt.labels,]) }            
    # Fisher's Exact test
    test =  cbind( genes.counts.set  , labels.counts, dim(GOA)[1]-labels.counts, labels.counts.set)
    # computes P[X > x] . Probabilities are given as log(p)
    pvals = phyper(test[,1]-1, test[,2], test[,3], test[,4], lower.tail=F)   
    # FDR Correction
    pvals.adj = p.adjust( pvals, method="BH")
    sigs = pvals.adj < ( 0.05/length(pvals.adj) )

    # Downloadable
    results = data.frame(   GO_term     = voc[v.g,1],
                            description = as.character(voc[v.g,2]),
                            N_sample    = test[v.f,1],
                            N_univ      = test[v.f,2],
                            pvals       = pvals[v.f],
                            adj_pvals   = pvals.adj[v.f],
                            sig         = sigs[v.f],
                            stringsAsFactors = FALSE
                            )


    # sort by pvals 
    results=results[order(results$pvals), ]

    return ( remove.factors(results) )
}

getMoreGenes <- function(gene_list, species , N=10, flag1="prio") {
# gene list needs to be in the proper format.
# Find genes not in current set of genes with coexpression greater than thrs for at least one gene in set
# Input:
# Output:   new_genes

    temp = toupper(gene_list)

    # load the data
    filename        = paste0("/ftp/data/networks/",species,"_", flag1, "AggNet.hdf5")

    # read in gene names 
    allGenes        = toupper(h5read(filename,name = "col")  )
    
    genesInSet      = intersect(temp, allGenes)
    genesNotInSet   = setdiff( allGenes,temp  )

    # if (length(genesInSet) == 0) { shinyalert("No input genes were found - unable to add any genes. ",type = "error",closeOnClickOutside= TRUE ) ; return(gene_list ) }
    ind1            = match(genesInSet    , allGenes)   # use these as the rows
    ind2            = match(genesNotInSet , allGenes)   # use these as the cols


    # Get total node degree of genes not in set 
    coexpMatrix     = h5read(filename, name="agg" ,index = list(ind1 , 1:length(allGenes) )  )
    coexpMatrix     = coexpMatrix[,ind2]
    coexpMatrix     = matrix(coexpMatrix,nrow = length(ind1), ncol=length(ind2)  )
    rownames(coexpMatrix ) = allGenes[ind1]
    colnames(coexpMatrix ) = allGenes[ind2]



    internalDeg = colSums(coexpMatrix)
    # totalDeg    = read.delim(paste0("/ftp/data/nodeDegree/", species,"_", flag1,"Degree.csv") , sep = ",", stringsAsFactors=F)
    # rownames(totalDeg) = toupper(rownames(totalDeg))

    # totalDeg    = totalDeg[names(internalDeg),] 
    ind = rank(-internalDeg  ,ties.method="min")  # low is good
    # internalDeg/totalDeg
    new_genes = allGenes[ind2][ind <= N]


    h5closeAll()
    return(toupper(c(gene_list, new_genes )))

}




orthologMapper <- function(genes,speciesA, speciesB){
##  Maps genes in speciesA to orthologous genes in speciesB
#   Input:  genes           - genes of SpeciesA
#           speciesA        - name of reference species
#           speciesB        - name of target species
#   Output: temp            - orthologous genes in speciesB

    # get gene info tables
    geneInfoA = read.delim(paste0("/ftp/data/geneInfo/",speciesA,"_geneInfo.tab"), stringsAsFactors=F, sep="\t")
    geneInfoB = read.delim(paste0("/ftp/data/geneInfo/",speciesB,"_geneInfo.tab"), stringsAsFactors=F, sep="\t")

    # look up the OrthogeneID for the first species
    reforthogenes   = geneInfoA[geneInfoA$EnsemblID %in% intersect(genes , geneInfoA$EnsemblID )  , 4 ]
    reforthogenes   = reforthogenes[! is.na(reforthogenes)]
    reforthogenes   = reforthogenes[! duplicated(reforthogenes)]


    filename1 = paste0('/ftp/data/orthologData/',speciesA,'_',speciesB,'_ortho.csv')
    filename2 = paste0('/ftp/data/orthologData/',speciesB,'_',speciesA,'_ortho.csv')

    # find target ortholog genes matching to reforthogenes from ortholog matrix
    if(file.exists(filename1) ){
        odata   = read.delim(filename1,sep=',',stringsAsFactors=F)
    } else if(file.exists(filename2 )){
        odata   = read.delim(filename2,sep=',',stringsAsFactors=F)
    }


    if (length(intersect(as.vector(odata[,1]),reforthogenes))>0){
        ind = match(reforthogenes,odata[,1])
        tarorthogenes <- as.vector(odata[ ind,2])
    } else if (length(intersect(as.vector(odata[,2]),reforthogenes))>0) {
        tarorthogenes <- as.vector(odata[match(reforthogenes,odata[,2]),1])        
    } else{ shinyalert('No ortholog genes not found',type="error" ,closeOnClickOutside=TRUE)}

    tarorthogenes = tarorthogenes[!is.na(tarorthogenes)]
    ind         = match(tarorthogenes,geneInfoB$OrthoID)
    tarorthogenes   = geneInfoB$EnsemblID[ind]
    tarorthogenes   = tarorthogenes[!is.na(     tarorthogenes)]
    tarorthogenes   = tarorthogenes[!duplicated(tarorthogenes)]
    
    return(tarorthogenes)
}



GBA_driver <- function(network, labels , op = "AUROC" ,nfold = 3,label_min=10,label_max=1000) {
# takes in network as full matrix

    # which input genes have an annotaton?
    rownames(network ) = toupper(rownames( network))
    colnames(network ) = toupper(colnames( network))

    rownames(labels) = toupper(rownames(labels))
    colnames(labels) = toupper(colnames(labels))
    m <- match(rownames(network),rownames(labels) )
    f <- !is.na(m)
    g <- m[f]

    
    if (network[1,2] != network[2,1]){
        network.sub <- network[f, f] + t(network[f,f])
        diag(network.sub) = 1
    } else {network.sub <- network[f, f]}

    #   genes.labels = labels[g, sig_GOs]
    genes.labels <- filter_network_cols(labels, label_min, label_max) # filters the GO terms
    

    # genes.labels= genes.labels[g,]
    genes.labels = filter_network_cols(genes.labels[g,], 1, label_max )
    

    roc.sub <- neighbor_voting(genes.labels, network.sub, nfold, output =op)

    return(roc.sub)

}

getCompleteNetwork <- function(species, flag1="prio"){
    # store all of the co-expression data in a single matrix
    # WARNING: RAM intensive for some species

    filename = paste0("/ftp/data/networks/" ,species,"_",flag1, "AggNet.hdf5")

    agg = h5read(filename, name= "agg"  )
    col = h5read(filename, name= "col"  )
    row = h5read(filename, name= "row"  )

    rownames(agg) = row
    colnames(agg) = col

    return(agg)

}

################################################################################################################################
########################################################### Example ############################################################
################################################################################################################################

# possible species include
# human, mouse, arabidopsis, boar, chicken, cow, fruitfly, maize, rat, rice, roundworm, soybean, yeast, zebrafish

speciesA = "yeast"
speciesB = "roundworm"
gene_listA = read.delim("/ftp/data/sample_input/top_yeastGenes.csv",stringsAsFactors=F, header=F)

thresholdA = 0.99
thresholdB = 0.99


# convert these genes to ensembl
gene_listA = getEnsemblID(speciesA, gene_listA )[[1]]

# get more genes if desired
gene_listA = getMoreGenes(gene_listA, speciesA , N=100) 

# get corresponding ortholog genes.
gene_listB = orthologMapper(gene_listA,speciesA, speciesB)

# construct full networks for both ( no threshold )
netdatA = networkMaker(speciesA, gene_listA)
netdatB = networkMaker(speciesB, gene_listB)

# network in coordinate format - i , j , x 
netA_coo = netdatA[[1]]
netB_coo = netdatB[[1]]

# dense network as adjacency matrix (lower tri + diag removed)
netA_full = netdatA[[1]]
netB_full = netdatB[[1]]

# apply threshold to get subnetworks
subnetA = netA_coo[netA_coo[,3] > thresholdA ,]
subnetB = netB_coo[netB_coo[,3] > thresholdB ,]

# apply gene set enrichment on subnetwork
subGeneListA = unique(c(subnetA[,1 ] , subnetA[,2] )) 
subGeneListB = unique(c(subnetB[,1 ] , subnetB[,2] )) 

gseA = gene_set_enrichment(subGeneListA, speciesA)
gseB = gene_set_enrichment(subGeneListB, speciesB)

# assess network quality using EGAD
labelA = getGOA(speciesA)
labelB = getGOA(speciesB)
rocsA = GBA_driver(netA_full, labelA , op = "AUROC",nfold=3)  # op can be "AUROC" or "PR"
    print(mean(rocsA, na.rm = T))
rocsB = GBA_driver(netB_full, labelB , op = "AUROC",nfold=3)
    print(mean(rocsB, na.rm = T))

# get gene set score using EGAD 
############################## WARNING: RAM INTENSIVE FOR SOME SPECIES ##############################
labelA = matrix(1, length(gene_listA) ,1 )
    rownames(labelA) = gene_listA
labelB = matrix(1, length(gene_listB) ,1 )
    rownames(labelB) = gene_listB

aggA = getCompleteNetwork(speciesA)
aggB = getCompleteNetwork(speciesB)

setScoreA = GBA_driver( aggA, labelA ,op = "AUROC", nnfold=3)
setScoreB = GBA_driver( aggB, labelB ,op = "AUROC", nnfold=3)
######################################################################################################



