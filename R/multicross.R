#' @title Multisample generalization of Rosenbaum's crossmatch test
#' @description In this packcage, we present a framework inspired by Rosenbaum's crossmatch idea to tackle the nonparametric, multisample problem wherein one is concerned with testing the equality of K unknown multivariate probability distributions.
#' We implement two tests: the first is a multisample generalization of Rosenbaum's crossmatch (MCM), and the other further introduces a Malahnobis-type modification to the test (MMCM). 
#' @param data_list is list of multifeature matrices corresponding to the K different classes, so each element of the list is a matrix, for a total of K matrices. Each matrix contains observations as the rows and features as the columns
#' @param level is the level alpha for hypothesis testing
#' @return The p-value corresponding to rejection of the alternative, along with the decision of the hypothesis testing (Null being accepted versus rejected)
#' @export
#' @examples
#' # Simulation Example when the user wants to test whether K=3 multivariate distributions are equal:
#' X1 = MASS::mvrnorm(10,rep(0,4),diag(2,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' X2 = MASS::mvrnorm(10,rep(0,4),diag(1,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' X3 = MASS::mvrnorm(10,rep(0,4),diag(3,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' mcm(list(X1,X2,X3),0.05)
#' @import nbpMatching
#' @import MASS 
#' @import crossmatch
mcm <- function(data_list,level) ## Function 1
{
  nvec<-rep(0,length(data_list))
  apmat<-c()
  for(i in (1:length(data_list))) 
  {
    nvec[i]<-nrow(data_list[[i]])
    apmat<-rbind(apmat,data_list[[i]])
  }
  K<-length(nvec)
  nvec<-as.matrix(nvec)
  N=sum(nvec)
  nvecsq<-nvec%*%t(nvec)
  sumninj<-sum(nvecsq[upper.tri(nvecsq, diag = FALSE)])
  nnminus1<-nvec*(nvec-1)
  nnminus1vecsq<-nnminus1%*%t(nnminus1)
  sumnnminus1<-sum(nnminus1vecsq[upper.tri(nnminus1vecsq, diag=FALSE)])
  s1<-0
  if(K >=3)
  {
    for(i in 1:(K-2))
    {
      for(j in (i+1):(K-1))
      {
        for(k in (j+1):(K))
        {
          s1<-s1+((nvec[i])*(nvec[j])*(nvec[k])*(nvec[i]-1))+((nvec[i])*(nvec[j])*(nvec[k])*(nvec[j]-1))+((nvec[i])*(nvec[j])*(nvec[k])*(nvec[k]-1))
        }
      }
    }
  }
  s2<-0
  if(K>=4)
  {
    for(i in 1:(K-3))
    {
      for(j in (i+1):(K-2))
      {
        for(k in (j+1):(K-1))
        {
          for(l in (k+1):K)
          {
            s2<-s2+((nvec[i])*(nvec[j])*(nvec[k])*(nvec[l]))
          }
        }
      }
    }
  }
  nullmean<-sumninj/(N-1)
  nullvar<-(1/((N-1)*(N-3)))*(sumnnminus1 + (2*s1) + (6*s2)) - (((N-2)/(N*(N-1)^2))*sumninj^2) + ((sumninj/(N-1))*(1-(2*sumninj/(N^2 - N))))
  #print(c(nullmean,nullvar,(nvec[1]*nvec[2])/(N-1),(2*nvec[2]*(nvec[2]-1)*nvec[1]*(nvec[1]-1))/((N-3)*(N-1)^2)))
  smatch <- as.matrix(nonbimatch(distancematrix(as.matrix(stats::dist(apmat, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))))$matches)
  multcm<-0
  cs<-c(0,cumsum(nvec))
  for(k in 1:K)
  {
    for(j in (cs[k]+1):(cs[k+1]))
    {
      multcm<-multcm+((as.numeric(smatch[j,4])<=cs[k])||(as.numeric(smatch[j,4])>cs[k+1]))
    }
  }
  multcm<-multcm/2
  #print(multcm)
  multstat<-(multcm-nullmean)/sqrt(nullvar)
  lowerpval <- stats::pnorm(multstat)
  dec<-noquote('Accept')
  if(lowerpval<level)
  {
    dec<-noquote('Reject')
  }
  #return(c(multcm,nullmean,sqrt(nullvar)))
  return(noquote(c(lowerpval,dec)))
}

#' Creates the null covariance matrix for mmcm, corresponding to the scenario when all K distributions are the same
#' @param nvec is a vector containing the sizes of the K different classes
#' @return The inputs for the Multisample Mahalanobis Crossmatch Test
#' @import nbpMatching
#' @import MASS
#' @import crossmatch
mhcccreate<-function(nvec)
{
  k<-length(nvec)
  n<-sum(nvec)
  mu1<-rep(0,k)
  sig1<-matrix(0,k,k)
  for(i in 1:k)
  {
    mu1[i]<-(nvec[i]*(nvec[i]-1))/(2*(n-1))
    sig1[i,i]<-(((nvec[i]*(nvec[i]-1))/(2*(n-1)))*(1-((nvec[i]*(nvec[i]-1))/(2*(n-1)))))+((nvec[i]*(nvec[i]-1)*(nvec[i]-2)*(nvec[i]-3))/(4*(n-1)*(n-3)))
  }
  for(i in 1:k)
  {
    for(j in setdiff(1:k,i))
    {
      sig1[i,j]<-((nvec[i])*(nvec[j])*(nvec[i]-1)*(nvec[j]-1))/(2*((n-1)^2)*(n-3))
    }
  }
  mu<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      mu[i,j]<-(((nvec[i])*(nvec[j]))/(n-1))
    }
  }
  muv<-t(mu)[lower.tri(t(mu))]
  bigsig<-matrix(0,(((k^2)-k)/2),(((k^2)-k)/2))
  for(i in 1:(k-1))
  {
    for(j in (i+1):k)
    {
      bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((i-1)*k)-((i*(i-1))/2)+j-i)]<- ((nvec[i]*nvec[j]*(nvec[i]-1)*(nvec[j]-1))/((n-1)*(n-3)))+(((nvec[i]*nvec[j])/(n-1))*(1-((nvec[i]*nvec[j])/(n-1))))
    }
  }
  if(k>=3)
  {
    for(i in 1:(k-2))
    {
      for(j in (i+1):(k-1))
      {
        for(l in (j+1):k)
        {
          bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-((nvec[i]*(nvec[i]-1)*nvec[j]*nvec[l])/((n-1)*(n-3))) - ((((nvec[i])^2)*nvec[j]*nvec[l])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-((nvec[i]*(nvec[i]-1)*nvec[j]*nvec[l])/((n-1)*(n-3))) - ((((nvec[i])^2)*nvec[j]*nvec[l])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-((nvec[j]*(nvec[j]-1)*nvec[i]*nvec[l])/((n-1)*(n-3))) - ((((nvec[j])^2)*nvec[i]*nvec[l])/((n-1)^2))
          bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-((nvec[j]*(nvec[j]-1)*nvec[i]*nvec[l])/((n-1)*(n-3))) - ((((nvec[j])^2)*nvec[i]*nvec[l])/((n-1)^2))
          bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-((nvec[l]*(nvec[l]-1)*nvec[j]*nvec[i])/((n-1)*(n-3))) - ((((nvec[l])^2)*nvec[j]*nvec[i])/((n-1)^2))
          bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-((nvec[l]*(nvec[l]-1)*nvec[j]*nvec[i])/((n-1)*(n-3))) - ((((nvec[l])^2)*nvec[j]*nvec[i])/((n-1)^2))
        }
      }
    }
  }
  if(k>=4)
  {
    for(i in 1:(k-3))
    {
      for(j in (i+1):(k-2))
      {
        for(l in (j+1):(k-1))
        {
          for(m in (l+1):k)
          {
            bigsig[(((i-1)*k)-((i*(i-1))/2)+j-i),(((l-1)*k)-((l*(l-1))/2)+m-l)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((i-1)*k)-((i*(i-1))/2)+l-i),(((j-1)*k)-((j*(j-1))/2)+m-j)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((i-1)*k)-((i*(i-1))/2)+m-i),(((j-1)*k)-((j*(j-1))/2)+l-j)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((j-1)*k)-((j*(j-1))/2)+l-j),(((i-1)*k)-((i*(i-1))/2)+m-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((j-1)*k)-((j*(j-1))/2)+m-j),(((i-1)*k)-((i*(i-1))/2)+l-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
            bigsig[(((l-1)*k)-((l*(l-1))/2)+m-l),(((i-1)*k)-((i*(i-1))/2)+j-i)]<-(2*nvec[i]*nvec[j]*nvec[l]*nvec[m])/(((n-1)^2)*(n-3))
          }
        }
      }
    }
  }
  return(list(as.matrix(mu1),sig1,as.matrix(muv),bigsig))
}

#' Split a data frame or matrix into subsets based on a particular categorical variable
#' @param obj is a data frame or matrix to be split into subsets, divided by the categorical variable
#' @param by is a character-string that specifies the columns that need to be subsetted 
#' @return A list containing the subsetted data sets. The names of the list corresponds to the value of the subsetted list
split_mat <- function (obj, by) {
  # Get the set of possible values
  column.levels <-if (is.factor(obj[, by])) {
    levels(obj[, by])
  } else {
    unique(obj[, by])
  }
  # A list used to store each individual data.frame
  res <- list()
  # Iterate through all possible values and store each subset in a separate
  # entry in the list
  for (val in column.levels) {
    # Determine which rows match this value
    hits <- obj[, by] == val
    # Store data set temporarily in a local value
    data.set <- obj[hits, ]
    # Assign levels to the column. This adds levels to string data.
    levels(data.set[, by]) <- column.levels
    # Store data set in list
    res[[val]] <- data.set
  }
  
  # Return list
  res
}

#' This function takes an scRNAseq counts matrix as an input (cells in rows x genes in cols) and outputs a matrix of cells x 3 covariates (number of genes detected, sequencing depth and mitochondrial gene expression). This covariate matrix can then be used to match cells and perform stratified permutation.
#' @param datamat is cells x genes matrix from an scRNAseq experiment. Could be non-UMI or UMI counts.
#' @return A list of scRNAseq matrices subsetted by KEGG Pathways, so list[[n]] corresponds to the n^th pathway.
#' @export
#' @import Seurat
#' @importFrom  Matrix rowSums
matchpar <- function(datamat) {
  library_size <- rowSums(datamat)
  genes_detected <- apply(datamat, 1, function(x) length(x[x>0]))
  mito.genes = grep(pattern = "^MT-", x = colnames(x = datamat), value = TRUE)
  percent.mito <- Matrix::rowSums(datamat[,mito.genes])/Matrix::rowSums(datamat)
  d1 <- cbind(library_size,genes_detected)
  d2 <- cbind(d1,percent.mito)
  return(d2)
}

#' This function takes an scRNAseq counts matrix as an input (cells in rows x genes in cols) and outputs a list of cells x pathway matrices
#' @param datamat is cells x genes matrix from an scRNAseq experiment. Could be non-UMI or UMI counts.
#' @param orgsm specifies whether the species is mouse ("Mm"), human ("Hs) or C. elegans ("Ce").
#' @return A list of scRNAseq matrices subsetted by KEGG Pathways, so list[[n]] corresponds to the n^th pathway.
#' @export
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Ce.eg.db
#' @importFrom  AnnotationDbi mappedkeys
scPath <- function(datamat, orgsm) {
  if(orgsm=="Hs")  { 
    kegg <- org.Hs.egPATH2EG 
    x <- org.Hs.egSYMBOL
  }
  if(orgsm=="Mm")  { 
    kegg <- org.Mm.egPATH2EG
    x <- org.Mm.egSYMBOL
  }
  if(orgsm=="Ce")  { 
    kegg <- org.Ce.egPATH2EG
    x <- org.Ce.egSYMBOL
  }
  mapped <- AnnotationDbi::mappedkeys(kegg)
  kegg2 <- as.list(kegg[mapped])
  ## kegg2 is a list where each [] index is a pathway and [[]] indices are genes, all in numeric/ID form
  KEGGList <- matrix(as.numeric(unlist(kegg2))) #Extract a List of all genes (>16K) which are functionally annotated in KEGG
  sums <- c()
  for (i in 1:length(kegg2)) { sums[i] <- length(kegg2[[i]])}
  pathnames <- c()
  for (i in 1:length(kegg2)) {
    pathnames[i] <- list(rep((names(kegg2)[i]),sums[i]))
  }
  pathlist <- matrix(unlist(pathnames))
  
  KEGGMap <- as.data.frame(cbind(pathlist,KEGGList))
  colnames(KEGGMap) <- c("Pathway ID","Gene Entrez ID")
  
  mapped_genes <- AnnotationDbi::mappedkeys(x)
  x2 <- as.list(x[mapped_genes])
  entrezID <- names(x2)
  kegg_gene <- matrix(unlist(x2))
  RefID <- as.data.frame(cbind(entrezID,kegg_gene))
  ##RefID has two cols, first with Entrez IDs, second with the corresponding gene names
  match_genes <- RefID[match(KEGGMap$`Gene Entrez ID`,RefID$entrezID),2]
  KEGGMap2 <- cbind(KEGGMap,match_genes) #KEGGMap2 contains Pathway ID, Gene Entrez ID, Gene Name
  
  Gene_Path <- split_mat(KEGGMap2,"Pathway ID")
  #Gene_Path is a sub-setted list of KEGGMap 2 with length 229 (# of pathways)
  ## Now should be able to subset sc Matrix using Gene_Path[[i]]$match_genes
  data_bypath <- c()
  for (i in 1:length(Gene_Path)) {
    data_bypath[[i]] <- datamat[,intersect(colnames(datamat),Gene_Path[[i]]$match_genes)]
  }
  return(data_bypath)
}


#' Calculates the pairwise crosscounts for the K classes being examined
#' @param nvec is a vector containing the sizes of the K different classes
#' @param apmat is the data matrix containing pooled data from each of the K classess, on which optimal non-bipartite matching is performed.
#' @return The inputs for the Multisample Mahalanobis Crossmatch Tests
#' @import nbpMatching
#' @import MASS
#' @import crossmatch
mhccexecutelong<-function(nvec, apmat)
{
  k<-length(nvec)
  n<-sum(nvec)
  smatch <- as.matrix(nonbimatch(distancematrix(as.matrix(stats::dist(apmat, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))))$matches)
  multcm<-rep(0,k)
  cs<-c(0,cumsum(nvec))
  for(i in 1:k)
  {
    for(j in (cs[i]+1):(cs[i+1]))
    {
      multcm[i]<-multcm[i]+((as.numeric(smatch[j,4])>cs[i])&&(as.numeric(smatch[j,4])<=cs[i+1]))
    }
  }
  multcm<-multcm/2
  A<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      for(l in (cs[i]+1):(cs[i+1]))
      {
        A[i,j]<-A[i,j]+((as.numeric(smatch[l,4])>cs[j])&&(as.numeric(smatch[l,4])<=cs[j+1]))
      }
    }
  }
  av<-t(A)[lower.tri(t(A))]
  return(list(as.matrix(multcm),as.matrix(av)))
}

#' Use the Mahalnobis-type multisample test based on optimal matching to compare K different multivariate distributions
#' @param data_list is list of multifeature matrices corresponding to the K different classes, so each element of the list is a matrix, for a total of K matrices.
#' @param level is the cutoff value (alpha) for hypothesis testing
#' @return The p-value corresponding to rejection of the alternative, along with the decision of the hypothesis testing (Null being accepted versus rejected)
#' @export
#' @examples
#' # Simulation Example when the user wants to test whether K=3 multivariate distributions are equal:
#' X1 = MASS::mvrnorm(10,rep(0,4),diag(2,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' X2 = MASS::mvrnorm(10,rep(0,4),diag(1,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' X3 = MASS::mvrnorm(10,rep(0,4),diag(3,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' mmcm(list(X1,X2,X3), 0.05)
#' @import nbpMatching
#' @import MASS
#' @import crossmatch
mmcm<-function(data_list,level)
{
  nvec<-rep(0,length(data_list))
  apmat<-c()
  for(i in (1:length(data_list))) 
  {
    nvec[i]<-nrow(data_list[[i]])
    apmat<-rbind(apmat,data_list[[i]])
  }
  ll<-mhcccreate(nvec)
  mu1<-ll[[1]]
  sig1<-ll[[2]]
  muv<-ll[[3]]
  bigsig<-ll[[4]]
  n<-sum(nvec)
  k<-length(nvec)
  lll<-mhccexecutelong(nvec,apmat)
  multcm<-lll[[1]]
  av<-lll[[2]]
  stbig<-t(as.matrix(av)-as.matrix(muv))%*%solve(bigsig)%*%(as.matrix(av)-as.matrix(muv))
  upperpval2 <- 1- stats::pchisq(stbig,(((k^2)-k)/2))
  decmhcc2<-noquote('Accept')
  if(upperpval2<level)
  {
    decmhcc2<-noquote('Reject')
  }
  return(noquote(c(upperpval2,decmhcc2)))
}

#' Given two input matrices with the same number of observations but differrent number of variables, this function returns the largest canonical correlation between variables of matrix 1 (X) and those of matrix 2 (Y).
#' @param X is a data matrix with observations as rows and features/genes as columns.
#' @param Y is a data matrix with observations as rows (same observations as X) and a different set of features/genes as columns.
#' @return The largest canonical correlation between X and Y as described in <https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Canonical_Correlation.pdf>. 
#' @export
#' @import MASS
#' @importFrom  stats cor
multigene <- function(X,Y) {
  X <- X[,colSums(X != 0) != 0]
  Y <- Y[,colSums(Y != 0) != 0]
  Rxx = stats::cor(X,X, method="spearman")
  Ryy = stats::cor(Y,Y, method="spearman")
  Rxy = stats::cor(X,Y, method="spearman")
  Ryx = stats::cor(Y,X, method="spearman")
  C = MASS::ginv(Ryy) %*% Ryx %*% MASS::ginv(Rxx) %*% Rxy
  S = svd(C)
  corcan <- sqrt(S$d[1])
  return (corcan)
}

#' When the MCM/MMCM tests reject the null, class selection can help determine which of the K classes are the likely contributors for rejection
#' @param data_list is list of multifeature matrices corresponding to the K different classes, so each element of the list is a matrix, for a total of K matrices.
#' @param level is the cutoff value (alpha) for hypothesis testing
#' @return A table of pairwise comparisons among the K classes, to further probe which class influences the rejection of the null the most. No p-value adjustment is made to these reported p-values
#' @export
#' @examples
#' # Simulation Example when the user wants to test whether K=3 multivariate distributions are equal:
#' X1 = MASS::mvrnorm(10,rep(0,4),diag(2,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' X2 = MASS::mvrnorm(10,rep(0,4),diag(1,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' X3 = MASS::mvrnorm(10,rep(0,4),diag(3,4),tol=1e-6, empirical=FALSE, EISPACK=FALSE)
#' select_class(list(X1,X2,X3), 0.05)
#' @import nbpMatching
#' @import MASS
#' @import crossmatch
select_class<-function(data_list,level)
{
  nvec<-rep(0,length(data_list))
  apmat<-c()
  for(i in (1:length(data_list))) 
  {
    nvec[i]<-nrow(data_list[[i]])
    apmat<-rbind(apmat,data_list[[i]])
  }
  lll<-mhccexecutelong(nvec,apmat)
  av<-lll[[2]]
  k<-length(nvec)
  n<-sum(nvec)
  mu<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in i:k)
    {
      mu[i,j]<-(((nvec[i])*(nvec[j]))/(n-1))
    }
  }
  muv<-t(mu)[lower.tri(t(mu))]
  sig<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in 1:k)
    {
      sig[i,j]<-((nvec[i]*nvec[j]*(nvec[i]-1)*(nvec[j]-1))/((n-1)*(n-3)))+(((nvec[i]*nvec[j])/(n-1))*(1-((nvec[i]*nvec[j])/(n-1))))
    }
  }
  sigv<-t(sig)[lower.tri(t(sig))]
  dvec<-(av-muv)/sqrt(sigv)
  pval<-rep(0,length(dvec))
  decvec<-rep('Accept',length(pval))
  for(i in 1:length(pval))
  {
    pval[i] <- stats::pnorm(dvec[i])
    if(pval[i]<level)
    {
      decvec[i]<-noquote('Reject')
    }
  }
  mt<-c()
  for(i in 1:(k-1))
  {
    for(j in (i+1):k)
    {
      mt<-rbind(mt,c(i,j))
    }
  }
  return(noquote(cbind(mt,as.matrix(pval),as.matrix(decvec))))
}
