#' CreatScDNSobject
#'
#' @param counts a matrix represents gene count expression, with rows being genes and columns being cells
#' @param data a matrix represents gene normalized expression, with rows being genes and columns being cells.
#' @param Network a two-column data frame represents the biological network
#' @param GroupLabel a string vector representing the grouping information of the cells
#' @param k a numerical value represents the kth (def: 10) nearest neighbor point
#' @param n.grid An integer value representing the number of grid (def:60) used to calculate the joint probability density
#' @param n.coarse An integer value representing the number of grid (def:20) after coarse-graining
#' @param loop.size An integer value representing the network size (def:3000) allocated to each compute core
#' @param parallel.sz An integer representing the number of parallel cores (def:1)
#' @param verbose Logical value, whether progression messages should be printed in the terminal.
#' @param exclude.zero Logical value indicating whether to exclude 0 values when calculating the joint probability density (def: FALSE)
#' @param NoiseRemove Logical value, indicating whether to remove noise that may lead to extremely low probability density (def:TRUE)
#' @param CoarseGrain Logical value, whether to coarse-grain the grid(def: TRUE).
#' @param returnCDS A logical value indicating whether to return the conditional probability density
#' @param ds.method Method to calculate joint probability density, optional knn and kde (def:knn)
#' @param divergence Method for calculating divergence, optional jsd (Jensen-Shannon divergence) and kld (Kullback–Leibler divergence)(def:jsd)
#' @param h A parameter used for kernel density estimation, which only takes effect when the density estimation method is kde (def:1).
#' @param n.dropGene Integer, number of cells sampled for model building (def:3000)
#' @param n.randNet Integer, number of cells sampled for model building (def:3000)
#' @param sdBias Numeric value, bias coefficient used to penalize low degree genes(>=1) (def:1.1)
#' @param parllelModel parllelModel=c('foreach','bplapply')[1]
#' @param rb.jsd
#' @param Div_weight
#'
#' @return
#' @export
#'
#' @examples
CreatScDNSobject <- function(counts,
                             data = NULL,
                             Network = NULL,
                             GroupLabel = NULL,
                             k= 10,
                             n.grid = 60,
                             NoiseRemove = TRUE,
                             CoarseGrain = TRUE,
                             n.coarse=20,
                             ds.method = c('knn','kde')[1],
                             divergence = c('jsd','kld')[1],
                             noiseSd=0.01,
                             h = 1,
                             rb.jsd=F,
                             Div_weight=sqrt(as.vector(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)*t(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)))),
                             loop.size=3000,
                             parallel.sz = 3,
                             verbose = TRUE,
                             exclude.zero = FALSE,
                             returnCDS = TRUE,
                             n.dropGene = 3000,
                             n.randNet = 20000,
                             sdBias=1.1,
                             parllelModel=c('foreach','bplapply')[1]){
  #01
  if(is.null(Network)){
    data(scDNSBioNet)
    message('No network is provided, we load the internal network.')
    Network <- scDNSBioNet
  }
  if(is.null(data)){
    message('Data is not provided, we will normalize with NormalizeData functions.')
    data <- NormalizeData(counts)
  }
  #02
  Div.Parameters <- list(k= k,
                         n.grid = n.grid,
                         NoiseRemove = NoiseRemove,
                         CoarseGrain = CoarseGrain,
                         n.coarse = n.coarse,
                         ds.method = ds.method,
                         divergence = divergence,
                         h = h,
                         noiseSd=noiseSd,
                         Div_weight = Div_weight,
                         rb.jsd=rb.jsd,
                         loop.size=loop.size,
                         parallel.sz = parallel.sz,
                         verbose = verbose,
                         exclude.zero = exclude.zero,
                         parllelModel=parllelModel)

  NEA.Parameters <- list(do.impute = FALSE,
                         n.dropGene = n.dropGene,
                         n.randNet = n.randNet,
                         sdBias=sdBias)
  #03
  message('get dropout for each gene pairs.')
  Network = getDoubleDropout(Network = Network,counts = counts)
  # smooth OutLier value
  # message('Smooth OutLier value.')
  # data = SmoothOutLier_UP(data)
  #
  GeneVariability <- CalGeneInfAmount(NormalizeData(counts),
                   probs = 0.99,
                   rescale = FALSE,
                   plot = FALSE)
  gc()

  scDNSobjcet <- new(Class = "scDNS",
                     counts = counts,
                     data = data,
                     Network = Network,
                     GroupLabel = GroupLabel,
                     Zscore = data.frame(),
                     scZscore = matrix(),
                     GeneVariability= GeneVariability,
                     Div.Parameters = Div.Parameters,
                     NEA.Parameters = NEA.Parameters,
                     NEAModel=list(),
                     JDensity_A = matrix(),
                     JDensity_B = matrix(),
                     uniCase = character())
  scDNSobjcet


}



#' scDNS_1_CalDivs
#'
#' @param scDNSobjcet scDNSobjcet
#' @param k a numerical value represents the kth (def: 10) nearest neighbor point
#' @param n.grid An integer value representing the number of grid (def:60) used to calculate the joint probability density
#' @param n.coarse An integer value representing the number of grid (def:20) after coarse-graining
#' @param loop.size An integer value representing the network size (def:3000) allocated to each compute core
#' @param parallel.sz An integer representing the number of parallel cores (def:1)
#' @param verbose Logical value, whether progression messages should be printed in the terminal.
#' @param exclude.zero Logical value indicating whether to exclude 0 values when calculating the joint probability density (def: FALSE)
#' @param NoiseRemove Logical value, indicating whether to remove noise that may lead to extremely low probability density (def:TRUE)
#' @param CoarseGrain Logical value, whether to coarse-grain the grid(def: TRUE).
#' @param returnCDS A logical value indicating whether to return the conditional probability density
#' @param ds.method Method to calculate joint probability density, optional knn and kde (def:knn)
#' @param divergence Method for calculating divergence, optional jsd (Jensen-Shannon divergence) and kld (Kullback–Leibler divergence)(def:jsd)
#' @param h A parameter used for kernel density estimation, which only takes effect when the density estimation method is kde (def:1).
#' @param returnCDS Logical value, whether to return conditional density matrix (def: FALSE)
#' @param parllelModel
#'
#' @return
#' scDNSobjcet
#' @export
#'
#' @examples
scDNS_1_CalDivs <- function(scDNSobjcet,
                            k = NULL,
                            n.grid = NULL,
                            NoiseRemove = NULL,
                            CoarseGrain = NULL,
                            n.coarse = NULL,
                            ds.method = NULL,
                            divergence = NULL,
                            h = NULL,
                            noiseSd=NULL,
                            loop.size = NULL,
                            parallel.sz = NULL,
                            verbose = NULL,
                            exclude.zero = NULL,
                            returnCDS = NULL,
                            parllelModel = NULL){
  if(!is.null(k)){scDNSobjcet@Div.Parameters$k=k}
  if(!is.null(n.grid)){scDNSobjcet@Div.Parameters$n.grid=n.grid}
  if(!is.null(NoiseRemove)){scDNSobjcet@Div.Parameters$NoiseRemove=NoiseRemove}
  if(!is.null(CoarseGrain)){scDNSobjcet@Div.Parameters$CoarseGrain=CoarseGrain}
  if(!is.null(n.coarse)){scDNSobjcet@Div.Parameters$n.coarse=n.coarse}
  if(!is.null(h)){scDNSobjcet@Div.Parameters$h=h}
  if(!is.null(loop.size)){scDNSobjcet@Div.Parameters$loop.size=loop.size}
  if(!is.null(parallel.sz)){scDNSobjcet@Div.Parameters$parallel.sz=parallel.sz}
  if(!is.null(verbose)){scDNSobjcet@Div.Parameters$verbose=verbose}
  if(!is.null(exclude.zero)){scDNSobjcet@Div.Parameters$exclude.zero=exclude.zero}
  if(!is.null(returnCDS)){scDNSobjcet@Div.Parameters$returnCDS=returnCDS}
  if(!is.null(parllelModel)){scDNSobjcet@Div.Parameters$parllelModel=parllelModel}

  # Div.Parameters <- list(k= k,
  #                        n.grid = n.grid,
  #                        NoiseRemove = NoiseRemove,
  #                        CoarseGrain = CoarseGrain,
  #                        n.coarse = n.coarse,
  #                        ds.method = ds.method,
  #                        divergence = divergence,
  #                        h = h,
  #                        loop.size=loop.size,
  #                        parallel.sz = parallel.sz,
  #                        verbose = verbose,
  #                        exclude.zero = exclude.zero,
  #                        parllelModel = parllelModel)

  # scDNSobjcet@Div.Parameters <- Div.Parameters

  NetDivs <- getKLD_cKLDnetwork(scDNSobjcet@data,
                     Network = scDNSobjcet@Network,
                     GroupLabel = scDNSobjcet@GroupLabel,
                     k = scDNSobjcet@Div.Parameters$k,
                     n.grid = scDNSobjcet@Div.Parameters$n.grid,
                     NoiseRemove = scDNSobjcet@Div.Parameters$NoiseRemove,
                     CoarseGrain = scDNSobjcet@Div.Parameters$CoarseGrain,
                     n.coarse = scDNSobjcet@Div.Parameters$n.coarse,
                     ds.method = scDNSobjcet@Div.Parameters$ds.method,
                     divergence = scDNSobjcet@Div.Parameters$divergence,
                     h = scDNSobjcet@Div.Parameters$h,
                     noiseSd = scDNSobjcet@Div.Parameters$noiseSd,
                     loop.size=scDNSobjcet@Div.Parameters$loop.size,
                     parallel.sz = scDNSobjcet@Div.Parameters$parallel.sz,
                     verbose = scDNSobjcet@Div.Parameters$verbose,
                     exclude.zero = scDNSobjcet@Div.Parameters$exclude.zero,
                     returnCDS = FALSE,
                     Div_weight = scDNSobjcet@Div.Parameters$Div_weight,
                     parllelModel = scDNSobjcet@Div.Parameters$parllelModel)
  scDNSobjcet@Network <- NetDivs$Network
  scDNSobjcet@JDensity_A <- NetDivs$ContextA_DS
  scDNSobjcet@JDensity_B <- NetDivs$ContextB_DS
  scDNSobjcet@uniCase <- NetDivs$uniCase
  rm(NetDivs)
  invisible(gc())
  scDNSobjcet

}

#' scDNS_2_creatNEAModel
#'
#' @param scDNSobjcet scDNSobjcet
#' @param n.dropGene Integer, number of cells sampled for model building (def:3000)
#' @param n.randNet Integer, number of cells sampled for model building (def:3000)
#' @param sdBias Numeric value, bias coefficient used to penalize low degree genes(>=1) (def:1.1)
#'
#' @return
#' scDNSobjcet
#' @export
#'
#' @examples
scDNS_2_creatNEAModel <- function(scDNSobjcet,
                                  n.dropGene = NULL,
                                  n.randNet = NULL,
                                  sdBias=NULL){

  if(!is.null(n.dropGene)){scDNSobjcet@NEA.Parameters$n.dropGene=n.dropGene}
  if(!is.null(n.randNet)){scDNSobjcet@NEA.Parameters$n.randNet=n.randNet}
  if(!is.null(sdBias)){scDNSobjcet@NEA.Parameters$sdBias=sdBias}


  # NEA.Parameters <- list(do.impute = FALSE,
  #                        n.dropGene = n.dropGene,
  #                        n.randNet = n.randNet,
  #                        sdBias = sdBias)
  # scDNSobjcet@NEA.Parameters <- NEA.Parameters
  NEAModel <- creatNEAModel(counts=scDNSobjcet@counts,
                            ExpData=scDNSobjcet@data,
                            do.impute = scDNSobjcet@NEA.Parameters$do.impute,
                            n.dropGene = scDNSobjcet@NEA.Parameters$n.dropGene,
                            n.randNet = scDNSobjcet@NEA.Parameters$n.randNet,
                            k = scDNSobjcet@Div.Parameters$k,
                            Likelihood=scDNSobjcet@GeneVariability,
                            GroupLabel=scDNSobjcet@GroupLabel,
                            n.grid = scDNSobjcet@Div.Parameters$n.grid,
                            n.coarse=scDNSobjcet@Div.Parameters$n.coarse,
                            noiseSd = scDNSobjcet@Div.Parameters$noiseSd,
                            Div_weight = scDNSobjcet@Div.Parameters$Div_weight,
                            CoarseGrain=scDNSobjcet@Div.Parameters$CoarseGrain,
                            loop.size=scDNSobjcet@Div.Parameters$loop.size,
                            parallel.sz = scDNSobjcet@Div.Parameters$parallel.sz,
                            verbose = scDNSobjcet@Div.Parameters$verbose,
                            exclude.zero = scDNSobjcet@Div.Parameters$exclude.zero,
                            NoiseRemove = scDNSobjcet@Div.Parameters$NoiseRemove,
                            divergence = scDNSobjcet@Div.Parameters$divergence,
                            ds.method = scDNSobjcet@Div.Parameters$ds.method,
                            h=scDNSobjcet@Div.Parameters$h,
                            sdBias=scDNSobjcet@NEA.Parameters$sdBias,
                            rb.jsd = scDNSobjcet@Div.Parameters$rb.jsd,
                            parllelModel =  scDNSobjcet@Div.Parameters$parllelModel)
  scDNSobjcet@NEAModel <- NEAModel
  scDNSobjcet
}


#' scDNS_3_GeneZscore
#'
#' @param scDNSobjcet scDNSobjcet
#'
#' @return
#' scDNSobjcet
#' @export
#'
#' @examples
scDNS_3_GeneZscore <- function(scDNSobjcet){
  scDNSobjcet <- getZscore(EdgeScore = scDNSobjcet,
                            NEAModel = scDNSobjcet@NEAModel,
                            Likelihood = scDNSobjcet@GeneVariability )
  scDNSobjcet
}

#' scDNS_4_scContribution
#'
#' @param scDNSobjcet scDNSobjcet
#' @param topGene Integer value, representing the number of top significant genes(def:100).When sigGene is not NULL, topGene is invalid.
#' @param sigGene A vector of strings representing the set of genes of interest (def:NULL)
#' @param q.th Adjusted p-value threshold
#'
#' @return
#' scDNSobjcet
#' @export
#'
#' @examples
scDNS_4_scContribution <- function(scDNSobjcet,
                                   topGene=100,
                                   sigGene=NULL,
                                   q.th=0.01,...){

  if(is.null(sigGene)){
    Zscore <- scDNSobjcet@Zscore
    Zscore$Qvalue.Zsplus <- p.adjust(Zscore$Pvalues.ZsPlus,method = 'BH')
    sigGene <- Zscore$Gene[Zscore$Qvalue.Zsplus<q.th]
    if(length(sigGene)==0){
      message('The number of sigGene is 0. We used top genes (100)')
    }else{
      message('The number of sigGene is ',length(sigGene),'.')
    }
  }

  nx <-  sqrt(ncol(scDNSobjcet@JDensity_A))
  # scCon_res <- scContribution(EdgeScore = scDNSobjcet@Network,
  #                    Zscores = scDNSobjcet@Zscore,
  #                    ExpData = scDNSobjcet@data,
  #                    nx = nx,
  #                    rmZero = scDNSobjcet@Div.Parameters$exclude.zero,
  #                    topGene = topGene,
  #                    sigGene = sigGene,
  #                    NEAModel = scDNSobjcet@NEAModel)
  # scCon_res <- scContribution(scDNSobjcet,
  #                          nx=nx,
  #                          rmZero = scDNSobjcet@Div.Parameters$exclude.zero,
  #                          topGene=topGene,
  #                          sigGene=sigGene)
  # scCon_res <- scContribution_v2(scDNSobjcet,
  #                                nx=nx,
  #                                rmZero = scDNSobjcet@Div.Parameters$exclude.zero,
  #                                topGene=topGene,
  #                                sigGene=sigGene,
  #                                rb = scDNSobjcet@Div.Parameters$rb.jsd)
  # scCon_res <- scContribution_v3(scDNSobjcet,
  #                             nx=nx,
  #                             # rmZero = scDNSobjcet@Div.Parameters$exclude.zero,
  #                             topGene=topGene,
  #                             sigGene=sigGene,
  #                             rb = scDNSobjcet@Div.Parameters$rb.jsd,...)
  scCon_res <- scContribution_v5(scDNSobjcet,
                                 nx=nx,
                                 # rmZero = scDNSobjcet@Div.Parameters$exclude.zero,
                                 topGene=topGene,
                                 sigGene=sigGene,
                                 rb = scDNSobjcet@Div.Parameters$rb.jsd,...)

  scDNSobjcet@scZscore <- scCon_res
  scDNSobjcet
}


#' scDNS_5_cluster
#'
#' @param scDNSobj scDNS object
#' @param Sobj seurat object
#' @param biasToscDNS The weight of scZscore during clustering, a floating point number (default 1).
#' When equal to 1, scZscore and PCA's principal component 1 have the same L2 norm.
#' @param resolution The resolution parameter in Seurat's FindClusters function controls the number of clusters.
#' A higher resolution will lead to more clusters, while a lower resolution will lead to fewer clusters.
#' @param merge.red Single-cell transcriptomics dimensionality reduction object (default:pca) used for integration with scZScore.
#' @param red.NewName New scDNS dimensionality reduction name (default: scDNS)
#'
#' @return
#' seurat object
#' @export
#'
#' @examples
scDNS_5_cluster <- function(scDNSobj,
                          Sobj,biasToscDNS=1,resolution=0.5,merge.red='pca',red.NewName='scDNS'){
  scCC <- Zscore <- scDNSobj@scZscore
  scCC_assay <- SeuratObject::CreateAssayObject(data = scCC)
  Sobj[['scCC']] <- scCC_assay
  scDNSLoading <- matrixStats::colSums2(scCC)
  names(scDNSLoading) <- colnames(scCC)
  Sobj$scDNSLoading <- scDNSLoading[colnames(Sobj)]
  Sobj <- scDNS_Cluster(Sobj,
                        scCC = scCC,
                        reduction = merge.red,
                        red.name = red.NewName,
                        biasToscDNS = biasToscDNS,
                        resolution = resolution)
  Sobj
}





