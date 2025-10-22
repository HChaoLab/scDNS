#' creatNEAModel
#'
#' @param ExpData a matrix represents gene normalized expression, with rows being genes and columns being cells.
#' @param counts a matrix represents gene count expression, with rows being genes and columns being cells
#' @param GroupLabel a string vector representing the grouping information of the cells
#' @param k a numerical value represents the kth (def: 10) nearest neighbor point
#' @param n.grid An integer value representing the number of grid (def:60) used to calculate the joint probability density
#' @param n.coarse An integer value representing the number of grid (def:20) after coarse-graining
#' @param loop.size An integer value representing the network size (def:3000) allocated to each compute core
#' @param parallel.sz An integer representing the number of parallel cores (def:1)
#' @param verbose Logical value, whether progression messages should be printed in the terminal.
#' @param exclude.zero Logical value indicating whether to exclude 0 values when calculating the joint probability density (def: FALSE)
#' @param NoiseRemove Logical value, indicating whether to remove noise that may lead to extremely low probability density (def:TRUE)
#' @param CoarseGrain Logical value, whether to coarse-grain the grid(def:TRUE).
#' @param ds.method Method to calculate joint probability density, optional knn and kde (def:knn)
#' @param divergence Method for calculating divergence, optional jsd (Jensen-Shannon divergence) and kld (Kullback–Leibler divergence)(def:jsd)
#' @param h A parameter used for kernel density estimation, which only takes effect when the density estimation method is kde.
#' @param do.impute Logical,whether to perform gene expression imputation
#' @param n.dropGene Integer, number of cells sampled for model building (def:3000)
#' @param n.randNet Integer, number of cells sampled for model building (def:3000)
#' @param sdBias Numeric value, bias coefficient used to penalize low degree genes(>=1) (def:1.1)
#' @param Likelihood Named numeric vector,
#' @param parllelModel parllelModel=c('foreach','bplapply')[1]
#' @param Div_weight sqrt(as.vector(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)*t(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)))),
#' @param rb.jsd (def:FALSE)
#'
#' @return
#' @keywords internal
#'
#' @examples
creatNEAModel<-function(counts,
                        ExpData=NULL,
                        do.impute = TRUE,
                        n.dropGene=3000,
                        n.randNet=20000,
                        k=10,
                        Div_weight=sqrt(as.vector(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)*t(matrix(1:n.coarse,nrow = n.coarse,ncol = n.coarse)))),
                        Likelihood,
                        GroupLabel=NULL,
                        n.grid = 60,
                        n.coarse=20,
                        CoarseGrain=TRUE,
                        loop.size=3000,
                        parallel.sz = 10,
                        noiseSd = 0.01,
                        verbose = TRUE,
                        exclude.zero = FALSE,
                        NoiseRemove = TRUE,
                        divergence = 'jsd',
                        ds.method = 'knn',
                        h=3,
                        sdBias=1.1,
                        rb.jsd=FALSE,
                        parllelModel=c('foreach','bplapply')[1]){

  # require('sn')

  if(do.impute&is.null(ExpData)){
    SOD=CreateSeuratObject(counts = counts)
    SOD <- Rmagic::magic(SOD)
    ExpData = SOD@assays$MAGIC_RNA@data
  }
  uniCase = unique(GroupLabel)
  #
  # ExpData = SommthOutLier_UP(mt = ExpData,zero.rm = T,percentile = 0.98)
  # SommthOutLier_UP(mt = ExpData['HBQ1',,drop=F],zero.rm = T,percentile = 0.99)%>%hist()
  #
  # ExpNc = rowSums(counts>0)
  # ExpData = ExpData[ExpNc>=minExpFeauture,]
  # sn = Likelihood[rownames(ExpData)]
  # RandGene = sample(names(sn[sn>0.90]),pmin(n.dropGene,length(sn[sn>0.90])))
  RandGene = sample(rownames(ExpData),pmin(n.dropGene,nrow(ExpData)))
  ExpData = ExpData[RandGene,]
  ExpData[,GroupLabel==uniCase[1]] = apply(ExpData[,GroupLabel==uniCase[1]] , 1, sample)%>%t()
  ExpData[,GroupLabel==uniCase[2]] = apply(ExpData[,GroupLabel==uniCase[2]] , 1, sample)%>%t()
  counts = counts[RandGene,]
  # 计算dropout矩阵
  # dropoutMatrix = ncol(counts)-Matrix::tcrossprod(counts==0 ) # Number of cells expressing neither A nor B
  dropoutMatrix = Matrix::tcrossprod(counts!=0 ) # number of cells expressing both G1 and G2
  dropoutMatrix_log = round(log10(dropoutMatrix+1),1)
  uniNpoints = unique(as.vector(dropoutMatrix_log))
  uniNpoints = sort(uniNpoints)
  # 生成随机网络
  message('Rand network')
  CadiateNet = NULL
  CadidateGene = rownames(dropoutMatrix_log)
  for(i in 1:length(uniNpoints)){
    CadiateNet_temp = BiocGenerics::which(dropoutMatrix_log==uniNpoints[i],arr.ind = TRUE,useNames = T)
    CadiateNet_temp = CadiateNet_temp[sample(nrow(CadiateNet_temp),
                                             round(n.randNet/length(uniNpoints),0),
                                             replace = T),] # sampling network for each dropout degree
    CadiateNet_temp = data.frame(source=CadidateGene[CadiateNet_temp[,1]],
                                 target=CadidateGene[CadiateNet_temp[,2]])
    CadiateNet = rbind(CadiateNet,CadiateNet_temp)

  }
  CadiateNet = preFilterNet(CadiateNet,counts)

  CadiateNet$npoint = dropoutMatrix[sub2ind(match(CadiateNet[,1],rownames(dropoutMatrix)),
                                            match(CadiateNet[,2],rownames(dropoutMatrix)),
                                            nrow = nrow(dropoutMatrix),ncol = ncol(dropoutMatrix))]
  CadiateNet$npoint_log = log10(CadiateNet$npoint+1)
  CadiateNet$npoint_ds = round(CadiateNet$npoint_log,1)
  CadiateNet <- CadiateNet[as.character(CadiateNet$npoint_ds)%in%names(table(CadiateNet$npoint_ds))[table(CadiateNet$npoint_ds)>30],]

  # CadiateNet2$npoint = dropoutMatrix[sub2ind(match(CadiateNet2[,1],rownames(dropoutMatrix)),
  #                                           match(CadiateNet2[,2],rownames(dropoutMatrix)),
  #                                           nrow = nrow(dropoutMatrix),ncol = ncol(dropoutMatrix))]
  # CadiateNet2$npoint_log = log10(CadiateNet2$npoint+1)
  # CadiateNet2$npoint_ds = round(CadiateNet2$npoint_log,1)
  #
  # CadiateNet2$npoint2 = dropoutMatrix[sub2ind(match(CadiateNet2[,1],rownames(dropoutMatrix)),
  #                                            match(CadiateNet2[,2],rownames(dropoutMatrix)),
  #                                            nrow = nrow(dropoutMatrix),ncol = ncol(dropoutMatrix))]
  # CadiateNet2$npoint_log2 = log10(CadiateNet2$npoint2+1)
  # CadiateNet2$npoint_ds2 = round(CadiateNet2$npoint_log2,1)
  #
  # ggplot(CadiateNet2,aes(npoint_log,Div,color=npoint_log3))+geom_point()+scale_color_viridis_b()
  # pec = log10(rowSums2(ExpData!=0)+1)
  # names(pec) = rownames(ExpData)
  # CadiateNet2$npoint_log3 = rowMins(matrix(c(pec[CadiateNet2[,1]],pec[CadiateNet2[,2]]),ncol = 2))


  # CadiateNet = CadiateNet[sample(nrow(CadiateNet),2000),]
  # randLabel = sample(GroupLabel)
  randLabel = GroupLabel
  # preFilterNet(CadiateNet,ExpData)

  RawKLD_R = getKLD_cKLDnetwork(ExpData = ExpData,
                                Network = CadiateNet,
                                k=k,
                                GroupLabel=randLabel,
                                n.grid = n.grid,
                                n.coarse = n.coarse,
                                loop.size = nrow(CadiateNet)/parallel.sz,
                                parallel.sz = parallel.sz,
                                verbose = verbose,
                                exclude.zero = exclude.zero,
                                NoiseRemove = NoiseRemove,
                                CoarseGrain = CoarseGrain,
                                returnCDS = TRUE,
                                Div_weight = Div_weight,
                                divergence = divergence,
                                ds.method = ds.method,
                                h = h,
                                noiseSd=noiseSd,
                                parllelModel=parllelModel)
  # EdgeID=900
  # Pdata <- data.frame(G1=ExpData[CadiateNet2[EdgeID,1],],G2=ExpData[CadiateNet2[EdgeID,2],],label=randLabel)
  # # Pdata <- data.frame(G1=ExpData['AP3B2',],G2=ExpData['FTL',],label=randLabel)
  # # Pdata <- data.frame(G1=counts[CadiateNet2[7659,1],],G2=counts[CadiateNet2[7659,2],],label=randLabel)
  # ggplot(Pdata,aes(G1,G2,color=label))+geom_point()
  # DensityPlotDF_withPoint(Pdata$G1,Pdata$G2)
  # showDenstiy(RawKLD_R,EdgeID = EdgeID)
  # densityCompare2(RawKLD_R,EdgeID = EdgeID)
  # RawKLD_R$Network[EdgeID,]
  # #
  # RawKLD_R2 = getKLD_cKLDnetwork(ExpData = ExpData,
  #                                Network = CadiateNet,
  #                                k=10,
  #                                GroupLabel=randLabel,
  #                                n.grid = n.grid,
  #                                n.coarse = n.coarse,
  #                                loop.size = loop.size,
  #                                parallel.sz = parallel.sz,
  #                                verbose = verbose,
  #                                exclude.zero = FALSE,
  #                                NoiseRemove = NoiseRemove,
  #                                CoarseGrain = TRUE,
  #                                returnCDS = TRUE,
  #                                divergence = 'jsd',
  #                                ds.method = 'knn')
  # RawKLD_R2 = getKLD_cKLDnetwork_fromDS(RawDS_1 = RawKLD_R$GontextA_DS,
  #                                       RawDS_2 = RawKLD_R$GontextB_DS,
  #                                       Network = RawKLD_R$Network,
  #                                       GroupLabel = RawKLD_R$uniCase,
  #                                       divergence = 'jsd')
  # RawKLD_R$Network$id = 1:nrow(RawKLD_R$Network)
  # RawKLD_R$Network$npoint_log = CadiateNet$npoint_log
  # View(RawKLD_R$Network)
  # EdgeID = 87
  # showDenstiy(RawKLD_R,EdgeID = EdgeID)
  # Network = CadiateNet[838,]
  # iAT2_data$randLabel = randLabel
  # FeatureScatter(iAT2_data, feature1 =RawKLD_R$Network[EdgeID,1], feature2 = RawKLD_R$Network[EdgeID,2],group.by = 'randLabel')
  # filteredNetid = RawKLD_R$Network$npoint_log>log10(minExpFeauture)

  #

  # RandHeight.KLD = getCutHeightWithDensityFeature(DS1 = RawKLD_R$ContextA_cDS$DS[filteredNetid,],
  #                                                 DS2 = RawKLD_R$ContextB_cDS$DS[filteredNetid,])
  #
  # RandHeight.cKLD = getCutHeightWithDensityFeature(DS1 = RawKLD_R$ContextA_cDS$cDS_1[filteredNetid,],
  #                                                  DS2 = RawKLD_R$ContextB_cDS$cDS_1[filteredNetid,])+
  #   getCutHeightWithDensityFeature(DS1 = RawKLD_R$ContextA_cDS$cDS_2[filteredNetid,],
  #                                  DS2 = RawKLD_R$ContextB_cDS$cDS_2[filteredNetid,])
  # RandHeight.cKLD = RandHeight.cKLD/2

  # RandHeight = getCutHeightWithDensityCombine(DS1 = RawKLD_R$ContextA_cDS$DS[filteredNetid,],
  #                                             DS2 = RawKLD_R$ContextB_cDS$DS[filteredNetid,],
  #                                             cDS_C1_D1 = RawKLD_R$ContextA_cDS$cDS_1[filteredNetid,],
  #                                             cDS_C1_D2 = RawKLD_R$ContextA_cDS$cDS_2[filteredNetid,],
  #                                             cDS_C2_D1 = RawKLD_R$ContextB_cDS$cDS_1[filteredNetid,],
  #                                             cDS_C2_D2 = RawKLD_R$ContextB_cDS$cDS_2[filteredNetid,])

  # Pdata = data.frame(G1=ExpData[CadiateNet2[EdgeID,1],],G2=ExpData[CadiateNet2[EdgeID,2],],group=randLabel)
  # ggplot(Pdata,aes(G1,G2,color=group))+geom_point()
  # Pdata = Pdata[!(Pdata$G1==0&Pdata$G2==0),]
  # DensityPlotDF_withPoint(x = Pdata$G1,y = Pdata$G2)

  CadiateNet2 = RawKLD_R$Network#[filteredNetid,] # 去掉了点少于10^1.6的关系对
  CadiateNet2$id = 1:nrow(CadiateNet2)
  RawDistrbution.Div <- NULL
  RawDistrbution.cDiv <- NULL
  for(i in unique(CadiateNet2$npoint_ds)%>%sort()){
    idx = CadiateNet2$npoint_ds==i
    RawDistrbution.Div <- c(RawDistrbution.Div,list(fit_selm_distrubution(CadiateNet2$Div[idx])))
    RawDistrbution.cDiv <- c(RawDistrbution.cDiv,list(fit_selm_distrubution(c(CadiateNet2$cDiv_D1[idx],CadiateNet2$cDiv_D2[idx]))))
  }
  names(RawDistrbution.Div) <- unique(CadiateNet2$npoint_ds)%>%sort()
  names(RawDistrbution.cDiv) <- unique(CadiateNet2$npoint_ds)%>%sort()
  # Pchi.Div=toChiSquareX(CadiateNet2$Div,CadiateNet2$npoint_ds,RawDistrbution.Div)
  # Pchi.cDiv_D1=toChiSquareX(CadiateNet2$cDiv_D1,CadiateNet2$npoint_ds,RawDistrbution.cDiv)
  # Pchi.cDiv_D2=toChiSquareX(CadiateNet2$cDiv_D2,CadiateNet2$npoint_ds,RawDistrbution.cDiv)
  # Pchi.Div = addLabel2Colnames(Pchi.Div,label = 'Div.',before = TRUE)
  # Pchi.cDiv_D1 = addLabel2Colnames(Pchi.cDiv_D1,label = 'cDiv_D1.',before = TRUE)
  # Pchi.cDiv_D2 = addLabel2Colnames(Pchi.cDiv_D2,label = 'cDiv_D2.',before = TRUE)
  # Pchi = cbind(Pchi.Div,Pchi.cDiv_D1,Pchi.cDiv_D2)
  ad.NetRScore3 <- function(netScore,DivFit,cDivFit){
    Pchi.Div=toChiSquareX(netScore$Div,netScore$npoint_ds,DivFit)
    Pchi.cDiv_D1=toChiSquareX(netScore$cDiv_D1,netScore$npoint_ds,cDivFit)
    Pchi.cDiv_D2=toChiSquareX(netScore$cDiv_D2,netScore$npoint_ds,cDivFit)
    Pchi.Div = addLabel2Colnames(Pchi.Div,label = 'Div.',before = TRUE)
    Pchi.cDiv_D1 = addLabel2Colnames(Pchi.cDiv_D1,label = 'cDiv_D1.',before = TRUE)
    Pchi.cDiv_D2 = addLabel2Colnames(Pchi.cDiv_D2,label = 'cDiv_D2.',before = TRUE)
    Pchi = cbind(Pchi.Div,Pchi.cDiv_D1,Pchi.cDiv_D2)
    cbind(netScore,Pchi)
  }
  # fit_poly(CadiateNet2$npoint_ds,y = Pchi.Div$chiSquare)
  # fit_poly(CadiateNet2$npoint_ds,y = Pchi.cDiv_D1$chiSquare)
  toChiSquareX<-function(rs,bias,DistrbutionList){
    Pvalues <-  rs
    uniBias <- unique(bias)
    Num_Model <- names(DistrbutionList)
    for(i in 1:length(uniBias)){
      ind=getCloseseData(data = as.numeric(Num_Model),uniBias[i],returnIndex = T)
      Pvalues[bias==uniBias[i]] = DistrbutionList[[ind]]$cal_Pvalue(rs[bias==uniBias[i]],
                                                                    ParmetersInput = DistrbutionList[[ind]]$ParmetersInput)$Pvalues
    }

    Pvalues[Pvalues==0]=1e-312
    chiSquare = sqrt(-log10(Pvalues))
    data.frame(Pvalues=Pvalues,chiSquare=chiSquare)
  }

  # test degree
  CadiateNet_samll = NULL
  message('creat random net')
  if(length(names(Likelihood[rownames(ExpData)][Likelihood[rownames(ExpData)]>0.95]))>=100){
    sn <-  sample(names(Likelihood[rownames(ExpData)][Likelihood[rownames(ExpData)]>0.95]),100)
  }else{
    sn <- sample(names(Likelihood[rownames(ExpData)]),
                 pmin(100,length(names(Likelihood[rownames(ExpData)]))),
                 prob = Likelihood[rownames(ExpData)])
  }

  CadiateNet_samll=lapply(1:length(sn), function(x)sample(rownames(ExpData),pmin(300,nrow(ExpData))))
  names(CadiateNet_samll) = sn
  CadiateNet_samll = List2dataFrame(CadiateNet_samll)
  colnames(CadiateNet_samll) = c('source','target')
  CadiateNet_samll$npoint_log = dropoutMatrix_log[sub2ind(match(CadiateNet_samll[,1],rownames(dropoutMatrix_log)),
                                                          match(CadiateNet_samll[,2],rownames(dropoutMatrix_log)),
                                                          nrow = nrow(dropoutMatrix_log),ncol = ncol(dropoutMatrix_log))]
  # CadiateNet_samll = CadiateNet_samll[CadiateNet_samll$dropout_log>log10(minExpFeauture),]
  # CadiateNet_samll = CadiateNet_samll[sampleByGroup(1:nrow(CadiateNet_samll),infGroup = CadiateNet_samll$source,s.num = 200,replace = T),]

  RawKLD_R_samll = getKLD_cKLDnetwork(ExpData = ExpData,
                                      Network = CadiateNet_samll,
                                      k = k,
                                      GroupLabel = GroupLabel,
                                      n.grid = n.grid,
                                      n.coarse = n.coarse,
                                      loop.size = nrow(CadiateNet_samll)/parallel.sz,
                                      parallel.sz = parallel.sz,
                                      verbose = verbose,
                                      exclude.zero = exclude.zero,
                                      NoiseRemove = NoiseRemove,
                                      CoarseGrain = CoarseGrain,
                                      returnCDS = FALSE,
                                      Div_weight = Div_weight,
                                      divergence = divergence,
                                      ds.method = ds.method,
                                      h = h,
                                      parllelModel=parllelModel,
                                      rb.jsd=rb.jsd)
  CadiateNet_samll = RawKLD_R_samll$Network
  CadiateNet_samll$npoint_ds = round(CadiateNet_samll$npoint_log,1)
  CadiateNet_samll = ad.NetRScore3(netScore = CadiateNet_samll,RawDistrbution.Div,RawDistrbution.cDiv)
  ad.dropModel = list(AdModelList_1=RawDistrbution.Div,
                      AdModelList_2=RawDistrbution.cDiv)


  # thresholdx = sqrt(-log10(0.5))
  # CadiateNet_samll$Div.chiSquare[CadiateNet_samll$Div.chiSquare<thresholdx]=0
  # CadiateNet_samll$cDiv_D1.chiSquare[CadiateNet_samll$cDiv_D1.chiSquare<thresholdx]=0
  # CadiateNet_samll$cDiv_D2.chiSquare[CadiateNet_samll$cDiv_D2.chiSquare<thresholdx]=0
  # DegreeDist = table(c(CadiateNet_samll$source,CadiateNet_samll$target))
  # ScoreType = 'KLD.ad.chiSquare'
  # ConnectedGenes = names(DegreeDist[DegreeDist>=2])
  # ConnectedEdge = lapply(ConnectedGenes, function(x)which(CadiateNet_samll[,1:2]==x,arr.ind = T)[,1]%>%as.numeric())
  # ConnectedEdge = lapply(ConnectedEdge, function(x)t(combn(x,2)))
  # ConnectedEdge = do.call('rbind',ConnectedEdge)
  # CEdgeScore = data.frame(x=CadiateNet_samll[ConnectedEdge_id[,1],ScoreType],
  #                         y=CadiateNet_samll[ConnectedEdge_id[,2],ScoreType]) # 相关联边的得分
  # xr = range(CEdgeScore$x)
  # dscr = findInterval_Meandata(CEdgeScore$x,seq(xr[1],xr[2],(xr[2]-xr[1])/10))%>%unique()
  # CEdgeScore$x_dc = findInterval_Meandata(CEdgeScore$x,seq(xr[1],xr[2],(xr[2]-xr[1])/10))
  # CEdgeScore$y_dc = findInterval_Meandata(CEdgeScore$y,seq(xr[1],xr[2],(xr[2]-xr[1])/10))
  # for(i in 1:length(dscr)){
  #   i=1
  #   hist(CEdgeScore$y[CEdgeScore$x_dc==dscr[i]])
  #   i=i+1
  # }
  #
  # Ras_3 = NULL
  # for(i in 1:30){
  #   degreeD = c(seq(1,9,1),seq(10,280,10))
  #   Ras_temp = data.frame(degree=degreeD,
  #                         RawScore = sapply( degreeD, function(x)sum(sample(CadiateNet_samll$Div.chiSquare,x))))
  #   Ras_3 = rbind(Ras_3,Ras_temp)
  # }
  # seaM=sea_fit(Dsize = Ras_3$degree,RawScores = Ras_3$RawScore,coarse = F)
  # DegreeD = table(unlist(CadiateNet_samll[,1:2]))
  # DegreeD = sort(DegreeD,decreasing = T)
  # RandNet_x = CadiateNet_samll
  #
  # NetList=tapply(c(1:nrow(CadiateNet_samll),1:nrow(CadiateNet_samll)), c(CadiateNet_samll$source,CadiateNet_samll$target), list)
  # Ras = NULL
  # for(i in c(seq(1,9,1),seq(10,max(DegreeD),round(max(DegreeD)/50,0)))){
  #   Netid = lapply(NetList, function(x)if(length(x)>=i&i!=1){sample(x,i)}else if(length(x)>=i&i==1){x}else{NA})
  #   Netid = Netid[sapply(Netid, length)==i]
  #   CaidGenes = names(Netid)
  #   if(length(CaidGenes)>=30){
  #     CaidGenes = sample(CaidGenes,30)
  #     for(j in 1:30){
  #       FeatureMatrix = cbind(RawKLD_R_samll$ContextA_cDS$DS[Netid[[CaidGenes[j]]],,drop=F],
  #                             RawKLD_R_samll$ContextB_cDS$DS[Netid[[CaidGenes[j]]],,drop=F])
  #       Ras_temp= getEdgeDviersity(Net = RandNet_x[Netid[[CaidGenes[j]]],,drop=F],
  #                        FeatureMatrix = FeatureMatrix,
  #                        cutHeight = RandHeight,
  #                        Nodes = CaidGenes[j],
  #                        scoreType = 'KLD.ad.chiSquare',nx=n.coarse)
  #       # Ras_temp=data.frame(gene = CaidGenes[j],degree=i,RawScore=sum(RandNet_x[Netid[[CaidGenes[j]]],'cKLD_1.ad.chiSquare']))
  #       Ras = rbind(Ras,Ras_temp)
  #     }
  #
  #   }
  # }
  #
  #
  # seaMR = sea_fit(Dsize = Ras$Degree,RawScores = Ras$RawScore1,coarse = F)
  #
  CadiateNet_samll$LR = pmin(Likelihood[CadiateNet_samll$source],Likelihood[CadiateNet_samll$target])
  CadiateNet_samll = CadiateNet_samll[CadiateNet_samll$LR>0.99,]
  NetList=tapply(c(1:nrow(CadiateNet_samll),1:nrow(CadiateNet_samll)), c(CadiateNet_samll$source,CadiateNet_samll$target), list)
  DegreeD = table(unlist(CadiateNet_samll[,1:2]))
  degreeRange = floor(10^seq(log10(1),log10(max(DegreeD)),0.05))%>%unique()
  RandRawScore = NULL
  # degreeRange = c(seq(1,9,1),seq(10,max(DegreeD),round(max(DegreeD)/50,0)))
  for(i in degreeRange){
    # Netid = lapply(NetList, function(x)if(length(x)>=i&i!=1){sample(x,i)}else if(length(x)>=i&i==1){x}else{NA})
    Netid = lapply(NetList, function(x)if(length(x)>=i){sample(x,i)}else{NA})
    Netid = Netid[sapply(Netid, length)==i]
    CaidGenes = names(Netid)
    if(length(CaidGenes)>=50){

      if(i < 20){
        CaidGenes = sample(CaidGenes,50)
        xtime = 1:50
      }else{
        CaidGenes = sample(CaidGenes,pmin(length(CaidGenes),200))
        xtime = 1:length(CaidGenes)
      }

      for(j in xtime){
        TempNet = CadiateNet_samll[Netid[[CaidGenes[j]]],]
        cDiv.data = TempNet[,c('cDiv_D1.chiSquare','cDiv_D2.chiSquare')]%>%as.matrix()
        cDiv.ind = which(TempNet[,1:2]==CaidGenes[j],arr.ind = TRUE)


        Raw.Div = sum(TempNet$Div.chiSquare*TempNet$LR)
        Raw.cDiv = sum(cDiv.data[cDiv.ind]*TempNet$LR)
        Raw.Or = sum(pmax(TempNet$Div.chiSquare*TempNet$LR,cDiv.data[cDiv.ind]*TempNet$LR))
        Raw.Plus = sum(TempNet$Div.chiSquare*TempNet$LR+cDiv.data[cDiv.ind]*TempNet$LR)
        Raw.Mpl = sum(TempNet$Div.chiSquare*TempNet$LR*cDiv.data[cDiv.ind]*TempNet$LR)
        Ras_temp = data.frame(gene = CaidGenes[j],
                              degree = i,
                              Raw.Div = Raw.Div,
                              Raw.cDiv = Raw.cDiv,
                              Raw.Or = Raw.Or,
                              Raw.Plus = Raw.Plus,
                              Raw.Mpl = Raw.Mpl)
        RandRawScore = rbind(RandRawScore,Ras_temp)
      }

    }
  }
  # RandRawScore$degree = ceiling(RandRawScore$degree)
  # RandRawScore$degree = round(findInterval_Meandata(RandRawScore$degree,v = degreeRange),0)+1
  # RandRawScore$degree = RandRawScore$degree%>%round(digits = 0)
  # degreeFreq = table(RandRawScore$degree)
  # RandRawScore = RandRawScore[RandRawScore$degree%in%as.numeric(names(degreeFreq[degreeFreq>50])),]
  # RandRawScore = RandRawScore[RandRawScore$degree>1,]
  # fit_poly(RandRawScore$degree,RandRawScore$Raw.Div)$p
  # RandRawScore$degree = round(findInterval_Meandata(RandRawScore$degree,v = c(1,2,3,4,5,6,7,8,9,10,50,100,150,200,250)),0)
  # # KLD
  # SmmD=group_by(RandRawScore,degree)%>%summarise(sd=sd(Raw.Div),mean=mean(Raw.Div))%>%as.data.frame()
  # MeanF = fit_poly(log10(RandRawScore$degree),RandRawScore$Raw.Div,degree = 5)$p
  # SDF = fit_poly(SmmD$degree,SmmD$mean,degree = 5)
  RandRawScore$degree = ceiling(RandRawScore$degree*mean(CadiateNet_samll$LR))
  ZscoreFit.Div = sea_fit_fix_st(Dsize = ceiling(RandRawScore$degree), # 改成了st分布
                              RawScores = RandRawScore$Raw.Div,coarse = F,
                              mean_1 = mean(RandRawScore$Raw.Div[RandRawScore$degree==1]),
                              sd_1 = sd(RandRawScore$Raw.Div[RandRawScore$degree==1])*sdBias)
  # ZscoreFit.KLD = sea_fit(Dsize = log10(RandRawScore$degree),RawScores = log2(RandRawScore$Raw.KLD),coarse = F)
  # RandRawScore_temp = RandRawScore[RandRawScore$degree>10,]
  # ZscoreFit.Div = sea_fit(Dsize = log10(RandRawScore$degree),
  #                         RawScores = RandRawScore$Raw.Div,
  #                         coarse = F)
  #
  Zs = ZscoreFit.Div$cal_Zscore_Pvalue(RawScore =  RandRawScore$Raw.Div,
                                       Dsize = RandRawScore$degree,
                                       ParmetersInput = ZscoreFit.Div$ParmetersInput)
  RandRawScore$Zscores = Zs$Zscores
  RandRawScore$Pvalues = Zs$Pvalues
  RandRawScore$Zscores.Div = Zs$Zscores
  p.div<-DensityPlotDF_withPoint(log10(RandRawScore$degree),RandRawScore$Zscores)+
    xlab(expression(Log[10]('degree')))+ylab('Z-score')
  # DensityPlotDF_withPoint(log10(RandZScore$degree),log10(RandZScore$Pvalues))+xlim(c(1,2.5))
  # fit_poly(log10(RandZScore$degree),RandZScore$Zscores)$p|ZscoreFit.Div$CombinePs
  # fit_poly(RandZScore$degree,RandZScore$Zscores)$p|ZscoreFit.Div$CombinePs
  # randNois = fit_selm_distrubution(Zs$Zscores[Zs$Zscores>-3&Zs$Zscores<8])
  # ZscoreFit.Div$ZscoreD = randNois$reslutP
  # ZscoreFit.Div$CombinePs = CombinePlots(list(ZscoreFit.Div$SizeVsMean,ZscoreFit.Div$SizeVsStd,ZscoreFit.Div$ZscoreD),ncol = 3)
  # ZscoreFit.Div$ParmetersInput$ParmetersZsD = randNois$ParmetersInput

  # Zs = ZscoreFit.Div$cal_Zscore_Pvalue(RawScore = RandRawScore$Raw.cDiv,
  #                                      Dsize = RandRawScore$degree,
  #                                      ParmetersInput = ZscoreFit.Div$ParmetersInput)
  # Zs = cbind(RandRawScore,Zs)
  # fit_poly(log10(Zs$degree),Zs$Zscores)$p

  #ZscoreFit.Div$CombinePs
  # cDiv
  ZscoreFit.cDiv = sea_fit_fix_st(Dsize = RandRawScore$degree,
                               RawScores = RandRawScore$Raw.cDiv,
                               coarse = F,
                               mean_1 = mean(RandRawScore$Raw.cDiv[RandRawScore$degree==1]),
                               sd_1 = sd(RandRawScore$Raw.cDiv[RandRawScore$degree==1])*sdBias)
  Zs = ZscoreFit.cDiv$cal_Zscore_Pvalue(RawScore =  RandRawScore$Raw.cDiv,
                                        Dsize = RandRawScore$degree,
                                        ParmetersInput = ZscoreFit.cDiv$ParmetersInput)
  RandRawScore$Zscores = Zs$Zscores
  RandRawScore$Pvalues = Zs$Pvalues
  RandRawScore$Zscores.cDiv = Zs$Zscores #
  p.cDiv<-DensityPlotDF_withPoint(log10(RandRawScore$degree),RandRawScore$Zscores)+
    xlab(expression(Log[10]('degree')))+ylab('Z-score')

  # or
  ZscoreFit.OR = sea_fit_fix_st(Dsize = RandRawScore$degree,
                             RawScores = RandRawScore$Raw.Or,
                             coarse = F,
                             mean_1 = mean(RandRawScore$Raw.Or[RandRawScore$degree==1]),
                             sd_1 = sd(RandRawScore$Raw.Or[RandRawScore$degree==1])*sdBias)
  Zs = ZscoreFit.OR$cal_Zscore_Pvalue(RawScore =  RandRawScore$Raw.Or,
                                      Dsize = RandRawScore$degree,
                                      ParmetersInput = ZscoreFit.OR$ParmetersInput)
  RandRawScore$Zscores = Zs$Zscores
  RandRawScore$Pvalues = Zs$Pvalues
  p.or<-DensityPlotDF_withPoint(log10(RandRawScore$degree),RandRawScore$Zscores)+
    xlab(expression(Log[10]('degree')))+ylab('Z-score')
  #
  ZscoreFit.Plus = sea_fit_fix_st(Dsize = RandRawScore$degree,
                               RawScores = RandRawScore$Raw.Plus,
                               coarse = F,
                               mean_1 = mean(RandRawScore$Raw.Plus[RandRawScore$degree==1]),
                               sd_1 = sd(RandRawScore$Raw.Plus[RandRawScore$degree==1])*sdBias)
  # ZscoreFit.Plus = sea_fit_fix(Dsize = RandRawScore$degree,
  #                              RawScores = RandRawScore$Raw.Plus,
  #                              coarse = F,
  #                              mean_1 = NULL,
  #                              sd_1 =NULL)

  Zs = ZscoreFit.Plus$cal_Zscore_Pvalue(RawScore =  RandRawScore$Raw.Plus,
                                        Dsize = RandRawScore$degree,
                                        ParmetersInput = ZscoreFit.Plus$ParmetersInput)
  RandRawScore$Zscores = Zs$Zscores
  RandRawScore$Pvalues = Zs$Pvalues
  p.plus<-DensityPlotDF_withPoint(log10(RandRawScore$degree),RandRawScore$Zscores)+
    xlab(expression(Log[10]('degree')))+ylab('Z-score')
  #
  ZscoreFit.Mpl = sea_fit_fix_st(Dsize = RandRawScore$degree,
                              RawScores = RandRawScore$Raw.Mpl,
                              coarse = F,
                              mean_1 = mean(RandRawScore$Raw.Mpl[RandRawScore$degree==1]),
                              sd_1 = sd(RandRawScore$Raw.Mpl[RandRawScore$degree==1])*sdBias)
  Zs = ZscoreFit.Mpl$cal_Zscore_Pvalue(RawScore =  RandRawScore$Raw.Mpl,
                                       Dsize = RandRawScore$degree,
                                       ParmetersInput = ZscoreFit.Mpl$ParmetersInput)
  RandRawScore$Zscores = Zs$Zscores
  RandRawScore$Pvalues = Zs$Pvalues
  p.Mpl<-DensityPlotDF_withPoint(log10(RandRawScore$degree),RandRawScore$Zscores)+
    xlab(expression(Log[10]('degree')))+ylab('Z-score')
  #
  # RandHeight = getCutHeightWithDensityCombine(DS1 = RawKLD_R_samll$ContextA_cDS$DS,
  #                                             DS2 = RawKLD_R_samll$ContextB_cDS$DS,
  #                                             cDS_C1_D1 = RawKLD_R_samll$ContextA_cDS$cDS_1,
  #                                             cDS_C1_D2 = RawKLD_R_samll$ContextA_cDS$cDS_2,
  #                                             cDS_C2_D1 = RawKLD_R_samll$ContextB_cDS$cDS_1,
  #                                             cDS_C2_D2 = RawKLD_R_samll$ContextB_cDS$cDS_2)

  RandRawScore$Zscore.Plus <- RandRawScore$Zscores.Div+RandRawScore$Zscores.cDiv
  ZscoreFit.ZsPlus<- fit_selm_distrubution_ST(RandRawScore$Zscore.Plus)

  NEAModel = list(DegreeData= RawKLD_R_samll, # save data for subsequent testing.
                  ad.dropModel = ad.dropModel,
                  RandRawScore = RandRawScore,
                  ZscoreFit.Div = ZscoreFit.Div,
                  ZscoreFit.cDiv = ZscoreFit.cDiv,
                  ZscoreFit.OR = ZscoreFit.OR,
                  ZscoreFit.Plus = ZscoreFit.Plus,
                  ZscoreFit.Mpl = ZscoreFit.Mpl,
                  ZscoreFit.ZsPlus=ZscoreFit.ZsPlus)

  return(NEAModel)




}










# showDenstiy<-function(EdgeScore,Nodes,EdgeID=NULL,interpolate=FALSE,filp=FALSE,subEdgeID=1,titleSzie=10){
#   if(is.null(EdgeID)){
#     if(length(Nodes)==1){
#       EdgeID=EdgeScore$Network[,1]%in%Nodes|EdgeScore$Network[,2]%in%Nodes
#       if(sum(EdgeID)>1){
#         print(sum(EdgeID))
#         EdgeID = which(EdgeID)[subEdgeID]
#       }
#     }else{
#       EdgeID=EdgeScore$Network[,1]%in%Nodes&EdgeScore$Network[,2]%in%Nodes
#     }
#
#   }
#   if(sum(EdgeID+0)!=1&length(EdgeID)!=1){
#     stop('EdgeID is mutiple or dose not exit')
#   }
#   subnet = EdgeScore$Network[EdgeID,]
#   ContextA_R=EdgeScore$ContextA_cDS
#   ContextB_R=EdgeScore$ContextB_cDS
#   if(filp){
#     p1<-DensityPlot_raster(ContextA_R$DS[EdgeID,],interpolate = interpolate,transpose = T)+
#       xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
#       DensityPlot_raster(ContextA_R$cDS_1[EdgeID,],addSideBar = T,transpose = F,interpolate = interpolate)+
#       xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
#       DensityPlot_raster(ContextA_R$cDS_2[EdgeID,],addSideBar = T,transpose = T,interpolate = interpolate)+
#       xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+plot_annotation('Conetxt A')+FontSize(main=titleSzie)
#     p2<-DensityPlot_raster(ContextB_R$DS[EdgeID,],interpolate = interpolate,transpose = T)+
#       xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
#       DensityPlot_raster(ContextB_R$cDS_1[EdgeID,],addSideBar = T,transpose = F,interpolate = interpolate)+
#       xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
#       DensityPlot_raster(ContextB_R$cDS_2[EdgeID,],addSideBar = T,transpose = T,interpolate = interpolate)+
#       xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+
#       plot_annotation(title = 'Conetxt B')+FontSize(main=titleSzie)
#     p1/p2+plot_layout(guides='collect')
#   }else{
#     p1<-DensityPlot_raster(ContextA_R$DS[EdgeID,],interpolate = interpolate)+
#       xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
#       DensityPlot_raster(ContextA_R$cDS_1[EdgeID,],addSideBar = T,transpose = F,interpolate = interpolate)+
#       xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
#       DensityPlot_raster(ContextA_R$cDS_2[EdgeID,],addSideBar = T,transpose = T,interpolate = interpolate)+
#       xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation('Conetxt A')+FontSize(main=titleSzie)
#     p2<-DensityPlot_raster(ContextB_R$DS[EdgeID,],interpolate = interpolate)+
#       xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
#       DensityPlot_raster(ContextB_R$cDS_1[EdgeID,],addSideBar = T,transpose = F,interpolate = interpolate)+
#       xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
#       DensityPlot_raster(ContextB_R$cDS_2[EdgeID,],addSideBar = T,transpose = T,interpolate = interpolate)+
#       xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation('Conetxt B')+FontSize(main=titleSzie)
#     p1/p2+plot_layout(guides='collect')
#   }
#
# }




showDenstiy<-function(EdgeScore,Nodes,EdgeID=NULL,interpolate=FALSE,filp=FALSE,subEdgeID=1,titleSzie=10){
  if(is.null(EdgeID)){
    if(length(Nodes)==1){
      EdgeID=EdgeScore$Network[,1]%in%Nodes|EdgeScore$Network[,2]%in%Nodes
      if(sum(EdgeID)>1){
        print(sum(EdgeID))
        EdgeID = which(EdgeID)[subEdgeID]
      }
    }else{
      EdgeID=EdgeScore$Network[,1]%in%Nodes&EdgeScore$Network[,2]%in%Nodes
    }

  }
  if(sum(EdgeID+0)!=1&length(EdgeID)!=1){
    stop('EdgeID is mutiple or dose not exit')
  }
  subnet = EdgeScore$Network[EdgeID,]
  ContextA_R=EdgeScore$ContextA_cDS
  ContextB_R=EdgeScore$ContextB_cDS
  if(filp){
    p1<-DensityPlot_raster(ContextA_R$DS[EdgeID,],interpolate = interpolate,transpose = T)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_1[EdgeID,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_2[EdgeID,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+plot_annotation('Conetxt A')+FontSize(main=titleSzie)
    p2<-DensityPlot_raster(ContextB_R$DS[EdgeID,],interpolate = interpolate,transpose = T)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+NoLegend()+labs(subtitle='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_1[EdgeID,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(subtitle='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_2[EdgeID,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(subtitle='Conditional probability')+
      plot_annotation(title = 'Conetxt B')+FontSize(main=titleSzie)
    p1/p2+plot_layout(guides='collect')
  }else{
    p1<-DensityPlot_raster(ContextA_R$DS[EdgeID,],interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_1[EdgeID,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextA_R$cDS_2[EdgeID,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation('Conetxt A')+FontSize(main=titleSzie)
    p2<-DensityPlot_raster(ContextB_R$DS[EdgeID,],interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Joint probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_1[EdgeID,],addSideBar = T,transpose = F,interpolate = interpolate)+
      xlab(subnet[1,1])+ylab(subnet[1,2])+NoLegend()+labs(title='Conditional probability')+FontSize(main=titleSzie)|
      DensityPlot_raster(ContextB_R$cDS_2[EdgeID,],addSideBar = T,transpose = T,interpolate = interpolate)+
      xlab(subnet[1,2])+ylab(subnet[1,1])+labs(title='Conditional probability')+plot_annotation('Conetxt B')+FontSize(main=titleSzie)
    p1/p2+plot_layout(guides='collect')
  }

}
getBoxCoxlambda<-function(x){
  b <- boxcox(lm(x ~ 1))
  lambda <- b$x[which.max(b$y)]
  lambda
}


boxcox_t<-function(x){
  new_x_exact <- (x ^ lambda - 1) / lambda
}




combineEdgeScore<-function(EdgeScore,NEAModel){
  addLabel2Colnames <-function(x,label,sep = ''){
    colnames(x) = paste(colnames(x),label,sep = sep)
    x
  }
  AdModelList_1 = NEAModel$ad.dropModel$AdModelList_1
  AdModelList_2 = NEAModel$ad.dropModel$AdModelList_2
  RandHeight = NEAModel$RandHeight$RandHeight
  ZscoreFit.KLD = NEAModel_ttt$ZscoreFit.KLD
  ZscoreFit.cKLD = NEAModel_ttt$ZscoreFit.cKLD
  ZscoreFit.cKLDorKLD = NEAModel_ttt$ZscoreFit.cKLDorKLD
  net = EdgeScore$Network
  net = adjustNetScore(net = net,
                       AdModelList_1 = AdModelList_1,
                       AdModelList_2 = AdModelList_2)
  net$cKLDorKLD.ad.chiSquare_1 = pmax(net$KLD.ad.chiSquare,net$cKLD_1.ad.chiSquare)
  net$cKLDorKLD.ad.chiSquare_2 = pmax(net$KLD.ad.chiSquare,net$cKLD_2.ad.chiSquare)

  RawScore = getEdgeDviersity_ds_cds(net = net,
                                     FeatureMatrix_DS = cbind(EdgeScore$ContextA_cDS$DS,EdgeScore$ContextB_cDS$DS),
                                     FeatureMatrix_cDS_D1 = cbind(EdgeScore$ContextA_cDS$cDS_1,EdgeScore$ContextB_cDS$cDS_1),
                                     FeatureMatrix_cDS_D2 = cbind(EdgeScore$ContextA_cDS$cDS_2,EdgeScore$ContextB_cDS$cDS_2),
                                     nodes = NULL,
                                     cutHeight = RandHeight)

  XXX = getEdgeDviersity_ds_cds_raw(net = net,
                                    FeatureMatrix_DS = cbind(EdgeScore$GontextA_cDS$DS,EdgeScore$GontextB_cDS$DS),
                                    FeatureMatrix_cDS_D1 = cbind(EdgeScore$GontextA_cDS$cDS_1,EdgeScore$GontextB_cDS$cDS_1),
                                    FeatureMatrix_cDS_D2 = cbind(EdgeScore$GontextA_cDS$cDS_2,EdgeScore$GontextB_cDS$cDS_2),
                                    nodes = NULL,
                                    cutHeight = RandHeight)


  Zscore.KLD = ZscoreFit.KLD$cal_Zscore_Pvalue(RawScore = RawScore$Raw.KLD.p,
                                               Dsize = RawScore$NodeDiversity,ParmetersInput = ZscoreFit.KLD$ParmetersInput)
  Zscore.KLD = addLabel2Colnames(Zscore.KLD,'.KLD')
  Zscore.cKLD = ZscoreFit.cKLD$cal_Zscore_Pvalue(RawScore = RawScore$Raw.cKLD.p,
                                                 Dsize = RawScore$NodeDiversity,ParmetersInput = ZscoreFit.cKLD$ParmetersInput)
  Zscore.cKLD = addLabel2Colnames(Zscore.cKLD,'.cKLD')
  Zscore.OR = ZscoreFit.cKLDorKLD$cal_Zscore_Pvalue(RawScore = RawScore$Raw.OR.p,
                                                    Dsize = RawScore$NodeDiversity,ParmetersInput = ZscoreFit.cKLDorKLD$ParmetersInput)
  Zscore.OR = addLabel2Colnames(Zscore.OR,'.OR')

  # Zscore.KLD = ZscoreFit.KLD$cal_Zscore_Pvalue(RawScore = RawScore$Raw.KLD.f,
  #                                              Dsize = RawScore$Degree,ParmetersInput = ZscoreFit.KLD$ParmetersInput)
  # Zscore.KLD = addLabel2Colnames(Zscore.KLD,'.KLD')
  # Zscore.cKLD = ZscoreFit.cKLD$cal_Zscore_Pvalue(RawScore = RawScore$Raw.cKLD.f,
  #                                                Dsize = RawScore$Degree,ParmetersInput = ZscoreFit.cKLD$ParmetersInput)
  # Zscore.cKLD = addLabel2Colnames(Zscore.cKLD,'.cKLD')
  # Zscore.OR = ZscoreFit.cKLDorKLD$cal_Zscore_Pvalue(RawScore = RawScore$Raw.OR.f,
  #                                                   Dsize = RawScore$Degree,ParmetersInput = ZscoreFit.cKLDorKLD$ParmetersInput)
  # Zscore.OR = addLabel2Colnames(Zscore.OR,'.OR')
  Zscore = cbind(RawScore[,1:3],Zscore.KLD,Zscore.cKLD,Zscore.OR)
  Zscore
  # fit_poly(log10(Zscore$NodeDiversity),Zscore$Zscores.KLD)$p|fit_poly(log10(Zscore$NodeDiversity),Zscore$Zscores.cKLD)$p|
  #   fit_poly(log10(Zscore$NodeDiversity),Zscore$Zscores.OR)$p
}






getCutHeightWithDensityFeature<-function(DS1,DS2){
  CbindDS = cbind(DS1,DS2)
  sampeID = sample(nrow(CbindDS),size = 1000)
  RandDist = dist(t(scale(t(CbindDS[sampeID,]))))
  Randhc = hclust(RandDist)
  randCluster = data.frame(Type = 'RandEDGE',height=seq(0,max(RandDist),max(RandDist)/200),
                           ratio = sapply(seq(0,max(RandDist),max(RandDist)/200),
                                          function(x)cutree(Randhc,h=x)%>%unique()%>%length()/length(Randhc$order)))
  RandHeight = max(randCluster$height[randCluster$ratio>=0.9])
  # ggplot(randCluster,aes(height,ratio))+geom_line(color='green')+
  # geom_hline(yintercept = min(randCluster$ratio[randCluster$ratio>=0.9]),colour='grey50',lty=2)+
  # geom_vline(xintercept = RandHeight,colour='grey50',lty=2)+theme_cowplot_i()+xlab('cutting height')
  RandHeight
}



getCutHeightWithDensityCombine<-function(DS1,DS2,cDS_C1_D1,cDS_C1_D2,cDS_C2_D1,cDS_C2_D2){
  CbindDS = cbind(DS1,DS2)
  CbindcDS_D1 = cbind(cDS_C1_D1,cDS_C2_D1)
  CbindcDS_D2 = cbind(cDS_C1_D2,cDS_C2_D2)
  sampeID = sample(nrow(CbindDS),size = min(c(1000,nrow(DS1))))
  RandDist = dist(t(scale(t(CbindDS[sampeID,]))))+dist(t(scale(t(CbindcDS_D1[sampeID,]))))+dist(t(scale(t(CbindcDS_D2[sampeID,]))))
  RandDist = RandDist/3
  Randhc = hclust(RandDist)
  randCluster = data.frame(Type = 'RandEDGE',height=seq(0,max(RandDist),max(RandDist)/200),
                           ratio = sapply(seq(0,max(RandDist),max(RandDist)/200),
                                          function(x)cutree(Randhc,h=x)%>%unique()%>%length()/length(Randhc$order)))
  RandHeight = max(randCluster$height[randCluster$ratio>=0.9])
  p<-ggplot(randCluster,aes(height,ratio))+geom_line(color='green')+
    geom_hline(yintercept = min(randCluster$ratio[randCluster$ratio>=0.9]),colour='grey50',lty=2)+
    geom_vline(xintercept = RandHeight,colour='grey50',lty=2)+theme_cowplot_i()+xlab('cutting height')
  return(list(RandHeight=RandHeight,p=p))
}


ZscoreAD<-function(RawScore,
                   sdModel,
                   meanModel,
                   covarData){
  (RawScore-predict(meanModel,data.frame(x=covarData)))/predict(sdModel,data.frame(x=covarData))
}
toChiSquare <- function(s,ddModel){
  Pvalues=ddModel$cal_Pvalue(s,ddModel$ParmetersInput)$Pvalue
  Pvalues[Pvalues==0]=1e-312
  # -log(Pvalues)*2
  sqrt(-log10(Pvalues))
}
toFC<-function(x,cova,meanFit){
  log2(x/predict(meanFit$predictP,data.frame(x=cova)))
}
adjustNetScore<-function(net,AdModelList_1,AdModelList_2,AdModelList_3){
  #1
  fitMean.Div.r = AdModelList_1$fitMean.Div.r
  fitMean.cDiv.r = AdModelList_1$fitMean.cDiv.r

  net$Div.fc = toFC(net$Div,cova = net$npoint_log,meanFit = fitMean.Div.r)
  net$cDiv_D1.fc = toFC(net$cDiv_D1,cova = net$npoint_log,meanFit = fitMean.cDiv.r)
  net$cDiv_D2.fc = toFC(net$cDiv_D2,cova = net$npoint_log,meanFit = fitMean.cDiv.r)
  #2
  fitMean.Div.fc = AdModelList_2$fitMean.Div.fc
  fitsd.Div.fc = AdModelList_2$fitsd.Div.fc
  fitsd.cDiv.fc = AdModelList_2$fitsd.cDiv.fc
  fitMean.cDiv.fc = AdModelList_2$fitMean.cDiv.fc

  net$Div.ad = ZscoreAD(net$Div.fc,
                        sdModel = fitsd.Div.fc$predictP,
                        meanModel = fitMean.Div.fc$predictP,
                        net$npoint_log)
  net$cDiv_D1.ad = ZscoreAD(net$cDiv_D1.fc,
                            sdModel = fitsd.cDiv.fc$predictP,
                            meanModel = fitMean.cDiv.fc$predictP,
                            net$npoint_log)
  net$cDiv_D2.ad = ZscoreAD(net$cDiv_D2.fc,
                            sdModel = fitsd.cDiv.fc$predictP,
                            meanModel = fitMean.cDiv.fc$predictP,
                            net$npoint_log)
  #3
  Div_db = AdModelList_3$Div_db
  cDiv_db = AdModelList_3$cDiv_db
  net$Div.chiSquare = toChiSquare(net$Div.ad,ddModel = Div_db)
  net$cDiv_D1.chiSquare = toChiSquare(net$cDiv_D1.ad,ddModel = cDiv_db)
  net$cDiv_D2.chiSquare = toChiSquare(net$cDiv_D2.ad,ddModel = cDiv_db)
  net
}

adjustNetScore.chiSquare<-function(net,AdModelList_1){

  #3
  Div_db = AdModelList_1$Div_db
  cDiv_db = AdModelList_1$cDiv_db
  net$Div.chiSquare = toChiSquare(net$Div,ddModel = Div_db)
  net$cDiv_D1.chiSquare = toChiSquare(net$cDiv_D1,ddModel = cDiv_db)
  net$cDiv_D2.chiSquare = toChiSquare(net$cDiv_D2,ddModel = cDiv_db)
  net
}

#' @keywords internal
getEdgeDviersity <- function(Net,
                             FeatureMatrix=NULL,
                             FeatureMatrix_C1,
                             FeatureMatrix_C2,
                             Nodes=NULL,
                             cutHeight = 8.5,
                             scoreType){
  if(is.null(Nodes)){
    Nodes = unique(unlist(Net[,1:2]))
  }


  if(is.null(FeatureMatrix)){
    FeatureMatrix = cbind(FeatureMatrix_C1,FeatureMatrix_C2)
  }
  FeatureMatrix = t(scale(t(FeatureMatrix)))
  #
  nx = sqrt(ncol(FeatureMatrix)/2)
  flipIndex_1 = t(matrix(1:nx^2,nx))%>%as.vector()
  flipIndex_2 = flipIndex_1+nx^2

  NodeResult = NULL
  for(i in 1:length(Nodes)){
    Direct1=Net[,1]%in%Nodes[i]
    Direct2=Net[,2]%in%Nodes[i]
    tempNetR = rbind(Net[Direct1,],Net[Direct2,])
    if(nrow(tempNetR)>1){
      FeatureMatrix_temp = rbind(FeatureMatrix[Direct1,],FeatureMatrix[Direct2,c(flipIndex_1,flipIndex_2)])
      ds_dist = dist(FeatureMatrix_temp)
      SumRawScore2 = tapply(tempNetR[,scoreType], cutree(hclust(ds_dist),h=cutHeight), mean)%>%base::sum()
      NodeDiversity = cutree(hclust(ds_dist),h=cutHeight)%>%unique()%>%length()
      tempNodeR = data.frame(Gene=Nodes[i],Degree = nrow(tempNetR),
                             NodeDiversity=NodeDiversity,
                             RawScore1=sum(tempNetR[,scoreType]),
                             RawScore2=SumRawScore2)
    }else{
      tempNodeR = data.frame(Gene=Nodes[i],
                             Degree = nrow(tempNetR),
                             NodeDiversity=1,
                             RawScore1=tempNetR[,scoreType],
                             RawScore2=tempNetR[,scoreType])
    }
    NodeResult = rbind(NodeResult,tempNodeR)
    # print(i)
  }
  NodeResult
}
#' @keywords internal
getEdgeDviersity_cds <- function(Net,
                                 FeatureMatrix_D1,
                                 FeatureMatrix_D2,
                                 Nodes=NULL,
                                 cutHeight = 8.5,
                                 scoreType){
  if(is.null(Nodes)){
    Nodes = unique(unlist(Net[,1:2]))
  }

  FeatureMatrix_D1 = t(scale(t(FeatureMatrix_D1)))
  FeatureMatrix_D2 = t(scale(t(FeatureMatrix_D2)))
  #
  nx = sqrt(ncol(FeatureMatrix_D1)/2)
  flipIndex_1 = t(matrix(1:nx^2,nx))%>%as.vector()
  flipIndex_2 = flipIndex_1+nx^2

  NodeResult = NULL
  for(i in 1:length(Nodes)){
    Direct1=Net[,1]%in%Nodes[i]
    Direct2=Net[,2]%in%Nodes[i]
    tempNetR = rbind(Net[Direct1,],Net[Direct2,])
    if(nrow(tempNetR)>1){
      FeatureMatrix_temp = rbind(FeatureMatrix_D1[Direct1,],FeatureMatrix_D2[Direct2,c(flipIndex_1,flipIndex_2)])
      ds_dist = dist(FeatureMatrix_temp)
      SumRawScore2 = tapply(tempNetR[,scoreType], cutree(hclust(ds_dist),h=cutHeight), mean)%>%base::sum()
      NodeDiversity = cutree(hclust(ds_dist),h=cutHeight)%>%unique()%>%length()
      tempNodeR = data.frame(Gene=Nodes[i],Degree = nrow(tempNetR),
                             NodeDiversity=NodeDiversity,
                             RawScore1=sum(tempNetR[,scoreType]),
                             RawScore2=SumRawScore2)
    }else{
      tempNodeR = data.frame(Gene=Nodes[i],
                             Degree = nrow(tempNetR),
                             NodeDiversity=1,
                             RawScore1=tempNetR[,scoreType],
                             RawScore2=tempNetR[,scoreType])
    }
    NodeResult = rbind(NodeResult,tempNodeR)
    # print(i)
  }
  NodeResult
}

getEdgeDviersity_ds_cds <- function(net,
                                    FeatureMatrix_DS,
                                    FeatureMatrix_cDS_D1,
                                    FeatureMatrix_cDS_D2,
                                    nodes=NULL,
                                    cutHeight = 8.5){
  if(is.null(nodes)){
    nodes = unique(unlist(net[,1:2]))
  }
  nx = sqrt(ncol(FeatureMatrix_DS)/2)
  flipIndex_1 = t(matrix(1:nx^2,nx))%>%as.vector()
  flipIndex_2 = flipIndex_1+nx^2

  FeatureMatrix_DS = t(scale(t(FeatureMatrix_DS)))
  FeatureMatrix_cDS_D1 = t(scale(t(FeatureMatrix_cDS_D1)))
  FeatureMatrix_cDS_D2 = t(scale(t(FeatureMatrix_cDS_D2)))
  NodeResult = NULL
  pb <- progress_bar$new(total = length(nodes))
  for(i in 1:length(nodes)){
    Direct1=Net[,1]%in%nodes[i]
    Direct2=Net[,2]%in%nodes[i]
    tempNetR = rbind(net[Direct1,],net[Direct2,])
    cKLD.ind = which(tempNetR[,1:2]==nodes[i],arr.ind = FALSE)
    if(nrow(tempNetR)>1){
      FeatureMatrix_temp_cds = rbind(FeatureMatrix_cDS_D1[Direct1,],FeatureMatrix_cDS_D2[Direct2,c(flipIndex_1,flipIndex_2)])
      FeatureMatrix_temp_ds = rbind(FeatureMatrix_DS[Direct1,],FeatureMatrix_DS[Direct2,c(flipIndex_1,flipIndex_2)])
      ds_dist = (dist(FeatureMatrix_temp_cds)+dist(FeatureMatrix_temp_ds))/2
      edgegroup = cutree(hclust(ds_dist),h=cutHeight)
      Raw.KLD.f = base::sum(tempNetR[,'KLD.ad.chiSquare'])
      Raw.cKLD.f = base::sum(as.matrix(tempNetR[,c( "cKLD_1.ad.chiSquare" ,  "cKLD_2.ad.chiSquare")])[cKLD.ind])
      Raw.OR.f = base::sum(as.matrix(tempNetR[,c( "cKLDorKLD.ad.chiSquare_1", "cKLDorKLD.ad.chiSquare_2")])[cKLD.ind])
      Raw.KLD.p = tapply(tempNetR[,'KLD.ad.chiSquare'], edgegroup, function(x)sum(x)/sqrt(length(x)))%>%base::sum()
      Raw.cKLD.p = tapply(as.matrix(tempNetR[,c( "cKLD_1.ad.chiSquare" ,  "cKLD_2.ad.chiSquare")])[cKLD.ind], edgegroup,
                          function(x)sum(x)/sqrt(length(x)))%>%base::sum()
      Raw.OR.p = tapply(as.matrix(tempNetR[,c( "cKLDorKLD.ad.chiSquare_1", "cKLDorKLD.ad.chiSquare_2")])[cKLD.ind], edgegroup,
                        function(x)sum(x)/sqrt(length(x)))%>%base::sum()
      NodeDiversity = edgegroup%>%unique()%>%length()
      tempNodeR = data.frame(Gene=nodes[i],
                             Degree = nrow(tempNetR),
                             NodeDiversity=NodeDiversity,
                             Raw.KLD.f = Raw.KLD.f,
                             Raw.cKLD.f = Raw.cKLD.f,
                             Raw.OR.f = Raw.OR.f,
                             Raw.KLD.p = Raw.KLD.p,
                             Raw.cKLD.p = Raw.cKLD.p,
                             Raw.OR.p = Raw.OR.p)
    }else{
      Raw.KLD.f = tempNetR[,'KLD.ad.chiSquare']
      Raw.cKLD.f = as.matrix(tempNetR[,c( "cKLD_1.ad.chiSquare" ,  "cKLD_2.ad.chiSquare")])[cKLD.ind]
      Raw.OR.f = as.matrix(tempNetR[,c( "cKLDorKLD.ad.chiSquare_1", "cKLDorKLD.ad.chiSquare_2")])[cKLD.ind]
      tempNodeR = data.frame(Gene = nodes[i],
                             Degree = nrow(tempNetR),
                             NodeDiversity=1,
                             Raw.KLD.f = Raw.KLD.f,
                             Raw.cKLD.f = Raw.cKLD.f,
                             Raw.OR.f = Raw.OR.f,
                             Raw.KLD.p = Raw.KLD.f,
                             Raw.cKLD.p = Raw.cKLD.f,
                             Raw.OR.p = Raw.OR.f)
    }
    NodeResult = rbind(NodeResult,tempNodeR)
    pb$tick()
  }
  NodeResult
}




getEdgeDviersity_ds_cds_raw <- function(net,
                                        FeatureMatrix_DS,
                                        FeatureMatrix_cDS_D1,
                                        FeatureMatrix_cDS_D2,
                                        nodes=NULL,
                                        cutHeight = 8.5){
  if(is.null(nodes)){
    nodes = unique(unlist(net[,1:2]))
  }
  nx = sqrt(ncol(FeatureMatrix_DS)/2)
  flipIndex_1 = t(matrix(1:nx^2,nx))%>%as.vector()
  flipIndex_2 = flipIndex_1+nx^2

  FeatureMatrix_DS = t(scale(t(FeatureMatrix_DS)))
  FeatureMatrix_cDS_D1 = t(scale(t(FeatureMatrix_cDS_D1)))
  FeatureMatrix_cDS_D2 = t(scale(t(FeatureMatrix_cDS_D2)))
  NodeResult = NULL
  pb <- progress_bar$new(total = length(nodes))
  for(i in 1:length(nodes)){
    Direct1=net[,1]%in%nodes[i]
    Direct2=net[,2]%in%nodes[i]
    tempNetR = rbind(net[Direct1,],net[Direct2,])
    cDiv.ind = which(tempNetR[,1:2]==nodes[i],arr.ind = FALSE)
    if(nrow(tempNetR)>1){
      FeatureMatrix_temp_cds = rbind(FeatureMatrix_cDS_D1[Direct1,],FeatureMatrix_cDS_D2[Direct2,c(flipIndex_1,flipIndex_2)])
      FeatureMatrix_temp_ds = rbind(FeatureMatrix_DS[Direct1,],FeatureMatrix_DS[Direct2,c(flipIndex_1,flipIndex_2)])
      ds_dist = (dist(FeatureMatrix_temp_cds)+dist(FeatureMatrix_temp_ds))/2
      edgegroup = cutree(hclust(ds_dist),h=cutHeight)
      tempD=data.frame(Netid = c(which(Direct1),which(Direct2)),cDiv.ind=cDiv.ind,cluster=edgegroup)



    }else{
      tempD=data.frame(Netid = c(which(Direct1),which(Direct2)),cDiv.ind=cDiv.ind,cluster=1)
    }
    NodeResult = c(NodeResult,list(tempD))
    pb$tick()
  }
  names(NodeResult) = nodes
  NodeResult
}




getKLD_cKLDnetwork_fromDS<-function(RawDS_1,
                                    RawDS_2,
                                    Network,
                                    GroupLabel,
                                    Condition_label_1 = NULL,
                                    Condition_label_2 = NULL,
                                    divergence=c('kld','jsd'),
                                    n.grid = 60,
                                    n.coarse=20,
                                    NoiseRemove = TRUE,
                                    CoarseGrain=TRUE,
                                    returnCDS = TRUE){
  if(!is.null(Condition_label_1)&!is.null(Condition_label_2)){
    uniCase <- unique(GroupLabel)%>%sort()
  }else{
    uniCase = c(Condition_label_1,Condition_label_2)
  }

  message('context A')
  ContextA_R = get_cDS_MI_DREMI_fromRawDS(RawDS = RawDS_1,
                                          Network,
                                          n.grid,
                                          n.coarse,
                                          label=uniCase[1],
                                          CoarseGrain=CoarseGrain)
  ## @@@ context B----
  message('context B')
  ContextB_R = get_cDS_MI_DREMI_fromRawDS(RawDS = RawDS_2,
                                          Network,
                                          n.grid,
                                          n.coarse,
                                          label=uniCase[2],
                                          CoarseGrain=CoarseGrain)
  NetNcol = ncol(ContextA_R$Network)
  Network_MI_DREMI = cbind(ContextA_R$Network[,(NetNcol-2):NetNcol],ContextB_R$Network[,(NetNcol-2):NetNcol])
  MI_types = colnames(Network_MI_DREMI)
  Network_MI_DREMI = Network_MI_DREMI[,c(MI_types[str_detect(MI_types,'^MI')],sort(MI_types[str_detect(MI_types,'^DREMI')]))]
  # Network = cbind(Network,Network_MI_DREMI)
  #
  if(NoiseRemove){ # Whether to remove some rasters with extremely low density
    if(CoarseGrain){
      NoNoiseCDS = removeNoiseInCDS_towDS(DS1 = ContextA_R$Density,
                                          DS2 = ContextB_R$Density,
                                          cDSList_1 = ContextA_R$cDS,
                                          cDSList_2 = ContextB_R$cDS,
                                          nx = n.coarse,
                                          min.pt = 0.5,
                                          min.density =k/((k)^2*pi)*(n.grid^2/n.coarse^2))
      ContextA_R$cDS=NoNoiseCDS$cDSList_1
      ContextB_R$cDS=NoNoiseCDS$cDSList_2
      # 重新计算DREMI
      Network[,paste('DREMI_1',uniCase[1],sep = '_')] = MI_from_ked2d_v(x = ContextA_R$cDS$cDS_1,logbase = 2,nx = n.coarse)
      Network[,paste('DREMI_2',uniCase[1],sep = '_')] = MI_from_ked2d_v(x = ContextA_R$cDS$cDS_2,logbase = 2,nx = n.coarse)
      Network[,paste('DREMI_1',uniCase[2],sep = '_')] = MI_from_ked2d_v(x = ContextB_R$cDS$cDS_1,logbase = 2,nx = n.coarse)
      Network[,paste('DREMI_2',uniCase[2],sep = '_')] = MI_from_ked2d_v(x = ContextB_R$cDS$cDS_2,logbase = 2,nx = n.coarse)
    }else{
      NoNoiseCDS = removeNoiseInCDS_towDS(DS1 = ContextA_R$Density,
                                          DS2 = ContextB_R$Density,
                                          cDSList_1 = ContextA_R$cDS,
                                          cDSList_2 = ContextB_R$cDS,
                                          nx = n.grid,
                                          min.pt = 0.5,
                                          min.density =k/((k)^2*pi))
      ContextA_R$cDS=NoNoiseCDS$cDSList_1
      ContextB_R$cDS=NoNoiseCDS$cDSList_2
      # calculate DREMI
      Network[,paste('DREMI_1',uniCase[1],sep = '_')] = MI_from_ked2d_v(x = ContextA_R$cDS$cDS_1,logbase = 2,nx = n.grid)
      Network[,paste('DREMI_2',uniCase[1],sep = '_')] = MI_from_ked2d_v(x = ContextA_R$cDS$cDS_2,logbase = 2,nx = n.grid)
      Network[,paste('DREMI_1',uniCase[2],sep = '_')] = MI_from_ked2d_v(x = ContextB_R$cDS$cDS_1,logbase = 2,nx = n.grid)
      Network[,paste('DREMI_2',uniCase[2],sep = '_')] = MI_from_ked2d_v(x = ContextB_R$cDS$cDS_2,logbase = 2,nx = n.grid)
    }


  }

  # @@@ calculate divergence----
  if(divergence[1] =='kld'){
    Network[,'Div'] = KLD_batch_matrix(Ctl = ContextA_R$cDS$DS,
                                       Pert = ContextB_R$cDS$DS,
                                       symmetrical = TRUE)

    Network[,'cDiv_D1'] = KLD_batch_matrix(Ctl = ContextA_R$cDS$cDS_1,
                                           Pert = ContextB_R$cDS$cDS_1,
                                           symmetrical = TRUE)
    Network[,'cDiv_D2'] = KLD_batch_matrix(Ctl = ContextA_R$cDS$cDS_2,
                                           Pert = ContextB_R$cDS$cDS_2,
                                           symmetrical = TRUE)
  }else{
    Network[,'Div'] = JSD_batch_matrix (Ctl = ContextA_R$cDS$DS,
                                        Pert = ContextB_R$cDS$DS)

    Network[,'cDiv_D1'] = JSD_batch_matrix (Ctl = ContextA_R$cDS$cDS_1,
                                            Pert = ContextB_R$cDS$cDS_1)
    Network[,'cDiv_D2'] = JSD_batch_matrix (Ctl = ContextA_R$cDS$cDS_2,
                                            Pert = ContextB_R$cDS$cDS_2)
  }

  # # i=5 # test
  # p1<-DensityPlot_raster(ContextA_R$CGDensity[i,])|DensityPlot_raster(ContextA_R$cDS$cDS_1[i,],addSideBar = F,transpose = F)|
  #   DensityPlot_raster(ContextA_R$cDS$cDS_2[i,],addSideBar = F,transpose = F)
  # p2<-DensityPlot_raster(ContextB_R$CGDensity[i,])|DensityPlot_raster(ContextB_R$cDS$cDS_1[i,],addSideBar = F,transpose = F)|
  #   DensityPlot_raster(ContextB_R$cDS$cDS_2[i,],addSideBar = F,transpose = F)
  # CombinePlots(list(p1,p2),ncol = 1)

  if(returnCDS){
    res  = list(Network=Network,
                DivergenceMethod=divergence[1],
                ContextA_DS = ContextA_R$Density,
                ContextB_DS = ContextB_R$Density,
                ContextA_cDS = ContextA_R$cDS,
                ContextB_cDS = ContextB_R$cDS,
                uniCase=uniCase)

  }else{
    res  = list(Network=Network,
                DivergenceMethod=divergence[1],
                ContextA_DS = ContextA_R$Density,
                ContextB_DS = ContextB_R$Density,
                uniCase = uniCase)
  }
  return(res)
}


DensityPlot_raster<-function(z,
                             nrow=NULL,
                             ncolor=100,
                             addSideBar = TRUE,
                             transpose=FALSE,
                             colors= c('black','red','yellow','white'),
                             interpolate=FALSE,
                             pt.color='black',
                             pt.alpha=1,
                             pt.size=0.1){
  if(is.vector(z)){
    if(is.null(nrow)){
      z = matrix(z,nrow = sqrt(length(z)))
    }else{
      z = matrix(z,nrow = nrow)
    }
  }
  if(transpose){
    z = t(z)
  }
  grid = expand.grid(1:nrow(z),1:ncol(z))
  Pdata = data.frame(gridx = grid[,1],gridy=grid[,2],z = matrix(z,ncol = 1))
  #Pdata_line = data.frame(x=1:nrow(z),y=apply(z, 1,which.max)%>%as.vector(),stringsAsFactors = F)
  if(!addSideBar){
    ggplot() +
      geom_raster(mapping = aes(gridx , gridy,fill = z), data = Pdata,interpolate = interpolate)+
      scale_fill_gradientn(colours = colorRampPalette(colors, interpolate = c("linear"))(ncolor),
                           breaks=c(range(z)),labels=c('low','high'))+  labs(fill='density')+
      geom_point(data = Pdata,mapping = aes(gridx,gridy),color=pt.color,alpha=pt.alpha,size=pt.size)+
      theme_cowplot_i()+AddBox()+xlab('G1')+ylab('G2')+AddBox()#+guides(fill=guide_legend(title="density"))
  }else{
    rightSideBar = data.frame(x=0,y=1:ncol(z),z=colSums(z))
    topSideBar = data.frame(x=1:nrow(z),y=0,z=rowSums(z))
    ggplot() +
      geom_raster(mapping = aes(gridx , gridy,fill = z), data = Pdata,interpolate = interpolate)+
      scale_fill_gradientn(colours = colorRampPalette(colors, interpolate = c("linear"))(ncolor),
                           breaks=c(range(z)),labels=c('low','high'))+ labs(fill='density')+
      geom_point(data = Pdata,mapping = aes(gridx,gridy),color=pt.color,alpha=pt.alpha,size=pt.size)+#guides(fill=guide_legend(title="density"))+
      theme_cowplot_i()+NoAxes2(keep.axis.title =  T)+xlab('G1')+ylab('G2')+
      ggnewscale::new_scale_fill()+ # right
      geom_raster(data = rightSideBar,mapping = aes(x,y,fill=z),interpolate=T)+
      geom_rect(data = data.frame(x=0,y=0),mapping = aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y+0.5, ymax = nrow(z)+0.5),
                fill='transparent',color='grey50')+
      scale_fill_gradientn(colours = colorRampPalette(colors, interpolate = c("linear"))(ncolor))+guides(fill='none')+
      ggnewscale::new_scale_fill()+ # top
      geom_raster(data = topSideBar,mapping = aes(x,y,fill=z),interpolate=T)+
      geom_rect(data = data.frame(x=0,y=0),mapping = aes(xmin = x +0.5, xmax = nrow(z)+0.5, ymin = y-0.5, ymax = y+0.5),
                fill='transparent',color='grey50')+
      scale_fill_gradientn(colours = colorRampPalette(colors, interpolate = c("linear"))(ncolor))+guides(fill='none')+
      geom_rect(data = data.frame(x=0.5,y=0.5),mapping = aes(xmin = x , xmax = nrow(z)+0.5, ymin = y, ymax =  nrow(z)+0.5),
                fill='transparent',color='grey30')
  }


}



DensityPlot_raster <- function(z,
         nrow=NULL,
         ncolor=100,
         addSideBar = TRUE,
         transpose=FALSE,
         colors= c('black','red','yellow','white'),
         interpolate=FALSE,
         pt.color='black',
         pt.alpha=1,
         pt.size=0.1){
  if(is.vector(z)){
    if(is.null(nrow)){
      z = matrix(z,nrow = sqrt(length(z)))
    }else{
      z = matrix(z,nrow = nrow)
    }
  }
  if(transpose){
    z = t(z)
  }
  grid = expand.grid(1:nrow(z),1:ncol(z))
  Pdata = data.frame(gridx = grid[,1],gridy=grid[,2],z = matrix(z,ncol = 1))
  #Pdata_line = data.frame(x=1:nrow(z),y=apply(z, 1,which.max)%>%as.vector(),stringsAsFactors = F)
  if(!addSideBar){
    ggplot() +
      geom_raster(mapping = aes(gridx , gridy,fill = z), data = Pdata,interpolate = interpolate)+
      scale_fill_gradientn(colours = colorRampPalette(colors, interpolate = c("linear"))(ncolor),
                           breaks=c(range(z)),labels=c('low','high'))+  labs(fill='density')+
      geom_point(data = Pdata,mapping = aes(gridx,gridy),color=pt.color,alpha=pt.alpha,size=pt.size)+
      theme_cowplot_i()+AddBox()+xlab('G1')+ylab('G2')+AddBox()#+guides(fill=guide_legend(title="density"))
  }else{
    rightSideBar = data.frame(x=0,y=1:ncol(z),z=colSums(z))
    topSideBar = data.frame(x=1:nrow(z),y=0,z=rowSums(z))
    ggplot() +
      geom_raster(mapping = aes(gridx , gridy,fill = z), data = Pdata,interpolate = interpolate)+
      scale_fill_gradientn(colours = colorRampPalette(colors, interpolate = c("linear"))(ncolor),
                           breaks=c(range(z)),labels=c('low','high'))+ labs(fill='density')+
      geom_point(data = Pdata,mapping = aes(gridx,gridy),color=pt.color,alpha=pt.alpha,size=pt.size)+#guides(fill=guide_legend(title="density"))+
      theme_cowplot_i()+NoAxes2(keep.axis.title =  T)+xlab('G1')+ylab('G2')+
      ggnewscale::new_scale_fill()+ # right
      geom_raster(data = rightSideBar,mapping = aes(x,y,fill=z),interpolate=T)+
      geom_rect(data = data.frame(x=0,y=0),mapping = aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y+0.5, ymax = nrow(z)+0.5),
                fill='transparent',color='grey50')+
      scale_fill_gradientn(colours = colorRampPalette(colors, interpolate = c("linear"))(ncolor))+guides(fill='none')+
      ggnewscale::new_scale_fill()+ # top
      geom_raster(data = topSideBar,mapping = aes(x,y,fill=z),interpolate=T)+
      geom_rect(data = data.frame(x=0,y=0),mapping = aes(xmin = x +0.5, xmax = nrow(z)+0.5, ymin = y-0.5, ymax = y+0.5),
                fill='transparent',color='grey50')+
      scale_fill_gradientn(colours = colorRampPalette(colors, interpolate = c("linear"))(ncolor))+guides(fill='none')+
      geom_rect(data = data.frame(x=0.5,y=0.5),mapping = aes(xmin = x , xmax = nrow(z)+0.5, ymin = y, ymax =  nrow(z)+0.5),
                fill='transparent',color='grey30')
  }


}
#' sea_fit_fix
#'
#' @param Dsize
#' @param RawScores
#' @param Dstd
#' @param winLength
#' @param gap
#' @param coarse
#' @param sd_1
#' @param mean_1
#'
#' @return
#' @export
#'
#' @examples
sea_fit_fix<-function(Dsize,RawScores,Dstd = NULL,
                      winLength = 10,gap = 0.6,
                      coarse = FALSE,sd_1=NULL,mean_1=NULL){
  #

  # require(sn)
  if(is.null(sd_1)){
    sd_1=sd(RawScores/Dsize)
  }
  if(is.null(mean_1)){
    mean_1=mean(RawScores/Dsize)
  }
  Dsize = Dsize-1
  #
  if(is.null(Dstd)){
    if(min(table(Dsize))>20){
      winLength = min(table(Dsize))
      gap = 1
      coarse = FALSE
    }
    Meandata = MeanSD_By_slidingWindows(ranks = Dsize,scores = RawScores,winLength = winLength,gap = gap,coarse = coarse)
    colnames(Meandata) = c('MeaS','MeaM','StdM')
  }else{
    Meandata= data.frame(MeaS =Dsize, MeaM = RawScores,StdM = Dstd)
  }
  # define functions
  # calculate Zscore
  cal_Zscore<-function(RawScore,Dsize,Meanfit,Stdfit){
    Zscores = (RawScore- predict(Meanfit,data.frame(MeaS = Dsize)))/
      predict(Stdfit,data.frame(MeaS = Dsize))
    Zscores = as.numeric(Zscores)
    Zscores
  }
  # calculate Pvalue
  cal_Zscore_Pvalue<-function(RawScore,Dsize,ParmetersInput){
    Zscores = (RawScore- predict(ParmetersInput$Meanfit,data.frame(MeaS = Dsize)))/
      predict(ParmetersInput$Stdfit,data.frame(MeaS = Dsize))
    Zscores = as.numeric(Zscores)
    Pvalues = 1-sn::psn(Zscores,xi=ParmetersInput$ParmetersZsD[1],
                    omega=ParmetersInput$ParmetersZsD[2],
                    alpha=ParmetersInput$ParmetersZsD[3]) #计算
    Pvalues2 = sn::dsn(Zscores,xi=ParmetersInput$ParmetersZsD[1],
                   omega=ParmetersInput$ParmetersZsD[2],
                   alpha=ParmetersInput$ParmetersZsD[3]) #计算
    Pvalues[Pvalues==0] = Pvalues2[Pvalues==0]
    data.frame(Zscores = Zscores,Pvalues = Pvalues)
  }

  # define fit function
  # f1=function(x,a, b,c) {a*I(log10(x)^b)+c};
  fmean=function(x,a, b) {a*I(x^b)+mean_1};
  fsd=function(x,a, b) {a*I(x^b)+sd_1};
  # size ~ sd fitting
  inital_coef <- lm(Meandata$StdM ~ Meandata$MeaS)

  resultSDPa=minpack.lm::nlsLM(StdM~fsd(MeaS,a,b),
                               data=Meandata,start =
                                 list(a= as.numeric(coef(inital_coef)[2]),b=3.36))
  for (i in 1:150){# 多迭代几次
    resultSDPa=minpack.lm::nlsLM(StdM~fsd(MeaS,a,b), data=Meandata,start = coef(resultSDPa))
  }

  resultSDPa=minpack.lm::nlsLM(StdM~fsd(MeaS,a,b),
                               data=Meandata,
                               start = coef(resultSDPa)) # 如果未出现warning 说明差不多
  summary(resultSDPa)

  # size ~ mean fitting
  inital_coef <- lm(Meandata$MeaM ~ Meandata$MeaS)
  resultMePa = minpack.lm::nlsLM(MeaM~fmean(MeaS,a,b),
                                 data=Meandata,start =
                                   list(a= as.numeric(coef(inital_coef)[2]),b=3.36))

  for (i in 1:150){# Iterate a few more times
    resultMePa=minpack.lm::nlsLM(MeaM~fmean(MeaS,a,b), data=Meandata,start = coef(resultMePa))
  }

  resultMePa=minpack.lm::nlsLM(MeaM~fmean(MeaS,a,b),
                               data=Meandata,
                               start = coef(resultMePa)) # 如果未出现warning 说明差不多
  summary(resultMePa)

  ## plot fit @@@@@@
  Meandata$PredictStdM = predict(resultSDPa,Meandata$MeaS)%>%as.numeric()
  Meandata$PredictMeanM = predict(resultMePa,Meandata$MeaS)%>%as.numeric()

  p1<-ggplot(Meandata,aes(MeaS+1,StdM))+geom_point(size=5)+
    geom_line(aes(y=PredictStdM),size=2,color='red')+
    theme_cowplot_i()+AddBox()+ggtitle('Degree ~ Sd')+
    xlab('Degree')+ylab('Std') # Size ~ sd

  p2<-ggplot(Meandata,aes(MeaS+1,MeaM))+geom_point(size = 5)+
    geom_line(aes(y=PredictMeanM),size=2,color='red')+
    theme_cowplot_i()+AddBox()+ggtitle('Degree ~ Mean')+
    xlab('Degree')+ylab('Mean') # Size ~ mean
  #

  #  Z score distribution for random data
  Inputdata =data.frame(Dsize=Dsize,RawScore=RawScores)
  Inputdata$Zscore = (Inputdata$RawScore- predict(resultMePa,data.frame(MeaS = Inputdata$Dsize)))/
    predict(resultSDPa,data.frame(MeaS = Inputdata$Dsize))
  Inputdata$Zscore = as.numeric(Inputdata$Zscore) # important for plot
  #hist(Inputdata$Zscore)

  # Skewed distribution fitting
  ZsDistribution <- sn::selm(Zscore ~ 1, family = "SN",data = Inputdata)
  # coef(ZsDistribution)
  ParmetersZsD = ZsDistribution@param$dp %>%as.numeric() # 分布的系数
  #Pa = ParmetersZsD
  # Pvalue<- 1-sn::psn(seq(-1,1,0.1),
  #                xi=ParmetersZsD[1],
  #                omega=ParmetersZsD[2],
                 # alpha=ParmetersZsD[3]) #计算
  # plot
  p3<-ggplot(data = Inputdata, aes(Zscore)) + # geom_line()
    geom_histogram(mapping = aes(y = after_stat(density)),colour = "white", fill = "black",bins = 50) +
    stat_function(fun = sn::dsn, args = list(omega = ParmetersZsD[2], alpha = ParmetersZsD[3] , xi = ParmetersZsD[1]),col = "red",size = 2)+
    theme_cowplot_i() +AddBox()+ggtitle('Z scores')# y = ..density.. 不能在ggplot被调用只能在geom_histogram


  ParmetersInput = list(Meanfit = resultMePa,
                        Stdfit = resultSDPa,
                        ParmetersZsD = ParmetersZsD)
  p4 <- p1|p2|p3


  return(list(cal_Zscore_Pvalue=cal_Zscore_Pvalue,ParmetersInput= ParmetersInput,
              SizeVsMean = p2,SizeVsStd = p1,ZscoreD = p3,CombinePs = p4))

}

#' sea_fit_fix_st
#'
#' @param Dsize
#' @param RawScores
#' @param Dstd
#' @param winLength
#' @param gap
#' @param coarse
#' @param sd_1
#' @param mean_1
#'
#' @return
#' @export
#'
#' @examples
sea_fit_fix_st<-function(Dsize,RawScores,Dstd = NULL,
                      winLength = 10,gap = 0.6,
                      coarse = FALSE,sd_1=NULL,mean_1=NULL){
  #

  # require(sn)
  if(is.null(sd_1)){
    sd_1=sd(RawScores/Dsize)
  }
  if(is.null(mean_1)){
    mean_1=mean(RawScores/Dsize)
  }
  Dsize = Dsize-1
  #
  if(is.null(Dstd)){
    if(min(table(Dsize))>20){
      winLength = min(table(Dsize))
      gap = 1
      coarse = FALSE
    }
    Meandata = MeanSD_By_slidingWindows(ranks = Dsize,scores = RawScores,winLength = winLength,gap = gap,coarse = coarse)
    colnames(Meandata) = c('MeaS','MeaM','StdM')
  }else{
    Meandata= data.frame(MeaS =Dsize, MeaM = RawScores,StdM = Dstd)
  }
  # define functions
  # calculate Zscore
  cal_Zscore<-function(RawScore,Dsize,Meanfit,Stdfit){
    Zscores = (RawScore- predict(Meanfit,data.frame(MeaS = Dsize)))/
      predict(Stdfit,data.frame(MeaS = Dsize))
    Zscores = as.numeric(Zscores)
    Zscores
  }
  # calculate Pvalue
  cal_Zscore_Pvalue<-function(RawScore,Dsize,ParmetersInput){
    Zscores = (RawScore- predict(ParmetersInput$Meanfit,data.frame(MeaS = Dsize)))/
      predict(ParmetersInput$Stdfit,data.frame(MeaS = Dsize))
    Zscores = as.numeric(Zscores)
    Pvalues <-  1-sn::pst(Zscores,xi=ParmetersInput$ParmetersZsD[1],
                          omega=ParmetersInput$ParmetersZsD[2],
                          alpha=ParmetersInput$ParmetersZsD[3],
                          nu =ParmetersInput$ParmetersZsD[4] ) #计算
    Pvalues2 <-  sn::dst(Zscores,xi=ParmetersInput$ParmetersZsD[1],
                         omega=ParmetersInput$ParmetersZsD[2],
                         alpha=ParmetersInput$ParmetersZsD[3],
                         nu =ParmetersInput$ParmetersZsD[4] ) #计算
    Pvalues[Pvalues==0] = Pvalues2[Pvalues==0]
    data.frame(Zscores = Zscores,Pvalues = Pvalues)
  }


  # define fit function
  # f1=function(x,a, b,c) {a*I(log10(x)^b)+c};
  fmean=function(x,a, b) {a*I(x^b)+mean_1};
  fsd=function(x,a, b) {a*I(x^b)+sd_1};
  # size ~ sd fitting
  inital_coef <- lm(Meandata$StdM ~ Meandata$MeaS)

  resultSDPa=minpack.lm::nlsLM(StdM~fsd(MeaS,a,b),
                               data=Meandata,start =
                                 list(a= as.numeric(coef(inital_coef)[2]),b=1))
  for (i in 1:150){# 多迭代几次
    resultSDPa=minpack.lm::nlsLM(StdM~fsd(MeaS,a,b), data=Meandata,start = coef(resultSDPa))
  }

  resultSDPa=minpack.lm::nlsLM(StdM~fsd(MeaS,a,b),
                               data=Meandata,
                               start = coef(resultSDPa)) # 如果未出现warning 说明差不多
  summary(resultSDPa)

  # size ~ mean fitting
  inital_coef <- lm(Meandata$MeaM ~ Meandata$MeaS)
  resultMePa = minpack.lm::nlsLM(MeaM~fmean(MeaS,a,b),
                                 data=Meandata,start =
                                   list(a= as.numeric(coef(inital_coef)[2]),b=1))

  for (i in 1:150){# Iterate a few more times
    resultMePa=minpack.lm::nlsLM(MeaM~fmean(MeaS,a,b), data=Meandata,start = coef(resultMePa))
  }

  resultMePa=minpack.lm::nlsLM(MeaM~fmean(MeaS,a,b),
                               data=Meandata,
                               start = coef(resultMePa)) # 如果未出现warning 说明差不多
  summary(resultMePa)

  ## plot fit @@@@@@
  Meandata$PredictStdM = predict(resultSDPa,Meandata$MeaS)%>%as.numeric()
  Meandata$PredictMeanM = predict(resultMePa,Meandata$MeaS)%>%as.numeric()

  p1<-ggplot(Meandata,aes(MeaS+1,StdM))+geom_point(size=5)+
    geom_line(aes(y=PredictStdM),size=2,color='red')+
    theme_cowplot_i()+AddBox()+ggtitle('Degree ~ Sd')+
    xlab('Degree')+ylab('Std') # Size ~ sd

  p2<-ggplot(Meandata,aes(MeaS+1,MeaM))+geom_point(size = 5)+
    geom_line(aes(y=PredictMeanM),size=2,color='red')+
    theme_cowplot_i()+AddBox()+ggtitle('Degree ~ Mean')+
    xlab('Degree')+ylab('Mean') # Size ~ mean
  #

  #  Z score distribution for random data
  Inputdata =data.frame(Dsize=Dsize,RawScore=RawScores)
  Inputdata$Zscore = (Inputdata$RawScore- predict(resultMePa,data.frame(MeaS = Inputdata$Dsize)))/
    predict(resultSDPa,data.frame(MeaS = Inputdata$Dsize))
  Inputdata$Zscore = as.numeric(Inputdata$Zscore) # important for plot
  #hist(Inputdata$Zscore)

  # Skewed distribution fitting
  fit0 <- sn::selm(Zscore ~ 1, data = Inputdata, family="SN")
  dp0  <- coef(fit0, param.type="DP")
  # 给一个中等的自由度起始值 nu = 5
  start0 <- c(dp0, nu = 3)
  ZsDistribution <- sn::selm(Zscore ~ 1, family = "ST",
                             data = Inputdata,start = start0)
  # ZsDistribution <- sn::selm(Zscore ~ 1, family = "ST",data = Inputdata)
  # coef(ZsDistribution)
  ParmetersZsD = ZsDistribution@param$dp %>%as.numeric() # 分布的系数
  #Pa = ParmetersZsD
  # Pvalue<- 1-sn::psn(seq(-1,1,0.1),
  #                    xi=ParmetersZsD[1],
  #                    omega=ParmetersZsD[2],
  #                    alpha=ParmetersZsD[3]) #计算
  # plot
  p3<-ggplot(data = Inputdata, aes(Zscore)) + # geom_line()
    geom_histogram(mapping = aes(y = after_stat(density)),colour = "white", fill = "black",bins = 50) +
    stat_function(fun = sn::dst, args = list(omega = ParmetersZsD[2],
                                             alpha = ParmetersZsD[3] ,
                                             xi = ParmetersZsD[1],
                                             nu=ParmetersZsD[4]),col = "red",size = 2)+
    theme_cowplot_i() +AddBox()+ggtitle('Z scores')# y = ..density.. 不能在ggplot被调用只能在geom_histogram


  ParmetersInput = list(Meanfit = resultMePa,
                        Stdfit = resultSDPa,
                        ParmetersZsD = ParmetersZsD)
  p4 <- p1|p2|p3


  return(list(cal_Zscore_Pvalue=cal_Zscore_Pvalue,ParmetersInput= ParmetersInput,
              SizeVsMean = p2,SizeVsStd = p1,ZscoreD = p3,CombinePs = p4))

}

fit_selm_distrubution <- function(x){
  cal_Pvalue<-function(x,ParmetersInput){
    Pvalues <-  1-sn::psn(x,xi=ParmetersInput[1],
                    omega=ParmetersInput[2],
                    alpha=ParmetersInput[3]) #计算
    Pvalues2 <-  sn::dsn(x,xi=ParmetersInput[1],
                   omega=ParmetersInput[2],
                   alpha=ParmetersInput[3]) #计算
    Pvalues[Pvalues==0] = Pvalues2[Pvalues==0]
    data.frame(x = x,Pvalues = Pvalues)
  }
  Inputdata <- data.frame(Zscore=x)
  ZsDistribution <- sn::selm(Zscore ~ 1, family = "SN",
                         data = Inputdata)
  ParmetersZsD = ZsDistribution@param$dp %>%as.numeric() # 分布的系数
  p<-ggplot(data = Inputdata, aes(Zscore)) + # geom_line()
    geom_histogram(mapping = aes(y = after_stat(density)),colour = "white", fill = "black",bins = 50) +
    stat_function(fun = sn::dsn, args = list(omega = ParmetersZsD[2], alpha = ParmetersZsD[3] , xi = ParmetersZsD[1]),col = "red",size = 2)+
    theme_cowplot_i() +AddBox()+ggtitle('Z scores')# y = ..density.. 不能在ggplot被调用只能在geom_histogram
  return(list(cal_Pvalue=cal_Pvalue,
              ParmetersInput=ParmetersZsD,
              reslutP=p))

}

fit_selm_distrubution_ST <- function(x){
  cal_Pvalue<-function(x,ParmetersInput){
    Pvalues <-  1-sn::pst(x,xi=ParmetersInput[1],
                          omega=ParmetersInput[2],
                          alpha=ParmetersInput[3],
                          nu =ParmetersInput[4] ) #计算
    Pvalues2 <-  sn::dst(x,xi=ParmetersInput[1],
                         omega=ParmetersInput[2],
                         alpha=ParmetersInput[3],
                         nu =ParmetersInput[4] ) #计算
    Pvalues[Pvalues==0] = Pvalues2[Pvalues==0]
    data.frame(x = x,Pvalues = Pvalues)
  }
  Inputdata <- data.frame(Zscore=x)
  fit0 <- sn::selm(Zscore ~ 1, data = Inputdata, family="SN")
  dp0  <- coef(fit0, param.type="DP")
  # 给一个中等的自由度起始值 nu = 5
  start0 <- c(dp0, nu = 3)
  ZsDistribution <- sn::selm(Zscore ~ 1, family = "ST",
                             data = Inputdata,start = start0)
  ParmetersZsD = ZsDistribution@param$dp %>%as.numeric() # 分布的系数
  p<-ggplot(data = Inputdata, aes(Zscore)) + # geom_line()
    geom_histogram(mapping = aes(y = after_stat(density)),colour = "white", fill = "black",bins = 50) +
    stat_function(fun = sn::dst, args = list(omega = ParmetersZsD[2], alpha = ParmetersZsD[3] ,
                                             xi = ParmetersZsD[1],nu=ParmetersZsD[4]),col = "red",size = 2)+
    theme_cowplot_i() +AddBox()+ggtitle('Z scores')# y = ..density.. 不能在ggplot被调用只能在geom_histogram
  return(list(cal_Pvalue=cal_Pvalue,
              ParmetersInput=ParmetersZsD,
              reslutP=p))

}


MeanSD_By_slidingWindows <- function(ranks,
         scores,
         winLength=10,
         gap=0.5,coarse=TRUE){
  if(coarse){
    coarseD = aggregate(scores,by=list(ranks),FUN=mean)
    scores = coarseD[,2]
    ranks = coarseD[,1]
  }
  scores2 = scores
  ranks2 = ranks
  scores2[rank(ranks,ties.method =  "random")] = scores2
  ranks2[rank(ranks,ties.method =  "random")] = ranks2

  gapNum = ceiling(winLength*gap)
  slidTimes = floor(length(ranks)/gapNum)
  Reslut = NULL
  for(i in 1:slidTimes){
    if((gapNum*(i-1)+winLength)>length(scores2)){
      break
    }
    scores3 <- scores2[(gapNum*(i-1)+1):(gapNum*(i-1)+winLength)]
    indices <- !detect_outliers_iqr(scores3)$indices
    MeanD = mean(scores3[indices])
    StdD = sd(scores3[indices])
    sizeD=mean(ranks2[(gapNum*(i-1)+1):(gapNum*(i-1)+winLength)][indices])
    temp=data.frame(size = sizeD,MeanD = MeanD,StdD = StdD)
    Reslut = rbind(Reslut,temp)
  }
  Reslut

}

detect_outliers_iqr <- function(x, coef = 1.5) {
  # 只对数值部分计算，忽略 NA
  x_num <- x[!is.na(x)]

  Q1   <- quantile(x_num, 0.25, names = FALSE)
  Q3   <- quantile(x_num, 0.75, names = FALSE)
  IQRv <- Q3 - Q1

  lower <- Q1 - coef * IQRv
  upper <- Q3 + coef * IQRv

  # 找出所有离群点的位置（相对于原始向量）
  is_outlier <- (x < lower) | (x > upper)
  out_idx    <- is_outlier & !is.na(is_outlier)
  out_vals   <- x[out_idx]

  list(
    indices = out_idx,
    values  = out_vals,
    lower   = lower,
    upper   = upper
  )
}


