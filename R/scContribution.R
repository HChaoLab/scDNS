#' scContribution
#'
#' @param scDNSobject scDNSobject
#' @param nx An integer value representing the number of grid (def:20)
#' @param rmZero Logical value indicating whether to exclude 0 values when calculating the joint probability density (def: FALSE)
#' @param topGene  An integer value representing The number of top genes is used to calculate the contribution of cells (def:100).The parameters will be invalid when the sigGene parameter is not NULL
#' @param sigGene Genes of interest are used to calculate the contribution of individual cells(def: NULL).
#'
#' @return
#' @keywords internal
#'
#' @examples
scContribution<-function(scDNSobject,
                          nx=20,
                          rmZero = F,
                          topGene=100,
                          sigGene=NULL){

  Zscores <- scDNSobject@Zscore

  if(is.null(sigGene)){
    sigGene = Zscores$Gene[top_n_vector(Zscores$Zscores.Plus,topGene,returnLogical = T)]

  }

  toChiSquareMat <- function(rs,bias,DistrbutionList){
    Pvalues <-  rs
    uniBias <- unique(bias)
    Num_Model <- names(DistrbutionList)
    for(i in 1:length(uniBias)){
      ind=getCloseseData(data = as.numeric(Num_Model),uniBias[i],returnIndex = T)
      Pvalues[bias==uniBias[i]] = DistrbutionList[[ind]]$cal_Pvalue(x=rs[bias==uniBias[i]],
                                                                    ParmetersInput = DistrbutionList[[ind]]$ParmetersInput)$Pvalues
    }

    Pvalues[Pvalues==0]=1e-312
    chiSquare = sqrt(-log10(Pvalues))
    chiSquare
  }

  pseudoGridNum = nx*nx # Used to calculate adjusted Div and cDiv.
  Network = scDNSobject@Network
  Network$id = 1:nrow(Network) # used to label the exracted sub-nework
  Network_sub = getSubNetByNode(Network,sigGene) # Extract sub-networks that cover significantly perturbed genes.
  # Network_sub$id = 1:nrow(Network_sub)
  ExpData <- scDNSobject@data
  ExpData = ExpData[unlist(Network_sub[,1:2])%>%unique(),]
  ExpData1_20 = mappingMinMaxRow(ExpData,minD = 1,maxD = nx,rmZero = rmZero) # mapping data
  ExpData1_20 = round(ExpData1_20,0) #

  SourceID = ExpData1_20[Network_sub[,1],] # the id of source or target nodes
  TargetID = ExpData1_20[Network_sub[,2],]
  filteredNetid = Network_sub$id
  ContextA_R <- getConditionalDenstiy(DS = scDNSobject@JDensity_A[filteredNetid,
                                                                  , drop = F])
  ContextB_R <- getConditionalDenstiy(DS = scDNSobject@JDensity_B[filteredNetid,
                                                                  , drop = F])
  Div.mat = JSD_batch_matrix_rawMt(Ctl = ContextA_R$DS,
                                   Pert = ContextB_R$DS) # Local divergence.
  Div.mat <- abs(Div.mat)
  # Div.mat[,1]=0
  cDiv.mat_D1 = JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS_1,
                                       Pert = ContextB_R$cDS_1)
  cDiv.mat_D1 <- abs(cDiv.mat_D1)
  cDiv.mat_D2 = JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS_2,
                                       Pert = ContextB_R$cDS_2)
  cDiv.mat_D2 <- abs(cDiv.mat_D2)
  NetInd = sub2ind(SourceID,TargetID,nrow = nx,ncol = nx) # index in the grid net

  NetInd_t = t(NetInd+(1:nrow(NetInd)-1)*ncol(Div.mat)) # index for multiple edge
  #


  # C3Cell <- WhichCells_HC(iAT2_data,idents = '3',group.by = 'RNA_snn_res.0.5')
  # Network_sub_CD9 = getSubNetByNode(Network_sub,Node = 'CD9')
  # Network_sub_CD9 = Network_sub_CD9[Network_sub_CD9$Symbol.2%in%'CD9',]
  # DensityPlot_raster(colSums2(Div.mat[Network_sub_CD9$id,]/rowSums2(Div.mat[Network_sub_CD9$id,])))
  # Network_sub_CD9 = getSubNetByNode(Network_sub,Node = 'CD9')
  # Network_sub_CD9 = Network_sub_CD9[Network_sub_CD9$Symbol.1%in%'CD9',]
  # DensityPlot_raster(colSums2(Div.mat[Network_sub_CD9$id,]/rowSums2(Div.mat[Network_sub_CD9$id,])))
  #
  # CD9Sum <- colSums2(Div.mat[Network_sub_CD9$id,]/rowSums2(Div.mat[Network_sub_CD9$id,]))
  #
  # SourceID = ExpData1_20[Network_sub_CD9[,1],C3Cell] # the id of source or target nodes
  # TargetID = ExpData1_20[Network_sub_CD9[,2],C3Cell]
  #
  # NetInd = sub2ind(SourceID,TargetID,nrow = nx,ncol = nx) # index in the grid net
  # hist(CD9Sum[NetInd])

  # Div
  # Div.mat_t = ZscoreAD(Div.mat_t,
  #                       sdModel = NEAModel$ad.dropModel$AdModelList_2$fitsd.Div.fc$predictP,
  #                       meanModel = NEAModel$ad.dropModel$AdModelList_2$fitMean.Div.fc$predictP,
  #                      Network_sub$npoint_log)
  # Div.mat_t[Div.mat_t<0]=0
  # # min(Div.mat_t)
  # # Div.mat_t = abs(Div.mat_t)
  # Div.mat_t = Div.mat_t/abs(Network_sub$Div)
  # Div.mat_t = Div.mat_t*Network_sub$Div.chiSquare.LR
  #
  # Network_sub$LR = pmin(Network_sub$LR_taret,Network_sub$LR_source)
  #
  # zs1 = sapply(1:nrow(Div.mat_t),  function(x)cal_Zscore_Pvalue(RawScore = Div.mat_t,
  #                                                         Dsize = Network_sub$LR,
  #                                                         ParmetersInput = NEAModel$ZscoreFit.Plus$ParmetersInput)$Zscores)
  # Dsize <- matrix(Network_sub$LR,ncol = ncol(Div.mat_t),nrow = nrow(Div.mat_t))
  #
  # zs1 = (RawScore- predict(NEAModel$ZscoreFit.Plus$ParmetersInput$Meanfit,data.frame(MeaS = Dsize)))/
  #   predict(NEAModel$ZscoreFit.Plus$ParmetersInput$Stdfit,data.frame(MeaS = Dsize))
  #
  # zs1[zs1<0] = 0 # Zscores
  # Div.mat_t = zs1

  # Div
  Div.mat[ContextA_R$DS<ContextB_R$DS] = -Div.mat[ContextA_R$DS<ContextB_R$DS]
  Div.mat[,1] <- rowMeans2(Div.mat[,c(3,4,31,34,23,24)])
  Div.mat[,2] <- rowMeans2(Div.mat[,c(3,4,31,34,23,24)])
  Div.mat[,21] <- rowMeans2(Div.mat[,c(3,4,31,34,23,24)])
  Div.mat[,22] <- rowMeans2(Div.mat[,c(3,4,31,34,23,24)])
  Div.mat_t = t(Div.mat)
  Div.mat_t = matrix(Div.mat_t[NetInd_t],nrow = ncol(NetInd))
  Div.mat_t = t(Div.mat_t)  # row is network,col is cells,
  # Div.mat_t[Div.mat_t<0]=0


  # GeneInteraion(DSP_scDNSob,Nodes = 'DSP',subEdgeID = 1,fillColorData = Div.mat_t[1,])+
  #   scale_fill_gradientn(colours = rev(c("black", "red", "white", "purple")))|
  #   scDNS::DensityPlot_raster(Div.mat[1,],colors =  rev(c("black", "red", "white","purple" )))

  # Data2Color_Seq()
  #
  # Div.mat_t[,colnames(iAT2_data)%in%Mutcells] = -Div.mat_t[,colnames(iAT2_data)%in%Mutcells]
  # colnames(Div.mat_t) <- colnames(iAT2_data)
  # rownames(Div.mat_t) <- paste(Network_sub$Symbol.1,Network_sub$Symbol.2,sep = '-')
  # PPIname <- paste(Network_sub$Symbol.1,Network_sub$Symbol.2,sep = '-')
  # PPIscDNS <- CreateAssayObject(data = Div.mat_t)
  # iAT2_data[['PPIscDNS']] <- PPIscDNS
  # #
  # FeaturePlot(iAT2_data,features = PPIname)&scale_color_gradientn(colours =  c('darkred','white','darkblue'))&shortAxixSeurat()
  #

  # Div.mat_t = Div.mat_t/rowSums2(cDiv.mat_t_d2)
  Div.mat_t = Div.mat_t/Network_sub$Div
  Div.mat_t = Div.mat_t*Network_sub$Div.chiSquare.LR
  # cDiv
  # direction 1--
  cDiv.mat_D1[ContextA_R$cDS_1<ContextB_R$cDS_1] = -cDiv.mat_D1[ContextA_R$cDS_1<ContextB_R$cDS_1]
  cDiv.mat_t_d1 = t(cDiv.mat_D1)
  cDiv.mat_t_d1 = matrix(cDiv.mat_t_d1[NetInd_t],nrow = ncol(NetInd))
  cDiv.mat_t_d1 = t(cDiv.mat_t_d1)

  # cDiv.mat_t_d1[cDiv.mat_t_d1<0]=0
  # cDiv.mat_t_d1 = abs(cDiv.mat_t_d1)
  # cDiv.mat_t_d1 = cDiv.mat_t_d1/rowSums2(cDiv.mat_t_d1)
  cDiv.mat_t_d1 = cDiv.mat_t_d1/Network_sub$cDiv_D1
  cDiv.mat_t_d1 = cDiv.mat_t_d1*Network_sub$cDiv_D1.chiSquare.LR


  # direction 2
  cDiv.mat_D2[ContextA_R$cDS_2<ContextB_R$cDS_2] = -cDiv.mat_D2[ContextA_R$cDS_2<ContextB_R$cDS_2]
  cDiv.mat_t_d2 = t(cDiv.mat_D2)
  cDiv.mat_t_d2 = matrix(cDiv.mat_t_d2[NetInd_t],nrow = ncol(NetInd))
  cDiv.mat_t_d2 = t(cDiv.mat_t_d2)
  # cDiv.mat_t_d2 = toFC(cDiv.mat_t_d2*pesudoGridNum,
  #                      cova = Network_sub$npoint_log,
  #                      meanFit = NEAModel$ad.dropModel$AdModelList_1$fitMean.cDiv.r)
  # cDiv.mat_t_d2 = ZscoreAD(cDiv.mat_t_d2,
  #                      sdModel = NEAModel$ad.dropModel$AdModelList_2$fitsd.cDiv.fc$predictP,
  #                      meanModel = NEAModel$ad.dropModel$AdModelList_2$fitMean.cDiv.fc$predictP,
  #                      Network_sub$npoint_log)
  # cDiv.mat_t_d2[cDiv.mat_t_d2<0]=0
  # cDiv.mat_t_d2 = abs(cDiv.mat_t_d2)
  # cDiv.mat_t_d2 = cDiv.mat_t_d2/rowSums2(cDiv.mat_t_d2)
  cDiv.mat_t_d2 = cDiv.mat_t_d2/Network_sub$cDiv_D2
  cDiv.mat_t_d2 = cDiv.mat_t_d2*Network_sub$cDiv_D2.chiSquare.LR


  RS=matrix(0,nrow = unlist(Network_sub[,1:2])%>%unique()%>%length(),ncol =ncol(Div.mat_t) )
  rownames(RS) = unlist(Network_sub[,1:2])%>%unique()%>%sort()
  colnames(RS) = colnames(NetInd)

  raw.1=base::rowsum(Div.mat_t,Network_sub[,1])
  raw.2=base::rowsum(Div.mat_t,Network_sub[,2])

  RS[rownames(raw.1),]=raw.1
  RS[rownames(raw.2),]=RS[rownames(raw.2),]+raw.2
  #
  cRS=matrix(0,nrow = unlist(Network_sub[,1:2])%>%unique()%>%length(),ncol =ncol(Div.mat_t) )
  rownames(cRS) = unlist(Network_sub[,1:2])%>%unique()%>%sort()
  colnames(cRS) = colnames(NetInd)

  raw.1=base::rowsum(cDiv.mat_t_d1,Network_sub[,1])
  raw.2=base::rowsum(cDiv.mat_t_d2,Network_sub[,2])

  cRS[rownames(raw.1),]=raw.1
  cRS[rownames(raw.2),]=cRS[rownames(raw.2),]+raw.2
  plusRS = RS+cRS
  plusRS = plusRS[sigGene,,drop=F]
  dim(plusRS)
  ScaleplusRS = plusRS/MatrixGenerics::rowSums2(abs(plusRS)) # Probabilizing
  ScaleplusRS = ScaleplusRS*Zscores[sigGene,]$Zscores.Plus
  ScaleplusRS # A matrix where rows represent genes,
  # columns represent cells and
  # values represent the contribution of the cell to the functional perturbation of the gene.
  #
  # zs1 = NEAModel$ZscoreFit.Plus$cal_Zscore_Pvalue(ScaleplusRS[1,]*400*Zscores[sigGene[1],]$Zscores.Plus,
  #                                                     Dsize = Zscores[sigGene[1],]$Zscores.Plus,
  #                                                     ParmetersInput = NEAModel$ZscoreFit.Plus$ParmetersInput)



}


JSD_batch_matrix_rawMt <- function(Ctl,Pert,base=exp(1)){
  Ctl = Ctl+1e-200
  Pert = Pert+1e-200
  Ctl = Ctl/base::rowSums(Ctl)
  Pert = Pert/base::rowSums(Pert)
  M <- (Pert + Ctl) / 2

  jsd_mt <- (Ctl * log(Ctl / M,base = base)+Pert * log(Pert / M,base = base))/2
  jsd_mt
}



#' likelihoodFromMixEM_3dToResult_constrMean
#' mean.constr = c(NA,0,NA)
#' @param x
#' @param muStart
#'
#' @return
#' @export
#'
#' @examples
likelihoodFromMixEM_3dToResult_constrMean <- function(x,muStart=c(mean(x[x<0]),0,mean(x[x>0])),plot=TRUE){
  n.distrubtion=3 # Mixture of three normal distributions.
  if(is.null(muStart)){
    muStart=seq(min(x),max(x),(max(x)-min(x))/(n.distrubtion+1))
    muStart = muStart[2:(length(muStart)-1)]
  }

  fitModel <- mixtools::normalmixEM(x, lambda = 1/n.distrubtion, mu = muStart, sigma = sd(x)/3,mean.constr = c(NA,0,NA))
  if(plot){
    plot(fitModel, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
         main2="", xlab2="x")
  }

  Cal_likelihood<-function(x,fitModel){
    D1=pnorm(x,fitModel$mu[3],sd = fitModel$sigma[3]) # perturbed
    D2=pnorm(x,fitModel$mu[2],sd = fitModel$sigma[2],lower.tail = F) # commn
    D3=pnorm(x,fitModel$mu[1],sd = fitModel$sigma[2],lower.tail = F) # control
    LR_right = D1/(D1+D2+D3)

    D2_left=pnorm(x,fitModel$mu[2],sd = fitModel$sigma[2],lower.tail = T) # commn
    LR_left = D3/(D1+D2_left+D3)
    diff = abs(x-fitModel$mu[2])
    D2_1=pnorm(diff+fitModel$mu[2],fitModel$mu[2],sd = fitModel$sigma[2],lower.tail = F) # commn
    D2_2=pnorm(fitModel$mu[2]-diff,fitModel$mu[2],sd = fitModel$sigma[2],lower.tail = T) # commn
    LR_mid = (D2_1+D2_2)/(D1+D2_1+D2_2+D3)
    data.frame(LR_left = LR_left,LR_mid=LR_mid,LR_right=LR_right)
  }
  LRs = Cal_likelihood(x,fitModel)
  return(list(LRs = LRs,fitModel=fitModel,Cal_likelihood=Cal_likelihood))
}



mappingMinMaxRow_w_Max_Min <- function(ma,
                                       minD=1,
                                       maxD=60,
                                       Maxs=NULL,
                                       Mins=NULL){

  ma <- (ma - Mins)/(Maxs - Mins) * (maxD - minD) + minD
  ma[ma>maxD] <- maxD
  ma[ma<minD] <- minD
  ma
}



scContribution_v2 <- function (scDNSobject, nx = 20, rmZero = F, topGene = 100, sigGene = NULL,rb=F)
{
  # Zscores <- IRF4_B_scDNSob@Zscore
  # Zscores <- IRF4_B_scDNSob@Zscore
  Zscores <- scDNSobject@Zscore
  if (is.null(sigGene)) {
    sigGene = Zscores$Gene[top_n_vector(Zscores$Zscores.Plus,
                                        topGene, returnLogical = T)]
  }

  Network = scDNSobject@Network
  Network$id = 1:nrow(Network)
  Network_sub = getSubNetByNode(Network, sigGene)




  NetPdata<- getCellDensity(scDNSob =  scDNSobject,
                            ExpData = NULL,
                            Div_weight = scDNSobject@Div.Parameters$Div_weight,
                            EdgeID = Network_sub$id,
                            rb = rb)

  Div.mat_t <- NetPdata$Div.mat_t
  cDiv.mat_t_d1 <- NetPdata$cDiv.mat_t_d1
  cDiv.mat_t_d2 <- NetPdata$cDiv.mat_t_d2


  Div.mat_t = Div.mat_t/rowSums2(abs(Div.mat_t)) * Network_sub$Div.chiSquare.LR
  #

  cDiv.mat_t_d1 = cDiv.mat_t_d1/rowSums2(abs(cDiv.mat_t_d1)) * Network_sub$cDiv_D1.chiSquare.LR
  #
  cDiv.mat_t_d2 = cDiv.mat_t_d2/rowSums2(abs(cDiv.mat_t_d2)) * Network_sub$cDiv_D2.chiSquare.LR


  RS = matrix(0, nrow = unlist(Network_sub[, 1:2]) %>% unique() %>%
                length(), ncol = ncol(Div.mat_t))
  rownames(RS) = unlist(Network_sub[, 1:2]) %>% unique() %>%
    sort()
  colnames(RS) = colnames(NetInd)
  raw.1 = base::rowsum(Div.mat_t, Network_sub[, 1])
  raw.2 = base::rowsum(Div.mat_t, Network_sub[, 2])
  RS[rownames(raw.1), ] = raw.1
  RS[rownames(raw.2), ] = RS[rownames(raw.2), ] + raw.2
  cRS = matrix(0, nrow = unlist(Network_sub[, 1:2]) %>% unique() %>%
                 length(), ncol = ncol(Div.mat_t))
  rownames(cRS) = unlist(Network_sub[, 1:2]) %>% unique() %>%
    sort()
  colnames(cRS) = colnames(NetInd)
  raw.1 = base::rowsum(cDiv.mat_t_d1, Network_sub[, 1])
  raw.2 = base::rowsum(cDiv.mat_t_d2, Network_sub[, 2])
  cRS[rownames(raw.1), ] = raw.1
  cRS[rownames(raw.2), ] = cRS[rownames(raw.2), ] + raw.2
  cRS[rownames(raw.2), ] =  raw.2

  Zscores$Ratio <- sqrt(-log10(Zscores$Pvalues.Div)) /(sqrt(-log10(Zscores$Pvalues.Div)) +sqrt(-log10(Zscores$Pvalues.cDiv)) )
  # Zscores$Ratio <- Zscores$RawScore.Div/(Zscores$RawScore.Div+Zscores$RawScore.cDiv)
  # plusRS = RS*Zscores[rownames(RS),]$Ratio + cRS*(1-Zscores[rownames(cRS),]$Ratio)
  plusRS <- cRS
  plusRS = plusRS[sigGene, , drop = F]
  # dim(plusRS)
  ScaleplusRS = plusRS/Zscores[sigGene, ]$RawScore.Plus
  ScaleplusRS = ScaleplusRS * Zscores[sigGene, ]$Zscores.Plus
  ScaleplusRS
  # Pxdata <- data.frame(Type = IRF4_B_scDNSob@GroupLabel,CellType=IRF4_Bcells_Mut_SOB_sub$CellType,
  #                      scZcores =as.vector(ScaleplusRS),Gene=scDNSobject@data['XAF1',])
  # ggplot(Pxdata,aes(CellType,scZcores,fill=Type))+geom_violin(scale='width')+stat_summary(fun.y = mean,geom='point',
  # color='red',position = position_dodge(width=0.9))
}


getCellDensity <- function(scDNSob,
                           ExpData=NULL,
                           EdgeID,
                           Div_weight,
                           rb=FALSE){
  uniCase <- scDNSob@uniCase
  GroupLabel <- scDNSob@GroupLabel
  nx = scDNSob@Div.Parameters$n.coarse
  Network = scDNSob@Network
  Network$id = 1:nrow(scDNSob@Network)
  Network_sub = Network[EdgeID,]# getSubNetByNode(Network, sigGene)
  if(is.null(ExpData)){
    ExpData <- scDNSob@data
    # ExpData <- ExpData[unlist(Network_sub[, 1:2]) %>% unique(),  ]
    # ExpData1_20 <- mappingMinMaxRow(ExpData, minD = 1, maxD = scDNSob@Div.Parameters$n.coarse,
    #                                 rmZero = scDNSob@Div.Parameters$exclude.zero)
    # ExpData1_20 <- round(ExpData1_20, 0)
    Exp20 <- (ExpData - rowMins(ExpData))/(rowMaxs(ExpData) - rowMins(ExpData)) * (scDNSob@Div.Parameters$n.coarse - 1) + 1
    Exp20 <- round(Exp20,digits = 0)
  }else{
    Exp20 <- (ExpData - rowMins(ExpData))/(rowMaxs(ExpData) - rowMins(ExpData)) * (scDNSob@Div.Parameters$n.coarse - 1) + 1
    Exp20 <- round(Exp20,digits = 0)
  }


  SourceID = Exp20[Network_sub[, 1], ,drop=F]
  TargetID = Exp20[Network_sub[, 2], ,drop=F]
  filteredNetid = Network_sub$id
  # if(ExpendDen){
  #   ContextA_R <- getConditionalDenstiy(DS = ExpendDensityBycopyMaxNeiber(scDNSob@JDensity_A[filteredNetid,
  #                                                                                            , drop = F]))
  #   ContextB_R <- getConditionalDenstiy(DS = ExpendDensityBycopyMaxNeiber(scDNSob@JDensity_B[filteredNetid,
  #                                                                                            , drop = F]))
  # }else{
  #   ContextA_R <- getConditionalDenstiy(DS = scDNSob@JDensity_A[filteredNetid,
  #                                                                                            , drop = F])
  #   ContextB_R <- getConditionalDenstiy(DS =scDNSob@JDensity_B[filteredNetid,
  #                                                                                            , drop = F])
  # }
  ContextA_R <- getConditionalDenstiy(DS = scDNSob@JDensity_A[filteredNetid,
                                                              , drop = F])
  ContextB_R <- getConditionalDenstiy(DS =scDNSob@JDensity_B[filteredNetid,
                                                             , drop = F])
  # Div.mat = JSD_batch_matrix_rawMt(Ctl = ContextA_R$DS*sqrt(as.vector(matrix(1:20,nrow = 20,ncol = 20)*t(matrix(1:20,nrow = 20,ncol = 20)))),
  #                                  Pert = ContextB_R$DS*sqrt(as.vector(matrix(1:20,nrow = 20,ncol = 20)*t(matrix(1:20,nrow = 20,ncol = 20))))
  # Div_weight = sqrt(as.vector(matrix(1:20,nrow = 20,ncol = 20)*t(matrix(1:20,nrow = 20,ncol = 20))))
  if(rb){
    Div.mat = JSD_batch_matrix_rawMt_robust( ContextA_R$DS,
                                             ContextB_R$DS)
    Div.mat = t(t(Div.mat)*Div_weight)
  }else{
    Div.mat = JSD_batch_matrix_rawMt(Ctl = ContextA_R$DS,
                                     Pert = ContextB_R$DS)
    Div.mat = t(t(Div.mat)*Div_weight)
  }

  # Div.mat = JSD_batch_matrix_rawMt(Ctl = ContextA_R$DS,
  #                                  Pert = ContextB_R$DS)

  Div.mat <- abs(Div.mat)
  # Div.mat <- Div.mat
  if(rb){
    cDiv.mat_D1 = JSD_batch_matrix_rawMt_robust(rawDensity1 =  ContextA_R$cDS_1,
                                                rawDensity2 = ContextB_R$cDS_1)
    cDiv.mat_D1 = t(t(cDiv.mat_D1)*Div_weight)
  }else{
    cDiv.mat_D1 = JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS_1,
                                         Pert = ContextB_R$cDS_1)
    cDiv.mat_D1 =t(t(cDiv.mat_D1)*Div_weight)
  }

  cDiv.mat_D1 <- abs(cDiv.mat_D1)
  if(rb){
    cDiv.mat_D2 = JSD_batch_matrix_rawMt_robust(rawDensity1 =  ContextA_R$cDS_2,
                                                rawDensity2 = ContextB_R$cDS_2)
    cDiv.mat_D2 = t(t(cDiv.mat_D2)*Div_weight)
  }else{
    cDiv.mat_D2 = JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS_2,
                                         Pert = ContextB_R$cDS_2)
    cDiv.mat_D2 = t(t(cDiv.mat_D2)*Div_weight)
  }
  cDiv.mat_D2 <- abs(cDiv.mat_D2)
  #--------
  NetInd = sub2ind(SourceID, TargetID, nrow = nx, ncol = nx)
  NetInd_t = t(NetInd + (1:nrow(NetInd) - 1) * ncol(Div.mat))
  #-----------
  Div.mat <- Div.mat/rowSums2(Div.mat)
  Div.mat[ContextA_R$DS < ContextB_R$DS] = -Div.mat[ContextA_R$DS <
                                                      ContextB_R$DS]

  Div.mat_t = t(Div.mat)
  Div.mat_t = matrix(Div.mat_t[NetInd_t], nrow = ncol(NetInd))
  Div.mat_t = t(Div.mat_t)
  #
  # Div.mat_t[,GroupLabel==uniCase[1]] <- pmax(Div.mat_t[,GroupLabel==uniCase[1]],0)
  # Div.mat_t[,GroupLabel==uniCase[2]] <- pmin(Div.mat_t[,GroupLabel==uniCase[2]],0)

  # Div.mat_t = Div.mat_t/Network_sub$Div
  # Div.mat_t = Div.mat_t * Network_sub$Div.chiSquare.LR
  #
  cDiv.mat_D1 <- cDiv.mat_D1/rowSums2(cDiv.mat_D1)
  cDiv.mat_D1[ContextA_R$cDS_1 < ContextB_R$cDS_1] = -cDiv.mat_D1[ContextA_R$cDS_1 <
                                                                    ContextB_R$cDS_1]

  cDiv.mat_t_d1 = t(cDiv.mat_D1)
  cDiv.mat_t_d1 = matrix(cDiv.mat_t_d1[NetInd_t], nrow = ncol(NetInd))
  cDiv.mat_t_d1 = t(cDiv.mat_t_d1)
  # cDiv.mat_t_d1[,GroupLabel==uniCase[1]] <- pmax(cDiv.mat_t_d1[,GroupLabel==uniCase[1]],0)
  # cDiv.mat_t_d1[,GroupLabel==uniCase[2]] <- pmin(cDiv.mat_t_d1[,GroupLabel==uniCase[2]],0)

  # cDiv.mat_t_d1 = cDiv.mat_t_d1/Network_sub$cDiv_D1
  # cDiv.mat_t_d1 = cDiv.mat_t_d1 * Network_sub$cDiv_D1.chiSquare.LR
  #
  cDiv.mat_D2 <- cDiv.mat_D2/rowSums2(cDiv.mat_D2)
  cDiv.mat_D2[ContextA_R$cDS_2 < ContextB_R$cDS_2] = -cDiv.mat_D2[ContextA_R$cDS_2 <
                                                                    ContextB_R$cDS_2]

  cDiv.mat_t_d2 = t(cDiv.mat_D2)
  cDiv.mat_t_d2 = matrix(cDiv.mat_t_d2[NetInd_t], nrow = ncol(NetInd))
  cDiv.mat_t_d2 = t(cDiv.mat_t_d2)
  # cDiv.mat_t_d2[,GroupLabel==uniCase[1]] <- pmax(cDiv.mat_t_d2[,GroupLabel==uniCase[1]],0)
  # cDiv.mat_t_d2[,GroupLabel==uniCase[2]] <- pmin(cDiv.mat_t_d2[,GroupLabel==uniCase[2]],0)
  # cDiv.mat_t_d2 = cDiv.mat_t_d2/Network_sub$cDiv_D2
  # cDiv.mat_t_d2 = cDiv.mat_t_d2 * Network_sub$cDiv_D2.chiSquare.LR
  return(list(Div.mat_t=Div.mat_t,cDiv.mat_t_d1=cDiv.mat_t_d1,cDiv.mat_t_d2=cDiv.mat_t_d2,Network_sub=Network_sub))
}



scContribution_v3 <- function (scDNSobject, nx = 20, topGene = 100, sigGene = NULL,rb=F,model=c('cDiv','Div','Plus')[3])
{
  # Zscores <- IRF4_B_scDNSob@Zscore
  # Zscores <- IRF4_B_scDNSob@Zscore
  Zscores <- scDNSobject@Zscore
  if (is.null(sigGene)) {
    sigGene = Zscores$Gene[top_n_vector(Zscores$Zscores.ZsPlus,
                                        topGene, returnLogical = T,)]
  }

  Network = scDNSobject@Network
  Network$id = 1:nrow(Network)
  # mapping to Zscore
  Div.Sum <- replace2(c(Network[,1],Network[,2]),RawData = Zscores$Gene,RepData = Zscores$RawScore.Div)%>%as.numeric()%>%matrix(ncol = 2)
  Div.Sum_ratio <- 1/(Div.Sum/Network$Div.chiSquare.LR)
  Div.Zs <- replace2(c(Network[,1],Network[,2]),RawData = Zscores$Gene,RepData = log10(1/Zscores$Pvalues.Div))%>%as.numeric()%>%matrix(ncol = 2)
  Div.Zs_contribution <-  Div.Zs*Div.Sum_ratio
  colnames(Div.Zs_contribution) <- c('Div_source_ZsCon','Div_Target_ZsCon')

  cDiv.Sum <- replace2(c(Network[,1],Network[,2]),RawData = Zscores$Gene,RepData = Zscores$RawScore.cDiv)%>%as.numeric()%>%matrix(ncol = 2)
  cDiv.Sum_ratio <- Network[,c('cDiv_D1.chiSquare.LR' ,'cDiv_D2.chiSquare.LR')]/cDiv.Sum
  cDiv.Zs <- replace2(c(Network[,1],Network[,2]),RawData = Zscores$Gene,RepData = log10(1/Zscores$Pvalues.cDiv))%>%as.numeric()%>%matrix(ncol = 2)
  cDiv.Zs_contribution <-  cDiv.Zs*cDiv.Sum_ratio
  colnames(cDiv.Zs_contribution) <- c('cDiv_D1_ZsCon','cDiv_D2_ZsCon')
  Network <- cbind(Network,Div.Zs_contribution,cDiv.Zs_contribution)
  #
  Network_sub = getSubNetByNode(Network, sigGene)


  loopTimes <- ceiling(nrow(Network_sub)/1000)
  Div.mat_t <- NULL
  cDiv.mat_t_d1 <- NULL
  cDiv.mat_t_d2 <-NULL
  for(i in 1:loopTimes){
    NetPdata<- getCellDensity(scDNSob =  scDNSobject,
                              ExpData = NULL,
                              Div_weight = scDNSobject@Div.Parameters$Div_weight,
                              EdgeID = Network_sub$id[((i-1)*1000+1):min(c(i*1000,nrow(Network_sub)))],
                              rb = rb)
    print(dim(NetPdata$Div.mat_t))
    Div.mat_t <- rbind(Div.mat_t,NetPdata$Div.mat_t)
    cDiv.mat_t_d1 <- rbind(cDiv.mat_t_d1,NetPdata$cDiv.mat_t_d1)
    cDiv.mat_t_d2 <-  rbind(cDiv.mat_t_d2,NetPdata$cDiv.mat_t_d2)
    suppressMessages(gc())
  }



  # Div.mat_t = Div.mat_t/rowSums2(abs(Div.mat_t)) * Network_sub$Div.chiSquare.LR
  # Div.mat_t <- Div.mat_t* Network_sub$Div.chiSquare.LR
  Div.mat_t_source = Div.mat_t * Network_sub$Div_source_ZsCon
  Div.mat_t_target = Div.mat_t * Network_sub$Div_Target_ZsCon
  # Div.mat_t = Div.mat_t/rowSums2(abs(Div.mat_t)) * (log(Network_sub$Div.Pvalues)*-2)
  # #
  #
  # cDiv.mat_t_d1 = cDiv.mat_t_d1/rowSums2(abs(cDiv.mat_t_d1)) * Network_sub$cDiv_D1.chiSquare.LR
  # cDiv.mat_t_d1 = cDiv.mat_t_d1 * Network_sub$cDiv_D1.chiSquare.LR
  cDiv.mat_t_d1 = cDiv.mat_t_d1 * Network_sub$cDiv_D1_ZsCon
  # cDiv.mat_t_d1 = cDiv.mat_t_d1/rowSums2(abs(cDiv.mat_t_d1)) * (log(Network_sub$cDiv_D1.Pvalues)*-2)
  # #
  # cDiv.mat_t_d2 = cDiv.mat_t_d2/rowSums2(abs(cDiv.mat_t_d2)) * Network_sub$cDiv_D2.chiSquare.LR
  # cDiv.mat_t_d2 = cDiv.mat_t_d2 * Network_sub$cDiv_D2.chiSquare.LR
  cDiv.mat_t_d2 = cDiv.mat_t_d2 * Network_sub$cDiv_D2_ZsCon
  # cDiv.mat_t_d2 = cDiv.mat_t_d2/rowSums2(abs(cDiv.mat_t_d2)) * (log(Network_sub$cDiv_D2.Pvalues)*-2)


  RS = matrix(0, nrow = unlist(Network_sub[, 1:2]) %>% unique() %>%
                length(), ncol = ncol(Div.mat_t))
  rownames(RS) = unlist(Network_sub[, 1:2]) %>% unique() %>%
    sort()
  colnames(RS) = colnames(scDNSobject@data)
  # raw.1 = base::rowsum(Div.mat_t, Network_sub[, 1])
  # raw.2 = base::rowsum(Div.mat_t, Network_sub[, 2])
  raw.1 = base::rowsum(Div.mat_t_source, Network_sub[, 1])
  raw.2 = base::rowsum(Div.mat_t_target, Network_sub[, 2])
  RS[rownames(raw.1), ] = raw.1
  RS[rownames(raw.2), ] = RS[rownames(raw.2), ] + raw.2

  cRS = matrix(0, nrow = unlist(Network_sub[, 1:2]) %>% unique() %>%
                 length(), ncol = ncol(Div.mat_t))
  rownames(cRS) = unlist(Network_sub[, 1:2]) %>% unique() %>%
    sort()
  colnames(cRS) = colnames(scDNSobject@data)
  raw.1 = base::rowsum(cDiv.mat_t_d1, Network_sub[, 1])
  raw.2 = base::rowsum(cDiv.mat_t_d2, Network_sub[, 2])
  cRS[rownames(raw.1), ] = raw.1
  cRS[rownames(raw.2), ] = cRS[rownames(raw.2), ] + raw.2
  # cRS[rownames(raw.2), ] =  raw.2
  suppressMessages(gc())
  # Zscores$Ratio <- Zscores$Zscores.Div/(Zscores$Zscores.Div +Zscores$Zscores.cDiv)
  # Zscores$Ratio <- log10(1/Zscores$Pvalues.Div)/(log10(1/Zscores$Pvalues.Div)+log10(1/Zscores$Pvalues.cDiv))
  # Zscores$Ratio <- Zscores$RawScore.Div/(Zscores$RawScore.Div+Zscores$RawScore.cDiv)
  # plusRS = RS*Zscores[rownames(RS),]$Ratio + cRS*(1-Zscores[rownames(cRS),]$Ratio)
  if(model=='cDiv'){
    plusRS <- cRS
    plusRS = plusRS[sigGene, , drop = F]
    # dim(plusRS)
    ScaleplusRS = plusRS
    # ScaleplusRS = ScaleplusRS/Zscores[sigGene, ]$RawScore.cDiv
    # ScaleplusRS = ScaleplusRS * Zscores[sigGene, ]$Zscores.cDiv
    ScaleplusRS
  }else if(model=='Div'){
    plusRS <- RS
    plusRS = plusRS[sigGene, , drop = F]
    # dim(plusRS)
    ScaleplusRS = plusRS
    # ScaleplusRS = ScaleplusRS/Zscores[sigGene, ]$RawScore.Div
    # ScaleplusRS = ScaleplusRS * Zscores[sigGene, ]$Zscores.Div
    ScaleplusRS
  }else{
    # print(Zscores[sigGene,]$Ratio)
    RS <- RS[sigGene, , drop = F]
    # print(RS[1:10])
    cRS <- cRS[sigGene, , drop = F]
    # print(RS[1:10])
    plusRS <- RS+cRS
    # plusRS <- RS*Zscores[sigGene,]$Ratio+cRS*(1-Zscores[sigGene,]$Ratio)
    # plusRS
    # plusRS = plusRS[sigGene, , drop = F]
    # dim(plusRS)
    # ScaleplusRS = plusRS/Zscores[sigGene, ]$RawScore.Plus
    # ScaleplusRS = ScaleplusRS * Zscores[sigGene, ]$Zscores.Plus
    plusRS
  }

  # Pxdata <- data.frame(Type = IRF4_B_scDNSob@GroupLabel,CellType=IRF4_Bcells_Mut_SOB_sub$CellType,
  #                      scZcores =as.vector(ScaleplusRS),Gene=scDNSobject@data['XAF1',])
  # ggplot(Pxdata,aes(CellType,scZcores,fill=Type))+geom_violin(scale='width')+stat_summary(fun.y = mean,geom='point',
  # color='red',position = position_dodge(width=0.9))
}



#' scContribution_v4
#'
#' @param scDNSobject
#' @param nx
#' @param topGene
#' @param sigGene
#' @param rb
#' @param model
#'
#' @return
#' @export
#'
#' @examples
scContribution_v4 <- function (scDNSobject, nx = 20, topGene = 100, sigGene = NULL,rb=F,model=c('cDiv','Div','Plus')[3])
{
  # Zscores <- IRF4_B_scDNSob@Zscore
  # Zscores <- IRF4_B_scDNSob@Zscore
  Zscores <- scDNSobject@Zscore
  if (is.null(sigGene)) {
    sigGene = Zscores$Gene[top_n_vector(Zscores$Zscores.ZsPlus,
                                        topGene, returnLogical = T,)]
  }
  #(0)
  #
  Network = scDNSobject@Network
  Network$id = 1:nrow(Network)
  ncells <- ncol(scDNSobject@data)

  # (1)mapping to Zscore
  message('mapping to Zscore')
  Div.Sum <- replace2(c(Network[,1],Network[,2]),RawData = Zscores$Gene,RepData = Zscores$RawScore.Div)%>%as.numeric()%>%matrix(ncol = 2)
  Div.Sum_ratio <- 1/(Div.Sum/Network$Div.chiSquare.LR)
  Div.Zs <- replace2(c(Network[,1],Network[,2]),RawData = Zscores$Gene,RepData = log10(1/Zscores$Pvalues.Div))%>%as.numeric()%>%matrix(ncol = 2)
  Div.Zs_contribution <-  Div.Zs*Div.Sum_ratio
  colnames(Div.Zs_contribution) <- c('Div_source_ZsCon','Div_Target_ZsCon')

  cDiv.Sum <- replace2(c(Network[,1],Network[,2]),RawData = Zscores$Gene,RepData = Zscores$RawScore.cDiv)%>%as.numeric()%>%matrix(ncol = 2)
  cDiv.Sum_ratio <- Network[,c('cDiv_D1.chiSquare.LR' ,'cDiv_D2.chiSquare.LR')]/cDiv.Sum
  cDiv.Zs <- replace2(c(Network[,1],Network[,2]),RawData = Zscores$Gene,RepData = log10(1/Zscores$Pvalues.cDiv))%>%as.numeric()%>%matrix(ncol = 2)
  cDiv.Zs_contribution <-  cDiv.Zs*cDiv.Sum_ratio
  colnames(cDiv.Zs_contribution) <- c('cDiv_D1_ZsCon','cDiv_D2_ZsCon')
  Network <- cbind(Network,Div.Zs_contribution,cDiv.Zs_contribution)

  #
  Network_sub = getSubNetByNode(Network, sigGene)
  nNet <- nrow(Network_sub)
  #
  #(01.)
  message('discretized expression level')
  Exp20 <- scDNSobject@data[rownames(scDNSobject@data)%in%(unlist(Network_sub[,1:2])%>%unique()),]
  Exp20 <- (Exp20 - rowMins(Exp20))/(rowMaxs(Exp20) - rowMins(Exp20)) * (scDNSobject@Div.Parameters$n.coarse - 1) + 1
  Exp20 <- round(Exp20,digits = 0)

  suppressMessages(gc())
  # table(Network_sub[,1]%in%rownames(Exp20))
  #(3)
  message('scZcores')
  RS = matrix(0, nrow = unlist(Network_sub[, 1:2]) %>% unique() %>%
                length(), ncol = ncells)
  rownames(RS) = unlist(Network_sub[, 1:2]) %>% unique() %>%
    sort()
  colnames(RS) = colnames(scDNSobject@data)

  cRS <- RS

  loopTimes <- ceiling(nrow(Network_sub)/5000)

  message('number of edeges')
  print(nrow(Network_sub))
  for(i in 1:loopTimes){
    id <- ((i-1)*5000+1):min(c(i*5000,nrow(Network_sub)))
    Network_sub_sub <- Network_sub[id,]

    NetPdata <- getDiv(scDNSob =  scDNSobject,
                       Exp20  = Exp20,
                       Div_weight = scDNSobject@Div.Parameters$Div_weight,
                       Network_sub = Network_sub_sub,
                       rb = F)
    # return(list(Div.mat=Div.mat,cDiv.mat_D1=cDiv.mat_D1,cDiv.mat_D2=cDiv.mat_D2,Network_sub=Network_sub,NetInd_t=NetInd_t))
    NetInd_t <- NetPdata$NetInd_t # column cells
    Div.mat <- NetPdata$Div.mat # column cells,
    cDiv.mat_D1 <- NetPdata$cDiv.mat_D1
    cDiv.mat_D2 <- NetPdata$cDiv.mat_D2
    rm(NetPdata)
    suppressMessages(gc())
    #
    Div.mat = matrix(Div.mat[NetInd_t], nrow = ncells)
    Div.mat = t(Div.mat)
    #

    cDiv.mat_D1 = matrix(cDiv.mat_D1[NetInd_t], nrow =ncells)
    cDiv.mat_D1 = t(cDiv.mat_D1)

    #

    cDiv.mat_D2 = matrix(cDiv.mat_D2[NetInd_t], nrow = ncells)
    cDiv.mat_D2 = t(cDiv.mat_D2)

    #
    Div.mat_t_source = Div.mat * Network_sub_sub$Div_source_ZsCon
    Div.mat_t_target = Div.mat * Network_sub_sub$Div_Target_ZsCon

    cDiv.mat_D1 = cDiv.mat_D1 * Network_sub_sub$cDiv_D1_ZsCon
    cDiv.mat_D2 = cDiv.mat_D2 * Network_sub_sub$cDiv_D2_ZsCon

    #
    raw.1 = base::rowsum(Div.mat_t_source, Network_sub_sub[, 1])
    raw.2 = base::rowsum(Div.mat_t_target, Network_sub_sub[, 2])
    RS[rownames(raw.1), ] = raw.1
    RS[rownames(raw.2), ] = RS[rownames(raw.2), ] + raw.2
    #
    raw.1 = base::rowsum(cDiv.mat_D1, Network_sub_sub[, 1])
    raw.2 = base::rowsum(cDiv.mat_D2, Network_sub_sub[, 2])
    cRS[rownames(raw.1), ] = raw.1
    cRS[rownames(raw.2), ] = cRS[rownames(raw.2), ] + raw.2
    print(c(i,loopTimes))

  }

  if(model=='cDiv'){
    plusRS <- cRS
    ScaleplusRS = plusRS[sigGene, , drop = F]
    ScaleplusRS
  }else if(model=='Div'){
    plusRS <- RS
    ScaleplusRS = plusRS[sigGene, , drop = F]
    ScaleplusRS
  }else{
    RS <- RS[sigGene, , drop = F]
    cRS <- cRS[sigGene, , drop = F]
    plusRS <- RS+cRS
    suppressMessages(gc())
    plusRS
  }

}

scContribution_v5 <- function (scDNSobject, nx = 20, topGene = 100, sigGene = NULL,
                               rb = F, model = c("cDiv", "Div", "Plus")[3])
{
  Zscores <- scDNSobject@Zscore
  if (is.null(sigGene)) {
    sigGene = Zscores$Gene[top_n_vector(Zscores$Zscores.ZsPlus,
                                        topGene, returnLogical = T, )]
  }
  Network = scDNSobject@Network
  Network$id = 1:nrow(Network)
  ncells <- ncol(scDNSobject@data)
  message("mapping to Zscore")
  Div.Sum <- replace2(c(Network[, 1], Network[, 2]), RawData = Zscores$Gene,
                      RepData = Zscores$RawScore.Div) %>% as.numeric() %>%
    matrix(ncol = 2)
  Div.Sum_ratio <- 1/(Div.Sum/Network$Div.chiSquare.LR)
  Div.Zs <- replace2(c(Network[, 1], Network[, 2]), RawData = Zscores$Gene,
                     RepData = log10(1/pmax(Zscores$Pvalues.Div, 1e-300))) %>% as.numeric() %>% # 防止inf
    matrix(ncol = 2)
  Div.Zs_contribution <- Div.Zs * Div.Sum_ratio
  colnames(Div.Zs_contribution) <- c("Div_source_ZsCon", "Div_Target_ZsCon")
  cDiv.Sum <- replace2(c(Network[, 1], Network[, 2]), RawData = Zscores$Gene,
                       RepData = Zscores$RawScore.cDiv) %>% as.numeric() %>%
    matrix(ncol = 2)
  cDiv.Sum_ratio <- Network[, c("cDiv_D1.chiSquare.LR", "cDiv_D2.chiSquare.LR")]/cDiv.Sum
  cDiv.Zs <- replace2(c(Network[, 1], Network[, 2]), RawData = Zscores$Gene,
                      RepData = log10(1/pmax(Zscores$Pvalues.cDiv, 1e-300))) %>% as.numeric() %>% # 防止inf
    matrix(ncol = 2)
  cDiv.Zs_contribution <- cDiv.Zs * cDiv.Sum_ratio
  colnames(cDiv.Zs_contribution) <- c("cDiv_D1_ZsCon", "cDiv_D2_ZsCon")
  Network <- cbind(Network, Div.Zs_contribution, cDiv.Zs_contribution)
  Network_sub = getSubNetByNode(Network, sigGene)
  nNet <- nrow(Network_sub)
  message("discretized expression level")
  Exp20 <- scDNSobject@data[rownames(scDNSobject@data) %in%
                              (unlist(Network_sub[, 1:2]) %>% unique()), ]
  # Exp20 <- (Exp20 - rowMins(Exp20))/(rowMaxs(Exp20) - rowMins(Exp20)) *
  #   (scDNSobject@Div.Parameters$n.grid + 1) + 0
  # Exp20 <- sqrt(Exp20)
  Exp20 <- (Exp20 - rowMins(Exp20))/(rowMaxs(Exp20)-rowMins(Exp20)) *
    (scDNSobject@Div.Parameters$n.coarse - 1) + 1
  # Exp20 <- (Exp20 - rowMins(Exp20))/(rowMaxs(Exp20) - rowMins(Exp20)) *
  #   (scDNSobject@Div.Parameters$n.coarse - 1) + 1
  Exp20 <- round(Exp20, digits = 0)
  suppressMessages(gc())
  message("scZcores")
  RS = matrix(0, nrow = unlist(Network_sub[, 1:2]) %>% unique() %>%
                length(), ncol = ncells)
  rownames(RS) = unlist(Network_sub[, 1:2]) %>% unique() %>%
    sort()
  colnames(RS) = colnames(scDNSobject@data)
  cRS <- RS
  loopTimes <- ceiling(nrow(Network_sub)/5000)
  message("number of edeges")
  print(nrow(Network_sub))
  for (i in 1:loopTimes) {
    id <- ((i - 1) * 5000 + 1):min(c(i * 5000, nrow(Network_sub)))
    Network_sub_sub <- Network_sub[id, ]
    NetPdata <- getDiv(scDNSob = scDNSobject, Exp20 = Exp20,
                       Div_weight = scDNSobject@Div.Parameters$Div_weight,
                       Network_sub = Network_sub_sub, rb = F)
    NetInd_t <- NetPdata$NetInd_t
    Div.mat <- NetPdata$Div.mat
    cDiv.mat_D1 <- NetPdata$cDiv.mat_D1
    cDiv.mat_D2 <- NetPdata$cDiv.mat_D2
    rm(NetPdata)
    suppressMessages(gc())
    Div.mat = matrix(Div.mat[NetInd_t], nrow = ncells)
    Div.mat = t(Div.mat)
    cDiv.mat_D1 = matrix(cDiv.mat_D1[NetInd_t], nrow = ncells)
    cDiv.mat_D1 = t(cDiv.mat_D1)
    cDiv.mat_D2 = matrix(cDiv.mat_D2[NetInd_t], nrow = ncells)
    cDiv.mat_D2 = t(cDiv.mat_D2)
    Div.mat_t_source = Div.mat * Network_sub_sub$Div_source_ZsCon
    Div.mat_t_target = Div.mat * Network_sub_sub$Div_Target_ZsCon
    cDiv.mat_D1 = cDiv.mat_D1 * Network_sub_sub$cDiv_D1_ZsCon
    cDiv.mat_D2 = cDiv.mat_D2 * Network_sub_sub$cDiv_D2_ZsCon
    raw.1 = base::rowsum(Div.mat_t_source, Network_sub_sub[,
                                                           1])
    raw.2 = base::rowsum(Div.mat_t_target, Network_sub_sub[,
                                                           2])
    RS[rownames(raw.1), ] = raw.1
    RS[rownames(raw.2), ] = RS[rownames(raw.2), ] + raw.2
    raw.1 = base::rowsum(cDiv.mat_D1, Network_sub_sub[,
                                                      1])
    raw.2 = base::rowsum(cDiv.mat_D2, Network_sub_sub[,
                                                      2])
    cRS[rownames(raw.1), ] = raw.1
    cRS[rownames(raw.2), ] = cRS[rownames(raw.2), ] + raw.2
    print(c(i, loopTimes))
  }
  if (model == "cDiv") {
    plusRS <- cRS
    ScaleplusRS = plusRS[sigGene, , drop = F]
    ScaleplusRS
  }
  else if (model == "Div") {
    plusRS <- RS
    ScaleplusRS = plusRS[sigGene, , drop = F]
    ScaleplusRS
  }
  else {
    RS <- RS[sigGene, , drop = F]
    cRS <- cRS[sigGene, , drop = F]
    plusRS <- RS + cRS
    plusRS
  }
}

getDiv <- function(scDNSob,
                   Exp20=NULL,
                   Network_sub,
                   Div_weight,
                   rb=FALSE){
  uniCase <- scDNSob@uniCase
  GroupLabel <- scDNSob@GroupLabel
  nx = scDNSob@Div.Parameters$n.coarse
  # Network = scDNSob@Network
  # Network$id = 1:nrow(scDNSob@Network)
  # Network_sub = Network[EdgeID,]# getSubNetByNode(Network, sigGene)
  if(is.null(Exp20)){
    Exp20 <- scDNSob@data[rownames(scDNSob@data)%in%unlist(Network_sub[,1:2])%>%unique(),]
    # ExpData <- ExpData[unlist(Network_sub[, 1:2]) %>% unique(),  ]
    # ExpData1_20 <- mappingMinMaxRow(ExpData, minD = 1, maxD = scDNSob@Div.Parameters$n.coarse,
    #                                 rmZero = scDNSob@Div.Parameters$exclude.zero)
    # ExpData1_20 <- round(ExpData1_20, 0)
    Exp20 <- (Exp20 - rowMins(Exp20))/(rowMaxs(Exp20) - rowMins(Exp20)) * (scDNSob@Div.Parameters$n.coarse - 1) + 1
    Exp20 <- round(Exp20,digits = 0)
  }


  SourceID = Exp20[Network_sub[, 1], ,drop=F]
  TargetID = Exp20[Network_sub[, 2], ,drop=F]
  filteredNetid = Network_sub$id
  # if(ExpendDen){
  #   ContextA_R <- getConditionalDenstiy(DS = ExpendDensityBycopyMaxNeiber(scDNSob@JDensity_A[filteredNetid,
  #                                                                                            , drop = F]))
  #   ContextB_R <- getConditionalDenstiy(DS = ExpendDensityBycopyMaxNeiber(scDNSob@JDensity_B[filteredNetid,
  #                                                                                            , drop = F]))
  # }else{
  #   ContextA_R <- getConditionalDenstiy(DS = scDNSob@JDensity_A[filteredNetid,
  #                                                                                            , drop = F])
  #   ContextB_R <- getConditionalDenstiy(DS =scDNSob@JDensity_B[filteredNetid,
  #                                                                                            , drop = F])
  # }
  ContextA_R <- getConditionalDenstiy(DS = scDNSob@JDensity_A[filteredNetid,
                                                              , drop = F])
  ContextB_R <- getConditionalDenstiy(DS =scDNSob@JDensity_B[filteredNetid,
                                                             , drop = F])
  # Div.mat = JSD_batch_matrix_rawMt(Ctl = ContextA_R$DS*sqrt(as.vector(matrix(1:20,nrow = 20,ncol = 20)*t(matrix(1:20,nrow = 20,ncol = 20)))),
  #                                  Pert = ContextB_R$DS*sqrt(as.vector(matrix(1:20,nrow = 20,ncol = 20)*t(matrix(1:20,nrow = 20,ncol = 20))))
  # Div_weight = sqrt(as.vector(matrix(1:20,nrow = 20,ncol = 20)*t(matrix(1:20,nrow = 20,ncol = 20))))
  if(rb){
    Div.mat = JSD_batch_matrix_rawMt_robust( ContextA_R$DS,
                                             ContextB_R$DS)
    Div.mat = t(t(Div.mat)*Div_weight)
  }else{
    Div.mat = JSD_batch_matrix_rawMt(Ctl = ContextA_R$DS,
                                     Pert = ContextB_R$DS)
    Div.mat = t(t(Div.mat)*Div_weight)
  }

  # Div.mat = JSD_batch_matrix_rawMt(Ctl = ContextA_R$DS,
  #                                  Pert = ContextB_R$DS)

  Div.mat <- abs(Div.mat)
  # Div.mat <- Div.mat
  if(rb){
    cDiv.mat_D1 = JSD_batch_matrix_rawMt_robust(rawDensity1 =  ContextA_R$cDS_1,
                                                rawDensity2 = ContextB_R$cDS_1)
    cDiv.mat_D1 = t(t(cDiv.mat_D1)*Div_weight)
  }else{
    cDiv.mat_D1 = JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS_1,
                                         Pert = ContextB_R$cDS_1)
    cDiv.mat_D1 =t(t(cDiv.mat_D1)*Div_weight)
  }

  cDiv.mat_D1 <- abs(cDiv.mat_D1)
  if(rb){
    cDiv.mat_D2 = JSD_batch_matrix_rawMt_robust(rawDensity1 =  ContextA_R$cDS_2,
                                                rawDensity2 = ContextB_R$cDS_2)
    cDiv.mat_D2 = t(t(cDiv.mat_D2)*Div_weight)
  }else{
    cDiv.mat_D2 = JSD_batch_matrix_rawMt(Ctl = ContextA_R$cDS_2,
                                         Pert = ContextB_R$cDS_2)
    cDiv.mat_D2 = t(t(cDiv.mat_D2)*Div_weight)
  }
  cDiv.mat_D2 <- abs(cDiv.mat_D2)
  #--------
  NetInd = sub2ind(SourceID, TargetID, nrow = nx, ncol = nx)
  NetInd_t = t(NetInd + (1:nrow(NetInd) - 1) * ncol(Div.mat)) # row: cells ; column: network
  #-----------
  Div.mat <- Div.mat/rowSums2(Div.mat)
  Div.mat[ContextA_R$DS < ContextB_R$DS] = -Div.mat[ContextA_R$DS <
                                                      ContextB_R$DS]

  Div.mat = t(Div.mat)
  # Div.mat_t = matrix(Div.mat_t[NetInd_t], nrow = ncol(NetInd))
  # Div.mat_t = t(Div.mat_t)
  #
  #
  cDiv.mat_D1 <- cDiv.mat_D1/rowSums2(cDiv.mat_D1)
  cDiv.mat_D1[ContextA_R$cDS_1 < ContextB_R$cDS_1] = -cDiv.mat_D1[ContextA_R$cDS_1 <
                                                                    ContextB_R$cDS_1]

  cDiv.mat_D1 = t(cDiv.mat_D1)
  # cDiv.mat_t_d1 = matrix(cDiv.mat_t_d1[NetInd_t], nrow = ncol(NetInd))
  # cDiv.mat_t_d1 = t(cDiv.mat_t_d1)

  #
  cDiv.mat_D2 <- cDiv.mat_D2/rowSums2(cDiv.mat_D2)
  cDiv.mat_D2[ContextA_R$cDS_2 < ContextB_R$cDS_2] = -cDiv.mat_D2[ContextA_R$cDS_2 <
                                                                    ContextB_R$cDS_2]

  cDiv.mat_D2 = t(cDiv.mat_D2)
  # cDiv.mat_t_d2 = matrix(cDiv.mat_t_d2[NetInd_t], nrow = ncol(NetInd))
  # cDiv.mat_t_d2 = t(cDiv.mat_t_d2)

  return(list(Div.mat=Div.mat,cDiv.mat_D1=cDiv.mat_D1,cDiv.mat_D2=cDiv.mat_D2,Network_sub=Network_sub,NetInd_t=NetInd_t))
}
