# scDNS define class
scDNSobjClass <- function(){
  setClassUnion("matrixORdgCMatrix", c("matrix", "dgCMatrix"))
  setClassUnion("matrixORarray", c("matrix", "array"))
  setClass("scDNS", representation(
    # attribute
    counts =  "matrixORdgCMatrix",
    data =  "matrixORdgCMatrix",
    Network = 'data.frame',
    GroupLabel = 'character',
    GeneVariability= 'numeric',
    Div.Parameters = 'list',
    NEA.Parameters = 'list',
    NEAModel = 'list',
    JDensity_A = 'matrixORarray',
    JDensity_B = 'matrixORarray',
    uniCase = 'character',
    Zscore = 'data.frame',
    scZscore = 'matrixORdgCMatrix',
    Other = 'list'
  ))
}

scDNSobjClass()
# # 创建类的方法
setMethod("show", "scDNS",
          function(object) {
            cat("Gene expression matrix size (gene × cell):", dim(object@data), "\n")
            cat("The number of edges in the network:", nrow(object@Network), "\n")
            cat("Group information:", object@uniCase, "\n")
          }
)

seurat2scDNSObj <- function(sob,imputedAssay ='MAGIC_RNA',GroupBy=NULL,...){
  if(is.null(GroupBy)){
    stop('Please provide cell groups information.')
  }
  scDNSob <- CreatScDNSobject(counts = sob@assays$RNA@counts,
                                  data = sob@assays[[imputedAssay]]@data,
                                  Network = scDNSBioNet,
                                  GroupLabel = sob@meta.data[,GroupBy]%>%as.character(),...)
  scDNSob
}


setMethod("colnames", signature(x = "scDNS"), function(x) {
  if (!is.null(x@data)) {
    return(colnames(x@data))
  } else {
    warning("Data attribute not available.")
    return(NULL)
  }
})
setMethod("rownames", signature(x = "scDNS"), function(x) {

  if (!is.null(x@data)) {
    return(rownames(x@data))
  } else {
    warning("Data attribute not available.")
    return(NULL)
  }
})






