#' The NeuronChat Class
#'# Class definitions
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom data.table data.table
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'AnyFactor', members = c("factor", "list"))
setClassUnion(name = 'AnyDF', members = c("data.frame", "data.table"))
#' The key slots used in the NeuronChat object are described below; adapted from CellChat https://github.com/sqjin/CellChat
#'
#' @slot data.raw raw count data matrix
#' @slot data normalized data matrix for CellChat analysis (Genes should be in rows and cells in columns)
#' @slot data.signaling a data.frame only containing signaling genes as well as cell group labels
#' @slot net0 a list contianing the original communication strength matrix without filtered by p-value
#' @slot pvalue a list containing the p-value (calculated by permutation test) matrices corresponding to net0
#' @slot net a list of length K contianing the communication strength matrices (N1 X N2) filtered by p-value, where N1 is the number of sending cell groups, N2 is the number of receiving cell groups, and K is the number of ligand-target pairs. Each row of the communication strength matrix indicates the communication probability originating from the sender cell group to other cell groups
#' @slot net_analysis a list containing various network analysis results
#' @slot DB ligand-target interaction database used in the analysis
#' @slot LR a list of information related with ligand-target pairs
#' @slot meta data frame storing the information associated with each cell
#' @slot idents a factor defining the cell identity used for all analysis. It becomes a list for a merged NeuronChat object
#' @slot dr List of the reduced 2D coordinates, one per method, e.g., umap/tsne/dm
#' @slot options List of miscellaneous data, such as parameters used throughout analysis, and a indicator whether the NeuronChat object is a single or merged
#' @slot fc deprecated
#' @slot info information flow of each ligand-target pair
#' @slot ligand.abundance ligand.abundance of sending cell groups (row) for ligand-target pairs (column)
#' @slot target.abundance target.abundance of receiving cell groups (row) for ligand-target pairs (column)
#' @exportClass NeuronChat
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
NeuronChat <- methods::setClass("NeuronChat",
                              slots = c(data.raw = 'AnyMatrix',
                                        data = 'AnyMatrix',
                                        data.signaling = "AnyDF",
                                        net0 = "list",
                                        pvalue = "list",
                                        net = "list",
                                        net_analysis = "list",
                                        fc = 'numeric',
                                        info ='numeric',
                                        ligand.abundance = 'matrix',
                                        target.abundance = 'matrix',
                                        meta = "data.frame",
                                        idents = "AnyFactor",
                                        DB = "list",
                                        LR = "character",
                                        dr = "list",
                                        options = "list")
)

#' Create a new NeuronChat object from a data matrix
#' Create a new NeuronChat object from a data matrix; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object a normalized (NOT count) data matrix (genes by cells)
#' @param meta a data frame (rows are cells with rownames) consisting of cell information, which will be used for defining cell groups.
#' @param group.by a vector to indicate group annotations of cells
#' @param assay Assay to use when the input is a Seurat object.
#' @param do.sparse whether use sparse format
# #' @param data is deprecated. Use `object`
#'
#' @return
#' @export
#' @importFrom methods as new
#' @importFrom Matrix t
#' @examples
createNeuronChat <- function(object, DB=c('mouse','human'),meta = NULL, group.by = NULL, assay = NULL, do.sparse = T) {
  DB <- match.arg(DB)
  # data matrix as input
  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix"))) {
    message("Create a NeuronChat object from a data matrix")
    data <- object
    if (is.null(group.by)) {
      group.by <- "labels"
    }
  }
  # if (do.sparse) {
  #   data <- as(data, "dgCMatrix")
  # }
  # signaling data
  DB_filename <- paste('interactionDB_',DB,sep='')
  # data(list=DB_filename)
  interactionDB <- eval(parse(text = DB_filename))
  # interactionDB <- readRDS(DB_filename)

  interactionDB <- interactionDB[!duplicated(names(interactionDB))]
  interactionDB <- interactionDB[sort(names(interactionDB),method='radix')]
  signaling_gene <- c()
  for(j in 1:length(interactionDB)){
    signaling_gene_j <- c(interactionDB[[j]]$lig_contributor,interactionDB[[j]]$receptor_subunit)
    signaling_gene <- c(signaling_gene,signaling_gene_j)
  }
  signaling_gene <- unique(signaling_gene)
  signaling_gene <- signaling_gene[signaling_gene %in% rownames(data)]
  data.signaling <- as.data.frame(Matrix::t(data[signaling_gene,]))
  data.signaling <- data.signaling/max(data.signaling) # normalized to range [0,1]
  data.signaling$cell_subclass <- group.by

  if (!is.null(meta)) {
    if (inherits(x = meta, what = c("matrix", "Matrix"))) {
      meta <- as.data.frame(x = meta)
    }
    if (!is.data.frame(meta)) {
      stop("The input `meta` should be a data frame")
    }
    if (!identical(rownames(meta), colnames(data))) {
      cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
      warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
      rownames(meta) <- colnames(data)
    }
  } else {
    meta <- data.frame()
  }

  object <- methods::new(Class = "NeuronChat",
                         data = data,
                         meta = meta,
                         data.signaling = data.signaling,
                         DB = interactionDB,
                         LR = names(interactionDB),
                         idents = factor(group.by))
  object@options$mode <- "single"
  return(object)
}


#' Merge NeuronChat objects
#' Merge NeuronChat objects; adapted from CellChat https://github.com/sqjin/CellChat
#'
#' @param object.list  A list of multiple NeuronChat objects
#' @param add.names A vector containing the name of each dataset
#' @param merge.data whether merging the data for ALL genes. Default only merges the data of signaling genes
#' @param cell.prefix whether prefix cell names
#' @importFrom methods slot new
#'
#' @return
#' @export
#'
#' @examples
mergeNeuronChat<- function(object.list, add.names = NULL, merge.data = FALSE, cell.prefix = FALSE) {
  if (is.null(add.names)) {
    add.names <- paste("Dataset",1:length(object.list),sep = "_")
  }
  slot.name <- c("net", "idents")
  slot.combined <- vector("list", length(slot.name))
  names(slot.combined) <- slot.name
  for (i in 1:length(slot.name)) {
    object.slot <- vector("list", length(object.list))
    for (j in 1:length(object.list)) {
      object.slot[[j]] <- slot(object.list[[j]], slot.name[i])
    }
    slot.combined[[i]] <- object.slot
    names(slot.combined[[i]]) <- add.names
  }

  if (cell.prefix) {
    warning("Prefix cell names!")
    for (i in 1:length(object.list)) {colnames(object.list[[i]]@data) <- paste(colnames(object.list[[i]]@data), add.names[i], sep = "_")}
  } else {
    cell.names <- c()
    for (i in 1:length(object.list)) {
      cell.names <- c(cell.names, colnames(object.list[[i]]@data))
    }
    if (sum(duplicated(cell.names))) {
      stop("Duplicated cell names were detected across datasets!! Please set cell.prefix = TRUE")
    }
  }

  meta.use <- colnames(object.list[[1]]@meta)
  for (i in 2:length(object.list)) {
    meta.use <- meta.use[meta.use %in% colnames(object.list[[i]]@meta)]
  }

  dataset.name <- c()
  cell.names <- c()
  meta.joint <- data.frame()
  for (i in 1:length(object.list)) {
    dataset.name <- c(dataset.name, rep(add.names[i], length(colnames(object.list[[i]]@data))))
    cell.names <- c(cell.names, colnames(object.list[[i]]@data))
    meta.joint <- rbind(meta.joint, object.list[[i]]@meta[ , meta.use, drop = FALSE])
  }
  if (!identical(rownames(meta.joint), cell.names)) {
    cat("The cell barcodes in merged 'meta' is ", head(rownames(meta.joint)),'\n')
    warning("The cell barcodes in merged 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of merged 'mata'!")
    rownames(meta.joint) <- cell.names
  }

  #dataset.name <- data.frame(dataset.name = dataset.name, row.names = cell.names)
  meta.joint$datasets <- factor(dataset.name, levels = add.names)

  genes.use <- rownames(object.list[[1]]@data)
  for (i in 2:length(object.list)) {
    genes.use <- genes.use[genes.use %in% rownames(object.list[[i]]@data)]
  }
  data.joint <- c()
  for (i in 1:length(object.list)) {
    data.joint <- cbind(data.joint, object.list[[i]]@data[genes.use, ])
  }
  gene.signaling.joint = unique(unlist(lapply(object.list, function(x) rownames(x@data.signaling))))
  data.signaling.joint <-as.data.frame(t(data.joint[rownames(data.joint) %in% gene.signaling.joint, ]))

  idents.joint <- c()
  idents.levels <- c()
  for (i in 1:length(object.list)) {
    idents.joint <- c(idents.joint, as.character(object.list[[i]]@idents))
    idents.levels <- union(idents.levels, levels(object.list[[i]]@idents))
  }
  names(idents.joint) <- cell.names
  idents.joint <- factor(idents.joint, levels = idents.levels)
  slot.combined$idents$joint <- idents.joint

  if (merge.data) {
    message("Merge the following slots: 'data','data.signaling','net', 'netP','meta', 'idents', 'var.features', 'DB', and 'LR'.")
    merged.object <- methods::new(
      Class = "NeuronChat",
      data = data.joint,
      data.signaling = data.signaling.joint,
      net = slot.combined$net,
      meta = meta.joint,
      idents = slot.combined$idents,
      LR = object.list[[1]]@LR,
      DB = object.list[[1]]@DB)
  } else {
    message("Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.")
    merged.object <- methods::new(
      Class = "NeuronChat",
      data.signaling = data.signaling.joint,
      net = slot.combined$net,
      meta = meta.joint,
      idents = slot.combined$idents,
      LR = object.list[[1]]@LR,
      DB = object.list[[1]]@DB)
  }
  merged.object@options$mode <- "merged"
  return(merged.object)
}

#' Calculate the average gene expression by cell group
#'
#' @param df a data frame or data.table containing signaling gene expression (row: cell; column: genes; the last column is cell group labels)
#' @param gene_used gene symbols used to calculate average expression for cell groups
#' @import data.table
#'
#' @return
#' @export
#'
#' @examples
cal_expr_by_group <- function(df,gene_used,mean_method=NULL) {
  df_used <- df[,c(gene_used,'cell_subclass')]
  if(is.null(mean_method)){
  q = c(.25, .5, .75)
  # ind_used <- which(gene_used %in% colnames(target_df))
  # ind_used <- c(ind_used,dim(target_df)[2])
  # target_used <- target_df[,ind_used]
  q1 <- setDT(df_used)[, lapply(.SD, quantile,q[1]), keyby = cell_subclass]
  q2 <- setDT(df_used)[, lapply(.SD, quantile,q[2]), keyby = cell_subclass]
  q3 <- setDT(df_used)[, lapply(.SD, quantile,q[3]), keyby = cell_subclass]
  yy <- 0.25*q1[,gene_used,with=FALSE]+0.5*q2[,gene_used,with=FALSE]+0.25*q3[,gene_used,with=FALSE]
  rownames(yy) <- q1$cell_subclass
  } else if(mean_method=='mean') {
  yy <- setDT(df_used)[, lapply(.SD, mean), keyby = cell_subclass]
  xx <- yy[,gene_used,with=FALSE]
  rownames(xx) <- yy$cell_subclass
  yy <- xx
  }
  return(yy)
}
# scPalette <- function(n) {
#   colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
#                   '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
#   if (n <= length(colorSpace)) {
#     colors <- colorSpace[1:n]
#   } else {
#     colors <- grDevices::colorRampPalette(colorSpace)(n)
#   }
#   return(colors)
# }


#' Calculate the communication strength matrix for a single ligand-target pair (without permutation test)
#' @param df a data frame or data.table containing signaling gene expression (row: cell; column: genes; the last column is cell group labels)
#' @param sender sending cell groups
#' @param receiver receiving cell groups
#' @param lig_contributor_new a vector with updated gene symbols (removing those excluded from the expression data) that are regarded as ligand contributors
#' @param lig_contributor_group a vector indicating the grouping of ligand contributors
#' @param lig_contributor_coeff a vector with stoichiometry of each group of ligand contributors for calculating ligand abundance
#' @param receptor_subunit_new a vector with updated gene symbols (removing those excluded from the expression data) that are regarded as target subunits
#' @param receptor_subunit_group a vector indicating the grouping of target subunits
#' @param receptor_subunit_coeff a vector with stoichiometry of each group of target subunits for calculating target abundance
#' @param targets_up deprecated
#' @param targets_down deprecated
#' @param targets_nichenet deprecated
#' @param N deprecated
#' @param K Hill coefficient for CellChat modeling method
#'
#' @return
#' @export
#'
#' @examples
cal_prob_mtx_downstream <- function(df,sender,receiver,lig_contributor_new,lig_contributor_group,lig_contributor_coeff,receptor_subunit_new,receptor_subunit_group,receptor_subunit_coeff,targets_up,targets_down,targets_nichenet,N,K=0.5,method=NULL,mean_method=NULL){
  gene_used <- c(lig_contributor_new,receptor_subunit_new)
  expr_gene <- cal_expr_by_group(df,gene_used,mean_method)
  cell.rownames <- rownames(expr_gene)
  expr_gene <- as.data.frame(expr_gene)
  #lig_contributor_expr <- apply(expr_gene[1:length(lig_contributor_new)],MARGIN = 1, function(x) exp(mean(log(x))))
  cellgroup_number <- dim(expr_gene)[1]
  targets_up_avg <- rep(0,cellgroup_number)
  targets_down_avg <- rep(0,cellgroup_number)
  targets_nichenet_avg <- rep(0,cellgroup_number)
  if(is.null(method)){
    # ligand
    lig_contributor_expr_tmp <- matrix(0,cellgroup_number,length(lig_contributor_coeff))
    for(i in 1:length(lig_contributor_coeff)){
      ind_i <- which(lig_contributor_group==i)
      if(length(ind_i)==1){
        lig_contributor_expr_tmp[1:cellgroup_number,i] <- expr_gene[1:cellgroup_number,ind_i]
      } else if(length(ind_i)>1) {
        lig_contributor_expr_tmp[1:cellgroup_number,i] <- apply(expr_gene[1:cellgroup_number,ind_i],MARGIN = 1, function(x) mean(x))
      } else {
        lig_contributor_expr_tmp[1:cellgroup_number,i] <- 1
      }
    }
    if(length(lig_contributor_coeff)>1){
      rep_coeff <- rep(1:length(lig_contributor_coeff),lig_contributor_coeff)
      lig_contributor_expr <- apply(lig_contributor_expr_tmp[,rep_coeff], MARGIN = 1, function(x) exp(mean(log(x))))} else {
        lig_contributor_expr <- lig_contributor_expr_tmp
      }
    names(lig_contributor_expr) <- cell.rownames
    # receptor
    receptor_subunit_expr_tmp <- matrix(0,cellgroup_number,length(receptor_subunit_coeff))
    for(i in 1:length(receptor_subunit_coeff)){
      ind_i <- which(receptor_subunit_group==i)
      if(length(ind_i)==1){
        receptor_subunit_expr_tmp[1:cellgroup_number,i] <- expr_gene[1:cellgroup_number,length(lig_contributor_new)+ind_i]
      } else if(length(ind_i)>1) {
        receptor_subunit_expr_tmp[1:cellgroup_number,i] <- apply(expr_gene[1:cellgroup_number,length(lig_contributor_new)+ind_i],MARGIN = 1, function(x) mean(x))
      } else {
        receptor_subunit_expr_tmp[1:cellgroup_number,i] <- 1
      }
    }
    if(length(receptor_subunit_coeff)>1){
      rep_coeff <- rep(1:length(receptor_subunit_coeff),receptor_subunit_coeff)
      receptor_subunit_expr <- apply(receptor_subunit_expr_tmp[,rep_coeff], MARGIN = 1, function(x) exp(mean(log(x))))} else {
        receptor_subunit_expr <- receptor_subunit_expr_tmp
      }
    names(receptor_subunit_expr) <- cell.rownames
    # # targets_up
    # if(length(targets_up)>0){
    #   targets_up_expr <- as.data.frame(cal_expr_by_group(df,targets_up,mean_method))
    #   targets_up_max <- apply(targets_up_expr,2,FUN = max)
    #   targets_up_nonzero <- which(targets_up_max>0)
    #   if(length(targets_up_nonzero)>0){
    #     targets_up_avg <- apply(targets_up_expr[,targets_up_nonzero,drop=FALSE],1,FUN = mean)}
    # }
    # # targets_down
    # if(length(targets_down)>0){
    #   targets_down_expr <- as.data.frame(cal_expr_by_group(df,targets_down,mean_method ))
    #   targets_down_max <- apply(targets_down_expr,2,FUN = max)
    #   targets_down_nonzero <- which(targets_down_max>0)
    #   if(length(targets_down_nonzero)>0){
    #     targets_down_avg <- apply(targets_down_expr[,targets_down_nonzero,drop=FALSE],1,FUN = mean)}
    # }
    #
    # # targets_nichenet_avg <- rep(0,cellgroup_number)
    # if(length(targets_nichenet)>0){
    #   targets_nichenet <- targets_nichenet[1:N]
    #   targets_nichenet_expr <- cal_expr_by_group(df,names(targets_nichenet),mean_method)
    #   targets_nichenet_max <- apply(targets_nichenet_expr,2,FUN = max)
    #   # targets_nichenet_nonzero <- which(targets_nichenet_max>0)
    #   # if(length(targets_nichenet_nonzero)>0){
    #   targets_nichenet_avg <- as.matrix(targets_nichenet_expr[,names(targets_nichenet),with=FALSE]) %*% as.matrix(targets_nichenet/sum(targets_nichenet))#}
    # }
    # compute prob_mtx without downstream weight
    prob_mtx0 <- ((as.matrix(lig_contributor_expr))%*%t(as.matrix(receptor_subunit_expr)))
    colnames(prob_mtx0) <- cell.rownames
    rownames(prob_mtx0) <- cell.rownames
    # # downstream weight
    # exp_prob_mtx0 <- prob_mtx0^2/(prob_mtx0^2+1) #exp(-1/prob_mtx0);
    # targets_up_weight <- rep(1,cellgroup_number)%*%t(targets_up_avg);
    # # exp_targets_up_weight <- exp(-1/targets_up_weight)
    # # exp_targets_up_weight <- 0.01+exp_targets_up_weight*exp_prob_mtx0/(exp_prob_mtx0+exp_targets_up_weight+1e-16)
    # targets_up_weight <- 1+prob_mtx0*targets_up_weight^1/(targets_up_weight^1+prob_mtx0+1e-6)#*(length(targets_up_avg))^0.5
    # targets_down_weight <- rep(1,cellgroup_number)%*%t(targets_down_avg);
    # # exp_targets_down_weight <- exp(-targets_down_weight)
    # # exp_targets_down_weight <- exp_targets_down_weight*exp_prob_mtx0/(exp_prob_mtx0+exp_targets_down_weight+1e-16)
    # targets_down_weight <- 1-targets_down_weight^1/(targets_down_weight^1+prob_mtx0+1e-6)
    # targets_nichenet_weight <- rep(1,cellgroup_number)%*%t(targets_nichenet_avg);
    # targets_nichenet_weight <- 1+1*prob_mtx0*targets_nichenet_weight^1/(targets_nichenet_weight^1+prob_mtx0+1e-6)
    # #targets_nichenet_weight <- prob_mtx0*targets_nichenet_weight^1/(targets_nichenet_weight^1+prob_mtx0+1e-6)
    # # compute prob_mtx with downstream weight
    # # prob_mtx <- prob_mtx0^1*targets_up_weight*targets_down_weight
    # prob_mtx <- prob_mtx0^1*targets_nichenet_weight
    prob_mtx <- prob_mtx0
  } else if(method=='CellChat'){
    if(length(lig_contributor_new)>1){
      lig_contributor_expr <- apply(expr_gene[1:cellgroup_number,lig_contributor_new],MARGIN = 1, function(x) exp(mean(log(x))))
    } else {  lig_contributor_expr <- expr_gene[1:cellgroup_number,lig_contributor_new] }
    if(length(receptor_subunit_new)>1){
      receptor_subunit_expr <- apply(expr_gene[1:cellgroup_number,receptor_subunit_new], MARGIN = 1, function(x) exp(mean(log(x))))
    } else {   receptor_subunit_expr <- expr_gene[1:cellgroup_number,receptor_subunit_new]}
    names(lig_contributor_expr) <- cell.rownames
    names(receptor_subunit_expr) <- cell.rownames
    prob_mtx0 <- ((as.matrix(lig_contributor_expr))%*%t(as.matrix(receptor_subunit_expr)))
    hill <- function(x,K){
      return(x^2/(x^2+K^2))
    }
    prob_mtx0 <- hill((as.matrix(lig_contributor_expr)),K)%*%hill(t(as.matrix(receptor_subunit_expr)),K)
    colnames(prob_mtx0) <- cell.rownames
    rownames(prob_mtx0) <- cell.rownames
    prob_mtx <- prob_mtx0
  } else if(method=='CellPhoneDB'){
    if(length(lig_contributor_new)>1){
      lig_contributor_expr <- apply(expr_gene[1:cellgroup_number,lig_contributor_new],MARGIN = 1, function(x) min(x))
    } else {  lig_contributor_expr <- expr_gene[1:cellgroup_number,lig_contributor_new] }
    if(length(receptor_subunit_new)>1){
      receptor_subunit_expr <- apply(expr_gene[1:cellgroup_number,receptor_subunit_new], MARGIN = 1, function(x) min(x))
    } else {   receptor_subunit_expr <- expr_gene[1:cellgroup_number,receptor_subunit_new]}
    names(lig_contributor_expr) <- cell.rownames
    names(receptor_subunit_expr) <- cell.rownames
    prob_mtx0_lig <- as.matrix(lig_contributor_expr)%*%matrix(rep(1,length(cell.rownames)),nrow = 1,ncol = length(cell.rownames))
    prob_mtx0_rec <- matrix(rep(1,length(cell.rownames)),nrow = length(cell.rownames),ncol =1 )%*%t(as.matrix(receptor_subunit_expr))
    # prob_mtx0 <- ((as.matrix(lig_contributor_expr))%*%matrix(rep(1,length(cell.rownames)),nrow = 1,ncol = length(cell.rownames))+ matrix(rep(1,length(cell.rownames)),nrow = length(cell.rownames),ncol =1 )%*%t(as.matrix(receptor_subunit_expr)))
    prob_mtx0 <- (prob_mtx0_lig+ prob_mtx0_rec)/2
    prob_mtx0 <- prob_mtx0*(prob_mtx0_lig>0)*(prob_mtx0_rec>0)
    colnames(prob_mtx0) <- cell.rownames
    rownames(prob_mtx0) <- cell.rownames
    prob_mtx <- prob_mtx0
  }
  rownames(prob_mtx) <- cell.rownames
  colnames(prob_mtx) <- cell.rownames
  prob_mtx <- prob_mtx[sender,receiver]
  prob_mtx0 <- prob_mtx0[sender,receiver]
  # prob_mtx <- exp_prob_mtx0*targets_up_weight*targets_down_weight
  # prob_mtx <- exp_prob_mtx0*exp_targets_up_weight*exp_targets_down_weight
  my_list <- list("prob_mtx" = prob_mtx, "ligand_score" = lig_contributor_expr, "receptor_score" = receptor_subunit_expr, "cell_group"=rownames(expr_gene),
                  "prob_mtx0" = prob_mtx0, "targets_up_score" = targets_up_avg, "targets_down_score" = targets_down_avg, "targets_nichenet_score"=targets_nichenet_avg)
  return(my_list)
}


#' Calculate the communication strength matrix for a single ligand-target pair (with permutation test)
#' @param df a data frame or data.table containing signaling gene expression (row: cell; column: genes; the last column is cell group labels)
#' @param sender sending cell groups
#' @param receiver receiving cell groups
#' @param interaction_name name of the ligand-target interaction
#' @param lig_contributor a vector with gene symbols that are regarded as ligand contributors
#' @param lig_contributor_group a vector indicating the grouping of ligand contributors
#' @param lig_contributor_coeff a vector with stoichiometry of each group of ligand contributors for calculating ligand abundance
#' @param receptor_subunit_new a vector with updated gene symbols (removing those excluded from the expression data) that are regarded as target subunits
#' @param receptor_subunit_group a vector indicating the grouping of target subunits
#' @param receptor_subunit_coeff a vector with stoichiometry of each group of target subunits for calculating target abundance
#' @param M number of permutation tests
#' @param fdr cutoff value that determines the significance of detected link
#' @param targets_up deprecated
#' @param targets_down deprecated
#' @param targets_nichenet deprecated
#' @param N deprecated
#' @param K deprecated
#'
#' @return
#' @export
#'
#' @examples
neuron_chat_downstream <- function(df,sender,receiver,interaction_name,lig_contributor,lig_contributor_group,lig_contributor_coeff,
                                   receptor_subunit,receptor_subunit_group,receptor_subunit_coeff,targets_up,targets_down,targets_nichenet,N,M,fdr,K=0.5,method=NULL,mean_method=NULL){
  ind_lig <- which(lig_contributor %in% colnames(df));lig_contributor_new <- lig_contributor[ind_lig]
  ind_re <- which(receptor_subunit %in% colnames(df));receptor_subunit_new <- receptor_subunit[ind_re]
  lig_contributor_group <- lig_contributor_group[ind_lig]
  receptor_subunit_group <- receptor_subunit_group[ind_re]
  targets_up <- targets_up[which(targets_up %in% colnames(df))]
  targets_down <- targets_down[which(targets_down %in% colnames(df))]
  targets_nichenet <- targets_nichenet[which(names(targets_nichenet) %in% colnames(df))]
  my_list <- cal_prob_mtx_downstream(df,sender,receiver,lig_contributor_new,lig_contributor_group,lig_contributor_coeff, receptor_subunit_new,receptor_subunit_group,receptor_subunit_coeff,targets_up,targets_down,targets_nichenet,N,K,method,mean_method)
  prob_mtx <- my_list$prob_mtx
  # FC.lig <- log2(max(my_list$ligand_score+1e-3)/min(my_list$ligand_score+1e-3))
  # FC.rec <- log2(max(my_list$receptor_score+1e-3)/min(my_list$receptor_score+1e-3))
  FC.lig <- max(my_list$ligand_score) - min(my_list$ligand_score)
  FC.rec <- max(my_list$receptor_score) - min(my_list$receptor_score)
  FC <- max(FC.lig,FC.rec)
  #FC <- min(FC.lig,FC.rec)
  ## permutation
  if(M==0 | max(prob_mtx)==0){net <- prob_mtx;pvalue <- NA*net} else {
    #prob_mtx_permutation <- array(0,c(dim(prob_mtx)[1],dim(prob_mtx)[2],M))
    pvalue <- matrix(0, dim(prob_mtx)[1],dim(prob_mtx)[2]);
    # for(j in 1:M){
    #   df_j <- df;
    #   df_j$cell_subclass[df$cell_subclass %in% sender] <- sample(df$cell_subclass[df$cell_subclass %in% sender],length(df$cell_subclass[df$cell_subclass %in% sender]), FALSE)
    #   df_j$cell_subclass[df$cell_subclass %in% receiver] <- sample(df$cell_subclass[df$cell_subclass %in% receiver],length(df$cell_subclass[df$cell_subclass %in% receiver]), FALSE)
    #   my_list_temp <- cal_prob_mtx_downstream(df_j,sender,receiver,lig_contributor_new,lig_contributor_group,lig_contributor_coeff,receptor_subunit_new,receptor_subunit_group,receptor_subunit_coeff,targets_up,targets_down,targets_nichenet,N)
    #   prob_mtx_permutation[,,j] <- my_list_temp$prob_mtx
    #   pvalue <- pvalue+1*(prob_mtx_permutation[,,j]>prob_mtx)
    # }
    df_j <- df;
    prob_mtx_permutation <- sapply(seq_along(1:M), function(j) {
      df_j$cell_subclass[df$cell_subclass %in% sender] <- sample(df$cell_subclass[df$cell_subclass %in% sender],length(df$cell_subclass[df$cell_subclass %in% sender]), FALSE)
      df_j$cell_subclass[df$cell_subclass %in% receiver] <- sample(df$cell_subclass[df$cell_subclass %in% receiver],length(df$cell_subclass[df$cell_subclass %in% receiver]), FALSE)
      my_list_temp <- cal_prob_mtx_downstream(df_j,sender,receiver,lig_contributor_new,lig_contributor_group,lig_contributor_coeff,receptor_subunit_new,receptor_subunit_group,receptor_subunit_coeff,targets_up,targets_down,targets_nichenet,N,K,method,mean_method)
      prob_mtx_permutation_j <- (my_list_temp$prob_mtx > prob_mtx)*1
      prob_mtx_permutation_j
    },simplify = 'array')
    pvalue <- apply(prob_mtx_permutation,1:2,sum)
    pvalue <- pvalue/M
    pvalue_v <- c(pvalue)
    ## FDR correction Benjamini-Hochberg Procedure
    pvalue_v_sort <- sort.int(pvalue_v, decreasing = FALSE, index.return = TRUE)
    alpha_i <- fdr*(1:length(pvalue_v))/length(pvalue_v)
    k <- which(pvalue_v_sort$x<alpha_i)
    prob_mtx_sig <- 0*pvalue
    prob_mtx_sig[pvalue_v_sort$ix[k]]<-prob_mtx[pvalue_v_sort$ix[k]]
    # if(max(prob_mtx>0)){
    # net <- prob_mtx/max(prob_mtx)} else {net <- prob_mtx}
    net <- prob_mtx_sig;
    # net_old <-net
    # net[pvalue>0.05] <- 0
    # net.rownames <- rownames(prob_mtx)
    }
  # net.rownames <- sort(unique(target_df$cell_subclass))
  # rownames(net)<- net.rownames;colnames(net)<- net.rownames
  # plot
  # library(GGally);library(network);library(sna);library(ggplot2);library(reshape2);library(circlize)
  # # q90 <- quantile(net,0.95);net[net<q90] <-0
  # group = structure(rep(1,dim(net)[1]),names=net.rownames)
  # par(mfrow=c(1,1))
  # netVisual_circle_zw(net, color.use = scPalette(nrow(net)),group=group, title.name = interaction_name, sources.use = net.rownames, targets.use = net.rownames, remove.isolate = FALSE, top = 1,
  #                     weight.scale = FALSE, vertex.weight =max(colSums(net)/(max(colSums(net))+0.01),0.1), vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= 'black',
  #                     edge.weight.max = NULL, edge.width.max=5, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
  #                     edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2,
  #                     arrow.width=1,arrow.size = 0.8)
  # ## barplot for ligand score & receptor score
  # data <- data.frame(cell_type=my_list$cell_group,ligand_score=my_list[[2]],receptor_score=my_list[[3]])
  # g1<- ggplot(data, aes(x=cell_type, y=ligand_score)) + geom_bar(stat = "identity")
  # g2<- ggplot(data, aes(x=cell_type, y=receptor_score)) + geom_bar(stat = "identity")
  # print(g1/g2)
  # col_map=colorspace::diverge_hsv(10);
  # par(mfrow=c(1,1))
  # col_map = brewer.pal(8,"YlOrBr")
  # # col_map = colorRamp2(seq(min(net), max(net), length = 3), c( "#EEEEEE",'pink',"red"), space = "RGB")
  # column_ha = HeatmapAnnotation( Rec_score = anno_barplot(my_list[[3]],border='FALSE',height=unit(2.5,'cm'),gp = gpar(fill = "lightskyblue1")),
  #                                annotation_name_side='left')
  # row_ha = rowAnnotation(Lig_score = anno_barplot(my_list[[2]],border='FALSE',width=unit(2.5,'cm'),gp = gpar(fill = "lightskyblue1")))
  # left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill =2:(length(unique(big_group))+1)), labels = unique(big_group)))
  # bottom_Annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill =2:(length(unique(big_group))+1)), labels = unique(big_group)))
  # if(max(net)>0){
  #   h1 <- Heatmap(net, name = "Prob", left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
  #                 top_annotation = column_ha, right_annotation = row_ha,col = col_map,
  #                 cluster_rows = FALSE,cluster_columns=FALSE,
  #                 split= big_group[rownames(net)],column_split = big_group[rownames(net)],
  #                 row_names_side='left',column_names_side='bottom',
  #                 row_title='Sender',row_title_side='left',row_title_gp = gpar(fontsize = 16),
  #                 column_title='Receiver',column_title_side = "bottom",column_title_gp = gpar(fontsize = 16))#,
  #   #   width = unit(11,'cm'),
  #   #   height = unit(11,'cm'))
  #   # print(h1)
  #   draw(h1, column_title = interaction_name, column_title_gp = gpar(fontsize = 16))}
  list_return <- list(net=net,FC=FC,pvalue=pvalue,net0=my_list$prob_mtx, info=sum(net),ligand.abundance = my_list$ligand_score, target.abundance = my_list$receptor_score)
  return(list_return)
  # net; color.use = NULL;title.name = NULL; sources.use = NULL; targets.use = NULL; remove.isolate = FALSE; top = 1;
  # weight.scale = FALSE; vertex.weight = 20; vertex.weight.max = NULL; vertex.size.max = NULL; vertex.label.cex=1;vertex.label.color= "black";
  # edge.weight.max = NULL; edge.width.max=8; alpha.edge = 0.6; label.edge = FALSE;edge.label.color='black';edge.label.cex=0.8;
  # edge.curved=0.2;shape='circle';layout=in_circle(); margin=0.2; vertex.size = NULL;
  # arrow.width=1;arrow.size = 0.2
  # CellChat::netVisual_chord_cell_internal(net, color.use = scPalette(nrow(net)), sources.use = net.rownames, targets.use = net.rownames,
  #                                         group=group,remove.isolate = FALSE,scale = FALSE,
  #                                         title.name = interaction_name)
  # chordDiagram(net, group = group, grid.col = group,
  #              reduce=-1,scale=FALSE,link.visible = net > 0.05,
  #              small.gap=1, big.gap=10,
  #              directional = -1, direction.type = c("diffHeight","arrows"),
  #              link.arr.type = 'big.arrow')
}

#' Calculate the communication strength matrix for all ligand-target pairs
#' @param object a NeuronChat object
#' @param sender sending cell groups
#' @param receiver receiving cell groups
#' @param M number of permutation tests
#' @param fdr cutoff value that determines the significance of detected link
#' @param N deprecated
#' @param K deprecated
#' @param strict indicating whether the calculation of communication score for an ligand-target pair requires that input data include at least one gene for each contributing gene group. strict=1 means 'require',  strict=0 means not 'require'. Default value is 1
#'
#' @return
#' @export
#'
#' @examples
run_NeuronChat <- function(object=object,sender=NULL,receiver=NULL,N=0,M=100,fdr=0.05,K=0.5,method=NULL,mean_method=NULL,strict=1){
start.time <- Sys.time()
# interactionDB <- readRDS('interactionDB_interactor.rds')
# interactionDB <- interactionDB_interactor_full_synGO
interactionDB <- object@DB
# ligand_target_matrix_mgi <- readRDS('~/Documents/CellChat/Rcode_cellchat_VISp_ALM/ligand_target_matrix_mgi.rds')
net0_all <- list() # original probability mtx
net_all <- list() # probability mtx filtered by pvalue with cutoff 0.05
pvalue_all <- list()
FC_all <- rep(0,length(interactionDB))
info_all <- rep(0,length(interactionDB))
net.rownames <- sort(unique(object@data.signaling$cell_subclass),method='radix')
if(is.null(sender)){sender <- net.rownames};if(is.null(receiver)){receiver <- net.rownames}
ligand.abundance_all <- matrix(0,nrow=length(net.rownames),ncol=length(interactionDB))
target.abundance_all <- ligand.abundance_all
for(j in 1:length(interactionDB)){
  interaction_name <- names(interactionDB)[j]
  lig_contributor <- interactionDB[[j]]$lig_contributor
  receptor_subunit <- interactionDB[[j]]$receptor_subunit
  lig_contributor_group <- interactionDB[[j]]$lig_contributor_group
  lig_contributor_coeff <- interactionDB[[j]]$lig_contributor_coeff
  receptor_subunit_group <- interactionDB[[j]]$receptor_subunit_group
  receptor_subunit_coeff <- interactionDB[[j]]$receptor_subunit_coeff
  # targets_up <- interactionDB[[j]]$targets_up
  # targets_up <- interactionDB[[j]]$interactors
  # targets_down <- interactionDB[[j]]$targets_down
  # targets_up <-  character(0)# interactionDB[[j]]$interactors
  targets_up  <- character(0)
  targets_down <- character(0)
  #if(receptor_subunit %in% colnames(ligand_target_matrix_mgi)){
    # targets_nichenet <- sort(ligand_target_matrix_mgi[,receptor_subunit],decreasing = T)
    targets_nichenet <- character(0)
  #} else {targets_nichenet <- character(0)}
  lig_boolean_group <- (lig_contributor %in% names(object@data.signaling))*lig_contributor_group;
  rec_boolean_group <- (receptor_subunit %in% names(object@data.signaling))*receptor_subunit_group;
  if(strict==1){
    lig_boolean <- prod(unique(lig_contributor_group) %in% lig_boolean_group)
    rec_boolean <- prod(unique(receptor_subunit_group) %in% rec_boolean_group)
  } else {
    lig_boolean <- sum(unique(lig_contributor_group) %in% lig_boolean_group)>0
    rec_boolean <- sum(unique(receptor_subunit_group) %in% rec_boolean_group)>0
  }
  if(prod(lig_boolean,rec_boolean)==0){net_all[[j]]=matrix(0,nrow=length(sender),ncol=length(receiver),
                                                                                                dimnames=list(sender,receiver))
  #if(prod((c(lig_contributor,receptor_subunit) %in% names(object@data.signaling)))==0){net_all[[j]]=matrix(0,nrow=length(sender),ncol=length(receiver),
  #                                                                                               dimnames=list(sender,receiver))
  net0_all[[j]] <- net_all[[j]]
  pvalue_all[[j]] <- net_all[[j]]
  } else
  {tmp <- neuron_chat_downstream(object@data.signaling,sender,receiver,interaction_name,lig_contributor,lig_contributor_group,lig_contributor_coeff,receptor_subunit,receptor_subunit_group,receptor_subunit_coeff,targets_up,targets_down,targets_nichenet,N,M,fdr,K,method,mean_method)
  net_all[[j]] <- tmp$net
  net0_all[[j]] <- tmp$net0
  pvalue_all[[j]] <- tmp$pvalue
  FC_all[j] <- tmp$FC
  info_all[j] <- tmp$info
  ligand.abundance_all[1:length(net.rownames),j] <- tmp$ligand.abundance
  target.abundance_all[1:length(net.rownames),j] <- tmp$target.abundance
  }
}
rownames(ligand.abundance_all) <- names(tmp$ligand.abundance);colnames(ligand.abundance_all) <- names(interactionDB)
rownames(target.abundance_all) <- names(tmp$target.abundance);colnames(target.abundance_all) <- names(interactionDB)
names(net_all) <- names(interactionDB)
names(net0_all) <- names(interactionDB)
# tmp_list <- future_lapply(seq_along(1:length(interactionDB)),function(j){
#   interaction_name <- names(interactionDB)[j]
#   lig_contributor <- interactionDB[[j]]$lig_contributor
#   receptor_subunit <- interactionDB[[j]]$receptor_subunit
#   lig_contributor_group <- interactionDB[[j]]$lig_contributor_group
#   lig_contributor_coeff <- interactionDB[[j]]$lig_contributor_coeff
#   receptor_subunit_group <- interactionDB[[j]]$receptor_subunit_group
#   receptor_subunit_coeff <- interactionDB[[j]]$receptor_subunit_coeff
#   # targets_up <- interactionDB[[j]]$targets_up
#   # targets_down <- interactionDB[[j]]$targets_down
#   targets_up <-  character(0)# interactionDB[[j]]$interactors
#   targets_down <- character(0)
#   if(receptor_subunit %in% colnames(ligand_target_matrix_mgi)){
#     # targets_nichenet <- sort(ligand_target_matrix_mgi[,receptor_subunit],decreasing = T)
#     targets_nichenet <- character(0)
#   } else {targets_nichenet <- character(0)}
#   if(prod((c(lig_contributor,receptor_subunit) %in% names(object@data.signaling)))==0){tmp$net <- matrix(0,nrow=length(sender),ncol=length(receiver),
#                                                                                                  dimnames=list(sender,receiver))
#   tmp$net0 <- tmp$net
#   tmp$pvalue<- tmp$net0
#   tmp$FC <- FC_all[j]
#   tmp$info <- 0
#   tmp$ligand.abundance <- ligand.abundance_all[1:length(net.rownames),j]
#   tmp$target.abundance <- ligand.abundance_all[1:length(net.rownames),j]
#   } else
#   {tmp <- neuron_chat_downstream(object@data.signaling,sender,receiver,interaction_name,lig_contributor,lig_contributor_group,lig_contributor_coeff,receptor_subunit,receptor_subunit_group,receptor_subunit_coeff,targets_up,targets_down,targets_nichenet,N,M)
#   }
#   tmp
# },future.seed = TRUE)
# for(j in length(tmp_list)){
# net_all[[j]] <- tmp_list[[j]]$net
# net0_all[[j]] <- tmp_list[[j]]$net0
# pvalue_all[[j]] <- tmp_list[[j]]$pvalue
# FC_all[j] <- tmp_list[[j]]$FC
# info_all[j] <- tmp_list[[j]]$info
# ligand.abundance_all[1:length(net.rownames),j] <- tmp_list[[j]]$ligand.abundance
# target.abundance_all[1:length(net.rownames),j] <- tmp_list[[j]]$target.abundance}

end.time <- Sys.time();time.taken <- end.time - start.time;
print(time.taken)
#list_return <- list(net_all=net_all,net0_all=net0_all,FC_all=FC_all,pvalue_all=pvalue_all)
#return(list_return)
object@net0 <- net0_all
object@pvalue <- pvalue_all
object@net <- net_all
object@fc <- FC_all
object@info <- info_all
object@ligand.abundance <- ligand.abundance_all
object@target.abundance <- target.abundance_all
return(object)
}

#' Aggregation of communication networks over all ligand-target pairs
#' @param net_list a list containing communication strength matrices for all ligand-target pairs, e.g., 'net' slot of NeuronChat object after run 'run_NeuronChat'
#' @param method one of 'weight','count','weighted_count','weighted_count2' and 'weight_threshold'
#' @param cut_off threshold used for methods 'weighted_count2' and 'weight_threshold'
#' @return
#' @export
#'
#' @examples
net_aggregation <- function(net_list,method=c('weight','count','weighted_count','weighted_count2','weight_threshold'),cut_off=0.05){
  method <- match.arg(method)
  if(method=='weighted_count'){
  net_aggregated <- 0*net_list[[1]]
  for(jj in 1:length(net_list)){
    net_aggregated <- net_aggregated+sum(net_list[[jj]])*(net_list[[jj]]>0);
  }} else if(method=='count'){
    net_aggregated <- 0*net_list[[1]]
    for(jj in 1:length(net_list)){
      net_aggregated <- net_aggregated+1*(net_list[[jj]]>0);
    }
  } else if(method=='weighted_count2'){
    net_aggregated <- 0*net_list[[1]]
    for(jj in 1:length(net_list)){
      net_aggregated <- net_aggregated+sum(net_list[[jj]])/(1e-6+sum(net_list[[jj]]>quantile(net_list[[jj]],probs = cut_off)))*(net_list[[jj]]>quantile(net_list[[jj]],probs = cut_off));
    }
  }else if(method=='weight_threshold'){
    net_aggregated <- 0*net_list[[1]]
    for(jj in 1:length(net_list)){
      net_aggregated <- net_aggregated+net_list[[jj]]*(net_list[[jj]]>quantile(net_list[[jj]],probs = cut_off));
    }
  }
   else {
    net_aggregated <- 0*net_list[[1]]
    for(jj in 1:length(net_list)){
      net_aggregated <- net_aggregated+net_list[[jj]];
    }
  }
  return(net_aggregated)
}
