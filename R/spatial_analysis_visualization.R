#' Calculation of cell proximity enrichment score for multiple tissue slices
#'
#' Calculation of cell proximity enrichment score & associated p-values for multiple tissue slices
#'
#' @param meta A dataframe containing spatial locations and cell type annotations of cells
#' @param celltype_label variable name (character) in meta, to indicate the cell type annotations of cells
#' @param centroid_x variable name (character) in meta, to indicate the x coordinates of cells
#' @param centroid_y variable name (character) in meta, to indicate the y coordinates of cells
#' @param slice_id variable name (character) in meta, to indicate slice id of cells
#' @param thresh_dist distance threshold (a scalar) used to determine the spatial neighbors of cells
#' @param permutation_number number of permutations (a scalar) used to calculate expected frequency matrices and p-values
#' @importFrom reshape2 melt dcast
#' @importFrom stats p.adjust
#' @return A list contains orginal frequency matrix, permutation frequency matrices, dataframe containing cell proximity enrichment score & pvalue, CPscore matrix, pvalue (adjusted) matrix
#' @export
cell_proximity_enrichment_score_multiple <- function(meta,celltype_label='subclass',centroid_x='centroid_x', centroid_y='centroid_y', slice_id='slice_id',
                                                     thresh_dist=400,permutation_number=1000){
  # meta <- ;thresh_dist=400;permutation_number=10;celltype_label='subclass';centroid_x='centroid_x';centroid_y='centroid_y';slice_id='slice_id';
  celltype <- sort(unique(as.character(meta[,celltype_label])),method='radix')
  ## calculating for each slice
  cellproximity_list <- lapply(unique(meta[,slice_id]),function(x){
    meta_slice <- meta[meta[,slice_id]==x,]
    meta_slice <- meta_slice[(meta_slice[,celltype_label] %in% celltype),]
    tmp <- cell_proximity_enrichment_score(meta_slice,celltype_label=celltype_label,centroid_x=centroid_x, centroid_y=centroid_x, thresh_dist=400,permutation_number=permutation_number,celltype= celltype)
    return(list(freq_mtx_simu_1K=tmp$freq_mtx_simu_1K,freq_mtx_orig=tmp$freq_mtx_orig))
  })
  ## summing over multiple slices
  freq_mtx_orig_sum <- 0*cellproximity_list[[1]]$freq_mtx_orig
  freq_mtx_simu_1K_sum <- 0*cellproximity_list[[1]]$freq_mtx_simu_1K
  for(j in 1:length(cellproximity_list)){
    freq_mtx_orig_sum <- freq_mtx_orig_sum+ cellproximity_list[[j]]$freq_mtx_orig
    freq_mtx_simu_1K_sum <- freq_mtx_simu_1K_sum + cellproximity_list[[j]]$freq_mtx_simu_1K
  }
  ## calculating observed-over-expected ratio & p-value
  freq_mtx_binary <- apply(freq_mtx_simu_1K_sum,3,function(x){(x>freq_mtx_orig_sum)*1},simplify = F)
  freq_mtx_binary <- simplify2array(freq_mtx_binary)
  freq_mtx_expected <- apply(freq_mtx_simu_1K_sum,MARGIN = c(1,2),FUN=mean)
  p_value <- apply(freq_mtx_binary,1:2,sum)/permutation_number; p_value_mtx <- p_value; p_value[lower.tri(p_value)] <- 0
  p_value_adjusted_mtx <- stats::p.adjust(c(p_value_mtx),method = 'BH'); dim(p_value_adjusted_mtx) <-  dim(p_value_mtx)

  CPscore <- (freq_mtx_orig_sum+1)/(freq_mtx_expected+1);  log2CPscore_mtx <- log2(CPscore); CPscore[lower.tri(CPscore)] <- 0

  CPscore_long <- reshape2::melt(CPscore,value.name = "CPscore")
  p_value_long <- reshape2::melt(p_value,value.name = "pvalue")
  CPs_df <- merge(CPscore_long,p_value_long);
  CPs_df <- CPs_df[CPs_df$CPscore>0,]
  CPs_df$log2CPscore <- log2(CPs_df$CPscore)
  CPs_df$group <- ifelse(CPs_df$log2CPscore>0,'enriched','depleted')
  CPs_df$p.adjusted <- stats::p.adjust(CPs_df$pvalue,method = 'BH')
  CPs_df$sig_0.05 <-''; CPs_df$sig_0.05[CPs_df$p.adjusted <0.05] <- '*'
  CPs_df <- CPs_df[order(CPs_df$CPscore,decreasing = T),]
  CPs_df$cellpair <- paste(CPs_df$Var1,CPs_df$Var2,sep='--')

  return(list(freq_mtx_simu_1K=freq_mtx_simu_1K_sum,freq_mtx_orig=freq_mtx_orig_sum,CPScore_df=CPs_df,log2CPscore_mtx=log2CPscore_mtx,p_value_adjusted_mtx=p_value_adjusted_mtx))
}

#' Calculation of cell proximity enrichment score for single tissue slice
#'
#' Calculation of cell proximity enrichment score & associated p-values for single tissue slice
#'
#' @param meta A dataframe containing spatial locations and cell type annotations of cells
#' @param celltype_label variable name (character) in meta, to indicate the cell type annotations of cells
#' @param centroid_x variable name (character) in meta, to indicate the x coordinates of cells
#' @param centroid_y variable name (character) in meta, to indicate the y coordinates of cells
#' @param thresh_dist distance threshold (a scalar) used to determine the spatial neighbors of cells
#' @param permutation_number number of permutations (a scalar) used to calculate expected frequency matrices and p-values
#' @importFrom reshape2 melt dcast
#' @importFrom stats p.adjust
#' @return A list contains orginal frequency matrix, permutation frequency matrices, dataframe containing cell proximity enrichment score & pvalue, CPscore matrix, pvalue (adjusted) matrix
#' @export
cell_proximity_enrichment_score_single <- function(meta,celltype_label='subclass',centroid_x='centroid_x', centroid_y='centroid_y',
                                                   thresh_dist=400,permutation_number=1000){
  celltype <- sort(unique(as.character(meta[,celltype_label])),method='radix')
  results <- cell_proximity_enrichment_score(meta,celltype_label=celltype_label,centroid_x=centroid_x, centroid_y=centroid_y,
                                             thresh_dist=400,permutation_number=1000,celltype= celltype)
  return(results)
}

#' Calculation of cell proximity enrichment score
#'
#' Calculation of cell proximity enrichment score & associated p-values (internal function)
#'
#' @param loc_centroid_subset_10K A dataframe containing spatial locations and cell type annotations of cells
#' @param celltype_label variable name (character) in loc_centroid_subset_10K, to indicate the cell type annotations of cells
#' @param centroid_x variable name (character) in loc_centroid_subset_10K, to indicate the x coordinates of cells
#' @param centroid_y variable name (character) in loc_centroid_subset_10K, to indicate the y coordinates of cells
#' @param thresh_dist distance threshold (a scalar) used to determine the spatial neighbors of cells
#' @param permutation_number number of permutations (a scalar) used to calculate expected frequency matrices and p-values
#' @param celltype uniqued cell type names (a vector)
#' @importFrom reshape2 melt dcast
#' @importFrom stats p.adjust
#' @return A list contains orginal frequency matrix, permutation frequency matrices, dataframe containing cell proximity enrichment score & pvalue, CPscore matrix, pvalue (adjusted) matrix
#' @export
cell_proximity_enrichment_score <- function(loc_centroid_subset_10K,celltype_label='subclass',centroid_x='centroid_x', centroid_y='centroid_y',
                                            thresh_dist=400,permutation_number=1000,celltype= celltype){
  # thresh_dist=400;permutation_number=1000;celltype_label='subclass';centroid_x='centroid_x';centroid_y='centroid_y';
  rownames(loc_centroid_subset_10K) <- NULL
  n_celltype <- length(celltype)
  coord_mtx <- as.matrix(loc_centroid_subset_10K[,c(centroid_x,centroid_y)])
  Euclidean_dist <- dist(coord_mtx, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  Euclidean_mtx <- as.matrix(Euclidean_dist)
  Euclidean_mtx[lower.tri(Euclidean_mtx,diag = T)] <- 0
  Euclidean_df <-  reshape2::melt(Euclidean_mtx)
  Euclidean_df <- Euclidean_df[Euclidean_df$value>0,]
  Euclidean_df_filter <- Euclidean_df[Euclidean_df$value <=thresh_dist,]
  idx_from <- match(Euclidean_df_filter$Var1,rownames(loc_centroid_subset_10K))
  idx_to <- match(Euclidean_df_filter$Var2,rownames(loc_centroid_subset_10K))

  # original freq matrix(observed)
  freq_mtx_orig <- matrix(0, n_celltype,n_celltype); dimnames(freq_mtx_orig) <- list(celltype,celltype)
  Euclidean_df_filter$from <- loc_centroid_subset_10K[,celltype_label][idx_from]
  Euclidean_df_filter$to <- loc_centroid_subset_10K[,celltype_label][idx_to]
  Euclidean_df_filter_cast <- reshape2::dcast(Euclidean_df_filter,from + to ~ . , fun.aggregate=length,value.var='value'); colnames(Euclidean_df_filter_cast)[3] <- 'freq'
  freq_mtx <- reshape2::dcast(Euclidean_df_filter_cast,from~ to,value.var='freq');freq_mtx[is.na(freq_mtx)] <-0
  rownames(freq_mtx) <- freq_mtx[,1]; freq_mtx <- freq_mtx[,-1]
  freq_mtx_orig[rownames(freq_mtx), colnames(freq_mtx)] <- as.matrix(freq_mtx);
  freq_mtx_orig <- freq_mtx_orig + t(freq_mtx_orig)- diag(diag(freq_mtx_orig))

  ## permutation freq matrix (expected)
  sampling_mtx <- sapply(1:permutation_number,FUN=function(x){sample(1:dim(loc_centroid_subset_10K)[1],dim(loc_centroid_subset_10K)[1])},simplify=T)
  freq_mtx_simu_1K <- sapply(1:permutation_number,FUN=function(x){
    freq_mtx_simu <- matrix(0, n_celltype,n_celltype); dimnames(freq_mtx_simu) <- list(celltype,celltype)
    Euclidean_df_filter$from <- loc_centroid_subset_10K[,celltype_label][sampling_mtx[,x]][idx_from]
    Euclidean_df_filter$to <- loc_centroid_subset_10K[,celltype_label][sampling_mtx[,x]][idx_to]
    Euclidean_df_filter_cast <- reshape2::dcast(Euclidean_df_filter,from + to ~ . , fun.aggregate=length,value.var='value'); colnames(Euclidean_df_filter_cast)[3] <- 'freq'
    freq_mtx <- reshape2::dcast(Euclidean_df_filter_cast,from~ to,value.var='freq')
    freq_mtx[is.na(freq_mtx)] <-0
    rownames(freq_mtx) <- freq_mtx[,1]; freq_mtx <- freq_mtx[,-1]
    freq_mtx_simu[rownames(freq_mtx), colnames(freq_mtx)] <- as.matrix(freq_mtx)
    freq_mtx_simu <- freq_mtx_simu + t(freq_mtx_simu)- diag(diag(freq_mtx_simu))
    freq_mtx_simu
  },simplify = F)
  freq_mtx_simu_1K <- simplify2array(freq_mtx_simu_1K)

  ## calculating observed-over-expected ratio & p-value
  freq_mtx_expected <- apply(freq_mtx_simu_1K,MARGIN = c(1,2),FUN=mean)
  CPscore <- (freq_mtx_orig+1)/(freq_mtx_expected+1);  log2CPscore_mtx <- log2(CPscore); CPscore[lower.tri(CPscore)] <- 0

  freq_mtx_binary <- apply(freq_mtx_simu_1K,3,function(x){(x>freq_mtx_orig)*1},simplify = F)
  freq_mtx_binary <- simplify2array(freq_mtx_binary)
  p_value <- apply(freq_mtx_binary,1:2,sum)/permutation_number; p_value[lower.tri(p_value)] <- 0

  CPscore_long <- reshape2::melt(CPscore,value.name = "CPscore")
  p_value_long <- reshape2::melt(p_value,value.name = "pvalue")
  CPs_df <- merge(CPscore_long,p_value_long);
  CPs_df <- CPs_df[CPs_df$CPscore>0,]
  CPs_df$log2CPscore <- log2(CPs_df$CPscore)
  CPs_df$group <- ifelse(CPs_df$log2CPscore>0,'enriched','depleted')
  CPs_df$p.adjusted <- stats::p.adjust(CPs_df$pvalue,method = 'BH')
  CPs_df$sig_0.05 <-''; CPs_df$sig_0.05[CPs_df$p.adjusted <0.05] <- '*'
  CPs_df <- CPs_df[order(CPs_df$CPscore,decreasing = T),]
  CPs_df$cellpair <- paste(CPs_df$Var1,CPs_df$Var2,sep='--')

  p_value_adjusted_mtx_tmp <- reshape2::dcast(CPs_df[,c('Var1','Var2','p.adjusted')],Var1~Var2,value.var='p.adjusted')
  rownames(p_value_adjusted_mtx_tmp) <- p_value_adjusted_mtx_tmp[,1];p_value_adjusted_mtx_tmp <- p_value_adjusted_mtx_tmp[,-1]
  p_value_adjusted_mtx_tmp[is.na(p_value_adjusted_mtx_tmp)] <- 0; #CP_net <- ifelse(CP_net>0,1,0)
  p_value_adjusted_mtx_tmp <- as.matrix(p_value_adjusted_mtx_tmp)
  p_value_adjusted_mtx_tmp <- p_value_adjusted_mtx_tmp + t(p_value_adjusted_mtx_tmp) - diag(diag(p_value_adjusted_mtx_tmp))
  p_value_adjusted_mtx <- log2CPscore_mtx; p_value_adjusted_mtx[rownames(p_value_adjusted_mtx_tmp),colnames(p_value_adjusted_mtx_tmp)] <- p_value_adjusted_mtx_tmp
  ## return value
  return(list(freq_mtx_simu_1K=freq_mtx_simu_1K,freq_mtx_orig=freq_mtx_orig,CPScore_df=CPs_df,log2CPscore_mtx=log2CPscore_mtx,p_value_adjusted_mtx=p_value_adjusted_mtx))
}

#' bar plot of spatial proximity for pairwise cell types
#'
#' bar plot of spatial proximity for pairwise cell types
#' @param cell_proximity_df a dataframe stored in the output list from function "cell_proximity_enrichment_score_single" or "cell_proximity_enrichment_score_multiple", e.g., cell_proximity$CPScore_df
#' @param font.size font.size for ggplot
#' @import ggplot2
#' @importFrom scales hue_pal
#' @return a ggplot object
#' @export
barplot_proximity <- function(cell_proximity_df,font.size=20){
CPs_df <- cell_proximity_df
p4 <- ggplot() + theme_bw() + geom_bar(data=CPs_df[order(CPs_df$log2CPscore,decreasing = T),], aes(x = reorder(cellpair, log2CPscore)  , y = log2CPscore, fill = group), stat = 'identity', show.legend = F) +
  labs(y ="cell proximity enrichment score", x = "interacting cell pair")  + coord_flip() + scale_fill_manual(values=rev(scales::hue_pal()(2))) +
  theme(text=element_text(size=font.size), #change font size of all text
        axis.text=element_text(size=font.size,color = 'black'), #change font size of axis text
        axis.title=element_text(size=font.size), #change font size of axis titles
        legend.text=element_text(size=font.size), #change font size of legend text
        legend.title=element_text(size=font.size)) #change font size of legend title
return(p4)
}

#' Circle plot of cell-cell proximity network
#'
#' Circle plot of cell-cell proximity network (adapted from CellChat https://github.com/sqjin/CellChat)
#'
#' The width of edges represent the absolute value of cell proximity enrichment score;
#' red line means enriched cell proximity while grey line means depleted cell proximity
#'
#' @param net_ori A weighted matrix representing the connections
#' @param group a vector indicating which cell types (cell subclass, e.g., L2/3 IT) belong to which big groups (cell class, e.g., Glutamatergic)
#' @param color.use Colors represent different cell groups
#' @param title.name the name of the title
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param idents.use a vector giving the index or the name of cell groups of interest.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.edge Whether or not shows the label of edges
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently “circle”, “square”,
#' “csquare”, “rectangle”, “crectangle”, “vrectangle”, “pie” (see
#' vertex.shape.pie), ‘sphere’, and “none” are supported, and only by the
#' plot.igraph command. “none” does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param layout The layout specification. It must be a call to a layout
#' specification function.
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
# #' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
#' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_ in_circle
#' @importFrom grDevices recordPlot
#' @import CellChat
#' @return  an object of class "recordedplot"
#' @export
netVisual_proximity <-function(net_ori, color.use = NULL,group=NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                   weight.scale = TRUE, vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = 10, vertex.label.cex=1.5,vertex.label.color= "black",
                                   edge.weight.max = NULL, edge.width.max=10, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                   edge.curved=0,shape=NULL,layout=in_circle(), margin=0.2, vertex.size = NULL,
                                   arrow.width=0,arrow.size = 0){
  # color.use = NULL;title.name = NULL; sources.use = NULL; targets.use = NULL; remove.isolate = FALSE; top = 1;
  # weight.scale = TRUE; vertex.weight = 1; vertex.weight.max = NULL; vertex.size.max = 10; vertex.label.cex=1.5;vertex.label.color= "black";
  # edge.weight.max = NULL; edge.width.max=10; alpha.edge = 0.6; label.edge = FALSE;edge.label.color='black';edge.label.cex=1.2;
  # edge.curved=0;shape='circle';layout=in_circle(); margin=0.2; vertex.size = NULL;arrow.width=0;arrow.size = 0
  # arrow.width=1;arrow.size = 0.2
  net_ori[lower.tri(net_ori)] <- 0
  net.names <- unique(c(rownames(net_ori),colnames(net_ori)))
  net <- matrix(0,nrow=length(net.names),ncol=length(net.names))
  dimnames(net) <- list(net.names,net.names)
  net[rownames(net_ori),colnames(net_ori)] <- net_ori
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0


  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  g1 <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  if (is.null(group)) {
    g <-igraph::permute(g1,match(V(g1)$name,net.names))
    group <- structure(rep(1,length(V(g))),names=names(V(g)))
    group <- group[names(V(g))]
  } else {
    g <-igraph::permute(g1,match(V(g1)$name,names(group)[order(group, names(group))]))
  }
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = CellChat::scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  # color.use.label  <- c('red','black')
  color.use.label <- c('black',CellChat::scPalette(length(unique(group))))[1:length(unique(group))]
  color.use.label <- structure(color.use.label,names=unique(group))
  # color.use.label <- color.use.label[-length(color.use.label)]
  # if(length(color.use.label)==1){ color.use.label='black'}
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  if(length(vertex.weight)==1){ igraph::V(g)$size<-vertex.weight} else {igraph::V(g)$size<-vertex.weight[names(V(g))]}
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$shape <- c('circle')#, 'square', 'csquare', 'rectangle', 'crectangle', 'vrectangle', 'pie', 'raster','sphere')[group[igraph::V(g)]]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  #igraph::V(g)$frame.color <- border.color.use[igraph::V(g)]
  igraph::V(g)$label.color <-color.use.label[group[names(V(g))]]#unlist(lapply(V(g), FUN= function(x){color.use.label[group[x]]}))
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    # igraph::E(g)$width<- 0.3+abs(igraph::E(g)$weight)/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*(E(g)$weight-min(E(g)$weight))
  }

  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight>0, 'red','grey');
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color,alpha.edge)
  igraph::E(g)$lty <- ifelse(igraph::E(g)$weight>0, 1,1)

  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+4
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  # plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
  #      vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica",vertex.label="") # "sans"
  # ## https://gist.github.com/ajhmohr/5337a5c99b504e4a243fad96203fa74f
  # text.pos <- coords*(1.15+strwidth(names(V(g)),font=12)/2)
  # angle = ifelse(atan(-(coords[,1]/coords[,2]))*(180/pi) < 0,  90 + atan(-(coords[,1]/coords[,2]))*(180/pi), 270 + atan(-coords[,1]/coords[,2])*(180/pi))
  # for (i in 1:length(V(g))) {
  #   text(x=text.pos[i,1], y=text.pos[i,2], labels=V(g)$name[i], adj=NULL, pos=NULL, cex=1.2, col=igraph::V(g)$label.color[i], srt=angle[i], xpd=T)
  # }

  if (!is.null(title.name)) {
    text(0,1.6,title.name, cex = 1.5)
  }
  gg <- recordPlot()
  return(gg)
}

#' Spatial plot of cells with their annotations
#'
#' Spatial plot of cells with their annotations
#'
#' @param meta A dataframe containing spatial locations and cell type annotations of cells
#' @param celltype_label variable name (character) in meta, to indicate the cell type annotations of cells
#' @param centroid_x variable name (character) in meta, to indicate the x coordinates of cells
#' @param centroid_y variable name (character) in meta, to indicate the y coordinates of cells
#' @import ggplot2
#' @importFrom scales hue_pal
#' @return a ggplot object
#' @export
spatialMap <- function(meta,celltype_label='subclass',centroid_x='centroid_x', centroid_y='centroid_y', pt.size=5,pt.alpha=0.8,font.size=20,legend.symbol.size=5){
  attr_x <- celltype_label
  p2 <- ggplot2::ggplot() + theme_bw()
  p2 <- p2+ ggplot2::geom_point(data = meta[order(meta[,attr_x],decreasing = F),], aes_string(x = centroid_x, y = centroid_y,colour = attr_x),
                                show.legend = T, size = pt.size, alpha =pt.alpha) + scale_color_manual(values=colours) +
    labs(x ="x coordinates", y = "y coordinates",color=attr_x) +  guides(color = guide_legend(override.aes = list(size=legend.symbol.size))) +
    theme(text=element_text(size=font.size), #change font size of all text
          axis.text=element_text(size=font.size), #change font size of axis text
          axis.title=element_text(size=font.size), #change font size of axis titles
          plot.title=element_text(size=font.size,hjust = 0.5), #change font size of plot title
          legend.text=element_text(size=font.size), #change font size of legend text
          legend.title=element_text(size=font.size)) #change font size of legend title
return(p2)
}

#' Spatial distribution of cell types along sptial axis
#'
#' Spatial distribution of cell types along sptial axis
#'
#' @param meta A dataframe containing spatial locations and cell type annotations of cells
#' @param celltype_label variable name (character) in meta, to indicate the cell type annotations of cells
#' @param centroid variable name (character) in meta, to indicate the  x or y coordinates of cells
#' @import ggplot2
#' @importFrom scales hue_pal
#' @return a ggplot object
#' @export
spatialDistribution <- function(meta,celltype_label='subclass',centroid='centroid_y', curve.alpha=0.5,font.size=20,legend.symbol.size=5){
attr_x <- celltype_label
p3 <- ggplot(data=meta[order(meta[,attr_x],decreasing = T),], aes_string(x=centroid, group=attr_x, fill=attr_x)) +
  geom_density(adjust=1.5, alpha=curve.alpha) + scale_fill_manual(values=colours) +
  labs(y ="density", x = "y coordinate",fill=attr_x) +  guides(color = guide_legend(override.aes = list(size=legend.symbol.size))) + coord_flip() +
  theme(text=element_text(size=font.size), #change font size of all text
        axis.text=element_text(size=font.size,color = 'black'), #change font size of axis text
        axis.title=element_text(size=font.size), #change font size of axis titles
        legend.text=element_text(size=font.size), #change font size of legend text
        legend.title=element_text(size=font.size)) #change font size of legend title
return(p3)
}
