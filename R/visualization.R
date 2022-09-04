####################
####################
#################### Visualization of communication network

####################
#################### For Heatmap plot
#' Heatmap plot with ligand abundance and target abundance for single ligand-target pair
#' @param object a NeuronChat object
#' @param interaction_name the ligand-target interaction used to plot
#' @param sender.names sending cell groups used for plot
#' @param receiver.names receiving cell groups used for plot
#' @param group a vector indicating which cell types (cell subclass, e.g., L2/3 IT) belong to which big groups (cell class, e.g., Glutamatergic)
#' @import circlize
#' @import ComplexHeatmap
#' @importFrom CellChat scPalette
#' @return
#' @export
heatmap_single <- function(object, interaction_name,group=NULL, sender.names=NULL,receiver.names=NULL){
if(is.null(sender.names)){sender.names=rownames(object@net[[1]])}
if(is.null(receiver.names)){receiver.names=colnames(object@net[[1]])}
col_map = colorRamp2(c(0,max(object@net[[interaction_name]][sender.names,receiver.names])/2, max(object@net[[interaction_name]][sender.names,receiver.names])), c("blue", "white", "red"))
column_ha = HeatmapAnnotation(target = anno_barplot(object@target.abundance[receiver.names,interaction_name],border = F,height = unit(2,'cm'),gp = gpar(fill = "lightskyblue1")),annotation_name_side='left',annotation_name_rot = 0,annotation_label = 'target \n abundance')
row_ha = rowAnnotation(ligand = anno_barplot(object@ligand.abundance[sender.names,interaction_name],border = F,width = unit(2,'cm'),gp = gpar(fill = "lightskyblue1")),annotation_label = 'ligand \n abundance')
if(is.null(group)){
ComplexHeatmap::Heatmap(object@net[[interaction_name]][sender.names,receiver.names], name = "Commu. \n Prob.", top_annotation = column_ha, right_annotation = row_ha,
                              cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                              column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                              heatmap_legend_param = list(color_bar='continuous'))} else {
left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill =scPalette(length(unique(group[sender.names])))), labels = sort(unique(group[sender.names]))))
bottom_Annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill =scPalette(length(unique(group[receiver.names])))), labels = sort(unique(group[receiver.names]))))
ComplexHeatmap::Heatmap(object@net[[interaction_name]][sender.names,receiver.names], name = "Commu. \n Prob.",
                        top_annotation = column_ha, right_annotation = row_ha,
                        left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
                        split=as.character(group[sender.names]),column_split = as.character(group[receiver.names]),
                        cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                        column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                        heatmap_legend_param = list(color_bar='continuous'))
                              }
}


#' Heatmap plot for the aggregated communication
#' @param object a NeuronChat object
#' @param interaction_use ligand-target interaction indexes  used for aggregation ('all' means use all interation pairs)
#' @param group a vector indicating which cell types (cell subclass, e.g., L2/3 IT) belong to which big groups (cell class, e.g., Glutamatergic)
#' @param method method used for aggregation; see also function `net_aggregation`
#' @param cut_off threshold used for aggregation; see also function `net_aggregation`
#' @param sender.names sending cell groups used for plot
#' @param receiver.names receiving cell groups used for plot
#' @import circlize
#' @import ComplexHeatmap
#' @importFrom CellChat scPalette
#' @return
#' @export
heatmap_aggregated <- function(object, method=c('weight','count','weighted_count','weighted_count2','weight_threshold'),cut_off=0.05, interaction_use='all',group=NULL, sender.names=NULL,receiver.names=NULL){
  method <- match.arg(method)
  if(is.null(sender.names)){sender.names=rownames(object@net[[1]])}
  if(is.null(receiver.names)){receiver.names=colnames(object@net[[1]])}
  if(interaction_use=='all') {
  net_aggregated <- net_aggregation(object@net,method=method,cut_off=cut_off)} else {net_aggregated <- net_aggregation(object@net[interaction_use],method=method,cut_off=cut_off)}
  col_map = colorRamp2(c(0,max(net_aggregated[sender.names,receiver.names])/2, max(net_aggregated[sender.names,receiver.names])), c("blue", "white", "red"))
  if(is.null(group)){
    ComplexHeatmap::Heatmap(net_aggregated[sender.names,receiver.names], name = method,
                            cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                            column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                            heatmap_legend_param = list(color_bar='continuous'))} else {
                              left_Annotation = rowAnnotation(foo = anno_block(gp = gpar(fill =scPalette(length(unique(group[sender.names])))), labels = sort(unique(group[sender.names]))))
                              bottom_Annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill =scPalette(length(unique(group[receiver.names])))), labels = sort(unique(group[receiver.names]))))
                              ComplexHeatmap::Heatmap(net_aggregated[sender.names,receiver.names], name = "Commu. \n Prob.",
                                                      left_annotation =left_Annotation,bottom_annotation = bottom_Annotation,
                                                      cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                                                      split=as.character(group[sender.names]),column_split = as.character(group[receiver.names]),
                                                      column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                                                      heatmap_legend_param = list(color_bar='continuous'))
                            }
}
#' A set of plots illustrating the communication network as well as violin plots for related genes for single ligand-target pair
#'
#' A set of plots illustrating the communication network as well as violin plots for related genes for single ligand-target pair
#' @param object a NeuronChat object
#' @param interaction_name the ligand-target interaction used to plot
#' @param sender.names sending cell groups used for plot
#' @param receiver.names receiving cell groups used for plot
#' @import circlize
#' @import ComplexHeatmap
#' @import Seurat
#' @import SeuratObject
#' @return
#' @export
lig_tar_heatmap <- function(object, interaction_name,width.vector=c(0.3,0.34,0.31), sender.names=NULL,receiver.names=NULL){
if(is.null(sender.names)){sender.names=rownames(object@net[[1]])}
if(is.null(receiver.names)){receiver.names=colnames(object@net[[1]])}
## create Seurat object from NeuronChat object
  meta <- data.frame(cell_subclass=object@data.signaling$cell_subclass)
  rownames(meta) <- rownames(object@data.signaling)
  seurat_VISp <- Seurat::CreateSeuratObject(counts = object@data, meta.data = meta)
  seurat_VISp$groups <- meta$cell_subclass
  SeuratObject::Idents(seurat_VISp) <- meta$cell_subclass
  levels(seurat_VISp) <- sort(unique(c(sender.names,receiver.names)),method='radix')
## plot #1: heatmap
  col_map = circlize::colorRamp2(c(0,max(object@net[[interaction_name]][sender.names,receiver.names])/2, max(object@net[[interaction_name]][sender.names,receiver.names])), c("blue", "white", "red"))
  column_ha = HeatmapAnnotation(target = anno_barplot(object@target.abundance[receiver.names,interaction_name],border = F,height = unit(2,'cm'),gp = gpar(fill = "lightskyblue1")),annotation_name_side='left',annotation_name_rot = 0,annotation_label = 'target \n abundance')
  row_ha = rowAnnotation(ligand = anno_barplot(object@ligand.abundance[sender.names,interaction_name],border = F,width = unit(2,'cm'),gp = gpar(fill = "lightskyblue1")),annotation_label = 'ligand \n abundance')
  h1 <- ComplexHeatmap::Heatmap(object@net[[interaction_name]][sender.names,receiver.names], name = "Commu. \n Prob.", top_annotation = column_ha, right_annotation = row_ha,
                cluster_rows = F,cluster_columns = F,column_names_rot=45,row_names_side = 'left',col = col_map,
                column_title = 'Receiver',row_title = 'Sender',column_title_side = 'bottom',
                heatmap_legend_param = list(color_bar='continuous'),show_row_dend = FALSE)
## plot #2: expression level for ligand related genes
  if(length(object@DB[[interaction_name]]$lig_contributor)>1){
    p1 <- Seurat::VlnPlot(subset(seurat_VISp, idents = sender.names),features = object@DB[[interaction_name]]$lig_contributor,group.by = 'cell_subclass',pt.size = 0,stack = T,flip = T)
    p1$labels$title <- c('Ligand genes')} else {
      p1 <- Seurat::VlnPlot(subset(seurat_VISp, idents = sender.names),features = object@DB[[interaction_name]]$lig_contributor,group.by = 'cell_subclass',pt.size = 0)+ theme(legend.position = 'none')
      p1$labels$title <- paste('Ligand gene: ',p1$labels$title)
    }
  p1$theme$plot.title$hjust <- 0.5
## plot #3: expression level for ligand related genes
  p2 <- Seurat::VlnPlot(subset(seurat_VISp, idents = receiver.names),features = object@DB[[interaction_name]]$receptor_subunit,group.by = 'cell_subclass',pt.size = 0) + theme(legend.position = 'none')
  p2$labels$title <- paste('Target gene: ',p2$labels$title)
  gb_h1 <- grid.grabExpr(draw(h1))
  gb_p1 <- grid.grabExpr(print(p1))
  gb_p2 <- grid.grabExpr(print(p2))
  grid.newpage()
  pushViewport(viewport(x = 0, y = 1, width = width.vector[1]-0.02, height =1, just = c("left", "top"),xscale=c(0,1),yscale = c(0,1)))
  grid.draw(gb_h1);popViewport()
  pushViewport(viewport(x = width.vector[1], y = 1, width = width.vector[2]-0.02, height = 1,just = c("left", "top"),xscale=c(0,1),yscale = c(0,1)))
  grid.draw(gb_p1);popViewport()
  pushViewport(viewport(x = width.vector[1]+width.vector[2], y = 1, width = width.vector[3]-0.02, height = 1,just = c("left", "top"),xscale=c(0,1),yscale = c(0,1)))
  grid.draw(gb_p2);popViewport()
}

####################
#################### For circle plot
#' Circle plot of cell-cell communication network
#'
#' Circle plot of cell-cell communication network (adapted from CellChat https://github.com/sqjin/CellChat)
#'
#' The width of edges represent the strength of the communication.
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
netVisual_circle_neuron <-function(net_ori, color.use = NULL,group=NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                   weight.scale = TRUE, vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1.2,vertex.label.color= "black",
                                   edge.weight.max = NULL, edge.width.max=5, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                   edge.curved=0.2,shape=NULL,layout=in_circle(), margin=0.2, vertex.size = NULL,
                                   arrow.width=1,arrow.size = 0.8){
  color.use = NULL;title.name = NULL; sources.use = NULL; targets.use = NULL; remove.isolate = FALSE; top = 1;
  weight.scale = TRUE; vertex.weight = 1; vertex.weight.max = NULL; vertex.size.max = NULL; vertex.label.cex=1;vertex.label.color= "black";
  edge.weight.max = NULL; edge.width.max=5; alpha.edge = 0.6; label.edge = FALSE;edge.label.color='black';edge.label.cex=0.8;
  edge.curved=0.2;shape='circle';layout=in_circle(); margin=0.2; vertex.size = NULL;
  arrow.width=1;arrow.size = 0.2
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
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }

  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)

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
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


#' Circle plot comparing two cell-cell communication networks (for benchmarking)
#'
#' Circle plot comparing two cell-cell communication networks (for benchmarking) (adapted from CellChat https://github.com/sqjin/CellChat)
#'
#' @param net0 A weighted matrix representing the connections, regarded as the true communication network
#' @param net1 A weighted matrix representing the connections, regarded as the predicted communication network
#' @param group a vector indicating which cell types (cell subclass) belong to which big groups (cell class)
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
netVisual_circle_compare <-function(net0,net1,color.use = NULL,group=NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                    weight.scale = TRUE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1.2,vertex.label.color= "black",
                                    edge.weight.max = NULL, edge.width.max=5, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                    edge.curved=0.2,shape=NULL,layout=in_circle(), margin=0.2, vertex.size = NULL,
                                    arrow.width=1,arrow.size = 0.2){
  # net0 <- real_net; net1 <- predicted_net;color.use = NULL;group=NULL;title.name = NULL; sources.use = net.rownames; targets.use = net.rownames; remove.isolate = FALSE; top = 1;
  # weight.scale = FALSE; vertex.weight = 20; vertex.weight.max = NULL; vertex.size.max = NULL; vertex.label.cex=1.2;vertex.label.color= "black";
  # edge.weight.max = NULL; edge.width.max=5; alpha.edge = 0.6; label.edge = FALSE;edge.label.color='black';edge.label.cex=0.8;
  # edge.curved=0.2;shape='circle';layout=in_circle(); margin=0.2; vertex.size = NULL;
  # arrow.width=1;arrow.size = 0.8
  net <- ((net0+net1)>0)*1
  net.names <- unique(c(rownames(net),colnames(net)))
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
  color.use.label <- c('black',CellChat::scPalette(length(unique(group))))
  if(length(color.use.label)==1){ color.use.label='black'}
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$shape <- c('circle')#, 'square', 'csquare', 'rectangle', 'crectangle', 'vrectangle', 'pie', 'raster','sphere')[group[igraph::V(g)]]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  #igraph::V(g)$frame.color <- border.color.use[igraph::V(g)]
  igraph::V(g)$label.color <- color.use.label[group[names(V(g))]]
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }

  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$lty <- 1
  net <- net[names(V(g)),names(V(g))]
  net1 <- net1[names(V(g)),names(V(g))]
  net0 <- net0[names(V(g)),names(V(g))]
  link_missed <- which(net-net1>0,arr.ind = T)
  link_FP <-  which(net1-net0>0,arr.ind = T)
  if(length(link_missed)>0){
    miss.ind <- match(paste(link_missed[,1],link_missed[,2]),paste(edge.start[,1],edge.start[,2]))
    igraph::E(g)[miss.ind]$color <- grDevices::adjustcolor('black',0.7)
    igraph::E(g)[miss.ind]$width <- 4
  }
  if(length(link_FP)>0){
    fp.ind <- match(paste(link_FP[,1],link_FP[,2]),paste(edge.start[,1],edge.start[,2]))
    igraph::E(g)[fp.ind]$color <- grDevices::adjustcolor('black',0.7)
    igraph::E(g)[fp.ind]$lty <- 6
    igraph::E(g)[fp.ind]$width <- 4
  }
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
  if (!is.null(title.name)) {
    text(0,1,title.name, cex = 1.5)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


#' Chord diagram of cell-cell communication network for single or multiple interaction pairs
#'
#' @param object a NeuronChat object
#' @param interaction_use ligand-target interaction indexes  used for aggregation ('all' means use all interation pairs); for single interaction pair, set interaction_use as the interaction name or index, e.g., 'Glu_Gria2', or simply use CellChat function netVisual_chord_cell_internal(object@net[['Glu_Gria2']], group = group,lab.cex=1.3)
#' @param group a vector indicating which cell types (cell subclass, e.g., L2/3 IT) belong to which big groups (cell class, e.g., Glutamatergic)
#' @param method method used for aggregation; see also function `net_aggregation`
#' @param cut_off threshold used for aggregation; see also function `net_aggregation`
#' @param sender.names sending cell groups used for plot
#' @param receiver.names receiving cell groups used for plot
#' @importFrom CellChat netVisual_chord_cell_internal
#' @return
#' @export
netVisual_chord_neuron <- function(object, method=c('weight','count','weighted_count','weighted_count2','weight_threshold'),cut_off=0.05, interaction_use='all',group=NULL, sender.names=NULL,receiver.names=NULL,lab.cex=1.3){
  method <- match.arg(method)
  if(is.null(sender.names)){sender.names=rownames(object@net[[1]])}
  if(is.null(receiver.names)){receiver.names=colnames(object@net[[1]])}
  if(interaction_use=='all') {
    net_aggregated <- net_aggregation(object@net,method=method,cut_off=cut_off)}
  else {net_aggregated <- net_aggregation(object@net[interaction_use],method=method,cut_off=cut_off)}
  CellChat::netVisual_chord_cell_internal(net_aggregated, group = group,lab.cex=lab.cex)
}
