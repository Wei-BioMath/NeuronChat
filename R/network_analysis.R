####################
####################
#################### Analysis of Network
####################
####################
#' Compute signaling network similarity for any pair of datasets
#'
#' Compute signaling network similarity for any pair of datasets; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object A merged NeuronChat object
#' @param slot.name the slot name of object that is used to store communication strength matrices, i.e.,  'net'
#' @param type "functional","structural"
#' @param comparison a numerical vector giving the datasets for comparison
#' @param k the number of nearest neighbors
#' @param thresh the fraction (0 to 0.25) of interactions to be trimmed before computing network similarity
#' @importFrom methods slot
#'
#' @return
#' @export
#'
computeNetSimilarityPairwise_Neuron <- function(object, slot.name = 'net', type = c("functional"), comparison = NULL, k = NULL, thresh = NULL) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("Compute signaling network similarity for datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")
  net <- list()
  signalingAll <- c()
  object.net.nameAll <- c()
  # 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    object.net.name <- names(methods::slot(object, slot.name))[comparison[i]]
    object.net.nameAll <- c(object.net.nameAll, object.net.name)
    net[[i]] = simplify2array(object.net)
    signalingAll <- c(signalingAll, paste0(dimnames(net[[i]])[[3]], "--", object.net.name))
    # signalingAll <- c(signalingAll, dimnames(net[[i]])[[3]])
  }
  names(net) <- object.net.nameAll
  net.dim <- sapply(net, dim)[3,]
  nnet <- sum(net.dim)
  position <- cumsum(net.dim); position <- c(0,position)

  if (is.null(k)) {
    if (nnet <= 25) {
      k <- ceiling(sqrt(nnet))
    } else {
      k <- ceiling(sqrt(nnet)) + 1
    }

  }
  if (!is.null(thresh)) {
    for (i in 1:length(net)) {
      neti <- net[[i]]
      neti[neti < quantile(c(neti[neti != 0]), thresh)] <- 0
      net[[i]] <- neti
    }
  }
  if (type == "functional") {
    # compute the functional similarity
    S3 <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
      }
    }

    # define the similarity matrix
    S3[is.na(S3)] <- 0;  diag(S3) <- 1
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        D_signalings[i,j] <- CellChat::computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    S_signalings <- 1-D_signalings
  }
  # smooth the similarity matrix using SNN
  SNN <- CellChat::buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- signalingAll
  colnames(Similarity) <- rownames(Similarity)

  if (!is.list(methods::slot(object, 'net_analysis')$similarity[[type]]$matrix)) {
    methods::slot(object, 'net_analysis')$similarity[[type]]$matrix <- NULL
  }
  # methods::slot(object, slot.name)$similarity[[type]]$matrix <- Similarity
  methods::slot(object, 'net_analysis')$similarity[[type]]$matrix[[comparison.name]] <- Similarity
  return(object)
}


#' Compute signaling network similarity for any pair of signaling networks
#'
#' Compute signaling network similarity for any pair of signaling networks; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store communication strength matrices, i.e.,  'net'
#' @param type "functional"
#' @param k the number of nearest neighbors
#' @param thresh the fraction (0 to 0.25) of interactions to be trimmed before computing network similarity
#' @importFrom methods slot
#' @import CellChat
#' @return
#' @export
#'
#' @examples
computeNetSimilarity_Neuron <- function(object, slot.name = "net", type = c("functional","structural"), k = NULL, thresh = NULL) {
  type <- match.arg(type)
  prob = simplify2array(methods::slot(object, slot.name))
  if (is.null(k)) {
    if (dim(prob)[3] <= 25) {
      k <- ceiling(sqrt(dim(prob)[3]))
    } else {
      k <- ceiling(sqrt(dim(prob)[3])) + 1
    }

  }
  if (!is.null(thresh)) {
    prob[prob < quantile(c(prob[prob != 0]), thresh)] <- 0
  }
  if (type == "functional") {
    # compute the functional similarity
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    S2 <- D_signalings; S3 <- D_signalings;
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
      }
    }
    # define the similarity matrix
    S3[is.na(S3)] <- 0; S3 <- S3 + t(S3); diag(S3) <- 1
    # S_signalings <- S1 *S2
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        D_signalings[i,j] <- CellChat::computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    D_signalings <- D_signalings + t(D_signalings)
    S_signalings <- 1-D_signalings
  }

  # smooth the similarity matrix using SNN
  SNN <- CellChat::buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- dimnames(prob)[[3]]
  colnames(Similarity) <- dimnames(prob)[[3]]

  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
    methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
  }
  methods::slot(object, 'net_analysis')$similarity[[type]]$matrix[[comparison.name]] <- Similarity
  return(object)
}


#' Select the number of the patterns for running `identifyCommunicationPatterns`
#'
#' Select the number of the patterns for running `identifyCommunicationPatterns`; adapted from CellChat https://github.com/sqjin/CellChat
#' We infer the number of patterns based on two metrics that have been implemented in the NMF R package, including Cophenetic and Silhouette. Both metrics measure the stability for a particular number of patterns based on a hierarchical clustering of the consensus matrix. For a range of the number of patterns, a suitable number of patterns is the one at which Cophenetic and Silhouette values begin to drop suddenly.
#'
#' @param object NeuronChat object
#' @param slot.namethe slot name of object that is used to store communication strength matrices, i.e.,  'net'
#' @param pattern "outgoing" or "incoming"
#' @param k.range a range of the number of patterns
#' @param title.name title of plot
#' @param do.facet whether use facet plot showing the two measures
#' @param nrun number of runs when performing NMF
#' @param seed.use seed when performing NMF
#' @importFrom methods slot
# #' @importFrom NMF nmfEstimateRank
#' @import NMF
# #' @importFrom ggplot2 scale_color_brewer
#' @import ggplot2
#' @return a ggplot object
#' @export
#'
#' @examples
selectK_Neuron <- function(object, slot.name = "net", pattern = c("outgoing","incoming"), title.name = NULL, do.facet = TRUE, k.range = seq(2,10), nrun = 30, seed.use = 10) {
  #slot.name = "net";pattern = c("outgoing");title.name = NULL; do.facet = TRUE; k.range = seq(2,10); nrun = 30; seed.use = 10;
  pattern <- match.arg(pattern)
  prob <- simplify2array(methods::slot(object, slot.name))
  if (pattern == "outgoing") {
    data_sender <- apply(prob, c(1,3), sum)
    data_sender = sweep(data_sender, 2L, apply(data_sender, 2, function(x) max(c(x,1e-6), na.rm = TRUE)), '/', check.margin = FALSE)
    data0 = as.matrix(data_sender)
  } else if (pattern == "incoming") {
    data_receiver <- apply(prob, c(2,3), sum)
    data_receiver = sweep(data_receiver, 2L, apply(data_receiver, 2, function(x) max(c(x,1e-6), na.rm = TRUE)), '/', check.margin = FALSE)
    data0 = as.matrix(data_receiver)
  }
  options(warn = -1)
  data <- data0
  data <- data[rowSums(data)!=0,]

  if (is.null(title.name)) {
    title.name <- paste0(pattern, " signaling \n")
    # title.name <- paste0(pattern, " signaling \n (nrun = ", nrun, ", seed = ", seed.use, ")")
  }
  library('NMF')
  res <- NMF::nmfEstimateRank(data, range = k.range, method = 'lee', nrun=nrun, seed = seed.use)
  df1 <- data.frame(k = res$measures$rank, score = res$measures$cophenetic, Measure = "Cophenetic")
  df2 <- data.frame(k = res$measures$rank, score = res$measures$silhouette.consensus, Measure = "Silhouette")
  # df3 <- data.frame(k = res$measures$rank, score = res$measures$dispersion, Measure = "Dispersion")
  df <- rbind(df1, df2)
  #df <- rbind(df1, df2, df3)
  gg <- ggplot(df, aes(x = k, y = score, group = Measure, color = Measure)) + geom_line(size=1) +
    geom_point() +
    theme_classic() + labs(x = 'Number of patterns', y='Measure score') +
    labs(title = title.name) +  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(legend.position = "right") + theme(text = element_text(size = 10)) + scale_x_discrete(limits = (unique(df$k))) +
    scale_color_brewer(palette="Set2") + guides(color=guide_legend("Measure type"))
  if (do.facet) {
    gg <- gg + facet_wrap(~ Measure, scales='free')
  }
  gg
  return(gg)
}

#' Identification of major signals for specific cell groups and general communication patterns
#'
#' Identification of major signals for specific cell groups and general communication patterns; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store communication strength matrices, i.e.,  'net'
#' @param pattern "outgoing" or "incoming"
#' @param k the number of patterns
#' @param k.range a range of the number of patterns
#' @param heatmap.show whether showing heatmap
#' @param color.use the character vector defining the color of each cell group
#' @param color.heatmap a color name in brewer.pal
#' @param title.legend the title of legend in heatmap
#' @param width width of heatmap
#' @param height height of heatmap
#' @param font.size fontsize in heatmap
#' @importFrom methods slot
#' @importFrom NMF nmfEstimateRank nmf
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom stats setNames
#' @importFrom grid grid.grabExpr grid.newpage pushViewport grid.draw unit gpar viewport popViewport
#'
#' @return
#' @export
#'
#' @examples

identifyCommunicationPatterns_Neuron <- function(object, slot.name = "net", pattern = c("outgoing","incoming"), k = NULL, k.range = seq(2,10), heatmap.show = TRUE,
                                                 color.use = NULL, color.heatmap = "Spectral", title.legend = "Contributions",
                                                 width = 4, height = 6, font.size = 8,thresh_quantile=0) {
  # heatmap.show = TRUE;color.use = NULL;color.heatmap = "Spectral";title.legend = "Contributions";
  # width = 4; height = 6; font.size = 8
  pattern <- match.arg(pattern)
  prob <- simplify2array(methods::slot(object, slot.name))
  if (pattern == "outgoing") {
    data_sender <- apply(prob, c(1,3), sum)
    data_sender = sweep(data_sender, 2L, apply(data_sender, 2, function(x) max(c(x,1e-6), na.rm = TRUE)), '/', check.margin = FALSE)
    data0 = as.matrix(data_sender)
  } else if (pattern == "incoming") {
    data_receiver <- apply(prob, c(2,3), sum)
    data_receiver = sweep(data_receiver, 2L, apply(data_receiver, 2, function(x) max(c(x,1e-6), na.rm = TRUE)), '/', check.margin = FALSE)
    data0 = as.matrix(data_receiver)
  }
  options(warn = -1)
  data <- data0
  data <- data[rowSums(data)!=0,colSums(data)!=0]
  data <- data[,colSums(data) >= quantile(colSums(data),thresh_quantile)]
  data <- data[rowSums(data)!=0,]
  if (is.null(k)) {
    stop("Please run the function `selectK` for selecting a suitable k!")
  }

  outs_NMF <- NMF::nmf(data, rank = k, method = 'lee', seed = 'nndsvd')
  W <- scaleMat(outs_NMF@fit@W, 'r1')
  H <- scaleMat(outs_NMF@fit@H, 'c1')
  colnames(W) <- paste0("Pattern ", seq(1,ncol(W))); rownames(H) <- paste0("Pattern ", seq(1,nrow(H)));
  if (heatmap.show) {
    net <- W
    if (is.null(color.use)) {
      color.use <- CellChat::scPalette(length(rownames(net)))
    }
    color.heatmap = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(255)

    df<- data.frame(group = rownames(net)); rownames(df) <- rownames(net)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    row_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),which = "row",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))

    ht1 = Heatmap(net, col = color.heatmap, na_col = "white", name = "Contribution",
                  left_annotation = row_annotation,
                  cluster_rows = T,cluster_columns = F,clustering_method_rows = "average",
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  show_heatmap_legend = F,
                  column_title = "Cell patterns",column_title_gp = gpar(fontsize = 10)
    )


    net <- t(H)

    ht2 = Heatmap(net, col = color.heatmap, na_col = "white", name = "Contribution",
                  cluster_rows = T,cluster_columns = F,clustering_method_rows = "average",
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  column_title = "Communication patterns",column_title_gp = gpar(fontsize = 10),
                  heatmap_legend_param = list(title = title.legend, title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, at = c(round(min(net, na.rm = T), digits = 1), round(max(net, na.rm = T), digits = 1)),
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 6),grid_width = unit(2, "mm"))
    )

    gb_ht1 = grid.grabExpr(draw(ht1))
    gb_ht2 = grid.grabExpr(draw(ht2))
    #grid.newpage()
    pushViewport(viewport(x = 0.1, y = 0.1, width = 0.2, height = 0.5, just = c("left", "bottom")))
    grid.draw(gb_ht1)
    popViewport()

    pushViewport(viewport(x = 0.6, y = 0.1, width = 0.2, height = 0.5, just = c("left", "bottom")))
    grid.draw(gb_ht2)
    popViewport()

  }

  data_W <- as.data.frame(as.table(W)); colnames(data_W) <- c("CellGroup","Pattern","Contribution")
  data_H <- as.data.frame(as.table(H)); colnames(data_H) <- c("Pattern","Signaling","Contribution")

  res.pattern = list("cell" = data_W, "signaling" = data_H)
  methods::slot(object, 'net_analysis')$pattern[[pattern]] <- list(data = data0, pattern = res.pattern)
  return(object)
}


#' 2D visualization of the joint manifold learning of signaling networks from two datasets
#'
#' 2D visualization of the joint manifold learning of signaling networks from two datasets; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store network analysis results, i.e.,  'net_analysis'
#' @param type "functional"
#' @param comparison a numerical vector giving the datasets for comparison. Default are all datasets when object is a merged object
#' @param pathway.labeled a char vector giving the signaling names to show when labeling each point
#' @param top.label the fraction of signaling pathways to label
#' @param pathway.remove a character vector defining the signaling to remove
#' @param pathway.remove.show whether show the removed signaling names
#' @param color.use defining the color for each cell group
#' @param point.shape a numeric vector giving the point shapes. By default point.shape <- c(21, 0, 24, 23, 25, 10, 12), see available shapes at http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @import CellChat
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingPairwise_Neuron <- function(object, slot.name = "net_analysis", type = c("functional","structural"), comparison = NULL, color.use = NULL, point.shape = NULL, pathway.labeled = NULL, top.label = 1, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2.5, dot.alpha = 0.5,
                                               xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,do.label = T, show.legend = T, show.axes = T) {
  # slot.name = "net_analysis";type = c("functional");comparison = NULL;color.use = NULL; point.shape = NULL;pathway.labeled = NULL;top.label = 1;pathway.remove = NULL; pathway.remove.show = TRUE; dot.size = c(2, 6); label.size = 2.5; dot.alpha = 0.5;
  # xlabel = "Dim 1";ylabel = "Dim 2";title = NULL;do.label = T;show.legend = T;show.axes = T
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- names(methods::slot(object, 'net'))[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, 'net')[[comparison[i]]]
    prob[[i]] = simplify2array(object.net)
  }
  names(prob) <-object.names

  if (is.null(point.shape)) {
    point.shape <- 1:15 #c(21, 0, 24, 23, 25, 10, 12,1:9)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    pathway.remove <- sub("--.*", "", pathway.remove)
  }

  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(paste0(dimnames(probi)[[3]],"--",object.names[i]) %in% pathway.remove)
      #  pathway.remove.idx <- which(dimnames(probi)[[3]] %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }
  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"--",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll
  prob_sum <- prob_sum[rownames(Y)]

  group <- sub(".*--", "", names(prob_sum))
  labels = sub("--.*", "", names(prob_sum))

  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  # color dots (light inside color and dark border) based on clustering and no labels
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob.,fill = clusters, colour = clusters, shape = group)) +
    CellChat_theme_opts() +
    theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) +
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) #+ scale_alpha(group, range = c(0.1, 1))
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  gg <- gg + scale_shape_manual(values = point.shape[1:length(prob)])
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = clusters, alpha=group), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))
  }

  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}


#' Zoom into the 2D visualization of the joint manifold learning of signaling networks from two datasets
#'
#' Zoom into the 2D visualization of the joint manifold learning of signaling networks from two datasets; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store network analysis results, i.e.,  'net_analysis'
#' @param type "functional"
#' @param comparison a numerical vector giving the datasets for comparison. Default are all datasets when object is a merged object
#' @param pathway.remove a character vector defining the signaling to remove
#' @param color.use defining the color for each cell group
#' @param nCol number of columns in the plot
#' @param point.shape a numeric vector giving the point shapes. By default point.shape <- c(21, 0, 24, 23, 25, 10, 12), see available shapes at http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @import CellChat
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingPairwiseZoomIn_Neuron <- function(object, slot.name = "net_analysis", type = c("functional","structural"), comparison = NULL, color.use = NULL, nCol = 1, point.shape = NULL, pathway.remove = NULL, dot.size = c(2, 6), label.size = 2.8, dot.alpha = 0.5,
                                                     xlabel = NULL, ylabel = NULL, do.label = T, show.legend = F, show.axes = T) {
  # slot.name = "net_analysis"; type = "functional"; comparison = NULL;color.use = NULL; nCol = 1; point.shape = NULL; pathway.remove = NULL;dot.size = c(2, 6); label.size = 2.8; dot.alpha = 0.5;
  # xlabel = NULL; ylabel = NULL; do.label = T; show.legend = F; show.axes = T;
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- names(methods::slot(object, 'net'))[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, 'net')[[comparison[i]]]
    prob[[i]] = simplify2array(object.net)
  }
  names(prob) <-object.names

  if (is.null(point.shape)) {
    point.shape <- 1:15 #c(21, 0, 24, 23, 25, 10, 12,
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    pathway.remove <- sub("--.*", "", pathway.remove)
  }

  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(paste0(dimnames(probi)[[3]],"--",object.names[i]) %in% pathway.remove)
      # pathway.remove.idx <- which(dimnames(probi)[[3]] %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }

  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"--",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll
  prob_sum <- prob_sum[rownames(Y)]

  group <- sub(".*--", "", names(prob_sum))
  labels = sub("--.*", "", names(prob_sum))

  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))

  # df <- df[group %in% c('VISp','AI'),]
  # clusters <- clusters[group %in% c('VISp','AI')]
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }

  # zoom into each cluster and do labels
  ggAll <- vector("list", length(unique(clusters)))
  for (i in 1:length(unique(clusters))) {
    clusterID = i
    title <- paste0("Cluster ",  clusterID)
    df2 <- df[df$clusters %in% clusterID,]
    gg <- ggplot(data = df2, aes(x, y)) +
      geom_point(aes(size = Commun.Prob., shape = group),fill = alpha(color.use[clusterID], alpha = dot.alpha), colour = alpha(color.use[clusterID], alpha = 1)) +
      CellChat_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
      guides(colour = guide_legend(override.aes = list(size = 3)))+
      labs(title = title, x = xlabel, y = ylabel) +
      scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
      theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
    idx <- match(unique(df2$group), levels(df$group), nomatch = 0)
    gg <- gg + scale_shape_manual(values= point.shape[idx])
    if (do.label) {
      gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels), colour = color.use[clusterID], size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))
    }

    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }

    if (!show.axes) {
      gg <- gg + theme_void()
    }
    ggAll[[i]] <- gg
  }
  gg.combined <- cowplot::plot_grid(plotlist = ggAll, ncol = nCol)

  gg.combined

}

#' Rank signaling networks based on the information flow or the number of interactions
#'
#' Rank signaling networks based on the information flow or the number of interactions; adapted from CellChat https://github.com/sqjin/CellChat
#'
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store communication strength matrices, i.e.,  'net'
#' @param measure "weight" or "count". "weight": comparing the total interaction weights (strength); "count": comparing the number of interactions;
#' @param mode "single","comparison"
#' @param comparison a numerical vector giving the datasets for comparison; a single value means ranking for only one dataset and two values means ranking comparison for two datasets
#' @param color.use defining the color for each cell group
#' @param stacked whether plot the stacked bar plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param signaling a vector giving the signaling pathway to show
#' @param pairLR a vector giving the names of L-R pairs to show (e.g, pairLR = c("IL1A_IL1R1_IL1RAP","IL1B_IL1R1_IL1RAP"))
#' @param do.stat whether do a paired Wilcoxon test to determine whether there is significant difference between two datasets. Default = FALSE
#' @param cutoff.pvalue the cutoff of pvalue when doing Wilcoxon test; Default = 0.05
#' @param tol a tolerance when considering the relative contribution being equal between two datasets. contribution.relative between 1-tol and 1+tol will be considered as equal contribution
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @param do.flip whether flip the x-y axis
#' @param x.angle,y.angle,x.hjust,y.hjust parameters for rotating and spacing axis labels
#' @param axis.gap whetehr making gaps in y-axes
#' @param ylim,segments,tick_width,rel_heights parameters in the function gg.gap when making gaps in y-axes
#' e.g., ylim = c(0, 35), segments = list(c(11, 14),c(16, 28)), tick_width = c(5,2,5), rel_heights = c(0.8,0,0.1,0,0.1)
#' https://tobiasbusch.xyz/an-r-package-for-everything-ep2-gaps
#' @param show.raw whether show the raw information flow. Default = FALSE, showing the scaled information flow to provide compariable data scale; When stacked = TRUE, use raw information flow by default.
#' @param return.data whether return the data.frame consisting of the calculated information flow of each signaling pathway or L-R pair
#' @param x.rotation rotation of x-labels
#' @param title main title of the plot
#' @param bar.w the width of bar plot
#' @param font.size font size

#' @import ggplot2
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
rankNet_Neuron <- function(object, slot.name = "net", measure = c("weight","count"), mode = c("comparison", "single"), comparison = c(1,2), color.use = NULL, stacked = FALSE, sources.use = NULL, targets.use = NULL,  signaling = NULL, pairLR = NULL, do.stat = FALSE, cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, show.raw = FALSE, return.data = FALSE, x.rotation = 90, title = NULL, bar.w = 0.75, font.size = 8,
                       do.flip = TRUE, x.angle = NULL, y.angle = 0, x.hjust = 1,y.hjust = 1,
                       axis.gap = FALSE, ylim = NULL, segments = NULL, tick_width = NULL, rel_heights = c(0.9,0,0.1)) {
  # slot.name = "net";measure = c("weight");mode = c("comparison");comparison = c(3,15);color.use = NULL;stacked = FALSE;sources.use = NULL; targets.use = NULL;  signaling = NULL;pairLR = NULL;do.stat = FALSE; cutoff.pvalue = 0.05;tol = 0.05; thresh = 0.05; show.raw = FALSE;return.data = FALSE; x.rotation = 90;title = NULL; bar.w = 0.75; font.size = 12;
  # do.flip = TRUE; x.angle = NULL; y.angle = 0; x.hjust = 1; y.hjust = 1;
  # axis.gap = FALSE; ylim = NULL; segments = NULL; tick_width = NULL; rel_heights = c(0.9,0,0.1)
  # ylim=c(0,90); axis.gap = TRUE; segments=c(10, 50);tick_width = c(1,10);rel_heights = c(0.9,0,0.1)
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  options(warn = -1)
  object.names <- names(methods::slot(object, slot.name))
  if (measure == "weight") {
    ylabel = "Information flow"
  } else if (measure == "count") {
    ylabel = "Number of links"
  }
  if (mode == "single") {
    object1 <- simplify2array(methods::slot(object, slot.name))
    object.names <- 'Single'
    prob = object1
    #prob[object1 > thresh] <- 0
    if (measure == "count") {
      prob <- 1*(prob > 0)
    }
    if (!is.null(sources.use)) {
      if (is.character(sources.use)) {
        if (all(sources.use %in% dimnames(prob)[[1]])) {
          sources.use <- match(sources.use, dimnames(prob)[[1]])
        } else {
          stop("The input `sources.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), sources.use)
      prob[idx.t, , ] <- 0
    }
    if (!is.null(targets.use)) {
      if (is.character(targets.use)) {
        if (all(targets.use %in% dimnames(prob)[[1]])) {
          targets.use <- match(targets.use, dimnames(prob)[[2]])
        } else {
          stop("The input `targets.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), targets.use)
      prob[ ,idx.t, ] <- 0
    }
    if (sum(prob) == 0) {
      stop("No inferred communications for the input!")
    }

    pSum <- apply(prob, 3, sum)
    pSum.original <- pSum
    if (measure == "weight") {
      # pSum <- -1/log(pSum)
      # pSum[is.na(pSum)] <- 0
      # idx1 <- which(is.infinite(pSum) | pSum < 0)
      # values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
      # position <- sort(pSum.original[idx1], index.return = TRUE)$ix
      # pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    } else if (measure == "count") {
      pSum <- pSum.original
    }

    pair.name <- names(pSum)

    df<- data.frame(name = pair.name, contribution = pSum.original, contribution.scaled = pSum, group = object.names[comparison[1]])
    idx <- with(df, order(df$contribution))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    for (i in 1:length(pair.name)) {
      df.t <- df[df$name == pair.name[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name[i]), ]
      }
    }

    if ((slot.name == "net") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "net") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }

    gg <- ggplot(df, aes(x=name, y=contribution.scaled)) + geom_bar(stat="identity",width = bar.w) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size=10)) +
      xlab("") + ylab(ylabel) + coord_flip()#+
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }

  } else if (mode == "comparison") {
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    for (i in 1:length(comparison)) {
      object.list <- methods::slot(object, slot.name)[[comparison[i]]]
      prob <- simplify2array(object.list)
      #prob[object.list$pval > thresh] <- 0
      if (measure == "count") {
        prob <- 1*(prob > 0)
      }
      prob.list[[i]] <- prob
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% dimnames(prob)[[1]])) {
            sources.use <- match(sources.use, dimnames(prob)[[1]])
          } else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), sources.use)
        prob[idx.t, , ] <- 0
      }
      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% dimnames(prob)[[1]])) {
            targets.use <- match(targets.use, dimnames(prob)[[2]])
          } else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), targets.use)
        prob[ ,idx.t, ] <- 0
      }
      if (sum(prob) == 0) {
        stop("No inferred communications for the input!")
      }
      pSum.original[[i]] <- apply(prob, 3, sum)
      if (measure == "weight") {
        pSum[[i]] <- pSum.original[[i]]
        # pSum[[i]] <- -1/log(pSum.original[[i]])
        # pSum[[i]][is.na(pSum[[i]])] <- 0
        # idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)
        # pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
      } else if (measure == "count") {
        pSum[[i]] <- pSum.original[[i]]
      }
      pair.name[[i]] <- names(pSum.original[[i]])
      object.names.comparison <- c(object.names.comparison, object.names[comparison[i]])
    }
    # if (measure == "weight") {
    #   values.assign <- seq(max(unlist(pSum))*1.1, max(unlist(pSum))*1.5, length.out = length(unlist(idx)))
    #   position <- sort(pSum.original.all, index.return = TRUE)$ix
    #   for (i in 1:length(comparison)) {
    #     if (i == 1) {
    #       pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
    #     } else {
    #       pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i-1]))+1:length(unlist(idx[1:i])), position)]
    #     }
    #   }
    # }



    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
      df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0, group = object.names[comparison[i]], row.names = pair.name.all)
      df[[i]][pair.name[[i]],3] <- pSum[[i]]
      df[[i]][pair.name[[i]],2] <- pSum.original[[i]]
    }


    # contribution.relative <- as.numeric(format(df[[length(comparison)]]$contribution/abs(df[[1]]$contribution), digits=2))
    # #  contribution.relative <- as.numeric(format(df[[length(comparison)]]$contribution.scaled/abs(df[[1]]$contribution.scaled), digits=1))
    # contribution.relative2 <- as.numeric(format(df[[length(comparison)-1]]$contribution/abs(df[[1]]$contribution), digits=2))
    # contribution.relative[is.na(contribution.relative)] <- 0
    # for (i in 1:length(comparison)) {
    #   df[[i]]$contribution.relative <- contribution.relative
    #   df[[i]]$contribution.relative2 <- contribution.relative2
    # }
    # df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    # idx <- with(df[[1]], order(-contribution.relative,  -contribution.relative2, contribution, -contribution.data2))

    contribution.relative <- list()
    for (i in 1:(length(comparison)-1)) {
      contribution.relative[[i]] <- as.numeric(format(df[[length(comparison)-i+1]]$contribution/df[[1]]$contribution, digits=2))
      contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    names(contribution.relative) <- paste0("contribution.relative.", 1:length(contribution.relative))
    for (i in 1:length(comparison)) {
      for (j in 1:length(contribution.relative)) {
        df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
      }
    }
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    if (length(comparison) == 2) {
      idx <- with(df[[1]], order(-contribution.relative.1, contribution, -contribution.data2))
    } else if (length(comparison) == 3) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2,contribution, -contribution.data2))
    } else if (length(comparison) == 4) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, contribution, -contribution.data2))
    } else {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, -contribution.relative.4, contribution, -contribution.data2))
    }



    for (i in 1:length(comparison)) {
      df[[i]] <- df[[i]][idx, ]
      df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL

    df <- do.call(rbind, df)
    df$group <- factor(df$group, levels = object.names.comparison)

    if (is.null(color.use)) {
      color.use =  ggPalette(length(comparison))
    }

    # https://stackoverflow.com/questions/49448497/coord-flip-changes-ordering-of-bars-within-groups-in-grouped-bar-plot
    df$group <- factor(df$group, levels = rev(levels(df$group)))
    color.use <- rev(color.use)

    # perform statistical analysis
    # if (do.stat) {
    #   pvalues <- c()
    #   for (i in 1:length(pair.name.all)) {
    #     df.prob <- data.frame()
    #     for (j in 1:length(comparison)) {
    #       if (pair.name.all[i] %in% pair.name[[j]]) {
    #         df.prob <- rbind(df.prob, data.frame(prob = as.vector(prob.list[[j]][ , , pair.name.all[i]]), group = comparison[j]))
    #       } else {
    #         df.prob <- rbind(df.prob, data.frame(prob = as.vector(matrix(0, nrow = nrow(prob.list[[j]]), ncol = nrow(prob.list[[j]]))), group = comparison[j]))
    #       }
    #
    #     }
    #     df.prob$group <- factor(df.prob$group, levels = comparison)
    #     if (length(comparison) == 2) {
    #       pvalues[i] <- wilcox.test(prob ~ group, data = df.prob)$p.value
    #     } else {
    #       pvalues[i] <- kruskal.test(prob ~ group, data = df.prob)$p.value
    #     }
    #   }
    #   df$pvalues <- pvalues
    # }
    if (do.stat & length(comparison) == 2) {
      for (i in 1:length(pair.name.all)) {
        if (nrow(prob.list[[j]]) != nrow(prob.list[[1]])) {
          stop("Statistical test is not applicable to datasets with different cellular compositions! Please set `do.stat = FALSE`")
        }
        prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * nrow(prob.list[[1]]), ncol = length(comparison))
        for (j in 1:length(comparison)) {
          if (pair.name.all[i] %in% pair.name[[j]]) {
            prob.values[, j] <- as.vector(prob.list[[j]][ , , pair.name.all[i]])
          } else {
            prob.values[, j] <- NA
          }
        }
        prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE]
        if (nrow(prob.values) >3 & sum(is.na(prob.values)) == 0) {
          pvalues <- wilcox.test(prob.values[ ,1], prob.values[ ,2], paired = TRUE)$p.value
        } else {
          pvalues <- 0
        }
        pvalues[is.na(pvalues)] <- 0
        df$pvalues[df$name == pair.name.all[i]] <- pvalues
      }
    }


    if (length(comparison) == 2) {
      if (do.stat) {
        colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2], ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
      } else {
        colors.text <- ifelse(df$contribution.relative < 1-tol, color.use[2], ifelse(df$contribution.relative > 1+tol, color.use[1], "black"))
      }
    } else {
      message("The text on the y-axis will not be colored for the number of compared datasets larger than 3!")
      colors.text = NULL
    }

    for (i in 1:length(pair.name.all)) {
      df.t <- df[df$name == pair.name.all[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name.all[i]), ]
      }
    }

    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }

    if (stacked) {
      gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w, position ="fill") # +
      # xlab("") + ylab("Relative information flow") #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      #  scale_y_discrete(breaks=c("0","0.5","1")) +
      if (measure == "weight") {
        gg <- gg + xlab("") + ylab("Relative information flow")
      } else if (measure == "count") {
        gg <- gg + xlab("") + ylab("Relative number of interactions")
      }

      gg <- gg + geom_hline(yintercept = 0.5, linetype="dashed", color = "grey50", size=0.5)
    } else {
      if (show.raw) {
        gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      } else {
        gg <- ggplot(df, aes(x=name, y=contribution.scaled, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      }

      if (axis.gap) {
        gg <- gg + theme_bw() + theme(panel.grid = element_blank())
        gg.gap::gg.gap(gg,
                       ylim = ylim,
                       segments = segments,
                       tick_width = tick_width,
                       rel_heights = rel_heights)
      }
    }
    gg <- gg +  CellChat_theme_opts() + theme_classic()
    if (do.flip) {
      gg <- gg + coord_flip() + theme(axis.text.y = element_text(colour = colors.text))
      if (is.null(x.angle)) {
        x.angle = 0
      }

    } else {
      if (is.null(x.angle)) {
        x.angle = 45
      }
      gg <- gg + scale_x_discrete(limits = rev) + theme(axis.text.x = element_text(colour = rev(colors.text)))

    }

    gg <- gg + theme(axis.text=element_text(size=font.size), axis.title.y = element_text(size=font.size))
    gg <- gg + scale_fill_manual(name = "", values = color.use)
    gg <- gg + guides(fill = guide_legend(reverse = TRUE))
    gg <- gg + theme(axis.text.x = element_text(angle = x.angle, hjust=x.hjust),
                     axis.text.y = element_text(angle = y.angle, hjust=y.hjust))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }
  }

  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}

#' Comparing the number of inferred communication links between different datasets
#'
#' Comparing the number of inferred communication links between different datasets; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object A merged NeuronChat object
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param color.use defining the color for each group of datasets
#' @param group a vector giving the groups of different datasets to define colors of the bar plot. Default: only one group and a single color
#' @param group.levels the factor level in the defined group
#' @param group.facet Name of one metadata column defining faceting groups
#' @param group.facet.levels the factor level in the defined group.facet
#' @param n.row Number of rows in facet_grid()
#' @param color.alpha transparency
#' @param legend.title legend title
#' @param width bar width
#' @param title.name main title of the plot
#' @param digits integer indicating the number of decimal places (round) to be used when `measure` is `weight`.
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param remove.xtick whether remove xtick
#' @param size.text font size of the text
#' @param show.legend whether show the legend
#' @param x.lab.rot,angle.x,vjust.x,hjust.x adjusting parameters if rotating xtick.labels when x.lab.rot = TRUE
#' @import ggplot2
#' @return A ggplot object
#' @export
#'
compareInteractions_Neuron <- function(object, measure = c("count", "weight"), color.use = NULL, group = NULL, comparison=c(1,2),group.levels = NULL, group.facet = NULL, group.facet.levels = NULL, n.row = 1, color.alpha = 1, legend.title = NULL, width=0.6, title.name = NULL, digits = 3,
                                   xlabel = NULL, ylabel = NULL, remove.xtick = FALSE,
                                   show.legend = TRUE, x.lab.rot = FALSE, angle.x = 45, vjust.x = NULL, hjust.x = 1, size.text = 10) {
  measure <- match.arg(measure)
  if (measure == "count") {
    df <- as.data.frame(sapply(object@net[comparison], function(x) sum(simplify2array(x)>0)))
    if (is.null(ylabel)) {
      ylabel = "Number of inferred interactions"
    }
  } else if (measure == "weight") {
    df <- as.data.frame(sapply(object@net[comparison], function(x) sum(simplify2array(x))))
    df[,1] <- round(df[,1],digits)
    if (is.null(ylabel)) {
      ylabel = "Interaction strength"
    }
  }
  colnames(df) <- "count"

  df$dataset <- names(object@net)[comparison]
  if (is.null(group)) {
    group <- 1
  }
  df$group <- group
  df$dataset <- factor(df$dataset, levels = names(object@net))
  if (is.null(group.levels)) {
    df$group <- factor(df$group)
  } else {
    df$group <- factor(df$group, levels = group.levels)
  }

  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(group)))
  }
  #   theme_classic() #+ scale_x_discrete(limits = (levels(df$x)))
  if (!is.null(group.facet)) {
    if (all(group.facet %in% colnames(df))) {
      gg <- ggplot(df, aes(x=dataset, y=count, fill = group)) +
        geom_bar(stat="identity", width=width, position=position_dodge())
      gg <- gg + facet_wrap(group.facet, nrow = n.row)
    } else {
      df$group.facet <- group.facet
      if (is.null(group.facet.levels)) {
        df$group.facet <- factor(df$group.facet)
      } else {
        df$group.facet <- factor(df$group.facet, levels = group.facet.levels)
      }
      gg <- ggplot(df, aes(x=dataset, y=count, fill = group)) +
        geom_bar(stat="identity", width=width, position=position_dodge())
      gg <- gg + facet_wrap(~group.facet, nrow = n.row)
    }
  } else {
    gg <- ggplot(df, aes(x=dataset, y=count, fill = group)) +
      geom_bar(stat="identity", width=width, position=position_dodge())
  }
  gg <- gg + geom_text(aes(label=count), vjust=-0.3, size=3, position = position_dodge(0.9))
  gg <- gg + ylab(ylabel) + xlab(xlabel) + theme_classic() +
    labs(title = title.name) +  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = size.text), axis.text = element_text(colour="black"))
  gg <- gg + scale_fill_manual(values = alpha(color.use, alpha = color.alpha), drop = FALSE)
  #  gg <- gg + scale_color_manual(values = alpha(color.use, alpha = 1), drop = FALSE) + guides(colour = FALSE)
  if (remove.xtick) {
    gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  if (is.null(legend.title)) {
    gg <- gg + theme(legend.title = element_blank())
  } else {
    gg <- gg + guides(fill=guide_legend(legend.title))
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x, size=size.text))
  }
  gg
  return(gg)
}


#' River plot showing the associations of latent patterns with cell groups and ligand-receptor pairs or signaling pathways
#'
#' River plot showing the associations of latent patterns with cell groups and ligand-receptor pairs or signaling pathways; adapted from CellChat https://github.com/sqjin/CellChat
#' River (alluvial) plot shows the correspondence between the inferred latent patterns and cell groups as well as ligand-receptor pairs or signaling pathways.
#'
#' The thickness of the flow indicates the contribution of the cell group or signaling pathway to each latent pattern. The height of each pattern is proportional to the number of its associated cell groups or signaling pathways.
#'
#' Outgoing patterns reveal how the sender cells coordinate with each other as well as how they coordinate with certain signaling pathways to drive communication.
#'
#' Incoming patterns show how the target cells coordinate with each other as well as how they coordinate with certain signaling pathways to respond to incoming signaling.
#'
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store communication strength matrices, i.e.,  'net'
#' @param pattern "outgoing" or "incoming"
#' @param cutoff the threshold for filtering out weak links
#' @param sources.use a vector giving the index or the name of source cell groups of interest
#' @param targets.use a vector giving the index or the name of target cell groups of interest
#' @param signaling a character vector giving the name of signaling pathways of interest
#' @param color.use the character vector defining the color of each cell group
#' @param color.use.pattern the character vector defining the color of each pattern
#' @param color.use.signaling the character vector defining the color of each signaling
#' @param do.order whether reorder the cell groups or signaling according to their similarity
#' @param main.title the title of plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @importFrom methods slot
#' @importFrom stats cutree dist hclust
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @import ggalluvial
# #' @importFrom ggalluvial geom_stratum geom_flow to_lodes_form
#' @importFrom ggplot2 geom_text scale_x_discrete scale_fill_manual theme ggtitle
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @return
#' @export
#'
#' @examples
netAnalysis_river_Neuron <- function(object, slot.name = "net", pattern = c("outgoing","incoming"), cutoff.1 = 0.25, cutoff.2 = 0.5,top.n1=1e4,top.n2=1e4,
                                 sources.use = NULL, targets.use = NULL, signaling = NULL,
                                 color.use = NULL, color.use.pattern = NULL, color.use.signaling = "grey50",
                                 do.order = FALSE, main.title = NULL,
                                 font.size = 2.5, font.size.title = 12){
  # object <- VISp[[1]]; slot.name = "net"; pattern = c("outgoing");  cutoff.1 = 0.25; cutoff.2 = 0.5; top.n1=40; top.n2=75;
  # sources.use = NULL; targets.use = NULL; signaling = NULL;
  # color.use = NULL; color.use.pattern = NULL; color.use.signaling = "grey50";
  # do.order = FALSE; main.title = NULL;
  # font.size = 2.5; font.size.title = 12
  message("Please make sure you have load `library(ggalluvial)` when running this function")
  requireNamespace("ggalluvial")
  #  suppressMessages(require(ggalluvial))
  res.pattern <- methods::slot(object, 'net_analysis')$pattern[[pattern]]
  data1 = res.pattern$pattern$cell
  data2 = res.pattern$pattern$signaling
  if (is.null(color.use.pattern)) {
    nPatterns <- length(unique(data1$Pattern))
    if (pattern == "outgoing") {
      color.use.pattern = ggPalette(nPatterns*2)[seq(1,nPatterns*2, by = 2)]
    } else if (pattern == "incoming") {
      color.use.pattern = ggPalette(nPatterns*2)[seq(2,nPatterns*2, by = 2)]
    }
  }
  if (is.null(main.title)) {
    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
    }
  }

  if (is.null(data2)) {
    data1$Contribution[data1$Contribution <= sort(data1$Contribution,decreasing = T)[min(dim(data1)[1],top.n1)]] <- 0
    data1$Contribution[data1$Contribution < cutoff.1] <- 0
    plot.data <- data1
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    if (is.null(color.use)) {
      color.use <- CellChat::scPalette(nCellGroup)
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }

    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
      color.use <- color.use[order.name]
    }
    color.use.all <- c(color.use, color.use.pattern)
    gg <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                    stratum = stratum, alluvium = connection,
                                    fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "backward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) +
      scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10))+
      ggtitle(main.title)

  } else {
    data1$Contribution[data1$Contribution <= sort(data1$Contribution,decreasing = T)[min(dim(data1)[1],top.n1)]] <- 0
    data1$Contribution[data1$Contribution < cutoff.1] <- 0
    plot.data <- data1
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    cells.level = levels(object@idents)
    if (is.null(color.use)) {
      color.use <- CellChat::scPalette(length(cells.level))[cells.level %in% unique(plot.data$CellGroup)]
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      plot.data <- subset(plot.data, CellGroup %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      plot.data <- subset(plot.data, CellGroup %in% targets.use)
    }
    ## connect cell groups with patterns
    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
      color.use <- color.use[order.name]
    }
    color.use.all <- c(color.use, color.use.pattern)
    StatStratum <- ggalluvial::StatStratum
    gg1 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                     stratum = stratum, alluvium = connection,
                                     fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "backward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) +
      scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    ## connect patterns with signaling
    data2$Contribution[data2$Contribution <= sort(data2$Contribution,decreasing = T)[min(dim(data2)[1],top.n2)]] <- 0
    data2$Contribution[data2$Contribution < cutoff.2] <- 0
    plot.data <- data2
    nPatterns<-length(unique(plot.data$Pattern))
    nSignaling<-length(unique(plot.data$Signaling))
    if (length(color.use.signaling) == 1) {
      color.use.all <- c(color.use.pattern, rep(color.use.signaling, nSignaling))
    } else {
      color.use.all <- c(color.use.pattern, color.use.signaling)
    }

    if (!is.null(signaling)) {
      plot.data <- plot.data[plot.data$Signaling %in% signaling, ]
    }

    plot.data.long <- ggalluvial::to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["Signaling"]], plot.data[["Pattern"]]), sum)
      mat[is.na(mat)] <- 0; mat <- mat[-which(rowSums(mat) == 0), ]
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(colnames(mat),names(cluster)[order.name]))
    }

    gg2 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("Pattern", "Signaling")),y= Contribution,
                                     stratum = stratum, alluvium = connection,
                                     fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "forward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) + # 2.5
      scale_x_discrete(limits = c(),  labels=c("Patterns", "Signaling")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size= 10))+
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    ## connect cell groups with signaling
    # data1 = data1[data1$Contribution > 0,]
    # data2 = data2[data2$Contribution > 0,]

    # data3 = merge(data1, data2, by.x="Pattern", by.y="Pattern")
    # data3$Contribution <- data3$Contribution.x * data3$Contribution.y
    # data3 <- data3[,colnames(data3) %in% c("CellGroup","Signaling","Contribution")]

    # plot.data <- data3
    # nSignaling<-length(unique(plot.data$Signaling))
    # nCellGroup<-length(unique(plot.data$CellGroup))
    #
    # if (length(color.use.signaling) == 1) {
    #   color.use.signaling <- rep(color.use.signaling, nSignaling)
    # }
    #
    #
    # ## connect cell groups with patterns
    # plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    # if (do.order) {
    #   mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Signaling"]]), sum)
    #   d <- dist(as.matrix(mat))
    #   hc <- hclust(d, "ave")
    #   k <- length(unique(grep("Signaling", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
    #   cluster <- hc %>% cutree(k)
    #   order.name <- order(cluster)
    #   plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
    #   color.use <- color.use[order.name]
    # }
    # color.use.all <- c(color.use, color.use.signaling)

    # gg3 <- ggplot(plot.data.long, aes(x = factor(x, levels = c("CellGroup", "Signaling")),y=Contribution,
    #                                  stratum = stratum, alluvium = connection,
    #                                  fill = stratum, label = stratum)) +
    #   geom_flow(width = 1/3,aes.flow = "forward") +
    #   geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
    #   geom_text(stat = "stratum", size = 2.5) +
    #   scale_x_discrete(limits = c(),  labels=c("Cell groups", "Signaling")) +
    #   scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
    #   theme_bw()+
    #   theme(legend.position = "none",
    #         axis.title = element_blank(),
    #         axis.text.y= element_blank(),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor  = element_blank(),
    #         panel.border = element_blank(),
    #         axis.ticks = element_blank(),axis.text=element_text(size=10)) +
    #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


    gg <- cowplot::plot_grid(gg1, gg2,align = "h", nrow = 1)
    title <- cowplot::ggdraw() + cowplot::draw_label(main.title,size = font.size.title)
    gg <- cowplot::plot_grid(title, gg, ncol=1, rel_heights=c(0.1, 1))
  }
  return(gg)
}


#' 2D visualization of the learned manifold of signaling networks
#'
#' 2D visualization of the learned manifold of signaling networks; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store network analysis results, i.e.,  'net_analysis'
#' @param type "functional"
#' @param pathway.labeled a char vector giving the signaling names to show when labeling each point
#' @param top.label the fraction of signaling pathways to label
#' @param pathway.remove a character vector defining the signaling to remove
#' @param pathway.remove.show whether show the removed signaling names
#' @param color.use defining the color for each cell group
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embedding_Neuron <- function(object, slot.name = "net_analysis", type = c("functional","structural"), color.use = NULL, pathway.labeled = NULL, top.label = 1, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2, dot.alpha = 0.5,
                                   xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,
                                   font.size = 10, font.size.title = 12, do.label = T, show.legend = T, show.axes = T) {
  type <- match.arg(type)
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  Groups <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  prob <- simplify2array(methods::slot(object, 'net'))
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
  }

  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(dimnames(prob)[[3]] %in% pathway.remove)
    prob <- prob[ , , -pathway.remove.idx]
  }

  prob_sum <- apply(prob, 3, sum)
  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum), labels = as.character(unlist(dimnames(prob)[3])), Groups = as.factor(Groups))
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(Groups)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob.,fill = Groups, colour = Groups), shape = 21) +
    CellChat_theme_opts() +
    theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size= font.size.title, face="plain"))+
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE)
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  if (do.label) {
    if (is.null(pathway.labeled)) {
      if (top.label < 1) {
        if (length(comparison) == 2) {
          g.t <- rankSimilarity(object, slot.name = slot.name, type = type, comparison1 = comparison)
          pathway.labeled <- as.character(g.t$data$name[(nrow(g.t$data)-ceiling(top.label * nrow(g.t$data))+1):nrow(g.t$data) ])
          data.label <- df[df$labels %in% pathway.labeled, , drop = FALSE]
        }
      } else {
        data.label <- df
      }

    } else {
      data.label <- df[df$labels %in% pathway.labeled, , drop = FALSE]
    }
    gg <- gg + ggrepel::geom_text_repel(data = data.label, mapping = aes(label = labels, colour = Groups), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))

    # gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = Groups), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5)
  }

  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}


#' Zoom into the 2D visualization of the learned manifold learning of the signaling networks
#'
#' Zoom into the 2D visualization of the learned manifold learning of the signaling networks; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store network analysis results, i.e.,  'net_analysis'
#' @param type "functional"
#' @param pathway.remove a character vector defining the signaling to remove
#' @param color.use defining the color for each cell group
#' @param nCol the number of columns of the plot
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingZoomIn_Neuron <- function(object, slot.name = "net_analysis", type = c("functional","structural"), color.use = NULL, pathway.remove = NULL,  nCol = 1, dot.size = c(2, 6), label.size = 2.8, dot.alpha = 0.5,
                                         xlabel = NULL, ylabel = NULL, do.label = T, show.legend = F, show.axes = T) {
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  prob <- simplify2array(methods::slot(object, 'net'))
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
  }

  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(dimnames(prob)[[3]] %in% pathway.remove)
    prob <- prob[ , , -pathway.remove.idx]
  }

  prob_sum <- apply(prob, 3, sum)
  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum), labels = as.character(unlist(dimnames(prob)[3])), clusters = as.factor(clusters))

  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }

  # zoom into each cluster and do labels
  ggAll <- vector("list", length(unique(clusters)))
  for (i in 1:length(unique(clusters))) {
    clusterID = i
    title <- paste0("Group ",  clusterID)
    df2 <- df[df$clusters %in% clusterID,]
    gg <- ggplot(data = df2, aes(x, y)) +
      geom_point(aes(size = Commun.Prob.), shape = 21, colour = alpha(color.use[clusterID], alpha = 1), fill = alpha(color.use[clusterID], alpha = dot.alpha)) +
      CellChat_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
      labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size=12))+
      scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
      theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
    if (do.label) {
      gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels), colour = color.use[clusterID], size = label.size, segment.size = 0.2, segment.alpha = 0.5)
    }

    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }

    if (!show.axes) {
      gg <- gg + theme_void()
    }
    ggAll[[i]] <- gg
  }
  gg.combined <- cowplot::plot_grid(plotlist = ggAll, ncol = nCol)

  gg.combined

}

