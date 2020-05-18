library(SingleCellExperiment)
library(ggplot2)
library(ggthemes)
library(scater)

####################################################################
# Loading the intersected datasets from step 1
####################################################################

load('intersected_scmix.rda')

####################################################################
# Helper functions
####################################################################

pca_preproc <- function(df){
  return(scale(df, center = TRUE, scale = FALSE))
}

pca_preproc_scale <- function(df){
  return(scale(df, center = TRUE, scale = TRUE))
}

list2mat <- function(matlist){
  concatted <- matlist[[1]]
  for(i in 2:length(matlist)){
    concatted <- rbind(concatted,matlist[[i]])
  }
  return(concatted)
}

get_pct_var_exp_svd <- function(thissvd){
  denom <- dot(thissvd$d, thissvd$d)
  return(thissvd$d^2 / denom)
}

####################################################################
# Pre-processing data then running SVD
####################################################################

raw_counts_list <- list('dropseq' = t(counts(int_sce_sc_Dropseq_qc)), 'celseq' = t(counts(int_sce_sc_CELseq2_qc)), '10x' = t(counts(int_sce_sc_10x_qc)))

pca_preproc_raw_counts_list <- lapply(raw_counts_list, pca_preproc)
pca_preproc_scale_raw_counts_list <- lapply(raw_counts_list, pca_preproc_scale)

logcounts_list <- list('dropseq' = t(logcounts(int_sce_sc_Dropseq_qc)), 'celseq' = t(logcounts(int_sce_sc_CELseq2_qc)), '10x' = t(logcounts(int_sce_sc_10x_qc)))

pca_preproc_logcounts_list <- lapply(logcounts_list, pca_preproc)
pca_preproc_scale_logcounts_list <- lapply(logcounts_list, pca_preproc_scale)

raw_counts_all_svd <- svd(list2mat(raw_counts_list))
pca_preproc_raw_counts_all_svd <- svd(list2mat(pca_preproc_raw_counts_list))
pca_preproc_scale_raw_counts_all_svd <- svd(list2mat(pca_preproc_scale_raw_counts_list))

logcounts_all_svd <- svd(list2mat(logcounts_list))
pca_preproc_logcounts_all_svd <- svd(list2mat(pca_preproc_logcounts_list))
pca_preproc_scale_logcounts_all_svd <- svd(list2mat(pca_preproc_scale_logcounts_list))

####################################################################
# Fig 2B plots
####################################################################

# For adding colors
batch_vec <- rep(c('Dropseq','Celseq','10X'), times = c(225,274,902))
by_cell_line <- c(colData(int_sce_sc_Dropseq_qc)$cell_line_demuxlet, 
                  colData(int_sce_sc_CELseq2_qc)$cell_line_demuxlet,
                  colData(int_sce_sc_10x_qc)$cell_line_demuxlet)



# Function for plotting

ggplot_pcs_svd_ellipse <- function(thissvd, plot_title, colorvec, color_title, xpc, ellipse_vec){
  xvar <- as.name(paste('X',xpc,sep = ''))
  yvar <- as.name(paste('X',xpc + 1,sep = ''))
  
  df <- cbind(data.frame(thissvd$u),ellipse_vec,colorvec)
  ellvar <- as.name(unlist(colnames(df)[dim(df)[2]-1]))
  colvar <- as.name(unlist(colnames(df)[dim(df)[2]]))
  
  pct_var_exp <- get_pct_var_exp_svd(thissvd)
  xlab = paste('PC ',xpc,' (', round(pct_var_exp[xpc]*100),'% of variance explained)', sep = '' )
  ylab = paste('PC ',xpc + 1,' (', round(pct_var_exp[xpc + 1]*100),'% of variance explained)', sep = '' )
  
  ggplot(df,aes(x = eval(xvar), y = eval(yvar), colour = eval(colvar))) + theme_classic() + 
    geom_hline(yintercept=0, color = 'gray') + geom_vline(xintercept=0, color = 'gray') + 
    geom_point() + scale_color_few() + 
    stat_ellipse(aes(x = eval(xvar), y = eval(yvar), group = eval(ellvar)), type= 'norm', linetype = 2) +
    labs(x = xlab, 
         y = ylab, 
         title = plot_title, 
         color = color_title) + 
    theme(axis.text.x = element_text(size = rel(1.7), colour = 'black'), axis.text.y = element_text(size = rel(1.7), colour = 'black'), axis.title = element_text(size = rel(1.8), colour = 'black'))
  ggsave(paste(plot_title,'ellipse',xpc,'.png',sep = ''))
}

# plots

ggplot_pcs_svd_ellipse(raw_counts_all_svd, 'Unprocessed raw counts', batch_vec, 'Platform', 1, by_cell_line)
ggplot_pcs_svd_ellipse(logcounts_all_svd, 'Unprocessed logcounts', batch_vec, 'Platform', 1, by_cell_line)

ggplot_pcs_svd_ellipse(pca_preproc_raw_counts_all_svd, 'PCA preprocessing, raw counts', batch_vec, 'Platform', 1, by_cell_line)
ggplot_pcs_svd_ellipse(pca_preproc_logcounts_all_svd, 'PCA preprocessing, logcounts', batch_vec, 'Platform', 1, by_cell_line) 

ggplot_pcs_svd_ellipse(pca_preproc_scale_raw_counts_all_svd, 'PCA preprocessing (scaled), raw counts', batch_vec, 'Platform', 1, by_cell_line)
ggplot_pcs_svd_ellipse(pca_preproc_scale_logcounts_all_svd, 'PCA preprocessing (scaled), logcounts', batch_vec, 'Platform', 1, by_cell_line) 

####################################################################
# Fig 2C plots
####################################################################

# for color
seq_depth <- c(colData(int_sce_sc_Dropseq_qc)$total_count_per_cell,
               colData(int_sce_sc_CELseq2_qc)$total_count_per_cell,
               colData(int_sce_sc_10x_qc)$total_count_per_cell)

plot_embedding_gradientcol <- function(embedding, xpc = 1, ypc = xpc + 1, plot_title = paste0('PC',xpc,' by PC',ypc), color_vec, color_title, ellipse_vec = NULL, saveplot = TRUE, plotfn = paste(plot_title,xpc, sep = '_'), showplot = TRUE, returngg = FALSE){
  xvar <- as.name(paste('X',xpc,sep = ''))
  yvar <- as.name(paste('X',ypc,sep = ''))
  
  xlab = paste0('PC',xpc)
  ylab = paste0('PC',ypc)
  
  df <- cbind(data.frame(embedding),color_vec)
  colvar <- as.name(unlist(colnames(df)[dim(df)[2]]))
  
  if(!is.null(ellipse_vec)){
    df <- cbind(df,ellipse_vec)
    ellvar <- as.name(unlist(colnames(df)[dim(df)[2]]))
  }
  
  gg_obj <- ggplot(df, aes(x = eval(xvar), y = eval(yvar), color = eval(colvar))) + theme_classic() + 
    geom_hline(yintercept=0, color = 'gray') + geom_vline(xintercept=0, color = 'gray') + 
    geom_point() + scale_color_gradient(low="blue", high="red") +
    labs(x = xlab, 
         y = ylab, 
         title = plot_title, 
         color = color_title) + 
    theme(axis.text.x = element_text(size = rel(1.4), colour = 'black'), axis.text.y = element_text(size = rel(1.4), colour = 'black'), axis.title = element_text(size = rel(1.4), colour = 'black'))
  
  if(!is.null(ellipse_vec)){
    gg_obj <- gg_obj + stat_ellipse(aes(x = eval(xvar), y = eval(yvar), group = eval(ellvar)), type = 'norm', linetype = 2)
  }
  
  if(saveplot){
    ggsave(paste0(plotfn,'.png'))
  }
  
  if(showplot){
    plot(gg_obj)
  }
  
  if(returngg){
    return(gg_obj)
  }
}

plot_embedding_gradientcol(embedding = pca_preproc_logcounts_all_svd$u,plot_title = 'PCA (covariance) on logcounts',color_vec = seq_depth/1000, color_title = 'Sequencing depth\n1000s of reads', ellipse_vec = celltype_vec, saveplot = FALSE, returngg = TRUE)

plot_embedding_gradientcol(embedding = pca_preproc_scale_raw_counts_all_svd$u,plot_title = 'PCA (correlation) on logcounts',color_vec = seq_depth/1000, color_title = 'Sequencing depth\n1000s of reads', ellipse_vec = celltype_vec, saveplot = FALSE, returngg = TRUE)
