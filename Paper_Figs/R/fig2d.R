library(SingleCellExperiment)
library(ggplot2)
library(ggthemes)
library(scater)
library(mogsa)

####################################################################
# Loading the intersected datasets from step 1
####################################################################

load('intersected_scmix.rda')

####################################################################
# Pre-processing and running CCA
####################################################################

raw_counts_list <- list('dropseq' = t(counts(int_sce_sc_Dropseq_qc)), 'celseq' = t(counts(int_sce_sc_CELseq2_qc)), '10x' = t(counts(int_sce_sc_10x_qc)))

logcounts_list <- list('dropseq' = t(logcounts(int_sce_sc_Dropseq_qc)), 'celseq' = t(logcounts(int_sce_sc_CELseq2_qc)), '10x' = t(logcounts(int_sce_sc_10x_qc)))

raw_moa <- mbpca(raw_counts_list, ncomp = 200, method = 'blockScore', option = 'uniform')

logmoa <- mbpca(logcounts_list, ncomp = 200, method = 'blockScore', option = 'uniform')


####################################################################
# Plots
####################################################################

# For adding colors
batch_vec <- rep(c('Dropseq','Celseq','10X'), times = c(225,274,902))
by_cell_line <- c(colData(int_sce_sc_Dropseq_qc)$cell_line_demuxlet, 
                  colData(int_sce_sc_CELseq2_qc)$cell_line_demuxlet,
                  colData(int_sce_sc_10x_qc)$cell_line_demuxlet)

ggplot_pcs_moa_ellipse <- function(thismoa, plot_title, colorvec, color_title, xpc, ellipse_vec){
  xvar <- as.name(paste('X',xpc,sep = ''))
  yvar <- as.name(paste('X',xpc + 1,sep = ''))
  
  df <- cbind(data.frame(thismoa@loading),ellipse_vec,colorvec)
  ellvar <- as.name(unlist(colnames(df)[dim(df)[2]-1]))
  colvar <- as.name(unlist(colnames(df)[dim(df)[2]]))
  
  xlab = paste('PC ',xpc, sep = '' )
  ylab = paste('PC ',xpc + 1, sep = '' )
  
  ggplot(df,aes(x = eval(xvar), y = eval(yvar), colour = eval(colvar))) + theme_classic() + 
    geom_hline(yintercept=0, color = 'gray') + geom_vline(xintercept=0, color = 'gray') + 
    geom_point() + scale_color_few() + 
    stat_ellipse(aes(x = eval(xvar), y = eval(yvar), group = eval(ellvar)), type= 'norm', linetype = 2) +
    labs(x = xlab, 
         y = ylab, 
         title = plot_title, 
         color = color_title) + 
    theme(axis.text.x = element_text(size = rel(1.7), colour = 'black'), axis.text.y = element_text(size = rel(1.7), colour = 'black'), axis.title = element_text(size = rel(1.8), colour = 'black'))
  
  ggsave(paste(plot_title,'_ellipse',xpc,'.png',sep = ''))
}

ggplot_pcs_moa_ellipse(raw_moa,'CCA on raw counts',batch_vec,'Platform',1,by_cell_line)
ggplot_pcs_moa_ellipse(raw_moa,'CCA on raw counts',batch_vec,'Platform',2,by_cell_line)

ggplot_pcs_moa_ellipse(logmoa,'CCA on log counts',batch_vec,'Platform',1,by_cell_line)
ggplot_pcs_moa_ellipse(logmoa,'CCA on log counts',batch_vec,'Platform',2,by_cell_line)
