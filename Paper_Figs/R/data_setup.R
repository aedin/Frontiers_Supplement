library(SingleCellExperiment)

##############################################################################
# DATA SET UP AND PREP
##############################################################################

# This dataset can be downloaded from: https://github.com/LuyiTian/sc_mixology/blob/master/data/sincell_with_class.RData

PATH = 'pathname/'

load(paste0(PATH, 'sincell_with_class.RData'))
# 3 objectcs are included in this .RData: *sce_sc_10x_qc*, *sce_sc_CELseq2_qc*, *sce_sc_Dropseq_qc*

# adding in the platform name to the colData
colData(sce_sc_10x_qc)$Method <- rep('10X',dim(colData(sce_sc_10x_qc))[1])
colData(sce_sc_CELseq2_qc)$Method <- rep('CELseq2',dim(colData(sce_sc_CELseq2_qc))[1])
colData(sce_sc_Dropseq_qc)$Method <- rep('Dropseq',dim(colData(sce_sc_Dropseq_qc))[1])

# identifying the intersect of genes in all 3
scmix_genes.use <- intersect(rownames(sce_sc_10x_qc),rownames(sce_sc_CELseq2_qc))
scmix_genes.use <- intersect(rownames(sce_sc_Dropseq_qc),scmix_genes.use)

# subset for the intersection of genes among all 3
int_sce_sc_10x_qc <- sce_sc_10x_qc[scmix_genes.use,]
int_sce_sc_CELseq2_qc <- sce_sc_CELseq2_qc[scmix_genes.use,]
int_sce_sc_Dropseq_qc <- sce_sc_Dropseq_qc[scmix_genes.use,]

# verifying that the genes are all in the same order across the objects
identical(rownames(int_sce_sc_Dropseq_qc),rownames(int_sce_sc_10x_qc))
identical(rownames(int_sce_sc_Dropseq_qc),rownames(int_sce_sc_CELseq2_qc))

save(int_sce_sc_10x_qc, int_sce_sc_CELseq2_qc, int_sce_sc_Dropseq_qc, file = 'intersected_scmix.rda')
