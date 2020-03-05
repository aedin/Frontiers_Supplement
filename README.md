Demystifying dimensionality reduction in R/Bioconductor 
 
(Based on this template: https://github.com/Bioconductor/BioC2018/blob/master/resources/workshop-syllabus.md)
 
Description:
Principal component analysis (PCA) is a key step in many bioinformatics pipelines.. But what exactly is PCA, anyway, and how does it work? What is the difference between the various implementations of PCA in R? How is PCA output interpreted?
 
In this interactive session lecture + lab session, we will take a deep dive into the various implementations of singular value decomposition (SVD) and principal component analysis (PCA) to clarify the relationship between these methods, and to demonstrate the equivalencies and contrasts between these methods. We will also discuss interpretation of outputs, as well as some common pitfalls and sources of confusion in utilizing these methods.
 
Pre-reqs:
A basic understanding of R syntax would be helpful, but not required. No prior knowledge of PCA necessary.
 
Workshop participation:
We invite audience members to engage with questions and examples from their own workflows. R notebooks will also be available in advance to run code interactively with the workshop.
 
R/Bioconductor packages to be used:
stats (prcomp, princomp, svd), FactoMineR, ade4, irlba, ggplot
 
Time outline:
Set-up + package installation (5 min)
Introduction to matrix factorization and PCA [conceptual] (15 min)
Interactive demonstration of methods (25 min)
Potential pitfalls, interpreting outputs, and how to decide what’s right for your pipeline (15 min)
 
Workshop goals + objectives:
 
Upon completion of this workshop, we expect participants to have gained an understanding of how to apply PCA and other SVD-based methods in research.
 
Learning goals:
Understand how PCA works, the variations of PCA, and how it relates to SVD
Suggest appropriate use cases for these dimensionality reduction techniques
Select appropriate methods for use in bioinformatics pipelines
 
Learning objectives:
 
Describe the similarities and differences between the different implementations of PCA and SVD in R/Bioconductor
Perform PCA/SVD on real data
Creating plots to interpret PCA/SVD outputs, including diagnosis of problems like arch/horseshoe effect


 
