# JDC
Some proposed methods for identifying essential proteins have better results by using biological information. Gene expression data is generally used to identify essential proteins. However, gene expression data is prone to fluctuations, which may affect the accuracy of essential protein identification. Therefore, we propose an essential protein identification method to calculate the similarity of "active" and "inactive" state of gene expression in a cluster of the PPI network based on gene expression and the PPI network data. Our experiments show that our method can improve the accuracy in predicting essential proteins.

We propose a new measure, named JDC, based on the PPI network data and gene expression data. The JDC method offers a dynamic threshold method to binarize gene expression data. After that, it combines the degree centrality and Jaccard similarity index to calculate the JDC score for each protein in the PPI network. We respectively perform experiments on Yeast data and E.coli data and evaluate our method by using ROC analysis, modular analysis, jackknife analysis, overlapping analysis, top analysis, and accuracy analysis. The results show that the performance of JDC is better than DC, IC, EC, SC, BC, CC, NC, PeC, and WDC. We compare JDC with both NF-PIN and TS-PIN methods, which predict essential proteins from active PPI networks constructed with dynamic gene expression.

The main ideas behind JDC are as follows: (1) Essential proteins are generally densely connected clusters in the PPI network. (2) Binarizing gene expression data can screen out fluctuations in gene expression profiles. (3) The essentiality of the protein depends on the similarity of "active" and "inactive" state of gene expression in a cluster of the PPI network.

%Data Prepare
data=xlsread('Yeast_topology.xlsx');%topology
data1=xlsread('Yeast_PPI_network.xlsx');%PPI network
data3=xlsread('geneexpress.xlsx');%gene expression
