# ASFWP
Adaptive sample-feature weighted Power  k-means
一种自适应的样本-特征双加权的改进k-means聚类算法，该算法引入两个不同的正则化项，通过一种MM算法向k-means的解进行退火，在每次迭代中自适应的调整样本和特征权重直至达到收敛条件，并且保持了经典k-means算法O(NKP)的简单性。
An adaptive sample-feature double weighting improved K-means clustering algorithm. The algorithm introduces two different regularization terms and anneals the k-means solution through a MM algorithm. In each iteration, the sample and feature weights are adjusted adaptically until convergence conditions are reached. And keep the simplicity of the classical K-means algorithm O(NKP).
