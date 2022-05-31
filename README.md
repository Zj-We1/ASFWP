由于具有适用范围广泛和算法伸缩性强等优点,k-means长久以来一直是数据聚类领域使用最广泛的算法之一。然而，在如今的大数据时代，大量的噪声数据使得k-means的使用越来越具有局限性，当许多特征不相关或具有一定数量离群点时，k-means的表现不尽人意。针对这些问题，本文提出一种自适应样本-特征双路加权的k-means聚类算法，该算法引入两个不同的正则化项，通过向k-means的解进行退火，在每次迭代中自适应的调整样本和特征权重，有效避免了得到较差的局部最小值，同时保持了经典k-means算法的简单性。在设计的合成数据集中，该算法成功区分了人为添加的噪声特征和样本。此外，通过在真实数据集上的验证以及与其他先进算法的对比，证实了该算法优越的聚类性能。

k-means has long been one of the most widely used algorithms in the field of data clustering due to its wide application range and strong algorithm scalability. However, in today's era of big data, the use of k-means is increasingly limited due to a large amount of noise data. When many features are irrelevant or have a certain number of outliers, k-means performance is unsatisfactory. To solve these problems, this paper proposes an adaptive sample-feature double-path weighted k-means clustering algorithm. The algorithm introduces two different regularization terms and adaptively adjusts the sample and feature weights in each iteration by annealing the k-means solution, which effectively avoids getting poor local minima. Meanwhile, the simplicity of the classical k-means algorithm is maintained. In the designed synthetic data set, the algorithm successfully distinguishes artificially added noise features and noise samples. In addition, the algorithm's superior clustering performance is verified by the validation on real data sets and comparison with other advanced algorithms.

email:20191002792@cug.edu.cn
