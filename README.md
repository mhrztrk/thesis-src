# Source code of Thesis work

MARKOV RANDOM FIELD BASED ROAD NETWORK EXTRACTION FROM HIGH RESOLUTION SATELLITE IMAGES

Road Networks play an important role in various applications such as urban and rural planning,infrastructure planning, transportation management, vehicle navigation. Extraction of Roads from Remote Sensed satellite images for updating road database in geographical information systems (GIS) is generally done manually by a human operator. However, manual extraction of roads is time consuming and labor intensive process. In the existing literature, there are a great number of researches published for the purpose of automating the road extraction process. However, automated processes still yield some erroneous and incomplete results and human intervention is still required.

The aim of this research is to propose a framework for road network extraction from high spatial resolution multi-spectral imagery (MSI) to improve the accuracy of road extraction systems. The proposed framework begins with a spectral classification using One-class Support Vector Machines (SVM) and Gaussian Mixture Models (GMM) classifiers. Spectral Classification exploits the spectral signature of road surfaces to classify road pixels. Then, an iterative template matching filter is proposed to refine spectral classification results. K-medians clustering algorithm is employed to detect candidate road centerline points. Final road network formation is achieved by Markov Random Fields. The extracted road network is evaluated against a reference dataset using a set of quality metrics.


