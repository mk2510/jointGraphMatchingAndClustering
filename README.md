# Joint Graph Matching and Clustering

![Teaser Image Graphic](https://github.com/mk2510/jointGraphMatchingAndClustering/blob/main/PaperTeaserImageNew.pdf)

## Abstract

This paper proposes a new algorithm for simultaneous graph matching and clustering. 
For the first time in the literature, these two problems are solved jointly and synergetically without relying on any training data, which brings advantages for identifying similar arbitrary objects in compound 3D scenes and matching them. 
To achieve joint processing, 
For joint reasoning, we first rephrase graph matching as a rigid point set registration problem operating on spectral graph embeddings. 
Consequently, we utilise efficient convex semidefinite program relaxations for aligning points in Hilbert spaces and add coupling constraints to model the mutual dependency and exploit synergies between both tasks. 
We outperform state of the art in challenging cases with non-perfectly matching and noisy graphs, and we show successful applications on real compound scenes with multiple 3D elements. 

## Setup

The code is currently configured to run on Windows systems.
In order to let it run with a Unix system it is necessary to adapt the filepaths in the following files (....) by replacing the \ with a /.

Before running the first experiments a few setup steps are necessary.

1) run make to compile all C++ files, which are necessary for the [FGM](https://github.com/zhfe99/fgm) algorithm to run

2) The following libraries should be added to the main folder:
    
    a) [Yalmip](https://yalmip.github.io)
    
    b) [Mosek](https://www.mosek.com). Here the code is configured and tested for Mosek 9.2

    c) ANN libray - supplied - make sure you compile it

    d) [Graph toolbox](https://de.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph)

This code was tested on Matlab R2019a and excetuted on  Matlab R2019b on Linux. Just from this version are all calculated results.   

## Quick start

In the folder code there are 3 files. The first file is called `cmuHouse.m`. This file performs the experiments on the cmu House dataset.
The second file is called `syntheticData.m`. Executing this file will display all graphs, which show the results of our synthetic data experiments.
The third file is called `princetonShapes.m`. This will execute the experiments, of the helicopter and ship scenes.

## 
