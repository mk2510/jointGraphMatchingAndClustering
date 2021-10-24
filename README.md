# Joint Graph Matching and Clustering

![Teaser Image Graphic](https://github.com/mk2510/jointGraphMatchingAndClustering/blob/main/PaperTeaserImageNew.png)

## Abstract

This paper proposes a new algorithm for simultaneous graph matching and clustering. 
For the first time in the literature, these two problems are solved jointly and synergetically without relying on any training data, which brings advantages for identifying similar arbitrary objects in compound 3D scenes and matching them. 
To achieve joint processing, 
For joint reasoning, we first rephrase graph matching as a rigid point set registration problem operating on spectral graph embeddings. 
Consequently, we utilise efficient convex semidefinite program relaxations for aligning points in Hilbert spaces and add coupling constraints to model the mutual dependency and exploit synergies between both tasks. 
We outperform state of the art in challenging cases with non-perfectly matching and noisy graphs, and we show successful applications on real compound scenes with multiple 3D elements. 

## Setup

The code is currently configured to run on Windows systems.
In order to let it run with a Unix system it is necessary to adapt the filepaths in the following files `PMSDP_relaxations/optimisation/interleaving2.m` and `code/clustered_setup/code/optimisation/interleaving2_out.m`
by replacing the \ with a /.

If the precomputed results should shown, you can skip the setup part and head 
straight to the **Quick start** section.
Before running the first experiments a few setup steps are necessary.

1) run make to compile all C++ files, which are necessary for the [FGM](https://github.com/zhfe99/fgm) algorithm to run

2) run 'code/clustered_setup/fgm-master/make' to compile all C++ files, which are necessary for the [FGM](https://github.com/zhfe99/fgm) algorithm to run

2) The following libraries should be added to the main folder:
    
    a) [Yalmip](https://yalmip.github.io)
    
    b) [Mosek](https://www.mosek.com). Here the code is configured and tested for Mosek 9.2

    c) ANN libray - supplied and compiled -

    d) [Graph toolbox](https://de.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph)

This code was tested on Matlab R2019a and excetuted on  Matlab R2019b on Linux. 
The precalculated files are from the Linux version.

## Quick start

The three excetution files for the experiments are in the main directory. The first file is called `cmuHouse.m`. This file performs the experiments on the cmu House dataset.
The second file is called `syntheticData.m`. Executing this file will display all graphs, which show the results of our synthetic data experiments. To display the pre-calculated results from those experiments no library installation is needed.
The third file is called `princetonShapes.m`. This will execute the demonstration of the algorithm executed with shape data of the [Princetion shape dataset](https://shape.cs.princeton.edu/benchmark/).
This execution will generate the helicopter and ship scenes from the paper.

## License

Permission is hereby granted, free of charge, to any person or company obtaining a copy of this software and associated documentation files (the "Software") from the copyright holders to use the Software for any non-commercial purpose. Publication, redistribution and (re)selling of the software, of modifications, extensions, and derivates of it, and of other software containing portions of the licensed Software, are not permitted. The Copyright holder is permitted to publically disclose and advertise the use of the software by any licensee.

Packaging or distributing parts or whole of the provided software (including code, models and data) as is or as part of other software is prohibited. Commercial use of parts or whole of the provided software (including code, models and data) is strictly prohibited. Using the provided software for promotion of a commercial entity or product, or in any other manner which directly or indirectly results in commercial gains is strictly prohibited.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
