# chronospace
A statistical exploration of time-calibrated phylogenies and the relative impact of methodological decisions.

Calibrating phylogenies to absolute time is complex, time-consuming and involves taking many methodological decisions that are often difficult to justify. The effect of these choices (gene sampling strategy, type of clock, model of molecular evolution, prior on divergence times, set of fossil calibrations, etc.) are generally explored by running multiple analyses under different settings, and summarizing the sensitivity of results by showing tables with node ages or plotting multiple consensus trees. However, a lot of information is lost in this process, or at least difficult to visualize. For example, it is often unclear whether older/younger ages for a particular node are realized through shrinkage/expansion of the branches leading to it, or from those that descend from it. Similarly, patterns of correlation between ages, which can illustrate ways in which changes in one region of the tree affect others, are lost. Another current limitation is the lack of methods that allow the effect of each of these choices on inferred dates (and thus their relative importance) to be quantified.

chronospace provides a way of visualizing and exploring the variation in time-calibrated topologies obtained through Bayesian dating methods. It also relies on between-group PCAs (bgPCAs) to estimate the impact of different methodological choices on inferred ages, providing a measure akin to an effect size.

## Usage
The script consists of a set of functions that can be easily loaded using ```source()```. Accompanying this script are six posterior distributions of time-calibrated trees obtained using PhyloBayes (Lartillot et al. 2013), found in folder [example_files](https://github.com/mongiardino/chronospace/tree/main/example_files). These were all run using the same contrained topology (necessary to perform these analyses, as the topology is assumed fixed) using three different sets of 100 genes subsampled from a larger phylogenomic dataset based on their level of clock-likeness or phylogenetic signal, or otherwise at random. Each subsample was also run under two models of molecular evolution, the site-homogeneous GTR+G and the site-heterogeneous CAT+GTR+G (Lartillot & Philippe 2004). This setting allows us to evaluate whether inferred ages are more sensitive to the choice of loci or to the choice of model of molecular evolution. This 6 files need to be downloaded and placed within a folder that contains no other file. The working directory needs to be then changed to this folder.

The script contains two different functions. The first of these is ```extract_ages``` which obtains and organizes the dates contained in all tree files present in the working directory. Two further inputs are required: type - a list of vectors, one for each factor being tested (in this case two, loci choice and model of evolution), specifying the group to which the chronograms from each file will be assigned; and sample - the number of topologies to retain from each file. For example:

```R
type <- list(c('clock', 'clock', 'random', 'random', 'signal', 'signal'), 
            c('CATGTR', 'GTR', 'CATGTR', 'GTR', 'CATGTR', 'GTR'))
sample <- 500
all_dates <- extract_ages(type, sample)
```

After obtaining the dates associated with each node across all files, and cataloguing trees depending on the setting of the runs that produced them, the data can be passed on to ```bgPCA_ages```, which uses bgPCA to rotate the multidimensional chornospace to obtain the axes that maximizing the variance explained by each factor. Separate bgPCAs are used for each factor. The script outputs the proportion of total variance explained by the bgPCA axes, and saves a series of plots including a chronospace, a summary of the most sensitive nodes, and a set of trees showing the relative stretching/compression of branches along the bgPCA axes (see Fig. 1 below). For the most sensitive nodes, one can specify to obtain those nodes that vary between analyses by more than a given amount (expressed in millions of years) using the parameter ```amount_of_change```, or the most sensitive X number of nodes using parameter ```chosen_clades```. A tree also needs to be provided for plotting, as well as a number of standard deviations away from the origin of the bgPCA plot that is used to draw the compressed/expanded topologies.

```R
topology <- read.tree('clockCATGTR_ln_sample.datedist')[[1]]
#run selecting the 5 most sensitive nodes and using 1 standard deviation
bgPCA_ages(data_ages = all_dates, tree = topology, chosen_clades = 5, amount_of_change = NA, sdev = 1)
```

## Citations
Lartillot N, Philippe H. 2004. A Bayesian mixture model for across-site heterogeneities in the amino-acid replacement process. Mol. Biol. Evol. 21, 1095–1109.

Lartillot N, Rodrigue N, Stubbs D, Richer J. 2013. PhyloBayes MPI. Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615.
