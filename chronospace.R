#Chronospace: a tool to visualize and quantify the effect of methodological
#decisions on time-calibrated inference
#Writen by Nicolas Mongiardino Koch 05/2021

#This script requires as input a series of tree files in newick format,
#representing separate Bayesian time-calibrated inferences of the same dataset
#(can be either different runs of the same analysis or separate analyses). The
#topology is assumed constrained, such that chronograms differ only in the
#branch lengths (and therefore, in their inferred node ages). A list with one or
#more character vectors is needed to assign labels to these files.

#A between-group PCA (bgPCA) is then used to summarize major axes of change, and
#quantify how much they contribute to the overall variance in node ages across
#analyses. Several plots are produced, including a chronospace plot showing the
#distribution of chronograms in multivariate space, posterior distributions of
#nodes that change the most between different analyses, and an overall summary
#of the effect of different choices on branch lengths.

#Parameters needing input are marked with 'INPUT' and are all in the first
#section below entitled 'Parameters'

#More details can be found in the following publications:
...

#Install and load packages------------------------------------------------------
packages <- c('ape','phangorn','phytools',
              'MASS','Morpho',
              'stringr','dplyr','ggplot2','ggpubr','plotrix','RColorBrewer')
new_packages <- packages[!packages %in% installed.packages()[,'Package']]
if(length(new_packages)) { install.packages(new_packages) }

rm(list=ls())
library(ape)
library(phangorn)
library(phytools)

library(MASS)
library(Morpho)

library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotrix)
library(RColorBrewer)
`%not in%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))

#add lines between geological periods
add_time_lines = function() {
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  t.max <- max(obj$xx)
  ii <- which(times <= t.max)
  times <- times[ii]
  
  y <- c(rep(0, 2), rep(par()$usr[4], 2))
  if (obj$direction == "rightwards") {
    old.usr <- par()$usr
    h <- max(obj$xx)
    new.xlim <- c(h - par()$usr[1], h - par()$usr[2])
    par(usr = c(new.xlim, old.usr[3:4]))
  }
  for (i in 1:length(times)) {
    lines(x = rep(times[i], 2), y = c(0, par()$usr[4]), 
          lty = "dotted", col = "grey")
  }
}

times = c(2.588, 23.03, 66, 145, 201.3, 252.17, 298.9, 
          252.17, 358.9, 419.2, 443.8, 485.4, 541, 4600)

#Set up-------------------------------------------------------------------
setwd('C:/Users/mongi/Documents/chronospace/example_files')
sample = 500
type = list(c('clock','random','signal'), 
            c('CATGTR','GTR','CATGTR','GTR','CATGTR','GTR'))

#load files---------------------------------------------------------------------
extract_ages = function(type, sample) {
  files = list.files()
  
  if(!all(sapply(type, length) == max(sapply(type, length)))) {
    to_stretch = which(sapply(type, length) != max(sapply(type, length)))
    if(length(to_stretch) == 1) {
      type[[to_stretch]] = rep(type[[to_stretch]], 
                               each = max(sapply(type, length))/sapply(type, length)[to_stretch])
    } else {
      for(i in 1:length(to_stretch)) {
        type[[to_stretch[i]]] = rep(type[[to_stretch[i]]], 
                                    each = max(sapply(type, length))/sapply(type, length)[to_stretch[i]])
      }
    }
  }
  
  #check that it set-up correctly
  cat("Check that labels are assigned correctly to the input files.\n", 
      "If there is an error, modify the order of factors in 'type'\n", 
      "or the name of input files for the two to match.\n\n", sep = '')
  
  for(i in 1:length(files)) {
    to_print = paste0('file = ', files[i], ' | type = ', 
                      paste0(unname(as.vector(as.data.frame(type)[i,])), 
                             collapse = ' - '))
    cat(to_print, '\n')
  }
  
  for(i in 1:length(files)) {
    trees = read.tree(paste0(getwd(), '/', files[i]))
    if(!is.na(sample)) {
      trees = trees[sample(1:length(trees), sample)]
    }
    
    if(i == 1) {
      all_trees = trees
    } else {
      all_trees = c(all_trees, trees)
    }
  }
  
  tree = all_trees[[1]]
  clades = list()
  for(i in 1:tree$Nnode) {
    clades[i] = list(tree$tip.label[unlist(Descendants(tree, length(tree$tip.label)+i, type = 'tips'))])
  }
  
  ages = matrix(0, ncol = length(clades), nrow = length(all_trees))
  
  root = sapply(all_trees, function(x) max(nodeHeights(x)))
  ages[,1] = root
  
  for(i in 2:ncol(ages)) {
    node = sapply(all_trees, findMRCA, clades[[i]])
    if(length(unique(node)) == 1) {
      node = node[1]
      pos = which(all_trees[[1]]$edge[,2] == node)
      ages[,i] = root - sapply(all_trees, function(x) nodeHeights(x)[pos,2])
    } else {
      for(j in 1:length(all_trees)) {
        ages[j,i] = root[j] - nodeHeights(all_trees[[j]])[which(all_trees[[j]]$edge[,2] == node[j]),2]
      }
    }
  }
  
  data_ages = data.frame(ages)
  colnames(data_ages) = paste0('clade_', 1:ncol(data_ages))
  
  types_of_runs = data.frame(matrix(NA, nrow = nrow(ages), ncol = length(type)))
  for(i in 1:length(type)) {
    types_of_runs[,i] = rep(type[[i]], each = nrow(ages)/length(type[[i]]))
  }
  
  colnames(types_of_runs) = paste('factor', LETTERS[1:ncol(types_of_runs)], 
                                  sep = '_')
  
  data_ages = cbind(data_ages, types_of_runs)
  data_ages = data_ages %>% mutate_if(sapply(data_ages, is.character), as.factor)

  return(data_ages)
}

#Run analyses-------------------------------------------------------------------
#LDA
bgPCA_ages = function(data_ages, chosen_clades = NA, amount_of_change = 20, sdev = 1) {
  
  #split into ages and factors
  ages = data_ages[,which(grepl('clade', colnames(data_ages)))]
  groups = data_ages[,which(grepl('factor', colnames(data_ages)))]
  
  #perform MANOVA
  if(ncol(groups) == 1) manova = manova(as.matrix(ages) ~ groups)
  if(ncol(groups) == 2) manova = manova(as.matrix(ages) ~ groups[,1] + groups[,2])
  if(ncol(groups) == 3) manova = manova(as.matrix(ages) ~ groups[,1] + groups[,2] + groups[,3])
  
  significant_axes = as.numeric(which(summary.manova(manova)[4]$stats[,6] < 0.5))
  
  #report MANOVA results
  if(ncol(groups) == 1) {
    if(length(significant_axes) > 0) {
      cat('The factor tested significantly affects divergence times')
    } else {
      cat('The factor tested does not significantly affect divergence times')
    }
  } else {
    if(length(significant_axes) == 1) {
      cat('Only', colnames(groups)[significant_axes], 'significantly affects divergence times')
    } else {
      if(length(significant_axes > 1)) {
        cat(paste0(colnames(groups)[significant_axes], collapse = ' and '), 
            'significantly affect divergence times')
      } else {
        cat('None of the factors tested significantly affect divergence times')
      }
    }
  }
  
  #perform bgPCA
  for(i in 1:length(significant_axes)) {
    bgPCA = groupPCA(ages, groups[,significant_axes[i]])
    cat(paste0('Proportion of variance explained by ', 
               paste0('factor_', LETTERS[significant_axes[i]]), ' = ', 
               round((bgPCA$combinedVar[max(which(grepl('bg', rownames(bgPCA$combinedVar)))),3]*100), 2), 
               '%', '\n'))
    scores = bgPCA$Scores
    
    #set axes to either 1 or 2 depending on number of groups
    if(ncol(scores) >= 2) {
      num_functions = 2
    } else {
      num_functions = 1
    }
    
    #plot chronospace
    if(num_functions == 1) {
      to_plot = data.frame(coordinates = scores, groups = unname(groups[significant_axes[i]]))
      colors_random = brewer.pal(7, 'Set3')[-2][sample(1:6, 2)]
      
      plot(ggplot(to_plot, aes(x = coordinates, fill = groups)) + 
             geom_histogram(alpha = 0.5, position = 'identity', bins = 30) + theme_bw() + 
             scale_fill_manual(values = colors_random) + ylab('Count') + 
             theme(legend.title = element_blank()) +
             xlab(paste0('bgPCA axis 1 (', round((bgPCA$combinedVar[1,2]*100), 2), '% of variance)')))
      
    } else {
      to_plot = data.frame(coordinates = scores[,1:2], groups = unname(groups[significant_axes[i]]))
      colors_random = brewer.pal(7, 'Set3')[-2][sample(1:6, nrow(unique(groups[significant_axes[i]])))]
      
      plot(ggplot(to_plot, aes(x = coordinates.1, y = coordinates.2, color = groups)) + 
             geom_point(alpha = 0.5) + theme_bw() + scale_color_manual(values = colors_random) + 
             xlab(paste0('bgPCA axis 1 (', round((bgPCA$combinedVar[1,2]*100), 2), '% of variance)')) + 
             theme(legend.title = element_blank()) +
             ylab(paste0('bgPCA axis 2 (', round((bgPCA$combinedVar[2,2]*100), 2), '% of variance)')))
    }
    
    
    
  }
  

  
  #plot the nodes with the strongest differences between runs
  if(!is.na(amount_of_change)) {
    num_nodes = length(which((apply(bgPCA$groupmeans, 2, max)-apply(bgPCA$groupmeans, 2, min)) > amount_of_change))
    if(num_nodes > 20) num_nodes = 20
    if(num_nodes == 0) num_nodes = 5
  } else {
    if(is.na(chosen_clades)) {
      num_nodes = 5
    } else {
      num_nodes = chosen_clades
    }
    if(num_nosed > 20) num_nodes = 20
  }
  
  plots = vector(mode = "list", length = num_nodes)
  for(j in 1:num_nodes) {
    if(length(unique(groups)) == 2) {
      clade = which(sort((apply(bgPCA$groupmeans, 2, max)-apply(bgPCA$groupmeans, 2, min)), decreasing = T)[j] == 
                      (apply(bgPCA$groupmeans, 2, max)-apply(bgPCA$groupmeans, 2, min)))
    } else {
      clade = which(sort(apply(bgPCA$groupmeans, 2, sd), decreasing = T)[j] == apply(bgPCA$groupmeans, 2, sd))
    }
    node = mrca.phylo(tree, clades[[clade]])
    desc = Descendants(tree, node = node, type = 'children')
    for(k in 1:length(desc)) {
      if(as.numeric(desc[k]) <= length(tree$tip.label)) {
        desc[k] = tree$tip.label[as.numeric(desc[k])]
      } else {
        desc[k] = tree$tip.label[sample(unlist(Descendants(tree, node = as.numeric(desc[k]), type = 'tips')), 1)]
      }
    }
    
    to_plot = data.frame(age = ages[,clade], group = groups)
    plots[[j]] = ggplot(to_plot, aes(x = -age, color = groups)) + geom_density(alpha = 0.3, size = 2) + 
      theme_bw() + scale_color_manual(values = colors_random) + 
      theme(plot.title = element_text(size = 8)) + 
      scale_x_continuous(breaks = pretty(-to_plot$age), labels = abs(pretty(-to_plot$age))) + 
      xlab('Age of MRCA') + ylab('Density')

    if(length(unique(groups)) == 2) {
      plots[[j]] = plots[[j]] + ggtitle(paste0(desc[1], ' and ', desc[2], ' (difference = ', 
                                               round((max(bgPCA$groupmeans[,clade])-min(bgPCA$groupmeans[,clade])), 1), ' Ma)'))
    } else {
      plots[[j]] = plots[[j]] + ggtitle(paste0(desc[1], ' and ', desc[2], ' (max difference = ', 
                                               round((max(bgPCA$groupmeans[,clade])-min(bgPCA$groupmeans[,clade])), 1), ' Ma, mean difference = ', 
                                               round(mean(dist(bgPCA$groupmeans[,clade])), 1), ' Ma)'))
    }
  }
  
  plot(annotate_figure(ggarrange(plotlist = plots, common.legend = T, legend = 'bottom', ncol = ceiling(num_nodes/5), nrow = 5)))


  ####################################Check this works when multivariate
  mean = matrix(bgPCA$Grandmean, ncol = 1)
  
  for(i in 1:num_functions) {
    tree_for_plotting = all_trees[[1]]
    tree_for_plotting$edge.length = rep(0, length(tree_for_plotting$edge.length))
    
    assign(paste0('plus_sd_', i), showPC(sdev*sd(scores[,i]), bgPCA$groupPCs[,i], mean))
    assign(paste0('minus_sd_', i), showPC(-sdev*sd(scores[,i]), bgPCA$groupPCs[,i], mean))
    
    clade_size = unlist(lapply(clades, length))
    
    tree_for_plotting_mean = tree_for_plotting
    tree_for_plotting_plus = tree_for_plotting
    tree_for_plotting_minus = tree_for_plotting
    for(j in 2:max(clade_size)) {
      which_clades = which(clade_size == j)
      if(length(which_clades) > 0) {
        for(k in 1:length(which_clades)) {
          node_to_change = getMRCA(tree_for_plotting, unlist(clades[which_clades[k]]))
          dif_minus = get(paste0('minus_sd_', i))[which_clades[k],]
          dif_mean = mean[which_clades[k],]
          dif_plus = get(paste0('plus_sd_', i))[which_clades[k],]
          
          if(j == 2) {
            branches_to_descendants = which(tree_for_plotting$edge[,1] == node_to_change)
            tree_for_plotting_minus$edge.length[branches_to_descendants] = dif_minus
            tree_for_plotting_mean$edge.length[branches_to_descendants] = dif_mean
            tree_for_plotting_plus$edge.length[branches_to_descendants] = dif_plus
          } else {
            nodes_of_descendants = tree_for_plotting$edge[,2][which(tree_for_plotting$edge[,1] == node_to_change)]
            
            if(any(nodes_of_descendants %in% 1:length(tree_for_plotting$tip.label))) {
              singletons = nodes_of_descendants[which(nodes_of_descendants %in% 1:length(tree_for_plotting$tip.label))]
              branches_to_singletons = which(tree_for_plotting$edge[,2] == singletons)
              tree_for_plotting_minus$edge.length[branches_to_singletons] = dif_minus
              tree_for_plotting_mean$edge.length[branches_to_singletons] = dif_mean
              tree_for_plotting_plus$edge.length[branches_to_singletons] = dif_plus
              nodes_of_descendants = nodes_of_descendants[-which(nodes_of_descendants == singletons)]
            }
            
            for(l in 1:length(nodes_of_descendants)) {
              tips = unlist(Descendants(tree_for_plotting, nodes_of_descendants[l], type = 'tips'))
              tree_for_plotting_minus$edge.length[which(tree_for_plotting_minus$edge[,2] == nodes_of_descendants[l])] = 
                dif_minus - dist.nodes(tree_for_plotting_minus)[tips[1], nodes_of_descendants[l]]
              tree_for_plotting_mean$edge.length[which(tree_for_plotting_mean$edge[,2] == nodes_of_descendants[l])] = 
                dif_mean - dist.nodes(tree_for_plotting_mean)[tips[1], nodes_of_descendants[l]]
              tree_for_plotting_plus$edge.length[which(tree_for_plotting_plus$edge[,2] == nodes_of_descendants[l])] = 
                dif_plus - dist.nodes(tree_for_plotting_plus)[tips[1], nodes_of_descendants[l]]
            }
          }
        }
      }
    }

    #work on this
    changes_plus = tree_for_plotting_plus$edge.length - tree_for_plotting_mean$edge.length
    changes_plus_plot = (changes_plus - min(changes_plus))/(max(changes_plus) - min(changes_plus))
    changes_minus = tree_for_plotting_minus$edge.length - tree_for_plotting_mean$edge.length
    changes_minus_plot = (changes_minus - min(changes_minus))/(max(changes_minus) - min(changes_minus))
    
    pal = colorRamp(c('blue', 'grey', 'red')) 
    palette = rgb(pal(c(changes_plus_plot, changes_minus_plot)), max = 255)
    
    par(mfrow=c(1,2))
    plot(tree_for_plotting_plus, show.tip.label = F, edge.width = 3, edge.color = palette[1:length(tree_for_plotting$edge.length)])
    add_time_lines()
    plot(tree_for_plotting_minus, show.tip.label = F, edge.width = 3, edge.color = palette[(length(tree_for_plotting$edge.length)+1):length(palette)])
    add_time_lines()
    
    obj = get("last_plot.phylo", envir = .PlotPhyloEnv)
    t.max = max(obj$xx)
    
    color.legend(t.max-(t.max*0.5), 60, t.max, 63, rect.col = rgb(pal(seq(0,1,0.1)), max = 255), cex = 0.8, 
                 legend = round(seq(min(c(changes_plus, changes_minus)), max(c(changes_plus, changes_minus)), 
                                    ((max(c(changes_plus, changes_minus))-min(c(changes_plus, changes_minus)))/10)), 1))
  }
}

bgPCA_clocks(ages, clock)
bgPCA_clocks(ages, model)
bgPCA_clocks(ages, genes)

#Plot chronograms----------------------------------------------------------------------------------------------
setwd('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Time calibration/final_runs')
folders = list.dirs(recursive = F, full.names = F)
for(i in 1:length(folders)) {
  files = list.files(paste0(getwd(), '/', folders[i]), pattern = 'chronogram')
  for(j in 1:length(files)) {
    if(i == 1 & j == 1) {
      trees = list(read.tree(paste0(paste0(getwd(), '/', folders[i], '/', files[j]))))
      class(trees) = "multiPhylo"
    } else {
      tree = list(read.tree(paste0(paste0(getwd(), '/', folders[i], '/', files[j]))))
      class(tree) = "multiPhylo"
      trees = c(trees, tree)
    }
  }
}

densiTree(trees, type = 'phylogram', width = 3)

#CorrTest-----------------------------------------------------------------------
source('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Time calibration/CorrTest.R')

tree = read.tree('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Phylogenomics_III/Analyses/ExaBayes/rooted_ExaBayes_ConsensusExtendedMajorityRuleNewick.echinoids.tre')
outgroup = tree$tip.label[1:12]

rate.CorrTest(tree, outgroup, sister.resample = 0, outputFile)

tree = read.tree('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Phylogenomics_III/Analyses/LG4X/unpartitioned_LG4X.raxml_rooted.bestTree')
rate.CorrTest(tree, outgroup, sister.resample = 0, outputFile)

tree = read.tree('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Phylogenomics_III/Analyses/Partitioned/rooted_echinoids3.txt.contree')
rate.CorrTest(tree, outgroup, sister.resample = 0, outputFile)

tree = read.tree('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Phylogenomics_III/Analyses/Phylobayes/rooted_echinoidscatgtr100.con.tre')
rate.CorrTest(tree, outgroup, sister.resample = 0, outputFile)

#Median heights-----------------------------------------------------------------
setwd('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Time calibration/final_runs/LN')
trees = read.tree('clockCATGTR1_ln_sample.datedist_combined.tre')

clades = list()
for(i in 1:trees[[1]]$Nnode) {
  clades[i] = list(trees[[1]]$tip.label[unlist(Descendants(trees[[1]], length(trees[[1]]$tip.label)+i, type = 'tips'))])
}

median_tree = trees[[1]]
median_tree$edge.length = 0

ages = matrix(0, ncol = length(clades), nrow = length(trees))

root = sapply(trees, function(x) max(nodeHeights(x)))
ages[,1] = root

for(i in 2:ncol(ages)) {
  node = sapply(trees, findMRCA, clades[[i]])
  if(length(unique(node)) == 1) {
    node = node[1]
    pos = which(median_tree$edge[,2] == node)
    ages[,i] = root - sapply(trees, function(x) nodeHeights(x)[pos,2])
  } else {
    for(j in 1:length(trees)) {
      ages[j,i] = root[j] - nodeHeights(trees[[j]])[which(trees[[j]]$edge[,2] == node[j]),2]
    }
  }
}

medians = round(apply(ages, 2, median), 2)
clade_size = unlist(lapply(clades, length))

for(j in 2:max(clade_size)) {
  which_clades = which(clade_size == j)
  if(length(which_clades) > 0) {
    for(k in 1:length(which_clades)) {
      node_to_change = getMRCA(median_tree, unlist(clades[which_clades[k]]))
      age = medians[which_clades[k]]
      
      median_tree$node.label[node_to_change-length(median_tree$tip.label)] = age
      
      if(j == 2) {
        branches_to_descendants = which(median_tree$edge[,1] == node_to_change)
        median_tree$edge.length[branches_to_descendants] = age
      } else {
        nodes_of_descendants = median_tree$edge[,2][which(median_tree$edge[,1] == node_to_change)]
        
        if(any(nodes_of_descendants %in% 1:length(median_tree$tip.label))) {
          singletons = nodes_of_descendants[which(nodes_of_descendants %in% 1:length(median_tree$tip.label))]
          branches_to_singletons = which(median_tree$edge[,2] == singletons)
          median_tree$edge.length[branches_to_singletons] = age
          nodes_of_descendants = nodes_of_descendants[-which(nodes_of_descendants == singletons)]
        }
        
        for(l in 1:length(nodes_of_descendants)) {
          tips = unlist(Descendants(median_tree, nodes_of_descendants[l], type = 'tips'))
          median_tree$edge.length[which(median_tree$edge[,2] == nodes_of_descendants[l])] = 
            age - dist.nodes(median_tree)[tips[1], nodes_of_descendants[l]]
        }
      }
    }
  }
}

write.tree(median_tree, file = 'usefulnessCATGTR_ln.consensus')

tree_as_character = write.tree(median_tree, file = '', digits = 4)

ages = apply(ages, 2, sort)

mins = round(ages[ceiling(nrow(ages)*0.025),], 2)
maxs = round(ages[floor(nrow(ages)*0.975),], 2)

for(i in 1:length(medians)) {
  comma = F
  split = unlist(strsplit(tree_as_character, medians[i]))
  if(length(split) > 2) {
    split[1] = paste0(split[1:(length(split)-1)], collapse = ' ')
    split[1] = gsub(' ', medians[i], split[1])
    split[2] = split[length(split)]
    split = split[1:2]
  }
  if(split[2] != ';') {
    brlength = unlist(strsplit(split[2], ')'))[1]
    if(grepl(',', brlength)) {
      brlength = unlist(strsplit(split[2], ','))[1]
      comma = T
    }
    if(comma) {
      rest = paste0(unlist(strsplit(split[2], ','))[2:length(unlist(strsplit(split[2], ',')))], collapse = ',')
      tree_as_character = paste0(split[1], medians[i], paste0('[age_95%HPD={', mins[i], ',', maxs[i], '}]'), brlength, ',', rest)
    } else {
      rest = paste0(unlist(strsplit(split[2], ')'))[2:length(unlist(strsplit(split[2], ')')))], collapse = ')')
      tree_as_character = paste0(split[1], medians[i], paste0('[age_95%HPD={', mins[i], ',', maxs[i], '}]'), brlength, ')', rest)
    }
  }
}

library(MCMCtreeR)
ages_list = split(ages, rep(1:ncol(ages), each = nrow(ages)))
for(k in 1:length(ages_list)) {
  names(ages_list)[k] = getMRCA(median_tree, unlist(clades[k]))-length(median_tree$tip.label)
}

MCMC.tree.plot(median_tree, method = 'user', node.ages = ages_list)

#Ages of selected nodes---------------------------------------------------------
setwd('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Time calibration/final_runs')
load(file = 'data_ages.Rda')
load(file = 'data_ages_cons.Rda')

data_ages = data_ages %>% mutate(type = paste0(clock, genes, model), 
                                 type_simple = paste0(clock, genes))

colors = c('#fee6ce', '#fdd0a2', '#fdae6b', '#fd8d3c', '#f16913', 
           '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6')

colors = c('#fde0dd', '#fcc5c0', '#fa9fb5', '#f768a1', '#dd3497', 
           '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6')

#Echinoidea
ggplot(data_ages, aes(x = -clade_13, fill = type_simple)) + geom_density(position="stack") + 
  scale_x_continuous(breaks = pretty(-data_ages$clade_13), 
                     labels = abs(pretty(-data_ages$clade_13))) + theme_bw() + 
  scale_fill_manual(values = colors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1)

length(which(data_ages$clade_13 < 298.9 & data_ages$clade_13 > 272.95))/nrow(data_ages)
length(which(data_ages$clade_13 > 298.9))/nrow(data_ages)

#Euechinoidea
ggplot(data_ages, aes(x = -clade_18, fill = type_simple)) + geom_density(position="stack") + 
  scale_x_continuous(breaks = pretty(-data_ages$clade_18), 
                     labels = abs(pretty(-data_ages$clade_18))) + theme_bw() + 
  scale_fill_manual(values = colors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1)

length(which(data_ages$clade_18 > 251.902))/nrow(data_ages)

#Aulodonta
ggplot(data_ages, aes(x = -clade_19, fill = type_simple)) + geom_density(position="stack") + 
  scale_x_continuous(breaks = pretty(-data_ages$clade_19), 
                     labels = abs(pretty(-data_ages$clade_19))) + theme_bw() + 
  scale_fill_manual(values = colors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1)

#Irregularia
ggplot(data_ages, aes(x = -clade_44, fill = type_simple)) + geom_density(position="stack") + 
  scale_x_continuous(breaks = pretty(-data_ages$clade_44), 
                     labels = abs(pretty(-data_ages$clade_44))) + theme_bw() + 
  scale_fill_manual(values = colors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1)

#Clypeasteroida
ggplot(data_ages, aes(x = -clade_57, fill = type_simple)) + geom_density(position="stack") + 
  scale_x_continuous(breaks = pretty(-data_ages$clade_57), 
                     labels = abs(pretty(-data_ages$clade_57))) + theme_bw() + 
  scale_fill_manual(values = colors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1)

1-(length(which(data_ages$clade_57 < 66))/nrow(data_ages))

#Scutelloida
ggplot(data_ages, aes(x = -clade_48, fill = type_simple)) + geom_density(position="stack") + 
  scale_x_continuous(breaks = pretty(-data_ages$clade_48), 
                     labels = abs(pretty(-data_ages$clade_48))) + theme_bw() + 
  scale_fill_manual(values = colors) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1)

1-(length(which(data_ages$clade_48 < 66))/nrow(data_ages))


#survivors P-T
survivors = vector(length = nrow(data_ages))
which = as.list(vector(length = nrow(data_ages)))
for(i in 1:nrow(data_ages)) {
  survivors[i] = length(which(data_ages[i,13:(ncol(data_ages)-4)] > 251.902)) + 1
  which[[i]] = which(data_ages[i,13:(ncol(data_ages)-4)] > 251.902)
}

surviving_clades = sort(unique(unlist(unique(which))))
which2 = which

for(i in 1:length(which)) {
  surv = surviving_clades[which(surviving_clades %in% which[[i]])]
  surv2 = c()
  if(1 %in% surv) surv2 = c(surv2, 'Cidaroida', 'Euechinoidea')
  if(2 %in% surv) {
    surv2 = c(surv2, 'Cidaroidea', 'Histocidaroidea')
    surv2 = surv2[-which(surv2 == 'Cidaroida')]
  }
  if(6 %in% surv) {
    surv2 = c(surv2, 'Aulodonta', 'Carinacea')
    surv2 = surv2[-which(surv2 == 'Euechinoidea')]
  }
  if(7 %in% surv) {
    surv2 = c(surv2, 'Echinothuriacea', 'Diadematacea')
    surv2 = surv2[-which(surv2 == 'Aulodonta')]
  }
  if(11 %in% surv) {
    surv2 = c(surv2, 'Pedinoida + Aspidodiadematoida', 'Echinothurioida')
    surv2 = surv2[-which(surv2 == 'Echinothuriacea')]
  }
  which2[[i]] = surv2
}

unique(which2)

survivors = data.frame(survivors)
table(survivors)/nrow(survivors)
mean(survivors$survivors)

ggplot(data.frame(survivors)) + geom_bar(mapping = aes(x = survivors, y = ..prop..), stat = "count") + 
  scale_y_continuous(labels = scales::percent_format()) + theme_bw() + coord_flip()

#LTT plots----------------------------------------------------------------------
echinoid_consensus = lapply(all_consensus, drop.tip, tip = which(all_consensus[[1]]$tip.label %not in% clades[[13]]))
class(echinoid_consensus) = "multiPhylo"
mltt.plot(echinoid_trees, log = 'y', legend = F, dcol = F)

#Plot logL traces---------------------------------------------------------------
setwd('C:/Users/mongi/Dropbox/Work/Labo/Echinoidea/Time calibration/final_runs')
folders = list.dirs(recursive = F, full.names = F)
for(i in 1:length(folders)) {
  if(i == 1) {
    traces = paste0(folders[i], '/', list.files(paste0(getwd(), '/', folders[i]), pattern = 'trace'))
  } else {
    traces = c(traces, paste0(folders[i], '/', list.files(paste0(getwd(), '/', folders[i]), pattern = 'trace')))
  }
}

for(i in 1:length(traces)) {
  name = gsub('.trace', '', traces[i])
  if(grepl('_ln', name)) name = gsub('_ln', '', name)
  name = gsub('/', '_', sapply(strsplit(name, '_'), `[`, 1))
  if(grepl('CAT', name)) {
    name = gsub('CAT', '_CAT', name)
  } else {
    name = gsub('GTR', '_GTR', name)
  }
  
  trace = read.table(paste0(getwd(), '/', traces[i]), header = F, sep = '\t')[,c(1,4)]
  name = rep(name, nrow(trace))
  
  if(i == 1) {
    all_data = name %>% as.tibble() %>% mutate(generation = trace[,1], 
                                               trace = trace[,2], 
                                               clock = sapply(strsplit(name, '_'), '[', 1), 
                                               genes = sapply(strsplit(name, '_'), '[', 2), 
                                               model = str_sub(sapply(strsplit(name, '_'), '[', 3), 1, -2), 
                                               run = str_sub(sapply(strsplit(name, '_'), '[', 3), -1, -1))
  } else {
    data = name %>% as.tibble() %>% mutate(generation = trace[,1], 
                                           trace = trace[,2], 
                                           clock = sapply(strsplit(name, '_'), '[', 1), 
                                           genes = sapply(strsplit(name, '_'), '[', 2), 
                                           model = str_sub(sapply(strsplit(name, '_'), '[', 3), 1, -2), 
                                           run = str_sub(sapply(strsplit(name, '_'), '[', 3), -1, -1))
    all_data = bind_rows(all_data, data)
  }
}

all_data = all_data %>% arrange(genes) %>% mutate(color = paste(clock, model, run, sep = '-'))

ggplot(subset(subset(all_data, model == 'CATGTR'), generation > 100), aes(x = generation, y = trace, color = color)) + 
  geom_line() + 
  facet_wrap(vars(clock, genes), nrow = 2, scales = 'free') + theme_bw() + 
  scale_color_manual(values = c('#e08214', '#fdb863', '#b2abd2', '#d8daeb')) + 
  geom_vline(xintercept = 10000, linetype = 'dashed')
