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

#More details can be found in the following publication:
...

#Plotting functions-------------------------------------------------------------
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

#Run bgPCA-------------------------------------------------------------------
bgPCA_ages = function(data_ages, tree = NA, chosen_clades = NA, amount_of_change = NA, sdev = 1) {
  
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
    num_functions = 1
    if(ncol(scores) >= 2) num_functions = 2
    
    #plot chronospace
    if(num_functions == 1) {
      to_plot = data.frame(coordinates = scores, groups = unname(groups[significant_axes[i]]))
      colors_random = brewer.pal(7, 'Set3')[-2][sample(1:6, 2)]
      
      chronospace = ggplot(to_plot, aes(x = coordinates, fill = groups)) + 
        geom_histogram(alpha = 0.5, position = 'identity', bins = 30) + theme_bw() + 
        scale_fill_manual(values = colors_random) + ylab('Count') + 
        theme(legend.title = element_blank()) +
        xlab(paste0('bgPCA axis 1 (', round((bgPCA$combinedVar[1,2]*100), 2), '% of variance)'))
      
    } else {
      to_plot = data.frame(coordinates = scores[,1:2], groups = unname(groups[significant_axes[i]]))
      colors_random = brewer.pal(7, 'Set3')[-2][sample(1:6, nrow(unique(groups[significant_axes[i]])))]
      
      chronospace = ggplot(to_plot, aes(x = coordinates.1, y = coordinates.2, color = groups)) + 
        geom_point(alpha = 0.5) + theme_bw() + scale_color_manual(values = colors_random) + 
        xlab(paste0('bgPCA axis 1 (', round((bgPCA$combinedVar[1,2]*100), 2), '% of variance)')) + 
        theme(legend.title = element_blank()) + 
        ylab(paste0('bgPCA axis 2 (', round((bgPCA$combinedVar[2,2]*100), 2), '% of variance)'))
    
    }
    
    plot(chronospace)
    if(length(significant_axes) == 1) {
      ggsave('chronospace.pdf', plot = chronospace, width = 10, height = 6, units = 'in')
    } else {
      ggsave(paste0('chronospace_', colnames(groups)[significant_axes[i]], '.pdf'), 
             plot = chronospace, width = 10, height = 6, units = 'in')
    }
    
    #plot the posterior distribution of nodes with the strongest differences
    #between runs
    
    #first decide how many nodes will be plotted
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
      if(num_nodes > 20) num_nodes = 20
    }
    
    plots = vector(mode = "list", length = num_nodes)
    for(j in 1:num_nodes) {
      clade = which(sort((apply(bgPCA$groupmeans, 2, max)-apply(bgPCA$groupmeans, 2, min)), decreasing = T)[j] == 
                        (apply(bgPCA$groupmeans, 2, max)-apply(bgPCA$groupmeans, 2, min)))
      node = mrca.phylo(tree, clades[[clade]])
      desc = Descendants(tree, node = node, type = 'children')
      
      for(k in 1:length(desc)) {
        if(as.numeric(desc[k]) <= length(tree$tip.label)) {
          desc[k] = tree$tip.label[as.numeric(desc[k])]
        } else {
          desc[k] = tree$tip.label[sample(unlist(Descendants(tree, node = as.numeric(desc[k]), type = 'tips')), 1)]
        }
      }
      
      to_plot = data.frame(age = ages[,clade], group = unname(groups[significant_axes[i]]))
      plots[[j]] = ggplot(to_plot, aes(x = -age, color = group)) + geom_density(alpha = 0.3, size = 2) + 
        theme_bw() + scale_color_manual(values = colors_random) + 
        theme(plot.title = element_text(size = 8)) + 
        scale_x_continuous(breaks = pretty(-to_plot$age), labels = abs(pretty(-to_plot$age))) + 
        xlab('Age of MRCA') + ylab('Density')
      
      if(length(unique(to_plot$group)) == 2) {
        plots[[j]] = plots[[j]] + ggtitle(paste0('MRCA of ', desc[1], ' and ', desc[2], 
                                                 ' (difference = ', 
                                                 round((max(bgPCA$groupmeans[,clade])-min(bgPCA$groupmeans[,clade])), 1), 
                                                 ' Ma)'))
      } else {
        plots[[j]] = plots[[j]] + ggtitle(paste0('MRCA of ', desc[1], ' and ', desc[2], 
                                                 ' (max difference = ', 
                                                 round((max(bgPCA$groupmeans[,clade])-min(bgPCA$groupmeans[,clade])), 1), 
                                                 ' Ma)'))
      }
    }
    
    most_affected = plot(annotate_figure(ggarrange(plotlist = plots, 
                                                   common.legend = T, legend = 'bottom', 
                                                   ncol = ceiling(num_nodes/5), nrow = 5)))
    plot(most_affected)
    
    width = (ceiling(num_nodes/5))*4 + 9
    if(length(significant_axes) == 1) {
      ggsave('nodes_most_affected.pdf', plot = most_affected, width = width, height = 16, units = 'in')
    } else {
      ggsave(paste0('nodes_most_affected.pdf_', colnames(groups)[significant_axes[i]], '.pdf'), 
             plot = most_affected, width = width, height = 16, units = 'in')
    }
    
    #compute changes in branches along the different bgPCA axes
    if(is.na(tree)) {
      cat('Plotting changes on branch lengths can only be shown if a tree is provided\n')
    } else {
      mean = matrix(bgPCA$Grandmean, ncol = 1)
      
      for(i in 1:num_functions) {
        tree$edge.length = rep(0, length(tree$edge.length))
        
        assign(paste0('plus_sd_', i), showPC(sdev*sd(scores[,i]), bgPCA$groupPCs[,i], mean))
        assign(paste0('minus_sd_', i), showPC(-sdev*sd(scores[,i]), bgPCA$groupPCs[,i], mean))
        
        clade_size = unlist(lapply(clades, length))
        
        tree_mean = tree
        tree_plus = tree
        tree_minus = tree
        for(j in 2:max(clade_size)) {
          which_clades = which(clade_size == j)
          if(length(which_clades) > 0) {
            for(k in 1:length(which_clades)) {
              node_to_change = getMRCA(tree, unlist(clades[which_clades[k]]))
              dif_minus = get(paste0('minus_sd_', i))[which_clades[k],]
              dif_mean = mean[which_clades[k],]
              dif_plus = get(paste0('plus_sd_', i))[which_clades[k],]
              
              if(j == 2) {
                branches_to_descendants = which(tree$edge[,1] == node_to_change)
                tree_minus$edge.length[branches_to_descendants] = dif_minus
                tree_mean$edge.length[branches_to_descendants] = dif_mean
                tree_plus$edge.length[branches_to_descendants] = dif_plus
              } else {
                nodes_of_descendants = tree$edge[,2][which(tree$edge[,1] == node_to_change)]
                
                if(any(nodes_of_descendants %in% 1:length(tree$tip.label))) {
                  singletons = nodes_of_descendants[which(nodes_of_descendants %in% 1:length(tree$tip.label))]
                  branches_to_singletons = which(tree$edge[,2] == singletons)
                  tree_minus$edge.length[branches_to_singletons] = dif_minus
                  tree_mean$edge.length[branches_to_singletons] = dif_mean
                  tree_plus$edge.length[branches_to_singletons] = dif_plus
                  nodes_of_descendants = nodes_of_descendants[-which(nodes_of_descendants == singletons)]
                }
                
                for(l in 1:length(nodes_of_descendants)) {
                  tips = unlist(Descendants(tree, nodes_of_descendants[l], type = 'tips'))
                  tree_minus$edge.length[which(tree_minus$edge[,2] == nodes_of_descendants[l])] = 
                    dif_minus - dist.nodes(tree_minus)[tips[1], nodes_of_descendants[l]]
                  tree_mean$edge.length[which(tree_mean$edge[,2] == nodes_of_descendants[l])] = 
                    dif_mean - dist.nodes(tree_mean)[tips[1], nodes_of_descendants[l]]
                  tree_plus$edge.length[which(tree_plus$edge[,2] == nodes_of_descendants[l])] = 
                    dif_plus - dist.nodes(tree_plus)[tips[1], nodes_of_descendants[l]]
                }
              }
            }
          }
        }
        
        changes_plus = tree_plus$edge.length - tree_mean$edge.length
        changes_plus_plot = (changes_plus - min(changes_plus))/(max(changes_plus) - min(changes_plus))
        changes_minus = tree_minus$edge.length - tree_mean$edge.length
        changes_minus_plot = (changes_minus - min(changes_minus))/(max(changes_minus) - min(changes_minus))
        
        pal = colorRamp(c('blue', 'grey', 'red')) 
        palette = rgb(pal(c(changes_plus_plot, changes_minus_plot)), max = 255)
        
        par(mfrow=c(1,2))
        plot(tree_plus, show.tip.label = F, edge.width = 3, edge.color = palette[1:length(tree$edge.length)])
        add_time_lines()
        plot(tree_minus, show.tip.label = F, edge.width = 3, edge.color = palette[(length(tree$edge.length)+1):length(palette)])
        add_time_lines()
        
        obj = get("last_plot.phylo", envir = .PlotPhyloEnv)
        t.max = max(obj$xx)
        
        color.legend(t.max-(t.max*0.5), 60, t.max, 63, rect.col = rgb(pal(seq(0,1,0.1)), max = 255), cex = 0.6, 
                     legend = round(seq(min(c(changes_plus, changes_minus)), max(c(changes_plus, changes_minus)), 
                                        ((max(c(changes_plus, changes_minus))-min(c(changes_plus, changes_minus)))/10)), 1))
        
        if(length(significant_axes) == 1) {
          dev.copy2pdf(file = paste0('branch_changes_', sdev, 'sd.pdf'), width = 16, height = 12)
        } else {
          dev.copy2pdf(file = paste0('branch_changes_', sdev, 'sd_', colnames(groups)[significant_axes[i]], '.pdf'), 
                       width = 16, height = 12)
        }
      }
    }
  }
}