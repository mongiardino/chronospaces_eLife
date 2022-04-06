#Chronospace: a tool to visualize and quantify the effect of methodological
#decisions on time-calibrated inference
#Writen by Nicolas Mongiardino Koch 05/2021

#This script requires as input a series of tree files in newick format,
#representing separate Bayesian time-calibrated inferences of the same dataset
#(can be either different runs of the same analysis or separate analyses). The
#topology is assumed constrained, such that chronograms differ only in the
#branch lengths (and therefore, in their inferred node ages). A list with one or
#more character vectors is needed to assign labels to these files.

#Code includes two main functions. The first of these, 'extract_ages', obtains
#node age for the different nodes in all of the trees provided and links each
#topology with a factor (i.e., information about the specific run). The second
#function, 'bgPCA_ages', runs a between-group principal component analysis
#(bgPCA) in order to obtain axes that summarize the the effects of
#methodological decisions, and quantify how much they contribute to the overall
#variance in node ages across analyses. Several plots are produced, including:
#1) A chronospace plot, showing the distribution of chronograms in multivariate space,
#2) The posterior distributions of nodes that change the most between analyses, and 
#3) a summary of the effect of different choices on all branch lengths.

#Install and load packages------------------------------------------------------
packages <- c('ape', 'phangorn', 'phytools',
              'stringr', 'dplyr', 
              'ggplot2', 'ggpubr', 'plotrix', 'RColorBrewer', 'ggtree', 'patchwork')
new_packages <- packages[!packages %in% installed.packages()[,'Package']]
if(length(new_packages)) { install.packages(new_packages) }

library(ape)
library(phangorn)
library(phytools)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotrix)
library(RColorBrewer)
library(ggtree)
library(patchwork)

#Plotting functions-------------------------------------------------------------
#add lines between geological periods THIS IS NOW OBSOLETE
add_time_lines <- function() {
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

times <- c(2.588, 23.03, 66, 145, 201.3, 252.17, 298.9, 
           252.17, 358.9, 419.2, 443.8, 485.4, 541, 4600)

#load files---------------------------------------------------------------------
extract_ages <- function(path = NA, type, sample) {
  
  #Obtain the names of all files in the specified directory (or the working
  #directory otherwise). These should all be newick tree files corresponding to
  #Bayesian posterior distributions of time-calibrated analyses with constrained
  #topology.
  if(is.na(path)) {
    files <- list.files()
  } else {
    files <- list.files(path = path)
  }
  
  #check that tree files and factors provided in 'types' match correctly
  cat("Check that labels are assigned correctly to the input files.\n", 
      "If there is an error, modify the order of factors in 'type'\n", 
      "or the name of input files for the two to match.\n\n", sep = '')
  
  for(i in 1:length(files)) {
    to_print <- paste0('file = ', files[i], ' | type = ', 
                      paste0(unname(as.vector(as.data.frame(type)[i,])), 
                             collapse = ' - '))
    cat(to_print, '\n')
  }
  
  #loop through the tree files, load them and subsample them to the number
  #specified in 'sample'
  for(i in 1:length(files)) {
    trees <- read.tree(paste0(getwd(), '/', files[i]))
    if(!is.na(sample)) {
      trees <- trees[sample(1:length(trees), sample)]
    }
    
    if(i == 1) {
      all_trees <- trees
    } else {
      all_trees <- c(all_trees, trees)
    }
  }
  
  #extract the number of nodes in the tree and the taxonomic composition of each
  #node
  tree <- all_trees[[1]]
  clades <- list()
  for(i in 1:tree$Nnode) {
    clades[i] <- list(tree$tip.label[unlist(Descendants(tree, length(tree$tip.label)+i, type = 'tips'))])
  }
  
  #build the matrix that will contain node ages (columns) for each tree (rows)
  ages <- matrix(0, ncol = length(clades), nrow = length(all_trees))
  
  #assign values from the oldest node (i.e., root)
  root <- sapply(all_trees, function(x) max(nodeHeights(x)))
  ages[,1] <- root
  
  #assign values to all other nodes
  for(i in 2:ncol(ages)) {
    #check node number associated with clade i
    node <- sapply(all_trees, findMRCA, clades[[i]])
    
    #if the clade is always associated with the same number this runs fast
    if(length(unique(node)) == 1) {
      node <- node[1]
      pos <- which(all_trees[[1]]$edge[,2] == node)
      ages[,i] <- root - sapply(all_trees, function(x) nodeHeights(x)[pos,2])
      
      #otherwise this takes a while, but the correct matching is confirmed such
      #that dates from the same clade are placed in the same column regardless
      #of its node number
    } else {
      for(j in 1:length(all_trees)) {
        ages[j,i] <- root[j] - 
          nodeHeights(all_trees[[j]])[which(all_trees[[j]]$edge[,2] == node[j]),2]
      }
    }
  }
  
  #change to data.frame, set column names and add types of runs as factors
  data_ages <- data.frame(ages)
  colnames(data_ages) <- paste0('clade_', 1:ncol(data_ages))
  
  types_of_runs <- data.frame(matrix(NA, nrow = nrow(ages), ncol = length(type)))
  for(i in 1:length(type)) {
    types_of_runs[,i] <- rep(type[[i]], each = nrow(ages)/length(type[[i]]))
  }
  
  colnames(types_of_runs) <- paste('factor', LETTERS[1:ncol(types_of_runs)], 
                                   sep = '_')
  
  data_ages <- cbind(data_ages, types_of_runs)
  data_ages <- data_ages %>% mutate_if(sapply(data_ages, is.character), as.factor)

  #export
  return(data_ages)
}


# internal between-group PCA function ---------------------------------------------
bgprcomp <- function(x, groups){
  
  grandmean <- colMeans(x)
  x_centered <- scale(x, scale = F, center = T)
  x_gmeans <- apply(X = x_centered, MARGIN = 2, FUN = tapply, groups, mean)
  
  V_g <- cov(x_gmeans)
  eig <- eigen(V_g)
  
  scores <- x_centered%*%eig$vectors
  scores <- cbind(scores[,1:(nlevels(groups) - 1)])
  rotation <- eig$vectors
  
  preds <- scores %*% t(rotation[,1:ncol(scores)])
  resids <- x - preds
  
  return(list(x = scores, residuals = resids, rotation = rotation, 
              values = eig$values, center = grandmean, gmeans = x_gmeans))
}


# internal reverse PCA function -----------------------------------------------------
revPCA<-function(scores, vectors, center){ t(t(scores%*%t(vectors))+center) }


#create chronospace-------------------------------------------------------------------
chronospace <- function(data_ages, factors, variation = "non-redundant")  {
  
  #set objects
  ages <- data_ages
  if(is.null(dim(factors))) groups <- data.frame(factor=factors) else groups <- factors
  
  #factors names
  facnames<-colnames(groups)
  
  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = ncol(groups))
  names(results) <- facnames
  
  #perform bgPCA using each factor separately
  for(i in 1:ncol(groups)) {
    
    #perform bgPCA between groups defined by factor i over original variation
    bgPCA1 <- bgprcomp(x = ages, groups = groups[,i])
    
    #compute percentage of variation explained
    totvar <- sum(apply(ages, 2, var))
    expvar <- sum(apply(bgPCA1$x, 2, var))
    perc_tot <- 100 * (expvar/totvar)
    
    #report proportion of total variation explained
    cat(paste0('Proportion of total variation in node ages explained by ', 
               facnames[i], ' = ', 
               round(perc_tot, digits=3), 
               '%', '\n'))
    
    if(ncol(groups)>1){
      #use bgPCA to compute an ordinaion that is residual to all factors but factor i
      bgPCA2.1 <- bgprcomp(x = ages, groups = groups[,-i])
      resids2.1 <- bgPCA2.1$residuals
      
      #perform bgPCA between groups defined by factor i over residual variation
      bgPCA2.2 <- bgprcomp(x = resids2.1, groups = groups[,i])
      expvar2.2 <- sum(apply(bgPCA2.2$x, 2, var))
      
      #compute percentage of non-redundant variation explained 
      perc_nonred <- 100 * (expvar2.2/totvar)
      
      #report proportion of non-redundant variation explained
      cat(paste0('Proportion of non-redundant variation in node ages explained by ', 
                 facnames[i], ' = ', 
                 round(perc_nonred, digits=3), 
                 '%', '\n'))
    } else {cat('(There is only one factor, non-redundant variation omitted)\n')}
    
    #select which bgPCA results are going to be used
    if(variation == "total" | ncol(groups)==1) bgPCA <- bgPCA1
    if(variation == "non-redundant" & ncol(groups)>1) bgPCA <- bgPCA2.2
    
    #store bgPCA results, along with total variation and groups of factor i
    bgPCA$totvar<-totvar
    bgPCA$groups<-groups[,i]
    bgPCA$ages<-ages
    results[[i]]<-bgPCA
  }
  
  return(invisible(results))
}


#plot chronospace-------------------------------------------------------------------
plot.chronospace<-function(obj, tree=NA, sdev=1, timemarks = NULL,
                           colors=1:5, factors=1:length(obj), axes=c(1,2), pt.alpha=0.5, pt.size=1.5, ell.width=1.2, dist.width=1, ct.size=5,
                           ellipses=TRUE, centroids=FALSE, distances=FALSE) {
  
  if(length(axes)!=2) axes<-c(1,2)
  
  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = length(obj))
  names(results) <- facnames <- names(obj)
  
  #get ordinations and Pc extremes for factor i
  for(i in 1:length(obj)){
    
    #create object for storing results of factor i, assing names
    results_i <- vector(mode = "list", length = 2)
    names(results_i) <- c("ordination", "PC_extremes")
    
    #extract information for factor i
    bgPCA <- obj[[i]]
    groups <- bgPCA$groups
    totvar <- bgPCA$totvar
    ages <- bgPCA$ages
    
    #set axes to either 1 (univariate plot) if the variable contains only two
    #groups, or 2 (bivariate plot) if it includes more groups
    num_functions <- 1
    if(ncol(bgPCA$x) >= 2) num_functions = 2
    
    #gather data for plotting
    to_plot <- data.frame(coordinates = bgPCA$x, groups = groups)
    
    #plot chronospace
    if(num_functions == 1) { #univariate
      
      chronospace <- ggplot(to_plot, aes(x = coordinates, fill = groups)) + 
        geom_histogram(alpha = 0.5, position = 'identity', bins = 30) + theme_bw() + 
        scale_fill_manual(values = colors) + ylab('Count') + 
        theme(legend.title = element_blank(), panel.grid = element_blank()) +
        xlab(paste0('bgPCA axis 1 (', round((100*apply(bgPCA$x,2,var)[1]/totvar), 2), '% of variance)'))
      
    } else { #bivariate
      #compute groups centroids from bgPCA scores
      cents<-apply(X = bgPCA$x, MARGIN = 2, FUN = tapply, groups, mean)
      cents_df<-data.frame(coordinates.1=cents[,1], coordinates.2=cents[,2], groups=rownames(cents))
      
      #compute groups centroids from original variables; calculate and standardize distances between centroids
      cents_original<-apply(X = ages, MARGIN = 2, FUN = tapply, groups, mean)
      dists<-as.matrix(dist(cents_original))
      dists_std<-dists/max(dists)
      
      #generate combinations
      combins<-combn(x = levels(groups), m = 2)
      
      #plot chronospace
      chronospace<-ggplot(to_plot, aes(x = coordinates.1, y = coordinates.2, color = groups)) + 
        geom_point(alpha = pt.alpha, size=pt.size, key_glyph = "point") + 
        theme_bw() + scale_color_manual(values = colors) + 
        theme(legend.title = element_blank(), panel.grid = element_blank()) + 
        xlab(paste0('bgPCA axis ', axes[1],  ' (', round((100*apply(bgPCA$x,2,var)[1]/totvar), 2), '% of variance)')) + 
        ylab(paste0('bgPCA axis ', axes[2],  ' (', round((100*apply(bgPCA$x,2,var)[2]/totvar), 2), '% of variance)'))
      
      if(ellipses){
        chronospace <- chronospace + 
          stat_ellipse(lwd=ell.width, key_glyph = "point")
      }
      
      if(distances){
        for(h in 1:ncol(combins)){
          rdf<-cents_df[combins[,h],]
          width<-(5*dists_std[combins[1,h], combins[2,h]])-2
          chronospace <- chronospace + geom_line(data=rdf, aes(x = coordinates.1, y = coordinates.2), color=gray.colors(n=10)[1], size=width*dist.width)
        }
      }
      
      if(centroids|distances){
        chronospace <- chronospace + 
          geom_point(shape=21, data=cents_df, color="black", fill = colors[1:nlevels(groups)], aes(x = coordinates.1, y = coordinates.2), size=ct.size)
      }
      
      chronospace <- chronospace + 
        guides(colour = guide_legend(override.aes = list(alpha=1, shape=21, color="black", fill = colors[1:nlevels(groups)], size=3.5)))
      
    }
    
    #save chronospace
    results_i$ordination <- chronospace
    
    #obtain clades from tree
    clades <- list()
    for(j in 1:tree$Nnode) {
      clades[j] <- list(tree$tip.label[unlist(Descendants(tree, length(tree$tip.label)+j, type = 'tips'))])
    }
    
    #Finally, compute changes in each branch captured by the bgPCA axes (needs a
    #tree!)
    if(is.na(tree)[1]) {
      cat('Plotting changes on branch lengths can only be shown if a tree is provided\n')
    } else {
      #ages implied by a position at the origin of the bgPCA plot
      mean <- matrix(colMeans(ages), ncol = 1)
      
      #create object for storing the extremes of the bgPC j
      PCextremes <- vector(mode = "list", length = num_functions)
      
      if(num_functions==1) ax<-1 else ax<-axes
      
      #loop through the bgPCA axes (depending on the number of groups in the
      #variable being tested)
      for(j in 1:num_functions) {
        
        #create a tree that contains topology but no branch lengths
        tree$edge.length <- rep(0, length(tree$edge.length))
        
        #ages implied by moving along this bgPCA axis 'sdev' number of standard
        #deviations to both sides
        assign(paste0('plus_sd_', j), revPCA(sdev*sd(bgPCA$x[,ax[j]]), bgPCA$rotation[,ax[j]], mean))
        assign(paste0('minus_sd_', j), revPCA(-sdev*sd(bgPCA$x[,ax[j]]), bgPCA$rotation[,ax[j]], mean))
        
        #check number of descendants stemming from each node
        clade_size <- unlist(lapply(clades, length))
        
        #setup trees that will have mean branch lengths, mean+sdev and mean-sdev
        #trees
        tree_mean <- tree_plus <- tree_minus <- tree
        
        #loop through clades from smallest to biggest (i.e., up the tree)
        for(k in 2:max(clade_size)) {
          #which nodes have the number of descendants
          which_clades <- which(clade_size == k)
          if(length(which_clades) > 0) {
            for(l in 1:length(which_clades)) {
              #which node are we talking about
              node_to_change <- getMRCA(tree, unlist(clades[which_clades[l]]))
              
              #get node ages for this node
              dif_minus <- get(paste0('minus_sd_', j))[,which_clades[l]]
              dif_mean <- mean[which_clades[l],]
              dif_plus <- get(paste0('plus_sd_', j))[,which_clades[l]]
              
              #if the clade is a cherry (i.e., 2 descendants)
              if(k == 2) {
                #get branches descending to both tips and assign them their
                #correct branches (which is == to the node age)
                branches_to_descendants <- which(tree$edge[,1] == node_to_change)
                tree_minus$edge.length[branches_to_descendants] <- dif_minus
                tree_mean$edge.length[branches_to_descendants] <- dif_mean
                tree_plus$edge.length[branches_to_descendants] <- dif_plus
                
                #if it is not a cherry
              } else {
                #get nodes of direct descendant
                nodes_of_descendants <- tree$edge[,2][which(tree$edge[,1] == node_to_change)]
                
                #if any descendant is a tip do as above, assign branch length ==
                #node age
                if(any(nodes_of_descendants %in% 1:length(tree$tip.label))) {
                  singletons <- nodes_of_descendants[which(nodes_of_descendants %in% 
                                                             1:length(tree$tip.label))]
                  branches_to_singletons <- which(tree$edge[,2] == singletons)
                  tree_minus$edge.length[branches_to_singletons] <- dif_minus
                  tree_mean$edge.length[branches_to_singletons] <- dif_mean
                  tree_plus$edge.length[branches_to_singletons] <- dif_plus
                  #remove it from descendants as its branch length is already
                  #set
                  nodes_of_descendants <- nodes_of_descendants[-which(nodes_of_descendants == 
                                                                        singletons)]
                }
                
                #for descendant clades do the following
                for(m in 1:length(nodes_of_descendants)) {
                  #obtain all descendants
                  tips <- unlist(Descendants(tree, nodes_of_descendants[m], 
                                             type = 'tips'))
                  
                  #remove from the age the age of the descendant node, which is
                  #already set up correctly as the loop goes from smaller to
                  #larger clades
                  tree_minus$edge.length[which(tree_minus$edge[,2] == nodes_of_descendants[m])] <- 
                    dif_minus - dist.nodes(tree_minus)[tips[1], nodes_of_descendants[m]]
                  tree_mean$edge.length[which(tree_mean$edge[,2] == nodes_of_descendants[m])] <- 
                    dif_mean - dist.nodes(tree_mean)[tips[1], nodes_of_descendants[m]]
                  tree_plus$edge.length[which(tree_plus$edge[,2] == nodes_of_descendants[m])] <- 
                    dif_plus - dist.nodes(tree_plus)[tips[1], nodes_of_descendants[m]]
                }
              }
            }
          }
        }
        
        #compute delta in branch lengths between the mean tree and the positive and negative extremes
        changes_plus <- tree_plus$edge.length - tree_mean$edge.length
        changes_minus <- tree_minus$edge.length - tree_mean$edge.length
        
        #if time marks have been specified, use them to  draw vertical lines in the corresponding tree
        if(!is.null(timemarks)){
          t.max <- max(nodeHeights(tree_minus))
          timemarks1.1 <- timemarks[which(timemarks <= t.max)]
          timemarks1.2 <- t.max - timemarks1.1
          
          t.max <- max(nodeHeights(tree_plus))
          timemarks2.1 <- timemarks[which(timemarks <= t.max)]
          timemarks2.2 <- t.max-timemarks2.1
        } else {
          timemarks1.2 <- timemarks2.2 <- NULL
        }
        
        #convert phylo trees into ggtrees, adding delta in branch length to the metadata
        tree_plus_gg <- ggtree(tree_plus, size = 1.5) %<+% 
          data.frame(node = tree_plus$edge[,2], delta = changes_plus)
        tree_minus_gg <- ggtree(tree_minus, size = 1.5) %<+% 
          data.frame(node = tree_minus$edge[,2], delta = changes_minus)
        
        #create graphics for each extreme of the bgPC j
        negative <- tree_minus_gg + aes(color=delta) + 
          scale_color_gradient2(limits = range(c(changes_minus, changes_plus)),
                                high = "red", low = "blue", mid = "gray", midpoint = 0) +
          ggtitle(paste0(facnames[i], " - bgPC", j, ", negative extreme")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_vline(xintercept = timemarks1.2, lty = 2, col = "gray")
        
        positive <- tree_plus_gg + aes(color = delta) + 
          scale_color_gradient2(limits = range(c(changes_minus, changes_plus)),
                                high = "red", low = "blue", mid = "gray", midpoint = 0) +
          ggtitle(paste0(facnames[i], " - bgPC", j, ", positive extreme")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_vline(xintercept = timemarks2.2, lty = 2, col = "gray")
        
        #combine both into a single graphic and store
        PCextremes[[j]] <- negative + positive + plot_layout(guides = "collect") & 
          theme(legend.position="bottom")
        
      }
      
      #assign list names and save
      names(PCextremes) <- paste0("bgPC", 1:j)
      results_i$PC_extremes <- PCextremes
      
    }
    
    #add to overall results list
    results[[i]] <- results_i
    
  }
  
  factors<-factors[!factors>length(results)]
  return(results[factors])
  
}

#get senstive nodes ----------------------------------------------------
sensitive_nodes <- function(obj, tree, amount_of_change, factors=1:length(obj), 
                             chosen_clades, colors=1:5){
  
  #create object for storing overall results, assign names
  results <- vector(mode = "list", length = length(obj))
  names(results) <- names(obj)
  
  #perform bgPCA on each variable
  for(i in 1:length(obj)) {
    
    #extract information for factor i
    bgPCA <- obj[[i]]
    groups <- bgPCA$groups
    ages <- bgPCA$ages
    
    #plot the posterior distribution of nodes with the strongest differences
    #between runs
    
    #first decide how many nodes will be plotted
    #if an minimum amount of change is specified, go with it
    if(!is.na(amount_of_change)) {
      num_nodes <- length(which((apply(bgPCA$gmeans, 2, max) -
                                  apply(bgPCA$gmeans, 2, min)) > amount_of_change))
      
      #reduce to a max of 20,or plot 5 if none changes by the specified amount
      if(num_nodes > 20) num_nodes <- 20
      if(num_nodes == 0) num_nodes <- 5
      
    } else { #if a minimum amount is not specified
      #if a number of clades is not specified, do 5
      if(is.na(chosen_clades)) {
        num_nodes <- 5
      } else {
        #else go with what the user chose, although cap at 20
        num_nodes <- chosen_clades
      }
      if(num_nodes > 20) num_nodes <- 20
    }
    
    #obtain clades from tree
    clades = list()
    for(j in 1:tree$Nnode) {
      clades[j] <- list(tree$tip.label[unlist(Descendants(tree, length(tree$tip.label)+j, 
                                                          type = 'tips'))])
    }
    
    #make room to save the individual plots
    plots <- vector(mode = "list", length = num_nodes)
    
    #loop through the nodes
    for(j in 1:num_nodes) {
      
      #sort clades starting by those that vary the most between analyses and
      #choose clade j
      clade <- which(sort((apply(bgPCA$gmeans, 2, max) - apply(bgPCA$gmeans, 2, min)), decreasing = T)[j] ==
                      (apply(bgPCA$gmeans, 2, max) - apply(bgPCA$gmeans, 2, min)))
      
      
      #obtain corresponding node number and the descendant taxa
      node <- mrca.phylo(tree, clades[[clade]])
      desc <- Descendants(tree, node = node, type = 'children')
      
      #obtain representative taxa from either side of the split
      for(k in 1:length(desc)) {
        #if it is a tip, extract the name
        if(as.numeric(desc[k]) <= length(tree$tip.label)) {
          desc[k] <- tree$tip.label[as.numeric(desc[k])]
          #else choose a random tip from the descendant clade
        } else {
          desc[k] <- tree$tip.label[sample(unlist(Descendants(tree, node = as.numeric(desc[k]), 
                                                              type = 'tips')), 1)]
        }
      }
      
      #make the plot
      to_plot <- data.frame(age = ages[,clade], group = groups)
      plots[[j]] <- ggplot(to_plot, aes(x = -age, color = group)) + 
        geom_density(alpha = 0.3, size = 2) +
        theme_bw() + scale_color_manual(values = colors) +
        theme(plot.title = element_text(size = 8)) +
        scale_x_continuous(breaks = pretty(-to_plot$age), labels = abs(pretty(-to_plot$age))) +
        xlab('Age of MRCA') + ylab('Density')
      
      if(length(unique(to_plot$group)) == 2) {
        plots[[j]] <- plots[[j]] + 
          ggtitle(paste0('MRCA of ', desc[1], ' and ', 
                         desc[2], ' (difference = ', 
                         round((max(bgPCA$gmeans[,clade]) - min(bgPCA$gmeans[,clade])), 1), 
                         ' Ma)'))
      } else {
        plots[[j]] <- plots[[j]] + 
          ggtitle(paste0('MRCA of ', desc[1], ' and ', 
                         desc[2], ' (max difference = ', 
                         round((max(bgPCA$gmeans[,clade])-min(bgPCA$gmeans[,clade])), 1), 
                         ' Ma)'))
      }
    }
    
    #plot and save, accounting for a varying number of columns depending on the
    #nodes plotted
    most_affected <- annotate_figure(ggarrange(plotlist = plots,
                                              common.legend = T, legend = 'bottom',
                                              ncol = ceiling(num_nodes/5), nrow = 5))
    #plot(most_affected)
    results[[i]] <- most_affected
  }
  
  factors<-factors[!factors>length(results)]
  return(results[factors])
}