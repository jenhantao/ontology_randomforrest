time.print <- function(...) {
  print(c(..., date()))
}


examine.tree <- function(j, rf.model, features, labels, feature.names, verbose=F) {
  # Examines a decision tree in a random forest.  Returns two things as a list
  #
  # 1) Takes a dataset and runs it down the tree.  This could been the
  # datset on which the tree was trained or a new test set.  It infers
  # the effective decision made at each branch from the ratios of
  # labels that go down the right and left branches of each decision.
  # Returns a matrix where each row represents a feature and the
  # columns count the frequencies of four different scenarios
  # (assuming an ontotype value of 0,1, or 2).  This makes a
  # subroutine call to function <traverse>.
  #
  # {0 | 1 2} split, higher "1" to "0" ratio on left branch than right branch 
  # {0 | 1 2} split, higher "1" t0 "0" ratio on right branch than left branch 
  # {0 1 | 2} split, higher "1" to "0" ratio on left branch than right branch 
  # {0 1 | 2} split, higher "1" to "0" ratio on right branch than left branch 
  #
  #
  # 2) Calculates the co-occurrence of every pair of features in the
  # tree.  Returns this as a matrix.
  #
  # Input
  # j : index of the tree in the random forest to inspect
  # rf.model : the randomForest object
  # features : the features matrix
  # labels : ordered factor of labels. First factor is "0", 2nd is "1"
  # feature.names : vector of feature names
  # verbose : verbose output

  
  if (verbose) time.print('Getting features for tree', j)

  suppressPackageStartupMessages(library(randomForest))
  
  # Convert tree to igraph object
  if (verbose) time.print('Creating igraph object from tree')
  g <- rf.tree.to.igraph(getTree(rf.model, j, labelVar=T))
  # Add in disconnected nodes for all features that did not occur in tree
  g <- g + vertices(name = setdiff(feature.names, V(g)$name[ ! is.na(V(g)$name) ] ) )

  # Infer decision reules
  if (verbose) time.print('Inferring decision rules')
  suppressPackageStartupMessages(library(Matrix))
#  features <- Matrix(features)
  feature.usage <- lapply(1:length(feature.names), function(x) c(0,0,0,0))
  decisions <- traverse(getTree(rf.model, j), features, labels, verbose=verbose)

  for (k in 1:nrow(decisions)) {
    # Ratio of "1" to "0" on left side of split
    left.ratio <- decisions[k,4] / decisions[k,3]
    # Ratio of "1" to "0" on right side of split
    right.ratio <- decisions[k,6] / decisions[k,5]
    # A likelihood ratio (divide left ratio by right ratio). >1 means
    # that "1" points got skewed to the left
    likelihood.ratio <- left.ratio / right.ratio

    if (verbose) time.print('Likelihood ratio', likelihood.ratio)
    
    if (is.nan(likelihood.ratio)) likelihood.ratio <- 1
    
    i <- ifelse(decisions[k, 2] < 1,
                if (likelihood.ratio > 1) 1 else 2,
                if (likelihood.ratio > 1) 3 else 4)

    if (verbose) print('Updating feature usage')
    feature.usage[[ decisions[k,1] ]] [i] <- feature.usage[[ decisions[k,1] ]] [i] + 1
    if (verbose) print('Done updating feature usage')
  }
  
  feature.usage <- matrix(unlist(feature.usage), nrow=length(feature.names), ncol=4, byrow=T)
  colnames(feature.usage) <- c('0,12_left', '0,12_right', '01,2_label_0_left', '01,2_right')
  rownames(feature.usage) <- feature.names
  
  if (verbose) time.print('Calculating cooccurrence matrix')
  if (verbose) time.print('Number of features', length(feature.names))
  ## Calculate co-occurrence matrix
  if (verbose) time.print('Calculating shortest paths')
  cooccur <- shortest.paths(g, v=V(g)[name %in% feature.names], to=V(g)[name %in% feature.names], mode='out')
  if (verbose) time.print(c('Dim of cooccur', dim(cooccur)))
  # Collapse columns that have the same name
  if (verbose) time.print('Collapsing matrix')
  dim.names <- colnames(cooccur)
  name.groups <- split(1:ncol(cooccur), dim.names)
  if (verbose) time.print('Collapsing columns')
  cooccur <- sapply(name.groups,
                    function(group) {if (length(group)==1) cooccur[,group] else apply(cooccur[,group], 1, min)} )
  
  ## Slow implementation
  ## cooccur <- sapply(name.groups,
  ##                   function(group) {
  ##                     x <- if (length(group)==1) cooccur[,group] else apply(cooccur[,group], 1, min)
  ##                     sapply(name.groups,
  ##                            function(group2) {
  ##                              if (length(group2)==1) x[group2] else min(x[group2])
  ##                            } )
  ##                   })
  
  ## Slow implementation
  #cooccur <- sapply(name.groups, function(group) Reduce(pmin, data.frame(cooccur[, group])))                           

  rownames(cooccur) <- dim.names
  cooccur <- t(cooccur)
  #cooccur <- sapply(name.groups, function(group) Reduce(pmin, data.frame(cooccur[, group])))
  if (verbose) time.print('Collapsing rows')
  # Currently, each instance of a feature in a decision tree is
  # treated as a separate node and will appear as a separate entry in
  # the <cooccur> matrix.  Here, we collapse the matrix for an
  # aggregate cooccurence of each feature.
  cooccur <- sapply(name.groups,
                    function(group) if (length(group)==1) cooccur[,group] else apply(cooccur[,group], 1, min) )
  rownames(cooccur) <- names(name.groups)
  cooccur <- t(cooccur)
  stopifnot( nrow(cooccur)==length(feature.names) && ncol(cooccur)==length(feature.names) )

  ### IDEA: Instead of having to collapsing the cooccurrence matrix,
  ### just modify the igraph by adding zero-weight edges to nodes of
  ### the same term and re-calculate the shortest paths

  if (verbose) time.print('Exiting examine.tree')
    
  return(list(feature.usage=feature.usage, cooccur=cooccur))
}

traverse <- function(tr, features, labels, verbose=F) {
  # Traverses a decision tree in a random forest and partitions a
  # dataset down the different paths
  #
  # Inp
  # tr : A decision tree from the getTree function of randomForest
  # features : the features matrix
  # labels : ordered factor of labels.  First factor is "0", 2nd is "1"
  # verbose : verbose output
  #
  # Notes:
  # Columns of tr
  # left daughter right daughter split var split point status prediction
  #
  # Output:
  # Matrix where each row represents a decision node
  # 1st column is split variable
  # 2nd column is split point  
  # 3rd through 6th columns show the partitioning of points to the
  # left and right branches  
  # -- 3) "0" to the left
  # -- 4) "1" to the left
  # -- 5) "0" to the right
  # -- 6) "1" to the right
  # 7th column is the row number of the split, w.r.t the original tree <tr>
  
  
  if (verbose) time.print('Entering traverse')
    
  partition.list <- list(1:nrow(features))
  bfs.walk <- c(1)
  #labels <- as.factor(labels)
    
  for (i in 1:nrow(tr)) {
    node <- tr[bfs.walk[i], ]

    # Non-leaves
    if (node['status']==1) {
      partition <- partition.list[[bfs.walk[i]]]
      vals <- features[partition, node['split var']]
      bfs.walk <- c(bfs.walk, node['left daughter'], node['right daughter'])
      partition.list[[ node['left daughter'] ]] <- partition[vals <= node['split point']]
      partition.list[[ node['right daughter'] ]] <- partition[vals > node['split point']]
    }
  }

  if (verbose) time.print('Made partition list of length')
  
  non.leaves <- which(tr[,'status'] == 1)

  if (verbose) time.print(length(non.leaves), 'non-leaves')
  
  h <- t(sapply(non.leaves,
                function(i) unlist(lapply(list(partition.list[[ tr[i, 'left daughter'] ]], partition.list[[ tr[i, 'right daughter'] ]]),
                                          function(part)  sapply(levels(labels), function(a) sum(labels[part] == a)) ) )
                ))

  if (verbose) time.print('Dim of h', dim(h))

  decisions <- cbind(tr[,'split var'][ non.leaves ], tr[,'split point'][ non.leaves ], h, as.numeric(rownames(tr)[non.leaves]))

  if (verbose) time.print('Dim of decisions', dim(decisions))
  
  if (verbose) time.print('Exiting traverse')
  return(decisions)
}

rf.tree.to.igraph <- function(rf.tree) {
  # Convert a decision tree in a random forest into an igraph object

  library(igraph)
  
  no.leaves <- rf.tree[rf.tree[, 'status'] != -1, ]

  # Create edgelist
  edgeL <- lapply(rownames(no.leaves),
                  function(i) { i <- as.numeric(i) ; list(c(i, rf.tree[i, 'left daughter']), c(i, rf.tree[i, 'right daughter']))} )
  edgeL <- matrix(unlist(edgeL), ncol=2, byrow=T)
  g <- graph.edgelist(edgeL)

  # Assign node attributes to igraph object
  V(g)[as.numeric(rownames(rf.tree))]$label <- levels(rf.tree[,'split var'])[rf.tree[,'split var']]
  V(g)[as.numeric(rownames(rf.tree))]$name <- levels(rf.tree[,'split var'])[rf.tree[,'split var']]
  V(g)$label[ is.na(V(g)$label) ] <- 'leaf'

  # Assign edge attributes to igraph object 
  E(g)[g[from=edgeL[,1], to=edgeL[,2], edges=T]]$label <- paste(rep(c('<=', '>'), nrow(no.leaves)),
                                         no.leaves[, 'split point'][as.vector(rbind(1:nrow(no.leaves), 1:nrow(no.leaves)))],
                                         sep=' ')

  return(g)
}
