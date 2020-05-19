# NODE EDGE INVENTORY ####
getNodeEdgeInventory = function(clust_locs, clust_order=NULL, ccf_cols=ccf_cols) {
  nodes = list()
  edges = list()
  edges_confidence = list()
  if (nrow(clust_locs) > 1) {
    clust_locs_temp = clust_locs
    for (i in 1:nrow(clust_locs)) {
      # Take highest total CCF cluster as starting point - this could be problematic in multi-sample cases as it may pick the wrong root node.
      highest_ccf = which.max(clust_locs_temp[,ccf_cols])
      curr_clust = clust_locs_temp[highest_ccf,]
      clust_locs_temp = clust_locs_temp[-highest_ccf,]

      if (length(nodes)==0) {
        nodes[[1]] = curr_clust
        edges[[1]] = NA
      } else {
        possible_ancestors = c()
        possible_ancestors_conf = c()
        # Find place where this node can go
        for (j in 1:length(nodes)) {
          # use either a provided order classification or list of CCF values
          if (!is.null(clust_order)) {

            # Get coordinates in the data.frame for this pair, using column names
            curr_clust_col_index = which(as.character(curr_clust$cluster.no)==colnames(clust_order))
            # current cluster row index - take into account there is a column for the timepoint
            curr_clust_row_index = curr_clust_col_index-1
            # Node j column index
            node_col_index = which(as.character(nodes[[j]]$cluster.no)==colnames(clust_order))

            # then for every timepoint (sample) in this case, check the to be investigated cluster against node j
            node_fits_beneath = unlist(lapply(unique(clust_order$timepoint), function(timepoint) {
              # temp = subset(clust_order, clust_order$timepoint==timepoint)
              temp = clust_order[clust_order$timepoint==timepoint,]
              # the current cluster must fit below the node, therefore look for classifications containing LT (less than)
              "LT" %in% temp[curr_clust_row_index, node_col_index] | "EQ" %in% temp[curr_clust_row_index, node_col_index]
            }))

          } else {
            # Simply check the CCFs of the node against the current cluster
            node_fits_beneath = sapply(ccf_cols, function(k) { nodes[[j]][k] >= curr_clust[k] })
          }

          if (all(node_fits_beneath)) {
            possible_ancestors = c(possible_ancestors, j)
            possible_ancestors_conf = c(possible_ancestors_conf, j)
          }
        }

        # If there are no ancestors, set the value to NA
        if (length(possible_ancestors)==0) {
          possible_ancestors = NA
        }

        nodes[[length(nodes)+1]] = curr_clust
        edges[[length(edges)+1]] = possible_ancestors
        edges_confidence[[length(edges_confidence)+1]] = possible_ancestors_conf

        # if (is.na(edges[[length(edges)]])) {
        #   print("Could not place node on the tree as there are no suitable ancestors:")
        #   print(curr_clust)
        # }
      }
    }
  } else {
    nodes[[1]] = clust_locs[1,]
    edges[[1]] = NA
  }
  return(list(nodes=nodes, edges=edges, edges_confidence=edges_confidence))
}


# TREE BUILDING FUNCTIONS ####
createNode = function(tree, node_index, cluster, parent, ccf_cols) {
  new_node = data.frame(clusterid=cluster[1],
                        nodeid=node_index,
                        parent=parent,
                        level=tree[tree$nodeid==parent, "level"]+1,
                        cluster[ccf_cols],
                        noevent=cluster[which(colnames(mergedseg)=="noevent")]
  )
  return(new_node)
}

buildUniquePlaceTree = function(nodes, edges, ccf_cols) {
  tree = data.frame()
  for (i in 1:length(nodes)) {
    if (i==1) {
      tree = data.frame(clusterid=nodes[[1]][1], nodeid=1, parent=0, level=1, nodes[[1]][ccf_cols], noevent=nodes[[1]][which(colnames(mergedseg)=="noevent")])
    } else {
      if (!is.na(edges[[i]]) & length(edges[[i]])==1) {
        # tree = rbind(tree, data.frame(clusterid=nodes[[i]][1], nodeid=i, parent=edges[[i]], level=tree[tree$nodeid==edges[[i]], "level"]+1, nodes[[i]][ccf_cols]))
        # tree = rbind(tree, createNode(tree, i, nodes, edges, ccf_cols))
        tree = rbind(tree, createNode(tree, i, nodes[[i]], edges[[i]], ccf_cols))
      }
    }
  }

  return(tree)
}


