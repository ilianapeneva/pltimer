#' Process the enriched data into the format required for the ordering step with the Plackett-Luce model
#'
#' @param annotated_segments_file Full path to the file with the annotated CNA
#' @param merged_segments_dir Full path to the directory with the LOH, Gain and HD mergedsegs files
#' @param tumour_type String that represents the type of tumour we work with
#' @param driver_mutations A logical: if TRUE, we include the driver mutations to the model
#' @param drivers_file Full path to the file with the drivers mutation information; default NULL
#' @param mixture_model A logical: if TRUE, it will run a mixture model, and if FALSE, it will run either the PlackettLuce or the PLMIX implementation of the Plackett Luce
#' @param model A string: it could be either PlackettLuce or PLMIX; if it is PlackettLuce, then we don't include the unobserved events, if it is PLMIX we include the unobserved events after the last observed events
#' @param clonal_driver_mutations A vector: the list of driver mutations to be included in the model; should be of the format c('cMut_TP53','cMut_ARID1A'), for example; default = NULL
#' @param output_dir Full path to the directory where the output files and plots are stored
#' @return files with a matrix with true & ordered events, a matrix with the values/ranking of the events, all segments with any CNA, a plot with the overall timing of CNA

order_events_across_chort <- function(annotated_segments_file, tumour_type,merged_segments_dir, output_dir, driver_mutations, drivers_file = NULL, mixture_model, model, clonal_driver_mutations = NULL){
  # STEP 1: LOAD THE DATA ####
  # load losses data for WGD samples - need to pick up the segments with nMin1_A~1####
  allsegs = read.table(annotated_segments_file, header = T, stringsAsFactors = F)
  loss_segments = allsegs[which(allsegs$CNA%in%c('cLoss', 'sLoss') & allsegs$tumour_ploidy>2 & allsegs$nMin1_A>0),]

  # create column with 0 for the noevent, and with 0 for w.mean
  loss_segments$noevent <- rep(0,nrow(loss_segments))
  loss_segments$w.mean  <- rep(0,nrow(loss_segments))

  # load enriched LOH events
  loh_mergedsegs = read.table(paste0(merged_segments_dir, tumour_type,"_LOH_mergedsegs.txt"), header = T,stringsAsFactors = F)
  loh_event_name <- c()
  for (i in 1:nrow(loh_mergedsegs)){
    loh_event_name[i] <- paste0('chr',loh_mergedsegs$chr[i],':',loh_mergedsegs$startpos[i],'-',loh_mergedsegs$endpos[i])
  }

  # get the unique LOH events
  unique_loh_events <- unique(loh_event_name)
  loh_regions <- loh_mergedsegs[!duplicated(loh_mergedsegs[,2:4]),2:4]

  # get the unique noevents
  unique_noevents <- unique(loh_mergedsegs$noevent)

  # add the names of the events and noevent to the table with the chromosomes and the start/endpos
  loh_regions <- cbind(loh_regions, unique_loh_events, unique_noevents)

  # pick only the losses that overlap with the enriched LOH events
  loss_data_wgd_samples <- NULL
  for (i in 1:nrow(loh_regions)){
    loss_segments_chr <- loss_segments[which(loss_segments$chr==loh_regions$chr[i]),]
    if (nrow(loss_segments_chr)>0){
    for (j in 1:nrow(loss_segments_chr)) {
      if (loss_segments_chr$startpos[j]<=loh_regions[i,2] & loss_segments_chr$endpos[j]>=loh_regions[i,2]){
        loss_segments_chr$startpos[j] <- loh_regions[i,2] # change the startpos to match the LOH event
        loss_segments_chr$endpos[j]   <- loh_regions[i,3] # change the endpos to match the LOH event
        loss_segments_chr$noevent[j]  <- loh_regions[i,5] # add noevent
        loss_segments_chr$w.mean[j]   <- loss_segments_chr$frac1_A[j] # add w.mean
        loss_data_wgd_samples <- rbind(loss_data_wgd_samples,loss_segments_chr[j,])
      }
      if (loss_segments_chr$startpos[j]>=loh_regions[i,2] & loss_segments_chr$startpos[j]<=loh_regions[i,3]){
        loss_segments_chr$startpos[j] <- loh_regions[i,2] # change the startpos to match the LOH event
        loss_segments_chr$endpos[j]   <- loh_regions[i,3] # change the endpos to match the LOH event
        loss_segments_chr$noevent[j]  <- loh_regions[i,5] # add noevent
        loss_segments_chr$w.mean[j]   <- loss_segments_chr$frac1_A[j] # add w.mean
        loss_data_wgd_samples <- rbind(loss_data_wgd_samples,loss_segments_chr[j,])
      }
    }
  }
  }

  # get rid of frac1_A and SDfrac_A, which are the 7th and 8th columns
  loss_data_wgd_samples <- loss_data_wgd_samples[,-c(7,8)]

  # remove duplicates
  loss_data_wgd_samples <- loss_data_wgd_samples[!duplicated(loss_data_wgd_samples),]

  # rename the cLoss to cLOH and sLoss to sLOH
  loss_data_wgd_samples$CNA[which(loss_data_wgd_samples$CNA=='cLoss')] <- 'cLOH'
  loss_data_wgd_samples$CNA[which(loss_data_wgd_samples$CNA=='sLoss')] <- 'sLOH'

  # load other CNA data ####
  cna_type <- c('LOH','Gain', 'HD')

  mergedseg = data.frame() # this will contain all CNA types with new noevent numbers assigned
  MAX = 0                  # event counter
  for (cna in 1:length(cna_type)){
    if(file.exists(paste0(merged_segments_dir, tumour_type,"_",cna_type[cna],"_mergedsegs.txt"))){
      # pick the list of enriched events for this type of CNA
      CNA = read.table(paste0(merged_segments_dir, tumour_type,"_",cna_type[cna],"_mergedsegs.txt"), header = T, stringsAsFactors = F, sep = "\t")
      CNA$noevent = CNA$noevent + MAX          # increase the count of the CNA for the subsequent lists we add
      mergedseg = rbind(mergedseg,CNA)         # add this CNA data to the other CNA data
      MAX = max(CNA$noevent)
    }
  }

  # create no.chrs.bearing.mut column in mergedseg
  mergedseg$no.chrs.bearing.mut <- rep('not_applicable', nrow(mergedseg))


  # load driver mutations data ####
  if (driver_mutations==T){
    if (file.exists(drivers_file)){
     # the data should be of format Tumour_Name chr startpos endpos nMaj1_A nMin1_A tumour_ploidy CNA noevent w.mean no.chrs.bearing.mut
     drivers_data <- read.table(drivers_file, header = T, stringsAsFactors = F)

     # increase the noevent counter (1,2,3 in the drivers file) for the driver mutations
     drivers_data$noevent <- drivers_data$noevent + MAX
     colnames(drivers_data) <- colnames(mergedseg)
     drivers_data <- drivers_data[!duplicated(drivers_data[,1:2]),]

     # combine the mergedseg with the drivers data
     mergedseg_with_drivers <- rbind(mergedseg, drivers_data)
     mergedseg <- mergedseg_with_drivers
   }
  }

  # STEP 2: SAMPLE TREES TO GET THE SUBCLONAL EVENTS ORDERING ####

  # pick up all the unique tumour names
  tumour.names <- unique(mergedseg$Tumour_Name)

  # create list to save all the trees
  trees.for.all = list()

  for(k in 1:length(tumour.names)){
    # get all the segments that are for this patient
    clust_locs <- mergedseg[as.character(mergedseg$Tumour_Name)==as.character(tumour.names[k]),]
    clust_locs <- clust_locs[clust_locs$w.mean!=1,]  # pick only the subclonal segments
    clust_locs_dup <- clust_locs[duplicated(clust_locs[,c(2,3,4)]),]

    clust_locs_d <- unique(clust_locs_dup[,c(2,3,4)]) # get the unique events from the duplicated events list

    # check if there are any duplicated events
    if (nrow(clust_locs_d)!=0){

      for (y in 1:nrow(clust_locs_d)){
        chr.c <- clust_locs_d$chr[y]
        start.p <- clust_locs_d$startpos[y]
        end.p <- clust_locs_d$endpos[y]

        # pick up all the segments for this event - fix the counter of the event?????
        clust_locs_all_dup <- clust_locs[clust_locs$startpos==start.p & clust_locs$endpos==end.p & clust_locs$chr==chr.c,]
        no.events <- unique(clust_locs_all_dup$noevent)
        events.collapse <- paste(no.events, collapse=" ")
        for(r in 1: nrow(clust_locs)){
          if(clust_locs$noevent[r] %in% no.events){
            clust_locs$noevent[r] <- events.collapse
          }
        }
      }
    }

    # remove duplicates
    clust_locs <- clust_locs[!duplicated(clust_locs[,c(2,3,4)]),]

    # get nodes and edges using the Node-Egde function
    ccf_cols = which(colnames(mergedseg)=="w.mean") # column in with weighted frac1_A
    all <- getNodeEdgeInventory(clust_locs, clust_order=NULL, ccf_cols=ccf_cols)
    nodes <- all$nodes
    edges <- all$edges

    #### build tree(s) ###
    # this deals with the ordering of the subclonal events (the trees that are built consist only of the subclonal events)
    base_tree <- buildUniquePlaceTree(nodes, edges, ccf_cols)

    trees1 = list()
    trees1[[1]] = base_tree
    violates_pigeonhole = c()

    # Add the remaining nodes
    for (i in 1:length(nodes)) {
      if (! i %in% base_tree$nodeid) {
        # For each possible node, walk over all the trees we have so far to see where it fits
        nodesaved = F
        for (j in 1:length(trees1)) {
          tree = trees1[[j]]
          nodesaved_intree = F
          for (possible_parent in edges[[i]]) {

            if (possible_parent %in% tree$nodeid) {
              lvl_parent = tree[tree$nodeid==possible_parent, "level"]
              node_meets_pigeonhole = sapply(ccf_cols, function(k) { sum(tree[tree$level==(lvl_parent+1), 5]) + nodes[[i]][,k] })

              if (all(node_meets_pigeonhole <= tree[tree$level==lvl_parent,5])) {
                # The node fits on the tree, so add it under its possible_parent
                new_tree = rbind(tree, createNode(tree, i, nodes[[i]], possible_parent, ccf_cols))
                # This node is saved in this tree in the next step, make the recording first
                nodesaved = T

                if (!nodesaved_intree) {
                  # This node was not saved in this tree yet, so add it and save this tree
                  trees1[[j]] = new_tree
                  nodesaved_intree = T
                } else {
                  # This node was already saved in this tree and has another position it can go, therefore save another tree
                  trees1[[length(trees1)+1]] = new_tree
                }
              }
            }
          }
        }

        if (!nodesaved) {
          print(paste("Cannot place node", i, "on tree as it violates the pigeonhole principle"))
          violates_pigeonhole = c(violates_pigeonhole, i)
        }
      }
    }
    trees.for.all[[k]] <- trees1
  }

  # STEP 3: ORDER THE EVENTS ACROSS THE COHORT ####
  events <- c(unique(mergedseg$noevent),max(mergedseg$noevent)+1) # the final event represents WGD

  # create a matrix to store the ordering of the events for each patient; patients in rows, events in columns
  orderingmatrix <- matrix(rep(0), ncol = length(events), nrow = length(tumour.names))

  # create a matrix to match the events to their ordered state
  matrix.coresp <- cbind(sort(events), c(1:length(events)))
  colnames(matrix.coresp) <- c("true","ordered")
  matrix.coresp <- as.data.frame(matrix.coresp)


  # data storage for info from all iterations if a single evolutionary trajectory
  matrix2 <- NULL
  matrixclass <- NULL
  allbics <- NULL
  classification <- NULL


  # if multiple evolutionary trajectories
  patients_list <- list()
  ordering_matrices_list <- list()

  matrix2_bic <- NULL
  matrixclass_bic <- NULL
  allbics_bic <- NULL
  classification_bic <- NULL
  class.ids_bic <- NULL
  class.ids <- NULL


  # get 1000 overall orderings for all patients
  for(r in 1:1000){

    print(paste0('Current iteration: ', r))

    for(i in 1:length(tumour.names)){

      v <- NULL
      vectoroforderedeventssubclonal.initial <- NULL         #  initial vector of subclonal events
      vectoroforderedeventssubclonal <- NULL                 #  vector of subclonal events
      list.events.subclonal <- list()


      mergedsegforsample    <- mergedseg[as.character(mergedseg$Tumour_Name)==as.character(tumour.names[i]), ]     # pick all segments for the current patient
      mergedsegforsample    <- mergedsegforsample[!duplicated(mergedsegforsample),]                                # remove any duplicated events
      vector.of.total.events<- unique(mergedsegforsample$noevent)                                                  # the unique events for this patient

      #take subclonal events separately and order them according to trees (using the build trees output from the first step)
      mergedsegforsamplesubclonal <- mergedsegforsample[mergedsegforsample$w.mean<1,]
      vectorofeventssubclonal     <- mergedsegforsamplesubclonal$noevent

      if (length(vectorofeventssubclonal)>0){
        l <- length(trees.for.all[[i]])
        notree <- sample(1:l, size = 1, TRUE)
        current.tree <- trees.for.all[[i]][[notree]]  # sample one of the trees for the patient
        current.tree <- current.tree[order(current.tree$level,decreasing=FALSE),]

        for(n in 1:nrow(current.tree)){
          parent.id <- current.tree$parent[n]
          if(parent.id==0){
            current.tree$index[1]=100              # sets maximum value 100 for the earliest event
          }
          else if(parent.id>0){
            #for each node on the tree generate a random value that is SMALLER than that for its parent.In this way we make sure that all nodes are ordered after their parents , but make no assumption about the relative timing at branching points
            parent.row <- current.tree[current.tree$nodeid==parent.id,]
            max.value <- parent.row$index
            current.tree$index[n] <- runif(1, min = 0, max = max.value)
          }
        }
        current.tree <- current.tree[order(current.tree$index, decreasing = TRUE),]  # order the tree in decreasing order - i.e. parent at top

        #check whether there are multiple events that occured simultaneously and assigns them a random ordering
        vectoroforderedeventssubclonal.initial <- current.tree$noevent  # sort out the order of events using the built trees


        for(x in 1:length(vectoroforderedeventssubclonal.initial)){
          list.events.subclonal[[x]] <- strsplit(as.character(vectoroforderedeventssubclonal.initial[x]), " ")
        }

        list.ordered.events.subclonal <- list.events.subclonal
        for(n in 1:length(list.events.subclonal)){
          if(length(list.events.subclonal[[n]][[1]])==1){
            list.ordered.events.subclonal[[n]][[1]]=list.events.subclonal[[n]][[1]]
          }else if(length(list.events.subclonal[[n]][[1]])>1){
            list.ordered.events.subclonal[[n]][[1]] = sample(list.events.subclonal[[n]][[1]], size = length(list.events.subclonal[[n]][[1]]), replace = FALSE)
          }else{
            list.ordered.events.subclonal[[n]][[1]] = NULL
          }
        }

        vectoroforderedeventssubclonal <- unlist(list.ordered.events.subclonal)

      } else if (length(vectorofeventssubclonal)==0){
        vectoroforderedeventssubclonal <- NULL
      }

      # analysis of clonal events
      tumour.ploidy = unique(mergedsegforsample$tumour_ploidy)


      if (length(tumour.ploidy)>1){
        print("different ploidies error")     # check just in case; there shouldn't be different ploidies
      }

      if (tumour.ploidy==2){
        #analysis of clonal events in diploid samples - we take all events and then randomly order them (we don't know their actual timing)
        mergedseg.for.sample.clonal<- mergedsegforsample[mergedsegforsample$w.mean==1,]
        vectorofeventsclonal       <- mergedseg.for.sample.clonal$noevent

        if (length(vectorofeventsclonal)==1){
          vector.of.ordered.events.clonal <- vectorofeventsclonal
        } else if (length(vectorofeventsclonal)>1){
          vector.of.ordered.events.clonal <- sample(vectorofeventsclonal, size=length(vectorofeventsclonal), replace = FALSE)
        } else {
          vector.of.ordered.events.clonal <- NULL
        }

      } else if (tumour.ploidy==4){
        #analysis of clonal events in tetraploid samples - we need to take into account the wgd

        mergedsegforsampleWGD <- mergedsegforsample[mergedsegforsample$w.mean==1,]           # these are the clonal events that we want to time relative to WGD
        mergedsegforsampleWGD <- mergedsegforsampleWGD[!duplicated(mergedsegforsampleWGD), ]  # remove any duplicates

        # BEFORE WGD ####

        # HOM DEL events - we're not considering if it has happened before or after WGD but we're assuming for some reason it's postWGD - WHY?????
        mergedsegforsampleHD <- mergedsegforsampleWGD[mergedsegforsampleWGD$CNA=="cHD" | mergedsegforsampleWGD$CNA=="sHD",]
        vector.of.events.HD  <- mergedsegforsampleHD$noevent

        # LOH events - the LOHs that have minor copy number = 0 have occurred before WGD (prostate cancer paper)
        mergedsegforsampleWGD.LOH       <- mergedsegforsampleWGD[mergedsegforsampleWGD$CNA=="cLOH" | mergedsegforsampleWGD$CNA=="sLOH",]
        mergedsegforsamplebeforeWGD.LOH <- mergedsegforsampleWGD.LOH[mergedsegforsampleWGD.LOH$nMin1_A==0,]

        # Gain events - if the major copy number is twice or greater than ploidy, then it is before WGD
        mergedsegforsampleWGD.Gain       <- mergedsegforsampleWGD[mergedsegforsampleWGD$CNA=="cGain" | mergedsegforsampleWGD$CNA=="sGain",]
        mergedsegforsamplebeforeWGD.Gain <- mergedsegforsampleWGD.Gain[mergedsegforsampleWGD.Gain$nMaj1_A>3,]

        # Driver mutations - use no.chrs.bearing.mutation to time the mutations before/after WGD: before ( if no.chrs.bearing.mutation > 1)
        mergedsegforsampleWGD.drivers       <- mergedsegforsampleWGD[mergedsegforsampleWGD$CNA%in%clonal_driver_mutations, ]
        mergedsegforsamplebeforeWGD.drivers <- mergedsegforsampleWGD.drivers[mergedsegforsampleWGD.drivers$no.chrs.bearing.mut>1, ]

        # merge LOH, GAIN and driver pre-WGD events
        mergedsegforsamplebeforeWGD <- rbind(mergedsegforsamplebeforeWGD.LOH, mergedsegforsamplebeforeWGD.Gain, mergedsegforsamplebeforeWGD.drivers)
        mergedsegforsamplebeforeWGD <- mergedsegforsamplebeforeWGD[!duplicated(mergedsegforsamplebeforeWGD), ]

        # randomly order the events before WGD
        vector.of.events.before.WGD <- mergedsegforsamplebeforeWGD$noevent
        if (length(vector.of.events.before.WGD)==1){
          vector.of.ordered.events.before.WGD <- vector.of.events.before.WGD
        } else if (length(vector.of.events.before.WGD)>1){
          vector.of.ordered.events.before.WGD <- sample(vector.of.events.before.WGD, size = length(vector.of.events.before.WGD), replace = FALSE)
        } else {
          vector.of.ordered.events.before.WGD <- NULL   # if no clonal events
        }

        # AFTER WGD  ####
        # LOH events post-WGD
        mergedsegforsampleafterWGD.LOH <- mergedsegforsampleWGD.LOH[mergedsegforsampleWGD.LOH$nMin1_A==1,]

        # Gain events post-WGD
        mergedsegforsampleafterWGD.Gain <- mergedsegforsampleWGD.Gain[mergedsegforsampleWGD.Gain$nMaj1_A<=3,]

        # Driver mutations events post-WGD - all that are on one or fewer copies
        mergedsegforsampleafterWGD.drivers <- mergedsegforsampleWGD.drivers[mergedsegforsampleWGD.drivers$no.chrs.bearing.mut<=1,]  # use no.chrs.bearing.mut

        # merge LOH, Gain and drivers post-WGD events
        mergedsegforsampleafterWGD <- rbind(mergedsegforsampleafterWGD.LOH,  mergedsegforsampleafterWGD.Gain, mergedsegforsampleafterWGD.drivers)
        mergedsegforsampleafterWGD <- mergedsegforsampleafterWGD[!duplicated(mergedsegforsampleafterWGD),]   # remove duplicates

        # randomly order post-wgd events
        vector.of.events.after.WGD <- mergedsegforsampleafterWGD$noevent
        if (length(vector.of.events.after.WGD)==1){
          vector.of.ordered.events.after.WGD <- vector.of.events.after.WGD
        } else if (length(vector.of.events.after.WGD)>1){
          vector.of.ordered.events.after.WGD <- sample(vector.of.events.after.WGD, size = length(vector.of.events.after.WGD), replace = FALSE)
        } else {
          vector.of.ordered.events.after.WGD <- NULL   # if no clonal events
        }

        #we create the final order of clonal events in which events that took place before the WGD, put in random order, are placed before the events which took place after WGD
        vector.of.ordered.events.clonal <- c(vector.of.ordered.events.before.WGD,max(events),vector.of.ordered.events.after.WGD,vector.of.events.HD)
      }

      # if no clonal events
      if (length(vector.of.ordered.events.clonal)==0){
        vector.of.ordered.events.clonal = NULL
      }

      #create the final vector of events in which events that are clonal take place before subclonal events
      vectoroforderedevents <- as.numeric(c(vector.of.ordered.events.clonal, base::setdiff(as.numeric(vectoroforderedeventssubclonal),vector.of.ordered.events.clonal)))


      if (model=='PlackettLuce'){
        map_vector = data.frame(vectoroforderedevents, rank=1:length(vectoroforderedevents))     # set the rankings of the observed events based on the order of occurrence
        map_vector_unobserved = data.frame(vectoroforderedevents = setdiff(events,vectoroforderedevents), rank = 0)  # set the rankings of the unobserved events to 0
        map_vector = rbind(map_vector,map_vector_unobserved)

        map_vector = map_vector[order(map_vector$vectoroforderedevents),]
        orderingmatrix[i,] = map_vector$rank  # the PlackettLuce model takes the rankings of the events as input
      }

      if (model=='PLMIX'){
        # include the unobserved events (we assume they have occurred after the last observed event)
        unobserved_events = base::setdiff(events,vectoroforderedevents)
        random_unobserved_events = sample(unobserved_events, length(unobserved_events), replace = F)

        vectoroforderedevents <- as.numeric(c(vectoroforderedevents, random_unobserved_events))
        vectoroforderedevents <- as.vector(vectoroforderedevents)

        # add the ordering from this iteration to the ordering matrix
        vectoroforderedevents.i <- NULL


        for (b in 1:length(vectoroforderedevents)){
          vectoroforderedevents[b] <- as.numeric(vectoroforderedevents[b])
        }
        for (b in 1:length(vectoroforderedevents)){
          curr.ev <- as.numeric(vectoroforderedevents[b])

          vectoroforderedevents.i[b] <- as.numeric(matrix.coresp[matrix.coresp$true==curr.ev,2])
        }

        for(v in 1:length(vectoroforderedevents)){
          orderingmatrix[i,v] <- vectoroforderedevents.i[v]
        }
      }
    }

    if (mixture_model==F){
      if (model=='PlackettLuce'){
        PL_output <- PlackettLuce(orderingmatrix,verbose = F)
        ordering_vector <- PL_output$coefficients[1:length(events)]
        matrix2 <- rbind(matrix2,ordering_vector)                   # all values of the events' timings so far
      }

      if (model=='PLMIX'){
        set.seed(123456)
        K <- ncol(orderingmatrix)   # number of events
        G <- 1                      # we consider single evolutionary trajectory

        outputMAP <- mapPLMIX_multistart(pi_inv=orderingmatrix, K=K, G=G, n_start=1, n_iter=400*G)  # n_start is the number of different starting points


        ordering_vector <- outputMAP$mod$P_map                      # values of the relative ordering of the events
        matrix2 <- rbind(matrix2,ordering_vector)                   # all values of the events' timings so far
        bic <- outputMAP$mod$bic

        classification <- outputMAP$mod$class_map          # the separation of samples/patients into clusters
        matrixclass <- rbind(matrixclass,classification)   # the separation of samples/patients into clusters from all iterations so far

        allbics <- rbind(bic,allbics)                     # BIC info
      }
    } else {
      K <- ncol(orderingmatrix)   # number of events
      # mixture model
      set.seed(123456)
      # initialisation values
      ordered_matrix <- as.top_ordering(orderingmatrix, format_input='ordering', aggr=F, ties_method = 'random')
      MAP_1 <- mapPLMIX_multistart(pi_inv = ordered_matrix, K = K, G = 1,
                                   n_start=5, n_iter = 400*1)
      MAP_2 <- mapPLMIX_multistart(pi_inv = ordered_matrix, K = K, G = 2,
                                   n_start=5, n_iter = 400*2)
      MAP_3 <- mapPLMIX_multistart(pi_inv = ordered_matrix, K = K, G = 3,
                                   n_start=5, n_iter = 400*3)
      MAP_4 <- mapPLMIX_multistart(pi_inv = ordered_matrix, K = K, G = 4,
                                   n_start=5, n_iter = 400*4)
      MAP_5 <- mapPLMIX_multistart(pi_inv = ordered_matrix, K = K, G = 5,
                                   n_start=5, n_iter = 400*5)

      all_map <- list(MAP_1, MAP_2, MAP_3, MAP_4, MAP_5)

      # using bic to pick a model
      bics <- rbind(MAP_1$mod$bic, MAP_2$mod$bic, MAP_3$mod$bic, MAP_4$mod$bic, MAP_5$mod$bic)
      chosen_model_bic <- which.min(bics)

      # saving the classification for the bic model
      outputMAP_bic <- all_map[[chosen_model_bic]]
      vector_bic <- outputMAP_bic$mod$P_map  # the bic model
      matrix2_bic <- rbind(matrix2_bic,vector_bic)

      classification_bic <- outputMAP_bic$mod$class_map
      matrixclass_bic <- rbind(matrixclass_bic, classification_bic)

      class.ids_bic <- rbind(class.ids_bic, classification_bic)

      # save the name of the patients
      patients_list <- c(patients_list, list(rownames(ordered_matrix)))
      ordering_matrices_list <- c(ordering_matrices_list, list(orderingmatrix))
    }
  }

  # save some relevant files for plotting
  if (mixture_model==F){
    write.table(matrix2, file = paste0(output_dir,tumour_type,"_wgd_matrix2_with_unobserved_events"), sep="\t",quote=F)       # values of how early/late each event is
    write.table(matrix.coresp, file = paste0(output_dir,tumour_type,"_wgd_matrix_coresp_with_unobserved_events"), sep="\t", quote=F) # matrix matching the events
    write.table(mergedseg, file = paste0(output_dir, tumour_type,"_mergedseg.txt"),col.names=T,row.names = F,quote = F,sep="\t")       # save the list of all mergesegs
  } else {
    write.table(matrix2_bic, file = paste0(output_dir,tumour_type,"_wgd_matrix2_mixture_model_with_unobserved_events"), sep="\t",quote=F)       # values of how early/late each event is
    write.table(matrix.coresp, file = paste0(output_dir,tumour_type,"_wgd_matrix_coresp_with_unobserved_events"), sep="\t", quote=F) # matrix matching the events
    write.table(mergedseg, file = paste0(output_dir, tumour_type,"_mergedseg.txt"),col.names=T,row.names = F,quote = F,sep="\t")       # save the list of all mergesegs
    write.table(class.ids_bic, file = paste0(output_dir, tumour_type,"_mixture_model_classification_1000its.txt"), col.names = T, row.names= T, quote = F,sep="\t")   # this has the classificationof each sample from each iteration
    # need to use class.ids_bic and count how many times every 2 patients are in the same cluster and then use hierarchical clustering to cluster the sample
  }



  # STEP 4: PLOT THE TIMING OF CNA EVENTS ####
 if (model=='PLMIX'){

  #1
  mergedseg$ID = paste(paste0("chr",mergedseg$chr), mergedseg$startpos, mergedseg$endpos, sep = "_")
  event.id = unique(mergedseg[ ,c("ID","noevent","CNA")])
  event.id = rbind(data.frame(ID = "Whole_Genome_Duplication", noevent = max(events), CNA = "WGD"), event.id)

  #2
  # load matrix.coresp from previous step
  names(matrix.coresp) = c("noevent","event")
  matrix.coresp$event = paste("p", matrix.coresp$event, sep = "_")
  #merge
  coresp = merge(event.id, matrix.coresp, "noevent")

  #3
  values <- melt(matrix2, id.vars = NULL)
  values$Var1 = NULL
  names(values) = c("event","value")
  out <- merge(values, coresp, "event")
  out$CNA = substr(out$CNA,2,5)
  write.table(out, paste0(output_dir, tumour_type ,"_wgd_ordering_plot_file"), col.names = T, row.names = F,quote = F, sep = "\t")
  #plot
  timing_plot <- ggplot(out, aes(x = ID, y = -value, col = CNA, fill = CNA)) +
    geom_violin() +
    coord_flip() +
    labs(x="CNA event",y="Timing scale") +
    ggtitle(paste0(tumour_type, " ordering of CNA events"))

  # save plot
  save_plot(paste0(output_dir, tumour_type,"_CNA_timing_order.pdf"), plot = timing_plot)
 }

  if (model=='PlackettLuce'){
    #1
    mergedseg$ID = paste(paste0("chr",mergedseg$chr), mergedseg$startpos, mergedseg$endpos, sep = "_")
    event.id = unique(mergedseg[,c("ID","noevent","CNA")])
    event.id = rbind(data.frame(ID = "Whole_Genome_Duplication", noevent = (length(unique(mergedseg$ID))+1), CNA = "WGD", stringsAsFactors = F), event.id)
    # Additional filtering # event.id should have same nrow as events
    event.id$cna = substring(event.id$CNA,first = 2)
    event.id = unique(event.id[,c("ID","noevent","cna")])

    #2
    names(matrix.coresp) = c("noevent","event")
    write.table(matrix.coresp, file = paste0(output_dir,tumour_type,"_wgd_matrix_coresp_with_unobserved_events"), sep="\t", quote=F) # matrix matching the events

    #merge
    coresp = merge(event.id,matrix.coresp,"noevent")
    #3
    values <- melt(matrix2, id.vars = NULL)
    values$Var1 = NULL
    names(values) = c("event","value")
    out <- merge(values,coresp,"event")

    PLOTDF = data.frame()
    unoevent = unique(out$noevent)
    for (i in 1:length(unique(out$noevent))){
      event = out[out$noevent==unoevent[i],]
      event = event[order(event$value),]
      PLOTDF = rbind(PLOTDF,event)
    }

    PLOTDF$ID = factor(PLOTDF$ID,levels=mixedsort(unique(PLOTDF$ID)))

    plotDF = ggplot(PLOTDF, aes(x = ID, y = -log10(value),fill=cna)) + geom_violin(aes(col=cna)) + coord_flip() + labs(x="CNA event",y="Timing scale")+
      ggtitle(paste0(tumour_type, " ordering of genomic events - 95% CI")) + labs(fill="Event")
    g = ggplot_build(plotDF)
    fillcol = unique(g$data[[1]]["fill"])$fill
    save_plot(paste0(output_dir, tumour_type,"_CNA_timing_order.pdf"), plot = plotDF)

  }


}




