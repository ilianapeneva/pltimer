#' Collate the subclones data from Battenberg into a single file and add the ploidy of the sample.
#'
#' @param data_dir Full path to the directory with the subclones files
#' @param tumour_type String that represents the type of tumour we work with
#' @param output_dir Full path to the directory where the allsegs.txt file output from this function will be stored
#' @return \emph{tumour_type}_allsegs.txt file with information about the samplenames, copy number aberrations (startpos, endpos, copy number state of major/minor allele), ploidy

subclone_collation <- function(data_dir,tumour_type, output_dir) {
  samples = list.files(path = data_dir,pattern = 'subclones.txt')

  allsegs <- NULL # this will store the collated info from the subclones
  for (i in 1:length(samples)) {
    name = sub("_subclones.txt","", samples[i])                                                    # pick the tumour platekey
    cndata <- read.table(paste0(data_dir, samples[i]), header = T, sep = '\t', stringsAsFactors = F) # load the subclones file for this patient
    cndata <- cndata[!is.na(cndata$nMaj1_A),]                        # remove any segments with major allele copy number of state 1 being NA
    cndata <- cbind(name,cndata)                      # add the tumour platekey to each of the segments of the subclone file

    # compute copy number for each of the segments
    totalcn <- rep(0, nrow(cndata))
    for (k in 1:nrow(cndata)){
      if (cndata$frac1_A[k]==1) {
        totalcn[k] <- cndata$nMaj1_A[k] + cndata$nMin1_A[k] }
      else if (cndata$frac1_A[k]!="Inf"&cndata$frac1_A[k]!="-Inf") {
        totalcn[k] <- (cndata$nMaj1_A[k] + cndata$nMin1_A[k])*cndata$frac1_A[k] + (cndata$nMaj2_A[k] + cndata$nMin2_A[k])*cndata$frac2_A[k]
      }
      cndataweighted <- totalcn*(cndata$endpos - cndata$startpos)/(sum(as.numeric(cndata$endpos-cndata$startpos)))
    }

    # compute ploidy of the sample using the copy number - the sample can be either diploid or tetraploid
    ploidy <- sum(cndataweighted)
    roundedploidy <- round(ploidy/2)*2

    # pick all the info from this sample needed to annotate with the type of CNA
    allsegsi <- data.frame(cndata[,1:4], cndata[,9:15], roundedploidy)
    names(allsegsi) <- c("Tumour_Name","chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A",
                         "SDfrac_A", "tumour_ploidy")

    # add this to the info from the other samples
    allsegs <- rbind(allsegs, allsegsi)
    print(paste(samples[i],"done"))
  }
  # save the info
  write.table(allsegs, paste0(output_dir,tumour_type,"_allsegs.txt"), sep="\t", row.names = F, quote=F)
}


#' Collate the subclones data from Battenberg into a single file and add the ploidy of the sample.
#'
#' @param allsegs_file Full path to the directory with the allsegs file with the subclones information
#' @param tumour_type String that represents the type of tumour we work with
#' @param output_dir Full path to the directory with the annotated allsegments.txt output from this function
#' @return \emph{tumour_type}_allsegments.txt file with information about the samplenames, copy number aberrations (startpos, endpos, copy number state of major/minor allele), ploidy and type of CNA


CNA_annotation <- function(allsegs_file, tumour_type, output_dir){
  column_names <- c("Tumour_Name",  "chr", "startpos", "endpos","nMaj1_A","nMin1_A", "frac1_A", "SDfrac_A", "tumour_ploidy") # names of the columns we want to keep from allsegs file

  # load the data with the collated subclones info
  allsegs <- read.table(allsegs_file, header = T, sep = '\t', stringsAsFactors = F)
  total_segs = nrow(allsegs)


  allsegsa <- NULL    # this stores the segments from the subclones files
  CNA <- NULL         # this stores the type of the CNA (subclonal/clonal HD, LOH, Gain, Amp, Loss)

  for (i in 1:nrow(allsegs)) {
    # annotate clonal events first
    if (is.na(allsegs$nMaj2_A[i])) {
      if (allsegs$nMin1_A[i] == 0 & allsegs$nMaj1_A[i] == 0) {    # clonal homozygous deletion - both copies are lost
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cHD")
      }
      if (allsegs$nMin1_A[i] == 0 | allsegs$nMaj1_A[i] == 0) {    # clonal LOH - one of the copies is lost
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cLOH")
      }
      if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]/2) |
          (allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]/2)) {
        if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]*2) |
            (allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]*2)) {
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
          CNA <- c(CNA, "cAmp")                                   # clonal amplification - if the number of copies of one of the alleles is more than 2* tumour ploidy
        }
        else {
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])    # clonal gain - whatever there is left
          CNA <- c(CNA, "cGain")
        }
      }
      if ((allsegs$nMaj1_A[i] == allsegs$tumour_ploidy[i]/2) &
          (allsegs$nMin1_A[i] == allsegs$tumour_ploidy[i]/2)) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "NoCNV")                                    # no CNA if there isn't any change in the number of copies for the alleles
      }
      if ((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2) |
          (allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2)) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cLoss")                                   # clonal loss - if the number of the copies of one of the alleles is less than the expected number for the type of sample
      }
    }
    # annotate now the subclonal events
    if (!is.na(allsegs$nMaj2_A[i])) {
      if (allsegs$nMin1_A[i] == 0 & allsegs$nMaj1_A[i] == 0) {  # if both copies are lost - subclonal homozygous deletion
        CNA <- c(CNA, "sHD")
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
      }
      if ((allsegs$nMin1_A[i] == 0 | allsegs$nMaj1_A[i] == 0) &     # if all of the copies of one of the alleles are lost in the clonal/subclonal cluster and then in the other subclonal cluster then it must be clonal LOH
          (allsegs$nMin2_A[i] == 0 | allsegs$nMaj2_A[i] == 0)) {
        CNA <- c(CNA, "cLOH")
        tmp <- allsegs[i,c(1:4,8:12)]
        names(tmp) <- column_names
        allsegsa <- rbind(allsegsa, tmp[,column_names])
      }
      else if (allsegs$nMin1_A[i] == 0 | allsegs$nMaj1_A[i] == 0) {   #subclonal LOH  - if all of the copies of one of the alleles are lost
        CNA <- c(CNA, "sLOH")
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
      }
      if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]/2 |      # if # copies of one of the alleles in the subclones are higher than 2*tumour_ploidy - clonal amplification
           allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]/2) &
          (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]/2 |
           allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]/2)) {
        if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]*2 |
             allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]*2) &
            (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]*2 |
             allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]*2)) {
          CNA <- c(CNA, "cAmp")
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        }
        else {
          CNA <- c(CNA, "cGain")                          # if #copies of one of the alleles in the subclones are higher than 1/2 * tumour_ploiody (i.e. the expected # copies given the ploidy) but less than 2*tumour ploidy
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        }
      }
      if ((allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])|
          (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]/2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])) {
        if ((allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]*2 ) |   # & allsegs$nMin1_A[i] != allsegs$nMin2_A[i] - already included in the previous if
            (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]*2 )) {  # & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i] - already included in the previous if
          CNA <- c(CNA, "sAmp")                        # if # copies of the subclone are higher than expected for the ploidy and not equal to the 1st subclone
          tmp <- allsegs[i,c(1:4,8:12)]
          names(tmp) <- column_names
          allsegsa <- rbind(allsegsa, tmp[,column_names])
        }
        else {
          CNA <- c(CNA, "sGain")
          tmp <- allsegs[i,c(1:4,8:12)]
          names(tmp) <- column_names
          allsegsa <- rbind(allsegsa, tmp[,column_names])
        }
      }
      if ((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2 |           # if some of the copies are lost in both subclones - clonal loss
           allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2) &
          (allsegs$nMaj2_A[i] < allsegs$tumour_ploidy[i]/2 |
           allsegs$nMin2_A[i] < allsegs$tumour_ploidy[i]/2)) {
        CNA <- c(CNA, "cLoss")
        tmp <- allsegs[i,c(1:4,8:12)]
        names(tmp) <- column_names
        allsegsa <- rbind(allsegsa, tmp[,column_names])
      }
      if (((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i]) |          # if some of the copies are lost in one of the subclones
           (allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i]))) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "sLoss")
      }
    }
    if(i %% 100 ==0){
      print(paste('Segment ',i,"/",total_segs, ' annotated with the type of CNA'))
    }
  }

  # add the CNA type to the subclones info
  allsegsa <- cbind(allsegsa, CNA)
  # set frac1_A to 1 for all clonal events
  allsegsa$frac1_A[allsegsa$CNA == "cGain"|allsegsa$CNA == "cAmp"|allsegsa$CNA == "cLOH"|allsegsa$CNA == "cHD"|allsegsa$CNA == "cLoss"] <- 1
  write.table(allsegsa, file = paste0(output_dir, tumour_type,"_annotated_segments.txt"),sep = "\t",quote = F,row.names = F)
}


#' Prepare data for plotting of landscape of CNAs across the whole genome.
#'
#' @param annotated_segments_file Full path to the file with the aggregated subclones information, which has been annotated with the CNA type
#' @param tumour_type String that represents the type of tumour we work with
#' @param genome_build String that represents the genome build (could be hg19 or hg38)
#' @param output_dir Full path to the directory with the refsegs files output from this function
#' @return files with information about which positions of the genome are aberrated across the whole cohort. There will be a separate file for all/clonal/subclonal LOH/gain/HD.

prepare_data_for_landscape <- function(annotated_segments_file, tumour_type, genome_build, output_dir){
   if (genome_build=='hg19'){
  hg_genome = BSgenome.Hsapiens.UCSC.hg19
  chr_lengths = seqlengths(hg_genome)
  chr_df = data.frame(chr = as.numeric(1:23),length = as.numeric(chr_lengths[1:23])) # Chromosome lengths with the corresponding chr
  maxchr = chr_df$length
 } else if (genome_build=='hg38'){
  hg_genome = BSgenome.Hsapiens.UCSC.hg38
  chr_lengths = seqlengths(hg_genome)
  chr_df = data.frame(chr = as.numeric(1:23),length = as.numeric(chr_lengths[1:23])) # Chromosome lengths with the corresponding chr
  maxchr = chr_df$length
 } else {
  print('Genome build should be hg19 or hg38')
 }

  # the CNAs we're interested to collect data for the landscape plot
  CNAs <- c("LOH","HD", "Gain","cLOH","cHD", "cGain","sLOH","sHD", "sGain")

  # load the subclones data annotated with the type of CNA
  copy.number.data <- read.table(annotated_segments_file, header = T, stringsAsFactors = F)
  copy.number.data <- copy.number.data[copy.number.data$frac1_A>0 & copy.number.data$frac1_A<=1 & copy.number.data$tumour_ploidy>0,]

  # loop over the overall/clonal/subclonal LOH, gains, HDs and count how many segments (i.e. samples) have the aberration
  for (id in 1:9){
    print(paste("Aggregating data for ",CNAs[id]))

    tmprss.data <- copy.number.data                                # create a temporary copy of the data
    current_CNA = CNAs[id]                                         # the type of CNA we're considering at this iteration
    input.data <- tmprss.data[grep(current_CNA,tmprss.data$CNA),]  # pick the data for this CNA type
    chr = c(1:22,"X")

    segs_data <- NULL   # this will store the chr, starting position, ending position and counter of samples with LOH/gain/HD in this interval
    for (i in c(1:22,"X")) {
      # pick the data for this chromosome (this will have the CNA info for all patients)
      data.chr <- input.data[input.data$chr==i,]

      if (nrow(data.chr) > 0) {
        # create data.frames with the startpos and endpos
        sp <- data.frame("loc" = data.chr$startpos, "val" = +1)
        ep <- data.frame("loc" = data.chr$endpos,"val" = -1)

        # put together all start and end positions and order them in increasing order
        locs <- rbind(sp,ep)
        locs <- locs[order(locs$loc),]

        # use the starting/ending positions to create intervals and count how many patients have the interval aberrated
        locsschr <- data.frame("chr" = i, "sp" = 0, "ep" = locs[1,1], "val" = 0)    # sp = starting position, ep = ending position, val = number of samples with the segment
        for (k in 2:nrow(locs)) {
          locsschri <- c(i, locs[k-1,1], locs[k,1], as.numeric(locsschr[k-1, 4])+locs[k-1,2])
          locsschr <- rbind(locsschr, locsschri)
        }

        locsschr <- rbind(locsschr, c(i,locs[k,1],maxchr[which(chr%in%i)],0))
        segs_data <- rbind(segs_data, locsschr)
      }

      if (nrow(data.chr) == 0) {
        locsschr <- data.frame("chr" = i, "sp"= 0, "ep" = maxchr[which(chr%in%i)], "val" = 0)
        segs_data <- rbind(segs_data, locsschr)
      }
    }
    segs_data <- segs_data[segs_data$sp!=segs_data$ep,]

    # save the landscape data for this CNA type
    write.table(segs_data, file = paste0(output_dir, tumour_type,"_", current_CNA,"_refsegs.txt"), sep="\t", quote=F)
  }
}


#' Plot the CNA data across the genome (all/clonal/subclonal aberrations)
#'
#' @param annotated_segments_file Full path to the aggregated subclones file, which has been annotated with the CNA type
#' @param refsegs_dir Full path to the directory with the refsegs files
#' @param tumour_type String that represents the type of tumour we work with
#' @param genome_build String that represents the genome build (could be hg19 or hg38)
#' @param output_dir Full path to the directory where the landscape plots will be saved
#' @return a pdf file with the overall landscape of all/clonal/subclonal events

plot_CN_landscape <- function(annotated_segments_file, refsegs_dir,  tumour_type, genome_build, output_dir){
  copy.number.data <- read.table(annotated_segments_file, header = T, stringsAsFactors = F)
  copy.number.data <- copy.number.data[copy.number.data$frac1_A>0 & copy.number.data$frac1_A<=1 & copy.number.data$tumour_ploidy>0,]
  samples = unique(copy.number.data$Tumour_Name) # number of samples in the cohort
  
  
 if (genome_build=='hg19'){
  hg_genome = BSgenome.Hsapiens.UCSC.hg19
  chr_lengths = seqlengths(hg_genome)
  chr_df = data.frame(chr = as.numeric(1:23),length = as.numeric(chr_lengths[1:23])) # Chromosome lengths with the corresponding chr
  maxchr = chr_df$length
 } else if (genome_build=='hg38'){
  hg_genome = BSgenome.Hsapiens.UCSC.hg38
  chr_lengths = seqlengths(hg_genome)
  chr_df = data.frame(chr = as.numeric(1:23),length = as.numeric(chr_lengths[1:23])) # Chromosome lengths with the corresponding chr
  maxchr = chr_df$length
 } else {
  print('Genome build should be hg19 or hg38')
 }
 
  # get the end positions for the chromosomes on the whole genome (i.e chr1 endpos1, chr2 endpos2=(endpos1+length of chr 2)...)
  loci <- 0
  for (i in 1:23) {
    loci[i+1] <- loci[i] + maxchr[i]
  }

  chrg<- NULL
  for (i in 2:length(loci)) {
    chrg<- c(chrg, (loci[i]+loci[i-1])/2)
  }

  # set up some parameters for the plotting of all/clonal/subclonal events
  chr       <- c(1:22,"X")
  CNAs      <- c("LOH","HD", "Gain")
  sign      <- c(-1,-1,1)
  clonality <- c("All", "Clonal", "Subclonal")
  prop      <- c("","c","s")
  cols      <- c("dodgerblue","black","red")   # LOH in blue, HD in black, gain in red

  # prepare the layout for the plots
  m <- matrix(1:3, nrow = 3, ncol = 1)
  layout(m)
  ln <- length(samples)  # number of samples

  pdf(paste0(output_dir, tumour_type, '_lansdscape.pdf'))
  for (j in 1:3) {
    plot(NA, xlim = c(0,max(loci)), ylim = c(-1, 1), main = paste(tumour_type, "-",clonality[j],"copy number profile"),
         xaxt = "n", xlab = "", ylab = "Fraction of Tumours",bty = "n")
    for (current_CNA in CNAs) {
      print(paste("plotting data for ",tumour_type, current_CNA, clonality[j]))
      data_for_current_CNA <- read.table(paste0(refsegs_dir, tumour_type,"_", prop[j], current_CNA,"_refsegs.txt"), sep="\t", header=T)

      x <- data_for_current_CNA$sp[1] + loci[which(chr%in%data_for_current_CNA[1,1])]   # starting x coord = 0 + the position of the chr on the genome from the 1st interval this CNA is observed
      y <- 0                                                                          # number of events
      for (i in 1:(nrow(data_for_current_CNA)-1)) {
        if (data_for_current_CNA$ep[i]!=data_for_current_CNA$sp[i+1]) {   # that's true when we move to a new chromosome (there might be +1 for the 3rd x)
          x <- c(x, data_for_current_CNA$ep[i] + loci[which(chr%in%data_for_current_CNA[i,1])], data_for_current_CNA$sp[i+1] + loci[which(chr%in%data_for_current_CNA[i,1])])
          y <- c(y, 0, 0)
        } else{
          x <- c(x, data_for_current_CNA$sp[i] + loci[which(chr%in%data_for_current_CNA[i,1])], data_for_current_CNA$ep[i] + loci[which(chr%in%data_for_current_CNA[i,1])])
          y <- c(y, data_for_current_CNA$val[i], data_for_current_CNA$val[i])
        }
      }
      # ending the data we need to plot for this CNA
      x <- c(x,data_for_current_CNA$ep[i+1] + loci[which(chr%in%data_for_current_CNA[i+1,1])])
      y <- c(y,0)

      # plot the data
      polygon(x,sign[which(CNAs%in%current_CNA)]*y/ln,col = "white", border = cols[which(CNAs%in%current_CNA)])
      text(x = chrg, y = rep(1, 22), labels = chr, cex = 0.8)   # add labels for the chromosomes
      abline(v = loci, lty = 2)                                 # add vertical lines to separate the chromosomes
      legend(0,-0.6, lty=c(1,1), lwd = c(1.5,1.5), paste(clonality[j],CNAs),cex = 0.5, col = cols, bg = "white")

    }
  }
  dev.off()

}
