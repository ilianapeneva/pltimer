#' Simulate random LOH, gains and HD to identify enriched events
#'
#' @param annotated_segments_file Full path to the aggregated subclones file, which has been annotated with the CNA type
#' @param refsegs_dir Full path to the directory with the refsegs files
#' @param gain_output_dir Full path to the directory where the simulations of the gains will be saved
#' @param loh_output_dir Full path to the directory where the simulations of the losses of heterozygosity will be saved
#' @param hd_output_dir Full path to the directoey where the simulations of the homozygous deletions will be saved
#' @param tumour_type String that represents the type of tumour we work with
#' @param genome_build String that represents the genome build - could be hg19 or hg38
#' @param run Number of simulation
#' @return files with the random placing of the LOH/gain/HD events


identify_enriched_regions <- function(annotated_segments_file, refsegs_dir, gain_output_dir, loh_output_dir, hd_output_dir, tumour_type,genome_build, run){

  # GAIN ####
  # load the annotated segments data
  allsegsa <- read.table(annotated_segments_file, header = T, stringsAsFactors = F, sep = "\t")
  all.data <- allsegsa[allsegsa$frac1_A>0&allsegsa$frac1_A<=1&allsegsa$tumour_ploidy>0, ]


  gains.data = all.data[(all.data$CNA=="sGain"|all.data$CNA=="cGain"|all.data$CNA=="sAmp"|all.data$CNA=="cAmp"),]   # pick the segments with gains/amplifications
  gains.segment.counts = read.table(paste0(refsegs_dir,tumour_type,"_Gain_refsegs.txt"), header = T, stringsAsFactors = F)                  # pick the info about how many samples have gain in specific intervals

  colnames(gains.segment.counts) <- c("chr", "startpos", "endpos", "count")

  chr.names = c(1:22,"X")
  chr.data = list(0)              # store the gains.segment.counts data
  chr.seg.counts = array(0,23)    # number of segments per chromosome
  for(chr in 1:23){
    chr.data[[chr]] = gains.segment.counts[gains.segment.counts$chr==chr.names[chr],]
    chr.seg.counts[chr] = nrow(chr.data[[chr]])
  }
  cum.chr.seg.counts = c(0,cumsum(chr.seg.counts))   # sort out genome and chr info in relation with the length of the chromosome and genome
  
 if (genome_build=='hg19'){
   chr.lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23]
 } else if (genome_build=='hg38'){
   chr.lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:23]
 } else {
  print('genome build should be hg19 or hg38')
 }
  #chr.lengths <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895)
  genome.length = sum(as.numeric(chr.lengths))
  half.genome.length = round(genome.length/2)
  chr.boundaries = c(0,cumsum(as.numeric(chr.lengths)))

  # create a vector to store the expected number of times an event occurring by chance
  sim.counts = vector(mode = "numeric",length = nrow(gains.segment.counts))

  gains.data$gain.length = gains.data$endpos - gains.data$startpos + 1   # lengths of the gain segments

  sim.gains = NULL  # where we will store the simulated data
  for(sample in unique(gains.data$Tumour_Name)){
    sample.sim.gains = NULL
    sample.data = gains.data[gains.data$Tumour_Name==sample,]          # pick the gains/ampl data for this sample

    for(i in 1:nrow(sample.data)){
      #April 2014 - allow deletion to cross over chromosome boundary, otherwise telomeres are never hit
      genome.start.pos = sample(half.genome.length,1)
      if(sample(2,1)==2){
        genome.start.pos = genome.start.pos + half.genome.length
      }
      genome.end.pos = genome.start.pos+sample.data$gain.length[i]-1
      chr = sum(genome.start.pos>chr.boundaries)
      chr.end = sum(genome.end.pos>chr.boundaries)
      if (chr == chr.end){
        chr.start.pos = genome.start.pos - chr.boundaries[chr]
        chr.end.pos = genome.end.pos - chr.boundaries[chr]
        sample.sim.gains = rbind(sample.sim.gains,c(chr,chr.start.pos,chr.end.pos))
      } else{
        #if the aberration overlaps chromosome boundary, add 2 separate aberrations, one on each chromosome
        chr.start.pos = genome.start.pos - chr.boundaries[chr]
        chr.end.pos = chr.lengths[chr]
        sample.sim.gains = rbind(sample.sim.gains,c(chr,chr.start.pos,chr.end.pos))
        if ((chr.end-chr)>1) {
          for (j in (chr+1):(chr.end-1)) {
            chr.start.pos = 1
            chr.end.pos = chr.lengths[j]
            sample.sim.gains = rbind(sample.sim.gains,c(j,chr.start.pos,chr.end.pos))
          }
        }
        chr.start.pos = 1
        chr.end.pos = genome.end.pos - chr.boundaries[chr.end]
        if (chr.end==24){chr.end = 1} # loop from chr 23 to chr 1
        sample.sim.gains = rbind(sample.sim.gains,c(chr.end,chr.start.pos,chr.end.pos))
      }
    }
    sim.gains = rbind(sim.gains,sample.sim.gains)
  }

  sim.hits = vector(mode = "numeric",length = nrow(gains.segment.counts))
  for (c in 1:23){
    if (dim(gains.segment.counts[gains.segment.counts$chr==chr.names[c],])[1]) {
      chr = chr.names[c]
      chr.sim.data = data.frame(startpos = sim.gains[sim.gains[,1]==c,2],endpos = sim.gains[sim.gains[,1]==c,3])
      all.breakpoints = c(1,sort(unique(unlist(chr.sim.data))),chr.lengths[c])
      no.breakpoints = length(all.breakpoints)-1
      chr.segment.data = cbind(chr,all.breakpoints[1:no.breakpoints],all.breakpoints[-1],0)
      chr.segment.data[,4] = sapply(1:nrow(chr.segment.data),function(d,s,i){sum(d$startpos >= as.numeric(s[i,2]) & d$startpos < as.numeric(s[i,3])| d$endpos > as.numeric(s[i,2]) & d$endpos <= as.numeric(s[i,3]) | d$startpos < as.numeric(s[i,2]) & d$endpos > as.numeric(s[i,3]))},d=chr.sim.data,s=chr.segment.data)
      chr.segment.data[1:(no.breakpoints-1),3] = as.numeric(chr.segment.data[1:(no.breakpoints-1),3])-1
      colnames(chr.segment.data) = c("chr","startpos","endpos","count")
      chr.segment.data = data.frame(chr = chr.segment.data[,1],startpos = as.numeric(chr.segment.data[,2]),endpos = as.numeric(chr.segment.data[,3]),count = as.numeric(chr.segment.data[,4]))
      #segment.data = rbind(segment.data,chr.segment.data)
      no.segments = chr.seg.counts[c]
      for (s in 1:no.segments){
        sim.hits[cum.chr.seg.counts[c]+s] = sum(sapply(1:nrow(chr.segment.data),function(c2,d,l,i){max(min(c2$endpos[s],d$endpos[i])-max(c2$startpos[s],d$startpos[i])+1,0)*d$count[i]/l},d=chr.segment.data,c2=chr.data[[c]],l=chr.data[[c]]$endpos[s] - chr.data[[c]]$startpos[s]+1))
        if (sim.hits[cum.chr.seg.counts[c]+s]>=chr.data[[c]]$count[s]){
          sim.counts[cum.chr.seg.counts[c]+s] = sim.counts[cum.chr.seg.counts[c]+s] + 1
        }
      }
    }
  }
  write.table(sim.counts,paste(gain_output_dir,"sim_counts_",tumour_type,"_Gain_",run,".txt",sep=""),sep="\t",row.names=F,quote=F)


  # LOH ####


  loh.data = all.data[(all.data$CNA=="sLOH"|all.data$CNA=="cLOH"),]   # pick the segments with LOH
  loh.segment.counts = read.table(paste0(refsegs_dir,tumour_type,"_LOH_refsegs.txt"),header=T,stringsAsFactors=F)                  # pick the info about how many samples have loss in specific intervals

  colnames(loh.segment.counts) <- c("chr", "startpos", "endpos", "count")

  chr.names = c(1:22,"X")
  chr.data = list(0)              # store the LOH.segment.counts data
  chr.seg.counts = array(0,23)  # number of segments per chromosome
  for(chr in 1:23){
    chr.data[[chr]] = loh.segment.counts[loh.segment.counts$chr==chr.names[chr],]
    chr.seg.counts[chr] = nrow(chr.data[[chr]])
  }
  cum.chr.seg.counts = c(0,cumsum(chr.seg.counts))   # sort out some genome and chr info
  #chr.lengths<-c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895)
  genome.length = sum(as.numeric(chr.lengths))
  half.genome.length = round(genome.length/2)
  chr.boundaries = c(0,cumsum(as.numeric(chr.lengths)))

  # create a vector to store the expected number of times an event occurring by chance
  sim.counts = vector(mode = "numeric",length = nrow(loh.segment.counts))

  loh.data$loh.length = loh.data$endpos - loh.data$startpos + 1   # lengths of the loh segments

  sim.loh = NULL  # we will store the simulated LOH events here
  for(sample in unique(loh.data$Tumour_Name)){
    sample.sim.loh = NULL
    sample.data = loh.data[loh.data$Tumour_Name==sample,]          # pick the loh/loss data for this sample

    for(i in 1:nrow(sample.data)){
      #April 2014 - allow deletion to cross over chromosome boundary, otherwise telomeres are never hit
      genome.start.pos = sample(half.genome.length,1)
      if(sample(2,1)==2){
        genome.start.pos = genome.start.pos + half.genome.length
      }
      genome.end.pos = genome.start.pos+sample.data$loh.length[i]-1
      chr = sum(genome.start.pos>chr.boundaries)
      chr.end = sum(genome.end.pos>chr.boundaries)
      if (chr == chr.end){
        chr.start.pos = genome.start.pos - chr.boundaries[chr]
        chr.end.pos = genome.end.pos - chr.boundaries[chr]
        sample.sim.loh = rbind(sample.sim.loh,c(chr,chr.start.pos,chr.end.pos))
      } else {
        #if deletion overlaps chromosome boundary, add 2 separate deletions, one on each chromosome
        chr.start.pos = genome.start.pos - chr.boundaries[chr]
        chr.end.pos = chr.lengths[chr]
        sample.sim.loh = rbind(sample.sim.loh,c(chr,chr.start.pos,chr.end.pos))
        if ((chr.end-chr)>1) {
          for (j in (chr+1):(chr.end-1)) {
            chr.start.pos=1
            chr.end.pos=chr.lengths[j]
            sample.sim.loh = rbind(sample.sim.loh,c(j,chr.start.pos,chr.end.pos))
          }
        }
        chr.start.pos = 1
        chr.end.pos = genome.end.pos - chr.boundaries[chr.end]
        if (chr.end==24){chr.end=1} # loop from chr 23 to chr 1
        sample.sim.loh = rbind(sample.sim.loh,c(chr.end,chr.start.pos,chr.end.pos))
      }
    }
    sim.loh = rbind(sim.loh,sample.sim.loh)
  }

  sim.hits = vector(mode="numeric",length = nrow(loh.segment.counts))
  for(c in 1:23){
    if (dim(loh.segment.counts[loh.segment.counts$chr==chr.names[c],])[1]) {
      chr = chr.names[c]
      chr.sim.data = data.frame(startpos = sim.loh[sim.loh[,1]==c,2],endpos = sim.loh[sim.loh[,1]==c,3])
      all.breakpoints = c(1,sort(unique(unlist(chr.sim.data))),chr.lengths[c])
      no.breakpoints = length(all.breakpoints)-1
      chr.segment.data = cbind(chr,all.breakpoints[1:no.breakpoints],all.breakpoints[-1],0)
      chr.segment.data[,4] = sapply(1:nrow(chr.segment.data),function(d,s,i){sum(d$startpos >= as.numeric(s[i,2]) & d$startpos < as.numeric(s[i,3])| d$endpos > as.numeric(s[i,2]) & d$endpos <= as.numeric(s[i,3]) | d$startpos < as.numeric(s[i,2]) & d$endpos > as.numeric(s[i,3]))},d=chr.sim.data,s=chr.segment.data)
      chr.segment.data[1:(no.breakpoints-1),3] = as.numeric(chr.segment.data[1:(no.breakpoints-1),3])-1
      colnames(chr.segment.data) = c("chr","startpos","endpos","count")
      chr.segment.data = data.frame(chr=chr.segment.data[,1],startpos=as.numeric(chr.segment.data[,2]),endpos=as.numeric(chr.segment.data[,3]),count=as.numeric(chr.segment.data[,4]))
      #segment.data = rbind(segment.data,chr.segment.data)
      no.segments = chr.seg.counts[c]
      for(s in 1:no.segments){
        sim.hits[cum.chr.seg.counts[c]+s] = sum(sapply(1:nrow(chr.segment.data),function(c2,d,l,i){max(min(c2$endpos[s],d$endpos[i])-max(c2$startpos[s],d$startpos[i])+1,0)*d$count[i]/l},d=chr.segment.data,c2=chr.data[[c]],l=chr.data[[c]]$endpos[s] - chr.data[[c]]$startpos[s]+1))
        if(sim.hits[cum.chr.seg.counts[c]+s]>=chr.data[[c]]$count[s]){
          sim.counts[cum.chr.seg.counts[c]+s] = sim.counts[cum.chr.seg.counts[c]+s] + 1
        }
      }
    }
  }
  write.table(sim.counts,paste(loh_output_dir,"sim_counts_",tumour_type,"_LOH_",run,".txt",sep=""),sep="\t",row.names=F,quote=F)



  # HD ####

  hd.data = all.data[(all.data$CNA=="sHD"|all.data$CNA=="cHD"),]   # pick the segments with HD
  hd.segment.counts = read.table(paste0(refsegs_dir, tumour_type ,"_HD_refsegs.txt"),header=T,stringsAsFactors=F)                  # pick the info about how many samples have loss in specific intervals

  colnames(hd.segment.counts) <- c("chr", "startpos", "endpos", "count")

  chr.names = c(1:22,"X")
  chr.data = list(0)              # store the hd.segment.counts data
  chr.seg.counts = array(0,23)  # number of segments per chromosome
  for(chr in 1:23){
    chr.data[[chr]] = hd.segment.counts[hd.segment.counts$chr==chr.names[chr],]
    chr.seg.counts[chr] = nrow(chr.data[[chr]])
  }
  cum.chr.seg.counts = c(0,cumsum(chr.seg.counts))   # sort out some genome and chr info
  #chr.lengths < -c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895)
  genome.length = sum(as.numeric(chr.lengths))
  half.genome.length = round(genome.length/2)
  chr.boundaries = c(0,cumsum(as.numeric(chr.lengths)))

  # create a vector to store the expected number of times an event occurring by chance
  sim.counts = vector(mode = "numeric",length = nrow(hd.segment.counts))

  hd.data$loh.length = hd.data$endpos - hd.data$startpos + 1   # lengths of the hd segments

  sim.hd = NULL
  for(sample in unique(hd.data$Tumour_Name)){
    sample.sim.hd = NULL
    sample.data = hd.data[hd.data$Tumour_Name==sample,]          # pick the HD data for this sample

    for(i in 1:nrow(sample.data)){
      #April 2014 - allow deletion to cross over chromosome boundary, otherwise telomeres are never hit
      genome.start.pos = sample(half.genome.length,1)
      if(sample(2,1)==2){
        genome.start.pos = genome.start.pos + half.genome.length
      }
      genome.end.pos = genome.start.pos+sample.data$hd.length[i]-1
      chr = sum(genome.start.pos>chr.boundaries)
      chr.end = sum(genome.end.pos>chr.boundaries)
      if(chr == chr.end){
        chr.start.pos = genome.start.pos - chr.boundaries[chr]
        chr.end.pos = genome.end.pos - chr.boundaries[chr]
        sample.sim.hd = rbind(sample.sim.hd,c(chr,chr.start.pos,chr.end.pos))
      }else{
        #if deletion overlaps chromosome boundary, add 2 separate deletions, one on each chromosome
        chr.start.pos = genome.start.pos - chr.boundaries[chr]
        chr.end.pos = chr.lengths[chr]
        sample.sim.hd = rbind(sample.sim.hd,c(chr,chr.start.pos,chr.end.pos))
        if ((chr.end-chr)>1) {
          for (j in (chr+1):(chr.end-1)) {
            chr.start.pos=1
            chr.end.pos=chr.lengths[j]
            sample.sim.hd = rbind(sample.sim.hd,c(j,chr.start.pos,chr.end.pos))
          }
        }
        chr.start.pos = 1
        chr.end.pos = genome.end.pos - chr.boundaries[chr.end]
        if(chr.end==24){chr.end=1} # loop from chr 23 to chr 1
        sample.sim.hd = rbind(sample.sim.hd,c(chr.end,chr.start.pos,chr.end.pos))
      }
    }
    sim.hd = rbind(sim.hd,sample.sim.hd)
  }

  sim.hits = vector(mode="numeric",length = nrow(hd.segment.counts))
  for(c in 1:23){
    if (dim(hd.segment.counts[hd.segment.counts$chr==chr.names[c],])[1]) {
      chr = chr.names[c]
      chr.sim.data = data.frame(startpos = sim.hd[sim.hd[,1]==c,2],endpos = sim.hd[sim.hd[,1]==c,3])
      all.breakpoints = c(1,sort(unique(unlist(chr.sim.data))),chr.lengths[c])
      no.breakpoints = length(all.breakpoints)-1
      chr.segment.data = cbind(chr,all.breakpoints[1:no.breakpoints],all.breakpoints[-1],0)
      chr.segment.data[,4] = sapply(1:nrow(chr.segment.data),function(d,s,i){sum(d$startpos >= as.numeric(s[i,2]) & d$startpos < as.numeric(s[i,3])| d$endpos > as.numeric(s[i,2]) & d$endpos <= as.numeric(s[i,3]) | d$startpos < as.numeric(s[i,2]) & d$endpos > as.numeric(s[i,3]))},d=chr.sim.data,s=chr.segment.data)
      chr.segment.data[1:(no.breakpoints-1),3] = as.numeric(chr.segment.data[1:(no.breakpoints-1),3])-1
      colnames(chr.segment.data) = c("chr","startpos","endpos","count")
      chr.segment.data = data.frame(chr=chr.segment.data[,1],startpos=as.numeric(chr.segment.data[,2]),endpos=as.numeric(chr.segment.data[,3]),count=as.numeric(chr.segment.data[,4]))
      #segment.data = rbind(segment.data,chr.segment.data)
      no.segments = chr.seg.counts[c]
      for(s in 1:no.segments){
        sim.hits[cum.chr.seg.counts[c]+s] = sum(sapply(1:nrow(chr.segment.data),function(c2,d,l,i){max(min(c2$endpos[s],d$endpos[i])-max(c2$startpos[s],d$startpos[i])+1,0)*d$count[i]/l},d=chr.segment.data,c2=chr.data[[c]],l=chr.data[[c]]$endpos[s] - chr.data[[c]]$startpos[s]+1))
        if(sim.hits[cum.chr.seg.counts[c]+s]>=chr.data[[c]]$count[s]){
          sim.counts[cum.chr.seg.counts[c]+s] = sim.counts[cum.chr.seg.counts[c]+s] + 1
        }
      }
    }
  }
  write.table(sim.counts,paste(hd_output_dir,"sim_counts_",tumour_type,"_HD_",run,".txt",sep=""),sep="\t",row.names=F,quote=F)
}
