#' Remove artefacts from the enriched regions, eg. segments near telomere, centromere, HLA region, in a few samples
#'
#' @param annotated_segments_file Full path to the file with the annotated CNA
#' @param tumour_type String that represents the type of tumour we work with
#' @param output_dir Full path to the directory where the files with the enriched regions (and the p-values) will be output
#' @param pvalues_dir Full path to the directory with the FDR-corrected p-values and the enriched regions
#' @param genome_chr_coordinates Full path to the file with the chromosome coordinates of the hg19 or hg38 build genome
#' @param min_region A numeric representing the minimum length of enriched region; default = 1e4
#' @param minN A numeric representing the minimum number of events per segment to be used; default = 3
#' @return files with enriched regions (separate files for LOH, HD and gain), with any artefacts removed


remove_artefacts <- function(annotated_segments_file, tumour_type, output_dir, pvalues_dir, genome_chr_coordinates, min_region=10000, minN=3){

  # LOOP TO GENERATE OUTPUT FOR ALL THREE CNA TYPES
  # load the all segments data
  allsegments <- read.table(annotated_segments_file, sep = "\t", stringsAsFactors = F, header = T)
  Ntumour = length(unique(allsegments$Tumour_Name))   # how many samples in the cohort
  cna_type  = c("Gain","LOH","HD")

  for (cna in 1:length(cna_type)){
    # load the pvalues data for this type of CNA
    pvalues_data = read.table(paste0(pvalues_dir, tumour_type,"_samples_significant_pvals_FDR_",cna_type[cna],".txt"), header = T, stringsAsFactors = F)
    names(pvalues_data) = c("chr","sp","ep","val","no.perms.better","p.val","bonf","fdr")

    # pick only the important columns (chr, start and end position, number of segments with the region being enriched and the FDR-corrected p-value)
    pvalues_data <- data.frame("chr" = pvalues_data$chr, "sp" = pvalues_data$sp, "ep" = pvalues_data$ep, "val" = pvalues_data$val, "fdr" = pvalues_data$fdr)
    pvalues_data = pvalues_data[which(pvalues_data$fdr<0.05),]  # pick only the ones with p-value < 0.05 (they should all be < 0.05)


    if (nrow(pvalues_data)>0){
      #cumulative chromosome coordinate generation
      v = NULL
      for (i in c(1:22,"X")){
        data.chr <- pvalues_data[pvalues_data$chr==i,]
        norow <- nrow(data.chr)
        if (i=="X"){
          j=23
        } else{
          j=i
        }
        v[j] = norow    #number of segments for each chromomsome
      }
      v2 = cumsum(v)    # cumulative sum of number of segments/rows


      all.sig.data = data.frame()    # a data frame with all the ordered enriched regions from the chromosomes (chr startpos endpos) and eventindex (which is the number/index of the enriched region in this data frame)
      for (i in c(1:22,"X")){
        data.chr <- pvalues_data[pvalues_data$chr==i,]

        if (nrow(data.chr)>0){
          segnumber <- nrow(data.chr)    # number of segments for this chromosome
          tempvector = NULL
          sig.data = NULL
          no.chr=ifelse(i=="X",22,as.numeric(i)-1)
          for (j in 1:segnumber){
            k = ifelse(no.chr==0,j,v2[no.chr]+j)
            sig.data.segment <- data.frame("chr"=i, "startpos"= data.chr$sp[j], "endpos" = data.chr$ep[j], "eventindices" = k)    # put together all the enriched segments
            all.sig.data     <- rbind(all.sig.data,sig.data.segment)
            sig.data.segment = NULL
          }
        } else {print(paste("chr",i,"has no enriched segments"))}
      }
      all.sig.data$chr = as.character(all.sig.data$chr)


      enriched_data = data.frame()      # create enriched_data data frame where we are going to merge the FDR-enriched regions in each chr
      chrs = unique(all.sig.data$chr)   # all the chromosomes having enriched regions

      # load the chromosome coordinates for hg19 or hg38
      chr_loc = read.table(genome_chr_coordinates, header = T, stringsAsFactors = F)


      for (i in 1:length(chrs)){
        seg_count = nrow(enriched_data)                                 # number of included segments so far
        CHR = all.sig.data[which(all.sig.data$chr==chrs[i]),]           # pick all the data for that chr
        if (nrow(CHR)==1){
          enriched_data = rbind(enriched_data, data.frame(chr = chrs[i], startpos = CHR$startpos, endpos = CHR$endpos, length = CHR$endpos-CHR$startpos))
        } else{
          # check how close the next enriched region is: if it is withing 1e4 of the previous, merge it to the previous segment; if not, don't merge
          start = CHR$startpos[1]
          for (j in 2:nrow(CHR)){
            if (CHR$startpos[j] <= CHR$endpos[j-1]+1e4){                # change interval to appropriate level (e.g. 1e4 or 1e5) ?
              end = CHR$endpos[j]                                       # include the new row (j) in the merge
            }
            else {
              end = CHR$endpos[j-1]                                      # stop merge at the previous row (j-1)
              enriched_data = rbind(enriched_data, data.frame(chr=chrs[i], startpos = start, endpos = end, length = end-start))
              start = CHR$startpos[j]
            }
          }

          # check if the merging of the smaller enriched regions into larger has resulted in the creation of one large segment
          if (start==CHR$startpos[1] & end==CHR$endpos[nrow(CHR)]){
            enriched_data = rbind(enriched_data,data.frame(chr = chrs[i], startpos = start, endpos = end, length = end-start))                 # if all segments were adjacent to each other
          } else {
            enriched_data = rbind(enriched_data,data.frame(chr = chrs[i], startpos = start, endpos = CHR$endpos[j], length = CHR$endpos[j]-start)) # if final segment is distant from previous
          }
        }
        segs_per_chr = nrow(enriched_data) - seg_count
        print(paste(chrs[i],"has ",segs_per_chr,"enriched region(s)"))
      }
      enriched_data$chr = as.character(enriched_data$chr)

      #identify enriched regions which are at telomeres and REMOVE #
      chrom = as.vector(unique(enriched_data$chr))
      ENRICHED = data.frame()             # this data frame will store all the enriched regions which are not near centromere/telomere or in the HLA region

      for (i in 1:length(chrom)){
        indices_remove = NULL                          # this will store the indices of the segments to be removed from the enriched_data
        CHROM = enriched_data[which(enriched_data$chr==chrom[i]),]
        # not sure why we have a limit on how big the size of segment - don't we want to remove any segment that is near the telomeres??

        for (j in 1:nrow(CHROM)){
          # if end of the enriched segment is within 1Mb of the end of the chromosome and the length of the segment is less than 5Mb; remove this segment
          if (CHROM$endpos[j]>=chr_loc[ifelse(chrom[i]!="X",chrom[i],23),]$end-1e6 & CHROM$length[j]<=5e6){ # length increased to 5Mb from 2Mb
            indices_remove = append(indices_remove,j)
          }
          # if the start of the enriched region is near the beginning of the chromosome and the length is less than 2Mb; remove this segment
          if (CHROM$startpos[j]<1e6 & CHROM$length[j]<=2e6){
            indices_remove = append(indices_remove,j)
          }
        }
        if (!is.null(indices_remove)){
          print(paste("enriched segment", indices_remove, "removed from chromosome", chrom[i]))
          CHROM = CHROM[-indices_remove,]
        }
        if (nrow(CHROM)>0){
          ENRICHED = rbind(ENRICHED,CHROM)
        }
        else {print(paste("Chromosome",chrom[i],": only telomeric and removed"))}
      }

      #identify enriched regions which are at the HLA region and remove #
      #https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37 MHC region based on hg37
      if (nrow(ENRICHED)>0){
        hla = c(28477797,33448354)        # coordinates of the hla region
        # added the last condition because the statement was missing the case where the enriched segment was covering the whole HLA region
        HLA = ENRICHED[which(ENRICHED$chr==6 & ((ENRICHED$startpos>hla[1] & ENRICHED$endpos<hla[2])|(ENRICHED$endpos>hla[1] & ENRICHED$startpos<hla[1])|(ENRICHED$startpos<hla[2] & ENRICHED$endpos>hla[2]|ENRICHED$startpos<hla[1] & ENRICHED$endpos>hla[2]))),]

        ENRICHED=setdiff(ENRICHED,HLA)    # remove any segments that are found in HLA region
      }


      # REMOVE SINGLETONS/RARE DOUBLETS when N is high:
      rare = NULL
      if (nrow(ENRICHED)>0){
        for (i in 1:nrow(ENRICHED)){
          ENRICHEDi = ENRICHED[i,]

          # pick up from the pvalues_data how many samples have each of the enriched events and remove the event if insufficient number of samples (max(3,Ntumour/100)) have it
          idata = pvalues_data[which(pvalues_data$chr==ENRICHEDi$chr & pvalues_data$sp<=ENRICHEDi$endpos & pvalues_data$ep>=ENRICHEDi$startpos),]
          if (max(idata$val) < max(minN,Ntumour/100)){
            rare = append(rare,i)
          }
        }
      }


      if (!is.null(rare)){
        ENRICHEDclean = ENRICHED[-rare,]
      } else {ENRICHEDclean = ENRICHED}

      # for noevent:
      if (nrow(ENRICHEDclean)>0){
        ENRICHEDclean = ENRICHEDclean[which(ENRICHEDclean$length>min_region),]
        ENRICHEDclean$length = NULL
        ENRICHEDclean$noevent = 1:nrow(ENRICHEDclean)
        write.table(ENRICHEDclean, paste0(output_dir,tumour_type,"_enriched_regions_",cna_type[cna],"_noevent.txt"), sep="\t",quote=F,row.names = F)
      } else {print(paste(data_dir,tumour_type,":no file written out - no enriched region remained after correction"))}
    } else {print(paste(data_dir ,tumour_type,"_samples_significant_pvals_FDR_",cna_type[cna],".txt")," is empty")}
  }
}


#' Merge segments overlapping the enriched regions
#'
#' @param annotated_segments_file Full path to the file with the annotated segments file
#' @param tumour_type String that represents the type of tumour we work with
#' @param enriched_dir Full path to the directory with the files with the enriched regions
#' @param output_dir Full path to the directory where the merged enriched segments and the plots of the enriched regions are output
#' @param genome_chr_coordinates Full path to the file with the chromosome coordinates of the hg19 or hg38 build genome
#' @return files with merged enriched regions (separate files for LOH, HD and gain) and plots of the segments overlapping the enriched regions


merge_enriched_regions <- function(annotated_segments_file, tumour_type, enriched_dir, genome_chr_coordinates, output_dir){

  cna_type = c("Gain","LOH","HD")
  allsegments <- read.table(annotated_segments_file, sep = "\t", stringsAsFactors = F, header = T) # load the allsegments file

  for (cna in 1:length(cna_type)){
    CNAtype <- paste0(c("s","c"), cna_type[cna])
    cna.data <- allsegments[allsegments$CNA %in% CNAtype, ]    # pick up the subclonal and clonal events for this CNA

    #remove subclonal after clonal AND order by sample, chr and startpos
    cna.data$chr[cna.data$chr=='X'] <- 23                   # set chr X as 23
    cna.data$chr = as.numeric(cna.data$chr)
    cna.data = cna.data[!duplicated(cna.data[,1:4]),]       # remove any duplicated events

    if (file.exists(paste0(enriched_dir, tumour_type,"_enriched_regions_",cna_type[cna],"_noevent.txt"))){
      enriched = read.table(paste0(enriched_dir, tumour_type,"_enriched_regions_",cna_type[cna],"_noevent.txt"),header=T,stringsAsFactors = F,sep="\t")
      enriched$chr[enriched$chr=="X"] = 23
      enriched$chr <- as.numeric(enriched$chr)

      print("files read")


      setDT(cna.data)                                              # convert cna.data into data.table
      setDT(enriched)                                              # all the enriched segments
      setkey(cna.data,"chr","startpos","endpos")                   # sort the cna.data by chr, startpos,endpos
      setkey(enriched,"chr","startpos","endpos")                   # sort the enriched by chr, startpos, endpos

      merged = foverlaps(cna.data, enriched, type = "any", nomatch = 0)     # find the overlaps between the cna.data and the enriched regions (will return only the matches)
      # merged is in the format first four columns corresponding to/being the same as the enriched regions intervals and the number of the event; then a column with the tumour names and columns with interval for
      # that tumour name overlapping this enriched region (with i.startpos and i.endpos being the start/end of the segment) and then the copy number information for this segment

      # remove regions outside of ENRICHED REGIONS - this is removing any parts of the segment overlapping the enriched region that are outside the enriched region
      merged$i.startpos = ifelse(merged$i.startpos<merged$startpos, merged$startpos, merged$i.startpos)        # i.startpos and i.endpos - the start and end of the segment overlapping with the enriched region
      merged$i.endpos = ifelse(merged$i.endpos>merged$endpos, merged$endpos, merged$i.endpos)
      print("files merged")

      # add 'y' to plot the subclones segments in the enriched regions
      tumour_names = unique(merged$Tumour_Name)
      for (i in 1:length(tumour_names)){
        for (j in 1:nrow(merged)){
          if (merged$Tumour_Name[j]==tumour_names[i]){
            merged$y[j] = i
          }
        }
      }

      # all chromosomes that have enriched regions
      chrs = unique(merged$chr)

      # load the chromosome coordinates for hg19 or hg38
      chr_loc = read.table(genome_chr_coordinates, header = T, stringsAsFactors = F)

      for (k in 1:length(chrs)){
        chr_enriched_regions = merged[which(merged$chr==chrs[k]),]   # plot the segments in the enriched regions
        enriched_regions_plot = ggplot() +
          geom_segment(data = chr_enriched_regions, aes(x = i.startpos, y = y, xend = i.endpos, yend = y, col = Tumour_Name))+
          theme(legend.position = "NONE") +
          ggtitle(paste(tumour_type, cna_type[cna]," - Enriched regions chr",chrs[k]))+
          expand_limits(y = 0, x = 0)+
          geom_vline(xintercept = chr_loc[chr_loc$chr==chrs[k],]$end, col = "black", linetype = "longdash")  # plot centromere/telomere

        # save the plot for this enriched region
        save_plot(paste0(output_dir, tumour_type,"_enriched_",cna_type[cna],"_chr",chrs[k],"_cowplot.pdf"), plot = enriched_regions_plot)

      }
        write.table(merged, paste0(output_dir, tumour_type,"_segments_in_enriched_regions_",cna_type[cna],".txt"), sep = "\t", quote = F,row.names = F)

      } else {print(paste(output_dir , tumour_type,"_enriched_regions_",cna_type[cna]," _noevent.txt"),"does not exist")}
    }
  }




#' Check if any new breakpoints should be added to the enriched regions
#'
#' @param enriched_dir Full path to the directory with the files with the enriched regions
#' @param tumour_type String that represents the type of tumour we work with
#' @param output_dir Full path to the directory where the enriched segments with the new breakpoints and the plots of start and end positions of the enriched regions are output
#' @return files with merged enriched regions (separate files for LOH, HD and gain) and plots of the segments overlapping the enriched regions

multipcf_new_breakpoints <- function(enriched_dir, tumour_type, output_dir){

  cna_type=c("Gain","LOH","HD")

  for (cna in 1:length(cna_type)){
    if (file.exists(paste0(enriched_dir, tumour_type,"_segments_in_enriched_regions_",cna_type[cna],".txt"))){
      # load the enriched data for this CNA type
      enriched_regions = read.table(paste0(enriched_dir, tumour_type,"_segments_in_enriched_regions_",cna_type[cna],".txt"), header = T, stringsAsFactors = F)
      enriched_regions = enriched_regions[!duplicated(enriched_regions[,c("Tumour_Name","i.startpos","i.endpos")]),] # remove any subclonal after clonal e.g. 20% 2:1, 80% 2:2 -> 100% 2:1 gain


      # adding weighted CCF to all segments in same individual
      enriched_regions$diff = enriched_regions$i.endpos - enriched_regions$i.startpos
      enriched_regions = enriched_regions[which(enriched_regions$diff>1),]               # remove single bp segments
      chrs = unique(enriched_regions$chr)


      # VISUALISE SEGMENT START POINTS#
      title = paste(tumour_type,"- CNA",cna_type[cna],"patterns")
      all_enriched_segments = data.frame()
      for (i in 1:length(chrs)){
        CHR = enriched_regions[which(enriched_regions$chr==chrs[i]),]
        sample_names = unique(CHR$Tumour_Name)   # all the samples that have data for this chromosome
        for (j in 1:length(sample_names)){
          sample_data = CHR[which(CHR$Tumour_Name==sample_names[j]),]
          sample_data$w.mean = weighted.mean(sample_data$frac1_A,sample_data$diff)   # ccf for the event; diff is the length of the segment
          all_enriched_segments = rbind(all_enriched_segments, sample_data)
        }
      }

      enriched_segments_plot = ggplot() + geom_point(data = all_enriched_segments, aes(x = i.startpos, y = w.mean, col = Tumour_Name)) + theme(legend.position = "NONE") +
        facet_grid(chr~.) +
        labs(x = "Chromosome Position")
      save_plot(paste0(output_dir, tumour_type,"_CNA_",cna_type[cna],"_patterns.pdf"), plot = enriched_segments_plot)


      # pivot table to prepare for multiPCF  #
      #  the data frame input for multiPCF should be in the following format:
      # the rows of the data frame should be the probes
      # the columns: column 1 = numeric/character chromosome numbers
      #              column 2 = the numeric local probe coordinates
      #              subsequent columns = the copy number measurements for 2 or more samples
      # the header should contain the sample IDs

      print("multiPCF input generation")
      all_enriched_segments$ID = paste(all_enriched_segments$chr, all_enriched_segments$i.startpos,sep = "_")        # add ID to the segments overlapping the enriched regions which are the chromosome and the starting position of this segment

      dc = dcast(all_enriched_segments, ID~Tumour_Name, value.var="w.mean")   # restructure the data for multipcf; the IDs of the events are in the first column, and the patient platekeys are in the subsequent columns; the matrix is filled with the CCF of the event for the corresponding sample
      dc[is.na(dc)] = 0

      out = dc %>% separate(ID, into = c("chr","startpos"),sep = "\\_")                    # split the ID column into 2 - chromosome and starting position
      out = out %>% mutate_if(is.character,as.numeric) %>% arrange(chr,startpos)           # sort by chromosome and position

      # run multiPCF
      mp = multipcf(out, gamma = 0.1) # gamma 10 or 40 = same breakpoints
      mp = mp[order(mp$chrom,mp$arm,mp$start.pos),]

      # save the multipcf output as the new enriched regions output
      ENRICHED = data.frame(chr = mp$chrom, startpos = mp$start.pos, endpos = mp$end.pos, noevent = 1:nrow(mp))
      write.table(ENRICHED, paste(output_dir, tumour_type,"_enriched_regions_",cna_type[cna],"_noevent.txt"), sep = "\t", quote = F, row.names = F)
    } else {print(paste(output_dir, tumour_type,"_segments_in_enriched_regions_",cna_type[cna],".txt"),"does not exist")}
  }
}



#' Process the enriched data into the format required for the ordering step with the Plackett-Luce model
#'
#' @param annotated_segments_file Full path to the file with the annotated CNA
#' @param tumour_type String that represents the type of tumour we work with
#' @param enriched_dir Full path the directory with the enriched regions
#' @param genome_chr_coordinates Full path to the file with the chromosome coordinates of the hg19 or hg38 build genome
#' @param output_dir Full path to the directory where the files with the enriched regions will be output
#' @return files with enriched regions (separate files for LOH, HD and gain) in the format for the ordering script
prepare_ordering_data <- function(annotated_segments_file, tumour_type, enriched_dir, genome_chr_coordinates, output_dir){
  # load the allsegments data
  allsegments <- read.table(annotated_segments_file, sep = "\t", stringsAsFactors = F, header = T)
  cna_type    <- c("LOH", "Gain", "HD")


  for (cna in 1:length(cna_type)){
    CNAtype  <- paste0(c("s","c"), cna_type[cna])
    cna.data <- allsegments[allsegments$CNA %in% CNAtype, ]   # get all the clonal and subclonal events for this type of CNA

    # remove subclonal after clonal AND order by sample, chr and startpos cna.data$chr[is.na(cna.data$chr)]=23
    cna.data$chr = as.numeric(cna.data$chr)
    cna.data$chr[is.na(cna.data$chr)] = 23
    cna.data = cna.data[!duplicated(cna.data[,1:4]), ]

    if (file.exists(paste0(enriched_dir, tumour_type,"_enriched_regions_",cna_type[cna],"_noevent.txt"))){
      # load the enriched regions for this type of CNA
      enriched = read.table(paste(enriched_dir, tumour_type,"_enriched_regions_",cna_type[cna],"_noevent.txt"), header = T, stringsAsFactors = F, sep = "\t")
      enriched$chr[enriched$chr=="X"] = 23
      enriched$chr = as.numeric(enriched$chr)
      print("files read")

      # convert cna.data and enriched to data.tables and order by chromosome, startpos and endpos
      setDT(cna.data)
      setDT(enriched)
      setkey(cna.data,"chr","startpos","endpos")
      setkey(enriched,"chr","startpos","endpos")

      # find the overlap between the cna.data and the enriched
      merged = foverlaps(cna.data, enriched, type = "any", nomatch = 0)

      # remove regions outside of ENRICHED REGIONS
      merged$i.startpos = ifelse(merged$i.startpos < merged$startpos, merged$startpos, merged$i.startpos)
      merged$i.endpos = ifelse(merged$i.endpos > merged$endpos, merged$endpos, merged$i.endpos)
      print("files merged")

      # add 'y' for geom_segment() to plot the segments covering each of the enriched regions following the multipcf step
      tumour.names = unique(merged$Tumour_Name)
      for (i in 1:length(tumour.names)){
        for (j in 1:nrow(merged)){
          if (merged$Tumour_Name[j]==tumour.names[i]){
            merged$y[j]=i
          }
        }
      }


      # plot the segments in each enriched region for all chromosomes
      chrs = unique(merged$chr)
      # load the chromosome coordinates for hg19 or hg38
      chr_loc = read.table(genome_chr_coordinates, header = T,stringsAsFactors = F)

      for (i in 1:length(chrs)){
        chr_enriched_regions=merged[which(merged$chr==chrs[i]),]
        enriched_regions_plot = ggplot() + geom_segment(data=chr_enriched_regions, aes(x = i.startpos, y = y, xend = i.endpos, yend=y, col = Tumour_Name)) +
          theme(legend.position = "NONE") +
          ggtitle(paste(tumour_type, cna_type[cna]," - Enriched regions chr",chrs[i])) +
          expand_limits(y=0,x=0) +
          geom_vline(xintercept = chr_loc[chr_loc$chr==chrs[i],]$end, col = "black", linetype = "longdash")

        # save the plot
        save_plot(paste0(output_dir, tumour_type,"_enriched_", cna_type[cna],"_chr",chrs[i],"cowplot_postFDR4.pdf"), plot = enriched_regions_plot)
      }

      # Prepare the input for the ordering script
      merged_enriched_regions = merged
      merged_enriched_regions = merged_enriched_regions[!duplicated(merged_enriched_regions[,c("Tumour_Name","i.startpos","i.endpos")]), ] # remove any subclonal after clonal e.g. 20% 2:1, 80% 2:2 -> 100% 2:1 gain

      # adding weighted CCF to all segments in same individual
      merged_enriched_regions$diff = merged_enriched_regions$i.endpos - merged_enriched_regions$i.startpos
      merged_enriched_regions = merged_enriched_regions[which(merged_enriched_regions$diff>1),]      # remove single bp segments

      chrs = unique(merged_enriched_regions$chr)

      data_for_ordering = data.frame()
      for (i in 1:length(chrs)){
        CHR = merged_enriched_regions[which(merged_enriched_regions$chr==chrs[i]),]    # pick the data for this chromosome
        sample_names = unique(CHR$Tumour_Name)                                         # all the sample names with data for this chromosome
        for (j in 1:length(sample_names)){
          sample_data = CHR[which(CHR$Tumour_Name==sample_names[j]),]                               # pick the info for this sample
          sample_data_weighted = data.frame(sample_data[1,c(5,1:3,8:9,12:13,4)], w.mean = weighted.mean(sample_data$frac1_A,sample_data$diff)) # weigh the ccf
          data_for_ordering = rbind(data_for_ordering, sample_data_weighted)
        }
      }

      # save the table
      write.table(data_for_ordering, paste0(output_dir, tumour_type,"_",cna_type[cna],"_mergedsegs.txt"), sep= "\t", row.names = F, quote = F)
    } else {print(paste(output_dir, tumour_type,"_enriched_regions_",cna_type[cna],"_noevent.txt" ),"does not exist")}
  }

}
