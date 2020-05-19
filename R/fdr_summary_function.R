#' Identify enriched events using the p-values
#'
#' @param refsegs_dir Full path to the directory with the refsegs files
#' @param simulation_dir Full path to the directory where the simulations of the CNA events are saved
#' @param tumour_type String that represents the type of tumour we work with
#' @param CNA_type String that represents the type of CNA we summarise the simulations of (could be LOH, gain or HD)
#' @param output_dir Full path to the directory where the files with the enriched regions (and the p-values) will be output
#' @return files with enriched regions, the Bonferroni and FDR-corrected p-values


fdr_summary <- function(refsegs_dir, simulation_dir, tumour_type, CNA_type, output_dir){

  # load the refsegs for this CNA_type
  segment.counts = read.table(paste(refsegs_dir,tumour_type,"_", CNA_type, "_refsegs.txt", sep=""), header = T)

  # set up a counter to count the number of expected events
  segment.counts$no.perms.better = 0
  no.files.found = 0
  runs = list.files(path = simulation_dir)

  # loop over all the files
  for (i in 1:length(runs)){
    simulation_file = paste(simulation_dir,"sim_counts_",tumour_type,"_",CNA_type,"_",i,".txt",sep = "")

    if (file.exists(simulation_file)){
      count.data = unlist(read.table(simulation_file, sep="\t", header = T)[,1])
      segment.counts$no.perms.better = segment.counts$no.perms.better + count.data
      no.files.found = no.files.found + 1
    }else{
      print(paste("ERROR. file ", simulation_file," does not exist."))
    }
  }

  # H0: test a region is not enriched by comparing the expected number of times the event occuring by chance in the 1000 simulations and the number of times it has occurred
  segment.counts$p.val = segment.counts$no.perms.better/(no.files.found)
  segment.counts$bonf = p.adjust(segment.counts$p.val,method = "bonferroni")
  segment.counts$fdr = p.adjust(segment.counts$p.val,method = "fdr")         # we will work with FDR-corrected p-values because less stringent that Bonferroni

  # save all pvalues
  write.table(segment.counts, paste(output_dir, tumour_type, "_samples_pvals_",CNA_type,".txt",sep = ""), sep = "\t", quote = F)
  # save all segments with Bonferroni-adjusted p-value <= 0.05
  write.table(segment.counts[segment.counts$bonf <= 0.05,],
              paste(output_dir, tumour_type, "_samples_significant_pvals_Bonferroni_",CNA_type,".txt",sep = ""), sep = "\t", quote = F)
  # save all segments with FDR-corrected p-value <= 0.05
  write.table(segment.counts[segment.counts$fdr <= 0.05,],
              paste(output_dir, tumour_type, "_samples_significant_pvals_FDR_", CNA_type,".txt",sep = ""), sep = "\t", quote = F)

}

