# pltimer

Implementation of different versions of Plackett Luce model for timing of CNA and driver mutations. 

## Installation instructions 

The instructions below will install the latest version of pltimer. 

### Installing from Github 

To install pltimer, run the following in R:

```R
devtools::install_github("ilianapeneva/pltimer")
```

### Required reference file
pltimer requires only the file with the hg19 build chromosome coordinates. 


# Description of the steps 
pltimer takes as an input the subclones information output from Battenberg. In order to be able to include the driver mutations
information in the timing model, it needs to be converted in certain format, with the columns being 

** Tumour_Name chr startpos endpos nMaj1_A nMin1_A tumour_ploidy CNA noevent w.mean no.chrs.bearing.mut **

To get the copy number information (nMaj1_A and nMin1_A) and the multiplicity of the mutation (no.chrs.bearing.mut), 
the preprocessing step of DPClust needs to be run. 
To get the CCF (w.mean) and the clonality of the mutations, DPClust is required to be run. CNA represents the type of driver
mutation, e.g. cMut_TP53 (clonal TP53 mutation) and sMut_TP53 (subclonal TP53 mutation) and noevent is the number of the event
in the list of enriched events. 

## STEP 1: LOAD AND COLLATE THE SUBCLONES DATA 
First, the copy number information from the subclones files need to be loaded and annotated based on the type of CNA (LOH,
Loss, HD, Gain, Amplification, NoCNV). Then using the annotated subclonal segments, non-overlapping intervals are generated over
the whole genome, and the number of samples (Gain/LOH/HD/cGain/cLOH/cHD/sGain/sLOH/sHD) that have an aberration in each interval is counted.
This step outputs *_allsegs.txt* (file with the collated subclones data), *_annotated_segments.txt* (file with the annotated collated
subclones data), *_refsegs.txt* (9 files for the different types of CNA, with the number of samples over the different genomic
intervals), *_landscape.pdf* (a pdf of the overall/clonal/subclonal landscape of CNA events).

## STEP 2: IDENTIFY ENRICHED EVENTS
The enriched events in the cohort are identified by randomly placing events that have occurred in the cohort on the genome.
The null hypothesis that the event is not enriched is tested by comparing the number of times an event has occurred in the 
cohort and the expected number of times the event occurs by chance in the 1000 simulations. The simulations for Gain, LOH and HD
should ideally be in separate directories. The resulting p-values are corrected using Bonferroni correction and False discovery 
rate (FDR). The FDR-corrected p-values are used for subsequent analysis.

## STEP 3: PREPARE THE DATA ABOUT THE ENRICHED EVENTS FOR THE ORDERING
First, known artefacts, such as segments near centromeres, telomeres, and in the HLA region are removed. Any enriched regions
found in fewer than 3 patients are not included for the subsequent analysis. Any segment from the *_annotated_segments.txt* file 
overlapping with an enriched region is added and the starting and ending positions of the enriched regions are updated if
required. After that the multipcf algorithm is run to test if any new breakpoints need to be introduces in the enriched 
segments. Finally, the *_annotated_segments.txt* file is loaded in order to separate the segments covering the enriched 
regions into *LOH_mergedsegs.txt, Gain_mergedsegs.txt, HD_mergedsegs.txt*. Plots showing the segments covering each enriched
region before and after the augmentation and the multipcf run are plotted.

## STEP 4: GET THE TEMPORAL ORDERING OF EVENTS FOR THIS COHORT USING A SUITABLE VERSION OF THE PLACKETT-LUCE MODEL
First, the *mergesegs.txt* files are loaded and the driver mutations data could be added here. The losses from the *_annotated_segments.txt*
file are added to time the LOH more accurately in regards with WGD. For each sample in the cohort, all possible phylogenetic
trees are built. Then for each sample in the cohort, a possible ordering of the CNA events and driver mutaitons is obtained.
The clonal events are assumed to have occurred before the subclonal events, and if the ploidy of the sample is 2, we randomly
assign the clonal events to early or late, and place the subclonal events as having occurred after the last clonal event. If
the tumour ploidy is 4, we time the clonal events with regards to WGD by taking into account the copy number states of the 
major and minor alleles and the possible routes to arriving at these states.
If the **PLMIX** version of the Plackett-Luce model is used, the unobserved enriched events in a patient are placed after 
the last observed enriched event. If the **PlackettLuce** version of the Plackett-Luce model is used, the unobserved enriched
events in a patient are assumed not to have ever occurred.  Using either version of the Plackett-Luce model, the event 
orderings from all patients are combined to arrive at a global value quantifying how early or late each event had occurred in the
cohort. This procedure is repeated 1000 times to get a distribution of the global values for each event.

The possibility of existing multiple evolutionary trajectories within the cohort can be investigated using thee mixture model
setting of **PLMIX**. In order to get the patient groups, the cluster membership from each iteration needs to be saved, and 
then the number of times every 2 patients were placed in the same cluster is counted to obtain a matrix with the rows and 
columns being the patients and filled with the number of times the patients are in the same cluster. The patients are then
clustered by hierarchical clustering using the generated matrix. 

The distributions of the global values for the enriched events are finally plotted, which enables the inference of probable 
temporal ordering of the events. 
