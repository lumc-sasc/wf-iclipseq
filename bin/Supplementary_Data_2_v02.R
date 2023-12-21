# SUPPLEMENTARY DATA 2
# Busch et al, Methods, 2019

# This file provides the R code to postprocess the output of PureCLIP peak calling as described in Chapter 5.2.

# Running this code requires a bed file with 'crosslink sites' output by PureCLIP (PureCLIP.crosslink_sites_short.bed)
# and bw files with crosslink events (sampleX.strand.bw). The code to obtain these files is described in Supplementary Data 1.

###############################

# added/adjusted by me (Amarise):

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("You're missing an argument to run this R script!", call. = FALSE)
}

# The following files need to be specified to run the code:
plus_bw.file  <- args[1]  #  < /path/to/sampleX.plus.bw >
minus_bw.file <- args[2]  #  < /path/to/sampleX.minus.bw >
peak.file     <- args[3]  #  < /path/to/PureCLIP.crosslink_sites_short.bed >

# adjustment ends here

###############################

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

###############################
### Import crosslink events ###
###############################

### Import sampleX.strand.bw for both strands and convert into coverage

plus_cov = import.bw(plus_bw.file, as = "Rle")
minus_cov = abs(import.bw(minus_bw.file, as = "Rle"))



############################################
### Import crosslink sites from PureCLIP ###
############################################

peaks = import(peak.file)
length(peaks)


#######################################
### Filter for standard chromosomes ###
#######################################

peaks_sc = keepStandardChromosomes(peaks, pruning.mode="coarse")
length(peaks_sc)


#############################################
### Merge peaks with gap width up of 7 nt ###
#############################################

peaks_sc_m = reduce(peaks_sc, min.gapwidth = 8)
length(peaks_sc_m)

export(peaks_sc_m,  "./bed_files/sign_sites_max_gap_width_7.bed")


##########################################
### Remove peaks with width 1 and 2 nt ###
##########################################

peaks_sc_m_pw = peaks_sc_m[width(peaks_sc_m) > 2]
length(peaks_sc_m_pw)


######################################################################################
### Get too short peaks (width <= 9 nt), position at maximum ( 5' to 3' direction) ###
### and extend by 4 nt up and downstream                                           ###
######################################################################################

### Get too short peaks
peaks_sc_m_pw_short = peaks_sc_m_pw[width(peaks_sc_m_pw) <= 9]
length(peaks_sc_m_pw_short)

### Resize peaks on plus strand
peaks_sc_m_pw_short_plus = peaks_sc_m_pw_short[strand(peaks_sc_m_pw_short) == "+"]
max_pos_plus = sapply(plus_cov[peaks_sc_m_pw_short_plus], which.max)
peaks_sc_m_pw_short_plus_centered = GRanges(seqnames = seqnames(peaks_sc_m_pw_short_plus),
                                     ranges = IRanges(start = start(peaks_sc_m_pw_short_plus) + max_pos_plus - 1, width = 1),
                                     strand = "+")

### Resize peaks on minus strand
peaks_sc_m_pw_short_minus = peaks_sc_m_pw_short[strand(peaks_sc_m_pw_short) == "-"]
max_pos_minus = sapply(lapply(minus_cov[peaks_sc_m_pw_short_minus], rev), which.max)
peaks_sc_m_pw_short_minus_centered = GRanges(seqnames = seqnames(peaks_sc_m_pw_short_minus),
                                      ranges = IRanges(start = end(peaks_sc_m_pw_short_minus) - max_pos_minus + 1, width = 1),
                                      strand = "-")

### Combine resized peaks
peaks_sc_m_pw_sc = c(peaks_sc_m_pw_short_plus_centered + 4, peaks_sc_m_pw_short_minus_centered + 4)

### Check widths before and after resizing
length(peaks_sc_m_pw_short)
table(width(peaks_sc_m_pw_short))

length(peaks_sc_m_pw_sc)
table(width(peaks_sc_m_pw_sc))

### Exclude overlaps
stopifnot(nrow(as.matrix(findOverlaps(peaks_sc_m_pw_sc))) == length(peaks_sc_m_pw_sc))


### Cleanup
rm(list = setdiff(ls(), c("peaks_sc_m_pw", "peaks_sc_m_pw_sc", "plus_cov", "minus_cov")))


###############################################
### Fit multiple peaks in too long peaks    ###
### fits are centered on decreasing maxima  ###
### only non-overlappping peaks are allowed ###
###############################################

###---------------------
### Plus strand
###---------------------

### Get too long peaks
peaks_sc_m_pw_long = peaks_sc_m_pw[width(peaks_sc_m_pw) > 9]
length(peaks_sc_m_pw_long)

### Process peaks on plus strand
peaks_sc_m_pw_long_plus = peaks_sc_m_pw_long[strand(peaks_sc_m_pw_long) == "+"]

peaks_sc_m_pw_long_plus_curated = GRanges()

### Interatively place peaks
for (i in 1:length(peaks_sc_m_pw_long_plus)){
  peak_tmp = peaks_sc_m_pw_long_plus[i]
  
  peak_set =  GRanges(seqnames = seqnames(peak_tmp),
                      ranges   = IRanges(start = c((start(peak_tmp)-4):(end(peak_tmp)-4)),
                                         width = 9),
                      strand   = strand(peak_tmp))
  max_ordered = order(sapply(plus_cov[peak_set-4], runValue), decreasing = T)
  
  # Reorder peaks by max signal
  peaks_ord = peak_set[max_ordered]

  peaks_final = peaks_ord[1]
  peaks_cands = peaks_ord[-1]
  peaks_cands = peaks_cands[-as.matrix(findOverlaps(peaks_final,peaks_cands))[,2]]
 
  while (length(peaks_cands) > 0){
    
    peaks_final = c(peaks_cands[1], peaks_final)
    peaks_cands = peaks_cands[-1]
    peaks_cands = peaks_cands[-as.matrix(findOverlaps(peaks_final,peaks_cands))[,2]]}
  
  peaks_sc_m_pw_long_plus_curated = c(peaks_sc_m_pw_long_plus_curated, peaks_final)
  
  # Check that there is really no overlap
  stopifnot(nrow(as.matrix(findOverlaps(peaks_sc_m_pw_long_plus_curated))) == length(peaks_sc_m_pw_long_plus_curated))
  
  if (i %% 500 == 0)
    cat(paste("Iteration ", i , " is done (out of ", length(peaks_sc_m_pw_long_plus), ").\n", sep = ""))
}

### Cleanup
rm(list = setdiff(ls(), c("peaks_sc_m_pw", "peaks_sc_m_pw_sc", "peaks_sc_m_pw_long_plus_curated", "plus_cov", "minus_cov")))


###---------------------
### Minus strand
###---------------------

### Get too long peaks
peaks_sc_m_pw_long = peaks_sc_m_pw[width(peaks_sc_m_pw) > 9]
length(peaks_sc_m_pw_long)

### Process peaks on minus strand
peaks_sc_m_pw_long_minus = peaks_sc_m_pw_long[strand(peaks_sc_m_pw_long) == "-"]

peaks_sc_m_pw_long_minus_curated = GRanges()

for (i in 1:length(peaks_sc_m_pw_long_minus)){
  peak_tmp = peaks_sc_m_pw_long_minus[i]
  
  # Peaks in direction 5' to 3' 
  peak_set =  GRanges(seqnames = seqnames(peak_tmp),
                      ranges   = IRanges(end = c((end(peak_tmp)+4):(start(peak_tmp)+4)),
                                         width = 9),
                      strand   = strand(peak_tmp))
  max_ordered = order(sapply(minus_cov[peak_set-4], runValue), decreasing = T)
  
  # Reorder peaks by max signal
  peaks_ord = peak_set[max_ordered]
  
  peaks_final = peaks_ord[1]
  peaks_cands = peaks_ord[-1]
  peaks_cands = peaks_cands[-as.matrix(findOverlaps(peaks_final,peaks_cands))[,2]]
  
  while (length(peaks_cands) > 0){
    
    peaks_final = c(peaks_cands[1], peaks_final)
    peaks_cands = peaks_cands[-1]
    peaks_cands = peaks_cands[-as.matrix(findOverlaps(peaks_final,peaks_cands))[,2]]}
  
  peaks_sc_m_pw_long_minus_curated = c(peaks_sc_m_pw_long_minus_curated, peaks_final)
  
  # Check that there is really no overlap
  stopifnot(nrow(as.matrix(findOverlaps(peaks_sc_m_pw_long_minus_curated))) == length(peaks_sc_m_pw_long_minus_curated))
  
  if (i %% 500 == 0)
    cat(paste("Iteration ", i , " is done (out of ", length(peaks_sc_m_pw_long_minus), ").\n", sep = ""))
}

### Cleanup
rm(list = setdiff(ls(), c("peaks_sc_m_pw", "peaks_sc_m_pw_sc", "peaks_sc_m_pw_long_plus_curated", "peaks_sc_m_pw_long_minus_curated", "plus_cov", "minus_cov")))


##########################################
### Merge resized short and long peaks ###
##########################################

peaks_sl_curated = c(peaks_sc_m_pw_sc, peaks_sc_m_pw_long_plus_curated, peaks_sc_m_pw_long_minus_curated)
length(peaks_sl_curated)

### Exclude overlaps
stopifnot(nrow(as.matrix(findOverlaps(peaks_sl_curated))) == length(peaks_sl_curated))

### Cleanup
rm(list = setdiff(ls(), c("peaks_sl_curated", "plus_cov", "minus_cov")))


#############################################################################################
### Check for at least 3 positions with crosslinks (any crosslinks, not only significant) ###
#############################################################################################

### Plus strand
peaks_sl_curated_plus = peaks_sl_curated[strand(peaks_sl_curated) == "+"]
peaks_curated_plus = peaks_sl_curated_plus[as.vector(which(sapply(lapply(plus_cov[peaks_sl_curated_plus], function(x)  x > 0), sum) >= 3))]

### Minus strand
peaks_sl_curated_minus = peaks_sl_curated[strand(peaks_sl_curated) == "-"]
peaks_curated_minus = peaks_sl_curated_minus[as.vector(which(sapply(lapply(minus_cov[peaks_sl_curated_minus], function(x)  x > 0), sum) >= 3))]

### Combine
peaks_curated = c(peaks_curated_plus, peaks_curated_minus)

### Check numbers
length(peaks_sl_curated_plus)
length(peaks_curated_plus)
length(peaks_sl_curated_minus)
length(peaks_curated_minus)
length(peaks_curated)

### Exclude overlaps
stopifnot(nrow(as.matrix(findOverlaps(peaks_curated))) == length(peaks_curated))


############################
### Save final peak file ###
############################

export(peaks_curated, "peaks_curated.bed")

