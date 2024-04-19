### this script uses the BEDTools suite! make sure to install it, for example using Conda. ###
### this shows how to perform gene annotation for ONE sample ###

### make directory ###
mkdir sample_1

### get annotations for gene and intergenic and split them in 2 separate files, force strandedness ###
bedtools intersect -s -a sites_with_peakID/sample_1_T1_with_peak_id.cmd.bed -b gtf/gene.gtf -loj > sample_1/sample_1_peaks.bed ;
cat sample_1/sample_1_peaks.bed | awk 'BEGIN {OFS="\t"} { if ($10 == -1 && $11 == -1) { print $0 }}' > sample_1/sample_1_peaks_in_intergenic.bed ;
cat sample_1/sample_1_peaks.bed | awk 'BEGIN {OFS="\t"} { if ($10 != -1 && $11 != -1) { print $0 }}' > sample_1/sample_1_peaks_in_genes_ann.bed ;
# ^  used for gene section of Rmarkdown GTF_annotation

### now we turn this into a BED6 again ###
cut -f1-6 sample_1/sample_1_peaks_in_genes_ann.bed > sample_1/sample_1_peaks_in_genes.bed ;

### intersect CDS, 5' UTR and 3' UTR with peaks.bed files, enforce strandedness ###
bedtools intersect -s -a sample_1/sample_1_peaks_in_genes.bed -b gtf/merged_cds.bed > sample_1/sample_1_peaks_in_cds.bed ;
bedtools intersect -s -a sample_1/sample_1_peaks_in_genes.bed -b gtf/merged_5_prime.bed > sample_1/sample_1_peaks_in_5prime.bed ;
bedtools intersect -s -a sample_1/sample_1_peaks_in_genes.bed -b gtf/merged_3_prime.bed > sample_1/sample_1_peaks_in_3prime.bed ;

### get peaks intersecting with gtf regions ###
bedtools intersect -s -a sample_1/sample_1_peaks_in_genes.bed -b gtf/introns3.gtf > sample_1/sample_1_peaks_in_introns.bed ;

### concatenate these files in the order of the priority list: CDS, 5', 3', intron ###
awk 'BEGIN {OFS="\t"}{$7 = "CDS"; print} ' sample_1/sample_1_peaks_in_cds.bed > sample_1/f_sample_1_peaks_in_cds.bed ;
awk 'BEGIN {OFS="\t"}{$7 = "5_UTR"; print} ' sample_1/sample_1_peaks_in_5prime.bed > sample_1/f_sample_1_peaks_in_5prime.bed ;
awk 'BEGIN {OFS="\t"}{$7 = "3_UTR"; print}  ' sample_1/sample_1_peaks_in_3prime.bed > sample_1/f_sample_1_peaks_in_3prime.bed ; 
awk 'BEGIN {OFS="\t"}{$7 = "intron"; print} ' sample_1/sample_1_peaks_in_introns.bed > sample_1/f_sample_1_peaks_in_introns.bed ;
cut -f1-6 sample_1/sample_1_peaks_in_intergenic.bed | awk '{print $0 "\t" "intergenic"} ' > sample_1/f_sample_1_peaks_in_intergenic.bed ;
cat sample_1/f_sample_1_peaks_in_cds.bed sample_1/f_sample_1_peaks_in_5prime.bed sample_1/f_sample_1_peaks_in_3prime.bed sample_1/f_sample_1_peaks_in_introns.bed > sample_1_merged.bed

