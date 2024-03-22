### make directories ###
mkdir gof_5

### get annotations for gene and intergenic and split them in 2 separate files, force strandedness ###
bedtools intersect -s -a sites_with_peakID/gof_5_T1_with_peak_id.cmd.bed -b gtf/gene.gtf -loj > gof_5/gof_5_peaks.bed ;
cat gof_5/gof_5_peaks.bed | awk 'BEGIN {OFS="\t"} { if ($10 == -1 && $11 == -1) { print $0 }}' > gof_5/gof_5_peaks_in_intergenic.bed ;
cat gof_5/gof_5_peaks.bed | awk 'BEGIN {OFS="\t"} { if ($10 != -1 && $11 != -1) { print $0 }}' > gof_5/gof_5_peaks_in_genes_ann.bed ;
# ^  used for gene section of Rmarkdown GTF_annotation

### now we turn this into a BED6 again ###
cut -f1-6 gof_5/gof_5_peaks_in_genes_ann.bed > gof_5/gof_5_peaks_in_genes.bed ;

### intersect CDS, 5' UTR and 3' UTR with peaks.bed files, enforce strandedness ###
bedtools intersect -s -a gof_5/gof_5_peaks_in_genes.bed -b gtf/merged_cds.bed > gof_5/gof_5_peaks_in_cds.bed ;
bedtools intersect -s -a gof_5/gof_5_peaks_in_genes.bed -b gtf/merged_5_prime.bed > gof_5/gof_5_peaks_in_5prime.bed ;
bedtools intersect -s -a gof_5/gof_5_peaks_in_genes.bed -b gtf/merged_3_prime.bed > gof_5/gof_5_peaks_in_3prime.bed ;

### get peaks intersecting with gtf regions ###
bedtools intersect -s -a gof_5/gof_5_peaks_in_genes.bed -b gtf/introns3.gtf > gof_5/gof_5_peaks_in_introns.bed ;

### concatenate these files in the order of the priority list: CDS, 5', 3', intron ###
awk 'BEGIN {OFS="\t"}{$7 = "CDS"; print} ' gof_5/gof_5_peaks_in_cds.bed > gof_5/f_gof_5_peaks_in_cds.bed ;
awk 'BEGIN {OFS="\t"}{$7 = "5_UTR"; print} ' gof_5/gof_5_peaks_in_5prime.bed > gof_5/f_gof_5_peaks_in_5prime.bed ;
awk 'BEGIN {OFS="\t"}{$7 = "3_UTR"; print}  ' gof_5/gof_5_peaks_in_3prime.bed > gof_5/f_gof_5_peaks_in_3prime.bed ; 
awk 'BEGIN {OFS="\t"}{$7 = "intron"; print} ' gof_5/gof_5_peaks_in_introns.bed > gof_5/f_gof_5_peaks_in_introns.bed ;
cut -f1-6 gof_5/gof_5_peaks_in_intergenic.bed | awk '{print $0 "\t" "intergenic"} ' > gof_5/f_gof_5_peaks_in_intergenic.bed ;
cat gof_5/f_gof_5_peaks_in_cds.bed gof_5/f_gof_5_peaks_in_5prime.bed gof_5/f_gof_5_peaks_in_3prime.bed gof_5/f_gof_5_peaks_in_introns.bed > gof_5_merged.bed

