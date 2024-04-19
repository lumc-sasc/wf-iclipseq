### this script uses the BedTools suite! make sure to install it, for example using Conda.. ###

### make directory ###
mkdir gtf

### make separate gtf files for the features (intron later) ###
grep -P "\tgene\t" Homo_sapiens.GRCh38.110.sorted.gtf > gtf/gene.gtf  ;
grep -P "\tCDS\t" Homo_sapiens.GRCh38.110.sorted.gtf > gtf/cds.gtf      ;
grep -P "\tfive_prime_utr\t" Homo_sapiens.GRCh38.110.sorted.gtf > gtf/5_prime.gtf ;
grep -P "\tthree_prime_utr\t" Homo_sapiens.GRCh38.110.sorted.gtf > gtf/3_prime.gtf ;
grep -P "\texon\t" Homo_sapiens.GRCh38.110.sorted.gtf > gtf/exon.gtf ;

### sort the exon file for later use ###
bedtools sort -g Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai -i gtf/exon.gtf > sorted_exon.gtf ;

### merge overlapping CDS regions, overlapping 5 UTR and 3 UTR separately, KEEP strand information (BED6 format) ###
bedtools merge -s -i gtf/cds.gtf -c 3,6,7 -o distinct > gtf/merged_cds.bed ;
bedtools merge -s -i gtf/5_prime.gtf -c 3,6,7 -o distinct > gtf/merged_5_prime.bed ;
bedtools merge -s -i gtf/3_prime.gtf -c 3,6,7 -o distinct > gtf/merged_3_prime.bed ;

### make the intron.gtf (also contains intergenic), intersect intron with peaks_in_genes.bed files later (get_features.sh) ###
### bedtools complement does not keep strand information so we intersect both strands separately ###
grep -v '+$' gtf/sorted_exon.gtf | bedtools complement -L -i stdin -g Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai | sed 's/$/\t-/' > gtf/introns.gtf ;
grep '+$' gtf/sorted_exon.gtf | bedtools complement -i stdin -g Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai | sed 's/$/\t+/' >> gtf/introns.gtf ;

### adding some extra columns, we will be using introns3.gtf for the gene annotation ###
awk 'BEGIN {OFS="\t"}{$3=$3 OFS "intron"}1 ' gtf/introns.gtf > gtf/introns2.gtf ;
awk 'BEGIN {OFS="\t"}{$4=$4 OFS "."}1 ' gtf/introns2.gtf > gtf/introns3.gtf ;
