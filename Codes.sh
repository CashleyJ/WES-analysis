trim_galore -q 30 --fastqc 8 raw/5_A5_ATGGCGAT_1.fastq.gz raw/5_A5_ATGGCGAT_2.fastq.gz
bwa mem -t 18 "@RG\tID:5A5_ATGGCGAT\tSM:5A5_ATGGCGAT\tPL:illumina\tLB:5A5_ATGGCGAT\tPU:1" Homo_sapiens_assembly38.fasta 5_A5_ATGGCGAT_1_trimmed.fq.gz 5_A5_ATGGCGAT_2_trimmed.fq.gz | samtools view -F 4 -@ 18 -O bam -h -o 5A5_ATGGCGAT_hg38_bwa_map.bam
samtools collate -o 5A5_ATGGCGAT_hg38.namecollate.bam 5A5_ATGGCGAT_hg38_bwa_map.bam
samtools fixmate -m 5A5_ATGGCGAT_hg38.namecollate.bam 5A5_ATGGCGAT_hg38.fixmate.bam
samtools sort -o 5A5_ATGGCGAT_hg38.sorted.bam 5A5_ATGGCGAT_hg38.fixmate.bam
samtools markdup 5A5_ATGGCGAT_hg38.sorted.bam 5A5_ATGGCGAT_hg38.markdup.bam
java -jar gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar BaseRecalibrator -I /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.markdup.bam -R Homo_sapiens_assembly38.fasta --known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.recal.table
java -jar gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar ApplyBQSR -I /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.markdup.bam -R Homo_sapiens_assembly38.fasta --bqsr-recal-file /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.recal.table -O /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.recal.bam
echo "#7 Apply BQSR to produce recal (recalibrated) .bam finished at `date`"
java -jar gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller -I /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.recal.bam -R Homo_sapiens_assembly38.fasta --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf.gz -L S04380110_Padded.bed -ERC GVCF -O /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.sureselectintervalhg38.dbsnp.gvcf.gz 
java -jar gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller -I /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.recal.bam -R Homo_sapiens_assembly38.fasta --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf.gz -L S04380110_Padded.bed -O /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.sureselectintervalhg38.dbsnp.vcf.gz
java -jar gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller -I /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.recal.bam -R Homo_sapiens_assembly38.fasta -O /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.vcf.gz
java -jar gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller -I /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.recal.bam -R Homo_sapiens_assembly38.fasta -L S04380110_Padded.bed -O /mnt/d/Catherine_WES/5A5_ATGGCGAT_hg38.ss_v5.vcf.gz
echo "#8 Call variants using HaplotypeCaller to produce .vcf finished at `date`"

#Merge vcf files here** before merging
#unzip files vcf.gz files in windows normally or gzip -d 
#sort .vcf files using command bcftools sort each file 
#zip the sorted files commang bgzip sorted.vcf> sorted.vcf.gz
#then index the sorted.gz files using the command bcftools index file name
bcftools sort 3H11_CATACGGA_hg38.sureselectintervalhg38.dbsnp.vcf -o 3H11_CATACGGA_hg38.sureselectintervalhg38.dbsnp_sorted.vcf
bgzip < 3H11_CATACGGA_hg38.sureselectintervalhg38.dbsnp_sorted.vcf > 3H11_CATACGGA_hg38.sureselectintervalhg38.dbsnp.sorted.vcf.gz
bcftools index 3H11_CATACGGA_hg38.sureselectintervalhg38.dbsnp.sorted.vcf.gz
#merge the sorted.vcf.gz files using bcftools merge the sorted files in vcf.gz > merged.vcf
bgzip merged.vcf
bcftools index merged.vcf.gz

 conda install bedops 
 #downlaod .bed file for kit used
  bgzip < gencode.v45.bed.gz > gencode.v45.bed.gz.gz -@ 8
  tabix /mnt/d/Catherine_WES/gencode.v45.bed.gz.gz
  tabix vcf /mnt/d/Catherine_WES/1C1_ACTATCGC_hg38.sureselectintervalhg38A.dbsnp.gvcf
 bcftools annotate -a gencode.v45.bed.gz.gz -c CHROM,FROM,TO,INFO/gene_name -h <(echo '##INFO=<ID=gene_name,Number=1,Type=String,Description="INFO from .BED">') /mnt/d/Catherine_WES/gvcf/mergedvcf.vcf > /mnt/d/Catherine_WES/gvcf/mergedvcf.annotated.vcf
 
## bcftools annotate -a gencode.v45.bed.gz -c -h <(echo '##INFO=<ID=gene_name,Number=1,Description="INFO from .BED">') /mnt/d/Catherine_WES/1C1_ACTATCGC_hg38.sureselectintervalhg38.dbsnp.gvcf.gz > /mnt/d/Catherine_WES/1C1_ACTATCGC_hg38.sureselectintervalhg38A.dbsnp.gvcf

##For SNPeff
 java -jar snpEff.jar download GRCh38.99
  java -jar snpEff.jar GRCh38.99 ../mergedvcf.annotated.vcf > ../SNPeff.merged.vcf  
  grep 'UGT1A3' SNPeff.merged.vcf > SNPeff.merged.vcf.txt (it will output a text file)













grep 'CYP1A1' SNPeff.merged.vcf >> SNPeff.merged.vcf.txt
grep 'CYP1A2' SNPeff.merged.vcf >> SNPeff.merged.vcf.txt
  sed '/^##/d' 2E2_TAGCTGAG_hg38.ss_v5.vcf > 2E2_TAGCTGAG__hg38ssv5.txt
  awk 'OFS="\t" print {$1,$2,$3,$4,$5,$8,$10}' 1C1_ACTATCGC_hg38ssv5.txt > 1C1_ACTATCGC_hg38ssv5.2.txt
   awk 'OFS="\t" {if ($6=="0/1") {print $1":"$2,$3,$4,$5,$4$5} else {print $1":"$2,$3,$4,$5,$5$5}}' 1C1_ACTATCGC_hg38ssv5.6.txt > 1C1_ACTATCGC_hg38ssv5.7.txt
    awk 'OFS="\t" {if ($6=="0/1") {print $1":"$2,$3,$4,$5,$6,$4$5} else {print $1":"$2,$3,$4,$5,$5$5}}' 1C1_ACTATCGC_hg38ssv5.6.txt > 1C1_ACTATCGC_hg38ssv5.8.txt
awk -F ":" '{print $1":"$2}' 1C1_ACTATCGC_hg38ssv5.5.txt > 1C1_ACTATCGC_hg38ssv5.quiz.txt
$ awk -F ":" 'OFS="\t" {print $1,$2}' 1C1_ACTATCGC_hg38ssv5.5.txt > 1C1_ACTATCGC_hg38ssv5.quiz3.txt
 bcftools annotate -a gencode.v45.bed.gz -c CHROM,FROM,TO,INFO/gene_name -h <(echo '##INFO=<ID=gene_name,Number=1,Type=String,Description="INFO from .BED">') /mnt/d/Catherine_WES/1C1_ACTATCGC_hg38.sureselectintervalhg38C.dbsnp.gvcf.gz > /mnt/d/Catherine_WES/1C1_ACTATCGC_hg38.sureselectintervalhg38D.dbsnp.gvcf
 $ tabix -p vcf /mnt/d/Catherine_WES/1C1_ACTATCGC_hg38.sureselectintervalhg38BD.dbsnp.gvcf.gz
  bgzip < gencode.v45.bed.gz > gencode.v45.bed.gz.gz -@ 8
 


 vcftools --vcf mergedvcf.vcf --freq -out mergedfreq  
 vcftools --vcf mergedvcf.vcf --012 -out mergedvcfGT.vcf
 vcftools --vcf mergedvcf.vcf --extract-FORMAT-info GT -out mergedgenotype.vcf
