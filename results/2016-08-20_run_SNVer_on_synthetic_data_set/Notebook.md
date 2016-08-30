2016-08-20 run SNVer on synthetic data set
==============================

Purpose
------------
Call variants using SNVer pool.


Conclusions
-----------------
SNVer performs good when VAF is 0.3%, 1.0%, 10.0%, and 100.0% compared with RVD2.  
However, it failed to call any variants in the highest read depth when VAF is 0.1%.
   
Background
-----------------
Download from http://snver.sourceforge.net/

Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
	tar zxvf SNVer-x.x.x.tar.gz 


	[fzhang@redwood 2016-08-20_run_SNVer_on_synthetic_data_set]$ java -jar SNVerPool.jar -i bam/10/100 -r ../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa -c bam/10/100/files.txt -o vcf/10/vcf100_0.vcf
	Input bam directory is : bam/10/100
	Output result file is : vcf/10/vcf100_0.vcf.raw.vcf,    vcf/10/vcf100_0.vcf.indel.raw.vcf
	
	Collecting total number of reads...
	bam/10/100/20100916_c3_p1.07_CGT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/100/20100916_c3_p1.12_GTT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/100/20100916_c3_p1.14_TCT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/100/20100916_c3_p2.07_CGT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/100/20100916_c3_p2.12_GTT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/100/20100916_c3_p2.14_TCT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	
	Ready for pileup... 
	Progress: 0%=========>10%=========>20%=========>30%=========>40%=========>50%=========>60%=========>70%=========>80%=========>90%=========>100%
	
	The p-value cutoff (method: bonferroni=0.05) is 1.25E-4
	The p-value cutoff for indel (method: bonferroni=0.05) is 0.005
	
	Filtering SNVs based on p-value cutoff for output file: vcf/10/vcf100_0.vcf.filter.vcf ...
	Filtering INDELs based on p-value cutoff for output file: vcf/10/vcf100_0.vcf.indel.filter.vcf ...
	Time usage is 1309 seconds
	Done!


	
	[fzhang@redwood 2016-08-20_run_SNVer_on_synthetic_data_set]$ java -jar SNVerPool.jar -i bam/10/10 -r ../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa -c bam/10/10/files.txt -o vcf/10/vcf10_0.vcf
	Input bam directory is : bam/10/10
	Output result file is : vcf/10/vcf10_0.vcf.raw.vcf,     vcf/10/vcf10_0.vcf.indel.raw.vcf
	
	Collecting total number of reads...
	bam/10/10/20100916_c3_p1.02_ACT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/10/20100916_c3_p1.04_ATT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/10/20100916_c3_p1.05_CAT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/10/20100916_c3_p2.02_ACT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/10/20100916_c3_p2.04_ATT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/10/20100916_c3_p2.05_CAT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	
	Ready for pileup... 
	Progress: 0%=========>10%=========>20%=========>30%=========>40%=========>50%=========>60%=========>70%=========>80%=========>90%=========>100%
	
	The p-value cutoff (method: bonferroni=0.05) is 1.25E-4
	The p-value cutoff for indel (method: bonferroni=0.05) is 0.0125
	
	Filtering SNVs based on p-value cutoff for output file: vcf/10/vcf10_0.vcf.filter.vcf ...
	Filtering INDELs based on p-value cutoff for output file: vcf/10/vcf10_0.vcf.indel.filter.vcf ...
	Time usage is 799 seconds
	Done!




	[fzhang@redwood 2016-08-20_run_SNVer_on_synthetic_data_set]$ java -jar SNVerPool.jar -i bam/10/1 -r ../../../../../freeze/h1n1/Synthetic_BAM_files/plasmid.fa -c bam/10/1/files.txt -o vcf/10/vcf1_0.vcf
	Input bam directory is : bam/10/1
	Output result file is : vcf/10/vcf1_0.vcf.raw.vcf,      vcf/10/vcf1_0.vcf.indel.raw.vcf
	
	Collecting total number of reads...
	bam/10/1/20100916_c2_p1.07_CGT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/1/20100916_c2_p1.12_GTT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/1/20100916_c2_p1.14_TCT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/1/20100916_c2_p2.07_CGT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/1/20100916_c2_p2.12_GTT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	bam/10/1/20100916_c2_p2.14_TCT.reheadered.sorted.bam:HAPLOIDS[2]MQ[20]BQ[17]
	
	Ready for pileup... 
	Progress: 0%=========>10%=========>20%=========>30%=========>40%=========>50%=========>60%=========>70%=========>80%=========>90%=========>100%
	
	The p-value cutoff (method: bonferroni=0.05) is 1.25E-4
	The p-value cutoff for indel (method: bonferroni=0.05) is 0.025
	
	Filtering SNVs based on p-value cutoff for output file: vcf/10/vcf1_0.vcf.filter.vcf ...
	Filtering INDELs based on p-value cutoff for output file: vcf/10/vcf1_0.vcf.indel.filter.vcf ...
	Time usage is 4125 seconds
	Done!





Results
----------- 


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _______Fan Zhang_______     Date: ___________2016/08/22_________


Witnessed by: ________________________
