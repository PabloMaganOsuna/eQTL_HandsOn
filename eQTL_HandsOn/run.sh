wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a
20130502.genotypes.vcf{.gz,.gz.tbi} --directory-prefix input/unprocessed/1000g

#Task2:
print("Ejecutando Task2")

# Get GEUVADIS samples from the metadata
cut -f1 input/unprocessed/geuvadis/geuvadis.metadata.txt | sed '1d' | sort | uniq > tmp/geuvadis.samples.txt
# Subset the VCF (common samples, biallelic SNPs and indels, MAF >= 0.05, no duplicates)
bcftools view -v snps,indels -m 2 -M 2 -q 0.05:minor -S tmp/geuvadis.samples.txt -Ob input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -d all -Oz -o tmp/genotypes.chr22.vcf.gz
# Subset the VCF so that there are at least 10 individuals per genotype group and compress it (for indexing we require 'bgzip' compression)
filter.genotype.py -t 10 -g <(zcat tmp/genotypes.chr22.vcf.gz) | bgzip > input/processed/genotypes.chr22.vcf.gz
# Index the VCF
tabix -p vcf input/processed/genotypes.chr22.vcf.gz

#Q1:
################################################################BCFTOOLSOPTIONS###############################################################################
#FILE
#Files can be both VCF or BCF, uncompressed or BGZF-compressed. The file "-" is interpreted as standard input. Some tools may require tabix- or CSI-indexed files.
#-c, --collapse snps|indels|both|all|some|none|id
#Controls how to treat records with duplicate positions and defines compatible records across multiple input files. Here by "compatible" we mean records which should be considered as identical by the tools. For example, when performing line intersections, the desire may be to consider as identical all sites with matching positions (bcftools isec -c all), or only sites with matching variant type (bcftools isec -c snps  -c indels), or only sites with all alleles identical (bcftools isec -c none).

#none
#only records with identical REF and ALT alleles are compatible
#some
#only records where some subset of ALT alleles match are compatible
#all
#all records are compatible, regardless of whether the ALT alleles match or not. In the case of records with the same position, only the first #will be considered and appear on output.
#snps
#any SNP records are compatible, regardless of whether the ALT alleles match or not. For duplicate positions, only the first SNP record will be considered and appear on output.
#indels
#all indel records are compatible, regardless of whether the REF and ALT alleles match or not. For duplicate positions, only the first indel record will be considered and appear on output.
#both
#abbreviation of "-c indels  -c snps"
#id
#only records with identical ID column are compatible. Supported by bcftools merge only.
#-f, --apply-filters LIST
#Skip sites where FILTER column does not contain any of the strings listed in LIST. For example, to include only sites which have no filters set, use -f .,PASS.
#--no-version
#Do not append version and command line information to the output VCF header.
#-o, --output FILE
#When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default.
#-O, --output-type b|u|z|v
#Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.
#-r, --regions chr|chr:pos|chr:from-to|chr:from-[,…]
#Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R.
#-R, --regions-file FILE
#Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file #are: CHROM, POS, and, optionally, POS_TO, where positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and POS_TO (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the ".bed" or ".bed.gz" suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, "chr20" is not the same as "20". Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r.
#-s, --samples [^]LIST
#Comma-separated list of samples to include or exclude if prefixed with "^". The sample order is updated to reflect that given on the command line. Note that in general tags such as INFO/AC, INFO/AN, etc are not updated to correspond to the subset samples. bcftools view is the exception where some tags will be updated (unless the -I, --no-update option is used; see bcftools view documentation). To use updated tags for the subset in another command one can pipe from view into that command.
###########################################################################################################################################################################
#Q2:74656
zcat input/processed/genotypes.chr22.vcf.gz | grep -v "#" | wc -l 
#Q3:2504 i 445
bcftools stats input/unprocessed/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bcftools stats input/processed/genotypes.chr22.vcf.gz




#Task3
print("Ejecutando Task3")

docker run -v $PWD:$PWD -w $PWD -it dgarrimar/eqtlmapping #Accedemos al docker
release=12
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$release/gencode.v$release.annotation.gtf.gz
mv gencode.v$release.annotation.gtf.gz input/unprocessed/gencode/gencode.annotation.gtf.gz
PATH=$PATH:$PWD/bin
zcat input/unprocessed/gencode/gencode.annotation.gtf.gz | grep "gene_type \"protein_coding\"\|gene_type \"lincRNA\"" | gtf2bed.sh > tmp/gencode.annotation.bed

#Q1:v12
#Q2:GRCh37
#Q3:19940
#Q4:PATH=$PATH:$PWD/bin

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$3-$2,$6}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$2+1,$4,$5,$6}else{print $1,$3-1,$3,$4,$5,$6}}' tmp/gencode.annotation.bed > tmpfile; mv tmpfile tmp/gencode.annotation.bed
sed -i "s/^chr//" tmp/gencode.annotation.bed



#Task4
print("Ejecutando Task4")

join -1 4 -2 1 -t $'\t' <(sort -k4,4 tmp/gencode.annotation.bed) <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | sort -k1,1) > tmp/joint.tsv
awk '$2==22' tmp/joint.tsv > tmp/joint.chr22.tsv
paste <(awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$5,$6}' tmp/joint.chr22.tsv) <(cut -f1-6 --complement tmp/joint.chr22.tsv) | sort -k1,1V -k2,2n > tmp/joint.chr22.bed
cat <(zcat input/unprocessed/geuvadis/geuvadis.gene_expr_rpkm.tsv.gz | head -1 | sed "s/TargetID/#chr\tstart\tend\tgene\tlength\tstrand/") tmp/joint.chr22.bed > tmp/genes.chr22.rpkm.bed

#Q1:protein coding
#Q2: Para evitar errores en la normalizacion
#Q3: Al hacer los plots podras ver que son iguales
#Q4: A traves de plots, (GAUSS)

normalize.R -i tmp/genes.chr22.rpkm.bed -o tmp/genes.chr22.norm.bed
bgzip tmp/genes.chr22.norm.bed
tabix -p bed tmp/genes.chr22.norm.bed.gz
mv tmp/genes.chr22.norm.bed.gz* input/processed

#Task5
print("Ejecutando Task5")

#Normalitzades
check.norm.R -i input/processed/genes.chr22.norm.bed.gz -o result/plots/check.norm.pdf

#Sense Normalitzar
check.norm.R -i tmp/genes.chr22.rpkm.bed -o result/plots/check.NOnorm.pdf

#Q1:Podemos ver que las normalizadas.......

#Task6
print("Ejecutando Task6")

head -1 input/unprocessed/1000g/1000g.phase3_metadata.txt  | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'
head -1 input/unprocessed/geuvadis/geuvadis.metadata.txt | awk 'BEGIN{FS=OFS="\t"}{for (i=1;i<=NF;i++) {print i, $i}}'

#Q1: Source Name, Comment [ENA_SAMPLE], Characteristics [Organism], Term Source REF, Term Accession Number, Characteristics [Strain], Characteristics [population], Protocol REF, Technology Type

#Task7
print("Ejecutando Task7")

#Q1:
#Q2: Cointiene Principal Components y el otro los porcentajes.

QTLtools pca --bed input/processed/genes.chr22.norm.bed.gz --scale --center --out result/expression 
QTLtools pca --vcf input/processed/genotypes.chr22.vcf.gz --scale --center --maf 0.05 --distance 50000 --out result/genotypes
pcaPlot.R -i result/expression -o result/plots/expression.pca.pdf
pcaPlot.R -i result/genotypes -o result/plots/genotypes.pca.pdf

#Q3: Grafico de expresion, distribucion homogena. En el otro podemos ver dos grupos claramente diferenciados.
#Q4: Population
pcaPlot.R -i result/genotypes --metadata input/unprocessed/1000g/1000g.phase3_metadata.txt --color super_pop --out result/plots/genotypes.pca.super_pop.pdf
join -j 1 -t $'\t' <(sort -k1,1 input/unprocessed/1000g/1000g.phase3_metadata.txt) <(cut -f1,20 input/unprocessed/geuvadis/geuvadis.metadata.txt | sort -k1,1 | uniq) > tmp/metadata.txt
sed -i '1s/^/sampleID\tpop\tsuper_pop\tgender\tlab\n/' tmp/metadata.txt
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m tmp/metadata.txt --formula "~ (1|gender) + (1|pop) + (1|lab)" -o result/plots/vp.pdf

#Q5: Lab, population, gender.

#Task8
print("Ejecutando Task8")

#Q1: Todos estan por debajo del 10%. PEER4 y PEER5 ronda el 1%.
peer.R -i input/processed/genes.chr22.norm.bed.gz -p 10 -o tmp/peer.tsv
var_partition.R -i input/processed/genes.chr22.norm.bed.gz -m <(paste tmp/peer.tsv tmp/metadata.txt) -f "~ (1|pop) + (1|lab) + PEER1 + PEER2 + PEER3 + PEER4 + PEER5" -o result/plots/vp.peer.pdf
join -j 1 -t $'\t' tmp/metadata.txt tmp/peer.tsv  | Rscript -e 'write.table(t(read.table(file("stdin", open = "r", blocking = T), h = F)), file = "input/processed/covariates.tsv", quote = F, sep = "\t", col.names = F, row.names = F)'
gzip input/processed/covariates.tsv

#Task9 
print("Ejecutando Task9")

QTLtools cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --nominal 0.01 --out result/nominals.txt

#Q1:Si, podrian estar relacionados.
#Q2:You have (on the surface) a set of well-behaved p-values, That flat distribution along the bottom is all your null p-values, which are uniformly distributed between 0 and 1.
pvdist.R -i result/nominals.txt --col 12 -o result/plots/pvdist.pdf


plink --ld rs4819926 rs5748768 --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/ld2 
#Q3:rs4819926 rs5748768

 #Haplotype     Frequency    Expectation under LE
 #  ---------     ---------    --------------------
 #         GC      0.232584                0.054095
 #         AC      0                       0.178489
 #         GG      0                       0.178489
 #         AG      0.767416                0.588927

 #   In phase alleles are GC/AG

#Task10
print("Ejecutando Task10")

for j in $(seq 1 16); do
  echo "cis --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --permute 1000 --chunk $j 16 --out result/permutations_$j.txt"
done | xargs -P4 -n14 QTLtools
cat result/permutations_*.txt > result/permutations.txt; rm result/permutations_*.txt 

R #Abrimos R

p <- read.table("result/permutations.txt")                                                      # Read input file
pdf("result/plots/pv-correlation.pdf",  paper = 'a4r', width = 9, height = 6)                   # Open PDF device
plot(p[, 18], p[, 19], xlab = "pv (perm)", ylab = "pv (beta)")                                  # Plot p-values
abline(0, 1, col = "red")                                                                       # Add red line 1=1
plot(-log10(p[, 18]), -log10(p[, 19]), xlab = "-log10 pv (perm)", ylab = "-log10 pv (beta)")    # Repeat in -log10 space to check the behaviour of the small p-values.
abline(0, 1, col = "red")
dev.off()                                                                                       # Close device
quit("no")    

#Task11
print("Ejecutando Task11")

mtc.R -n result/nominals.txt -p result/permutations.txt --method 'bonferroni' --alpha 0.05 --out tmp/bonferroni.txt
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'fdr' --alpha 0.05 --out tmp/fdr.txt
mtc.R -n result/nominals.txt -p result/permutations.txt --method 'perm-fdr' --alpha 0.05 --out result/eqtls.tsv

wc -l result/nominals.txt
wc -l result/eqtls.tsv

#Q1:37941 y 12032


#Task12
print("Ejecutando Task12")

eQTLviewer.R -i <(head -n 10 result/eqtls.tsv) -g input/processed/genotypes.chr22.vcf.gz -e input/processed/genes.chr22.norm.bed.gz -o result/plots/eQTLs_head.pdf --verbose

#Task13
print("Ejecutando Task13")

rsync -av rsync://ftp.ensembl.org/ensembl/pub/grch37/release-86/regulation/homo_sapiens/AnnotatedFeatures.gff.gz input/unprocessed/ensembl
zcat input/unprocessed/ensembl/AnnotatedFeatures.gff.gz | awk 'BEGIN{FS=OFS="\t"}{print $1, $4-1, $5, $9}' | sed -r 's/Name=([^;]+);.*/\1/' | grep -v '^GL' | sort -V > tmp/ERB.bed

for feat in $(cut -f4 tmp/ERB.bed | sort | uniq); do 
 bedtools merge -i <(grep -Fw $feat tmp/ERB.bed) -c 4 -o distinct
done > input/processed/ERB.collapsed.bed
sed -i "s/^chr//" input/processed/ERB.collapsed.bed

for feat in $(cut -f4 input/processed/ERB.collapsed.bed | sort | uniq); do
 QTLtools fenrich --qtl <(sed '1d' result/eqtls.tsv | awk '{print $9, $10-1, $10, $8, $1, "."}') --tss tmp/gencode.annotation.bed  --bed <(grep -Fw $feat input/processed/ERB.collapsed.bed) --out tmp/enrich.txt > /dev/null; echo "$(cat tmp/enrich.txt) $feat" 
done | grep -Fwv inf | grep -Fwv nan > result/enrichments.txt
plot.enrich.R -i result/enrichments.txt -o result/plots/enrich.pdf

#Q1: H3K36me3,marcador de histonas.
#Q2:Menor probabilidad de resultado

#Task14
print("Ejecutando Task14")

#Q1:
#intron_variant: 47%
#upstream_gene_variant: 15%
#downstream_gene_variant: 14%
#non_coding_transcript_variant: 11%
#NMD_transcript_variant: 6%
#regulatory_region_variant: 2%
#intergenic_variant: 2%
#non_coding_transcript_exon_variant: 1%
#3_prime_UTR_variant: 1%
#Other


#Q2:23.
#splice_acceptor_variant,   non_coding_transcript_variant,stop_gained, splice_donor_variant,   frameshift_variant,  NMD_transcript_variant.
#Q3:3

sed '1d' result/eqtls.tsv | cut -f8 | sort | uniq > tmp/eqtls_snps.tsv
head -n100 tmp/eqtls_snps.tsv > 100eqtls.tsv


#Task15
print("Ejecutando Task15")

# Generate a list of sGenes
cut -f1 result/eqtls.tsv | sed '1d' | sed 's/\..\+//' | sort | uniq > tmp/egenes.txt

# We will use as background all the genes (PC and lincRNA) in chr22
awk '{if($1==22) print $4}' tmp/gencode.annotation.bed | sed 's/\..\+//' | sort | uniq > tmp/bg.txt

#Q1:response to lipopolysaccharide, response to molecule of bacterial origin, hormone metabolic process.

#Task16
print("Ejecutando Task16")

# Generate input files for QTLtools rtc
grep -Fwf <(cut -f1 result/eqtls.tsv ) result/permutations.txt > tmp/rtc_input
cut -f4,7 input/unprocessed/gwas/gwas.catalog.hg19.bed > tmp/gwas_trait

# Download the file 'hotspots_b37_hg19.bed' from QTLtools website
wget http://jungle.unige.ch/QTLtools_examples/hotspots_b37_hg19.bed --directory-prefix tmp

# Remove 'chr' from chromosome names (-i option to modify the file 'in place')
sed -i 's/^chr//' tmp/hotspots_b37_hg19.bed

# Run RTC
QTLtools rtc --vcf input/processed/genotypes.chr22.vcf.gz --bed input/processed/genes.chr22.norm.bed.gz --cov input/processed/covariates.tsv.gz --hotspot tmp/hotspots_b37_hg19.bed --gwas-cis tmp/gwas_trait tmp/rtc_input --out result/rtc.txt

#Q1: 39
awk '{if($20>0.9) print $20}' result/rtc.txt | wc -l 
#Q2: awk '$20=="1"' result/rtc.txt
#Q3:!!!!!!!!!!!!!!100% intergenic variant

gene=ENSG00000099968.13
compZscore.R --gene $gene --nominal result/nominals.txt -k 50 --output tmp/$gene.rs_z
plink --r square --snps $(cut -f1 tmp/$gene.rs_z) --vcf input/processed/genotypes.chr22.vcf.gz --out tmp/$gene

#Task17
print("Ejecutando Task17")

CAVIAR -z tmp/$gene.rs_z -l tmp/$gene.ld -o result/$gene
wc -l result/ENSG00000099968.13_set
wc -l result/ENSG00000099968.13.log
#rs8137591
#rs8141347
#rs7290691
#rs1978967
#Q1:  rs8137591	0.37092	0.74184, rs8141347	0.127845	0.25569, rs7290691	0.292331	0.584661, rs1978967	0.206424	0.412848
awk '$1=="rs8137591"' result/ENSG00000099968.13_post
awk '$1=="rs8141347"' result/ENSG00000099968.13_post
awk '$1=="rs7290691"' result/ENSG00000099968.13_post
awk '$1=="rs1978967"' result/ENSG00000099968.13_post
#Q2:
 awk '{print $1, $12, $13}' <(awk '$1=="ENSG00000099968.13"' result/nominals.txt)
#Q3:
#intron_variant: 38%
#non_coding_transcript_variant: 26%
#NMD_transcript_variant: 17%
#missense_variant: 11%
#upstream_gene_variant: 4%
#non_coding_transcript_exon_variant: 2%
#3_prime_UTR_variant: 2%
                   
#Task18
print("Ejecutando Task18")

gene=ENSG00000099968.13
cat <(echo "MarkerName P.value") <(grep $gene result/nominals.txt | cut -d " " -f8,12) > tmp/metal.$gene

#Task19
print("Ejecutando Task18")








































