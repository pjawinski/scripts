#!/bin/bash

# ======================================================================================================================
# === provide hgnc symbol or entrez id and get gene synonyms, gene biotype, gene description, hgnc id, and entrez id ===
# ======================================================================================================================
# January 6, 2023 | philippe.jawinski@hu-berlin.de

# check input arguments
if [ -z "$3" ]; then
    echo "Three arguments expected: genome build [hg19/hg38] + type of input [symbol/entrez] + path to input file"
    exit 1
fi
if [ "$1" != "hg19" ] && [ "$1" != "hg38" ]; then
	echo "Please specify genome build as first argument [hg19/hg39]"
    exit 1
fi
if [ "$2" != "symbol" ] && [ "$2" != "entrez" ]; then
	echo "Please specify the type of input as second argument [symbol/entrez]"
    exit 1
fi

# get variables
build=$1 # genome build [hg19|hg38]
input=$2 # type of input [symbol|entrez], 'symbol' for hgnc gene symbol; 'entrez' for entrez id
genelist=$3 # file that contains list of genes, e.g., genelist="mygenes.txt"
BASEDIR=$(dirname "$0")

# set weblink and file handler according to genome build
if [ "$build" == "hg19" ]; then
	refseq="$BASEDIR/GCF_000001405.25_GRCh37.p13_genomic.edit.gff.gz"
	weblink="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
	fileHandler="$BASEDIR/GCF_000001405.25_GRCh37.p13_genomic"
	ncbiTranslate=$(echo NC_000001.10$'\t'1$'\n'NC_000002.11$'\t'2$'\n'NC_000003.11$'\t'3$'\n'NC_000004.11$'\t'4$'\n'NC_000005.9$'\t'5$'\n'NC_000006.11$'\t'6$'\n'NC_000007.13$'\t'7$'\n'NC_000008.10$'\t'8$'\n'NC_000009.11$'\t'9$'\n'NC_000010.10$'\t'10$'\n'NC_000011.9$'\t'11$'\n'NC_000012.11$'\t'12$'\n'NC_000013.10$'\t'13$'\n'NC_000014.8$'\t'14$'\n'NC_000015.9$'\t'15$'\n'NC_000016.9$'\t'16$'\n'NC_000017.10$'\t'17$'\n'NC_000018.9$'\t'18$'\n'NC_000019.9$'\t'19$'\n'NC_000020.10$'\t'20$'\n'NC_000021.8$'\t'21$'\n'NC_000022.10$'\t'22$'\n'NC_000023.10$'\t'X$'\n'NC_000024.9$'\t'Y$'\n'NC_012920.1$'\t'MT) # # see https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_NCBI2UCSC.txt
elif  [ "$build" == "hg38" ]; then
	refseq="$BASEDIR/GCF_000001405.40_GRCh38.p14_genomic.edit.gff.gz"
	weblink="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
	fileHandler="$BASEDIR/GCF_000001405.40_GRCh38.p14_genomic"
	ncbiTranslate=$(echo NC_000001.11$'\t'1$'\n'NC_000002.12$'\t'2$'\n'NC_000003.12$'\t'3$'\n'NC_000004.12$'\t'4$'\n'NC_000005.10$'\t'5$'\n'NC_000006.12$'\t'6$'\n'NC_000007.14$'\t'7$'\n'NC_000008.11$'\t'8$'\n'NC_000009.12$'\t'9$'\n'NC_000010.11$'\t'10$'\n'NC_000011.10$'\t'11$'\n'NC_000012.12$'\t'12$'\n'NC_000013.11$'\t'13$'\n'NC_000014.9$'\t'14$'\n'NC_000015.10$'\t'15$'\n'NC_000016.10$'\t'16$'\n'NC_000017.11$'\t'17$'\n'NC_000018.10$'\t'18$'\n'NC_000019.10$'\t'19$'\n'NC_000020.11$'\t'20$'\n'NC_000021.9$'\t'21$'\n'NC_000022.11$'\t'22$'\n'NC_000023.11$'\t'X$'\n'NC_000024.10$'\t'Y$'\n'NC_012920.1$'\t'MT)
fi
	
# test if refseq file exists
if [ ! -f "$refseq" ] || [ ! -s "$refseq" ]; then
    echo $'\n'"Refseq file does not exist or is empty and will thus be downloaded from NCBI ftp server."

	# download refseq variants in .ggf format
	echo "  - downloading RefSeq file in GFF3 format."
	curl $weblink -s --output "${fileHandler}.gff.gz" 

	# get chromsome (translate), gene, biotype, description, synonyms from resources/GRCh37_latest_genomic.gff.gz
	echo "  - translating ncbi chromosome id and extracting gene name, gene biotype, gene description, and gene synonyms"
	awk -F'\t' 'NR==FNR { chr[$1]=$2; next }
		!($1 in chr) { next }
		$9!~/Name/ || $9!~/biotype/ { next }
		{ gene=$9; sub(/.*Name=/,"",gene); sub(/[;].*/,"",gene)}
		{ biotype=$9; sub(/.*gene_biotype=/,"",biotype); sub(/[;].*/,"",biotype) } $9!~/biotype/
		{ description=$9; sub(/.*description=/,"",description); sub(/[;].*/,"",description) } $9!~/description/ { description="" }
		{ synonym=$9; sub(/.*gene_synonym=/,"",synonym); sub(/[;].*/,"",synonym) } $9!~/synonym/ { synonym="" }
		{ chrom=chr[$1] }
		{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, chrom, gene, biotype, description, synonym }
		' OFS='\t' <(echo "${ncbiTranslate}") <(gzip -dc ${fileHandler}.gff.gz) > ${fileHandler}.edit.gff

	# sort by name, start, and stop coordinate and remove duplicate genes
	echo "  - sorting refseq entries by gene name, start, and stop coordinate and removing duplicates"
	cat ${fileHandler}.edit.gff | sort -k11,11 -k4,4n -k5,5n > ${BASEDIR}/tmp && \mv ${BASEDIR}/tmp ${fileHandler}.edit.gff

	# gzip
	echo "  - compressing file"
	chmod 770 ${fileHandler}.edit.gff
	gzip -f ${fileHandler}.edit.gff
	rm -f ${fileHandler}.gff.gz
fi

# remove ^M character and "GeneID:" (for entrez ids)
genelistVar=$(cat -v "$genelist" | sed -e "s/\^M//g" | sed -e "s/GeneID://g")

# remove .txt or .tsv ending
genelist=$(echo $genelist | sed 's/.txt$//g' | sed 's/.tsv$//g')

# if input = symbol
if [ "$input" = "symbol" ]; then

	# annotate biotype and description
	echo $'\n'"==================================================================="$'\n'"(1/3) Finding matches in RefSeq gene names and synonyms..."
	header=INPUT$'\t'NAME$'\t'CHR$'\t'START$'\t'STOP$'\t'SYNONYMS$'\t'BIOTYPE$'\t'DESCRIPTION
	awk -F'\t' -v header="$header" '
		BEGIN { print header }
		NR==FNR { 
			synonyms=$14; if(length(synonyms)==0) { synonyms = "NA" };
			gene[$11]=$11"\t"$10"\t"$4"\t"$5"\t"synonyms"\t"$12"\t"$13; next }
		$1 in gene { print $1, gene[$1]; next }
		{ print $1,"NA","NA","NA","NA","NA","NA","NA" }
		' OFS='\t' <(gzip -dc "$refseq") <(echo "$genelistVar") > "${genelist}.refseq.${build}.txt"

	# no match? get refseq rows that contain id as synonym
	missings=$(awk -F'\t' '$2 == "NA" { print $1}' "${genelist}.refseq.${build}.txt")
	awk -F'\t' -v header="$header" 'NR==FNR { missings[$1]; next }
		{ synonyms=$14; if(length(synonyms)==0) { synonyms = "NA" };
		  nSynonyms=split(synonyms,syns,","); for(k=1;k<=nSynonyms;++k) {
			if(syns[k] !~ /^[0-9]+$/ && syns[k] in missings) { 
				print syns[k], $11, $10, $4, $5, synonyms, $12, $13; delete missings[syns[k]] } } }
		' OFS='\t' <(echo "$missings") <(gzip -dc $refseq) > "${genelist}.synonyms.refseq.${build}.txt"

	# merge with annotation list
	nFindings=$(awk 'END { print NR}' "${genelist}.synonyms.refseq.${build}.txt")
	if [ "$nFindings" != 0 ]; then
		header=INPUT$'\t'NAME$'\t'CHR$'\t'START$'\t'STOP$'\t'SYNONYMS$'\t'BIOTYPE$'\t'DESCRIPTION
		awk -F'\t' -v header="$header" 'BEGIN { print header }
			NR==FNR { id[$1]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8; next }
			FNR==1 { next }
			$1 in id { print id[$1]; next }
			{ print $0 }' OFS='\t' ${genelist}.synonyms.refseq.${build}.txt "${genelist}.refseq.${build}.txt" > "${genelist}.refseq.${build}.tmp.txt"
		\mv "${genelist}.refseq.${build}.tmp.txt" "${genelist}.refseq.${build}.txt"
	fi

# if input = entrez
elif [ "$input" = "entrez" ]; then

# annotate biotype and desciption
echo $'\n'"==================================================================="$'\n'"(1/3) Finding matches in RefSeq entrez..."
header=INPUT$'\t'NAME$'\t'CHR$'\t'START$'\t'STOP$'\t'SYNONYMS$'\t'BIOTYPE$'\t'DESCRIPTION
awk -F'\t' -v header="$header" 'BEGIN { print header }
	NR==FNR && length($12) != 0 && $9 ~ /GeneID:/ { 
		synonyms=$14; if(length(synonyms)==0) { synonyms = "NA" };                           
		entrez=$9; gsub(/.*GeneID:/,"",entrez); gsub(/,.*/,"",entrez); gsub(/;.*/,"",entrez);
		gene[entrez]=$11"\t"$10"\t"$4"\t"$5"\t"synonyms"\t"$12"\t"$13; next }
	NR==FNR { next }
	$1 in gene { print $1, gene[$1]; next }
	{ print $1,"NA","NA","NA","NA","NA","NA","NA"}' OFS='\t' <(gzip -dc "$refseq") <(echo "$genelistVar") > "${genelist}.refseq.${build}.txt"

# end  if
fi

# add HGNC
echo "(2/3) Adding HGNC identifier..."
header=INPUT$'\t'NAME$'\t'CHR$'\t'START$'\t'STOP$'\t'SYNONYMS$'\t'BIOTYPE$'\t'DESCRIPTION$'\t'HGNC
awk -F'\t' -v header="$header" 'BEGIN { print header }
	NR==FNR { info=$9; gsub(/.*HGNC:HGNC:/,"HGNC:",info); gsub(/,.*/,"",info); { if(info ~ /^HGNC:/) { gene[$11]=info } }; next}
	FNR==1 { next } $2 in gene { print $0, gene[$2]; next }
	{ print $0, "NA"}' OFS='\t' <(gzip -dc $refseq) "${genelist}.refseq.${build}.txt" > "${genelist}.refseq.${build}.tmp.txt"
\mv "${genelist}.refseq.${build}.tmp.txt" "${genelist}.refseq.${build}.txt"

# add ENTREZ
echo "(3/3) Adding ENTREZ identifier..."
header=INPUT$'\t'NAME$'\t'CHR$'\t'START$'\t'STOP$'\t'SYNONYMS$'\t'BIOTYPE$'\t'DESCRIPTION$'\t'HGNC$'\t'ENTREZ
awk -F'\t' -v header="$header" 'BEGIN { print header }
	NR==FNR { info=$9; gsub(/.*GeneID:/,"GeneID:",info); gsub(/,.*/,"",info); { if(info ~ /^GeneID:/) { gene[$11]=info } }; next }
	FNR==1 { next } $2 in gene { print $0, gene[$2]; next }
	{ print $0, "NA"}' OFS='\t' <(gzip -dc $refseq) "${genelist}.refseq.${build}.txt" > "${genelist}.refseq.${build}.tmp.txt"
\mv "${genelist}.refseq.${build}.tmp.txt" "${genelist}.refseq.${build}.txt"

# summary
awk -F'\t' 'BEGIN { count_match=0; count_unique=0; count_hgnc=0; count_entrez=0 }
	NR==1 { next }
	$2!="NA" { count_match++ }
	$2!="NA" && !seen[$2]++ { count_unique++ }
	$9!="NA" { count_hgnc++ }
	$10!="NA" { count_entrez++ }
	END { print "Input id count:         "NR-1"\nRefSeq matches:         "count_match"\nUnique RefSeq genes:    "count_unique"\nHGNC ids assigned:      "count_hgnc"\nENTREZ ids assigned:    "count_entrez }' "${genelist}.refseq.${build}.txt"

# clean up
rm -f "${genelist}.synonyms.refseq.${build}.txt"
echo "Output written to ${genelist}.refseq.${build}.txt"$'\n'"==================================================================="$'\n'
