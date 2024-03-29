# Scripts for sharing      
This page has been created as unstructured repository for script sharing.

## Getting started
In order to download the content of this repository, open your terminal window (Linux/Mac) and type in the following commands:
```
git clone https://github.com/pjawinski/scripts.git
cd scripts
```

## RefSeq annotations 
- [biotype.sh](biotype.sh) takes gene symbols or entrez ids as input and returns RefSeq information (gene name, chromosome, start and stop position, gene synonyms, gene biotype, gene description, hgnc ids, and entrez ids)
- RefSeq files in GFF3 format (hg19 release 105.20220307, hg38 release 110) will be downloaded automatically from NCBI
- biotype.sh expects three arguments: genome build [hg19/hg38], type of input [symbol/entrez], and path to input file
- the input file (mygenes.txt) is expected to contain one gene symbol / entrez id per line

Example for running biotype.sh:
```
printf '%s\n' 'MAPT' 'APOE' 'DRD2' > mygenes.txt
./biotype.sh hg19 symbol mygenes.txt
```
- the output file will be located in the same folder as the input file with ending .refseq.[hg19/hg38].txt

<br>

## PGS height
- [pgsheight.sh](pgsheight.sh) and [pgsheight.R](pgsheight.R) can be used create an animated plot of height vs. polygenic score of height
- feel free to adapt these code snippets for your own purposes
- the animated plot has been posted on twitter: [follow this link](https://twitter.com/PJawinski/status/1638137015421276163)
- see below for a static plot
<img src="pgsheight.png" width="600" unselectable="on">
