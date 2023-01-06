# Scripts for sharing      
This page has been created as unstructured repository for script sharing.

## Getting started
In order to download the content of this repository, open your terminal window (Linux/Mac) and type in the following commands:
```
git clone https://github.com/pjawinski/scripts.git
cd scripts
```

## Scripts
### biotype.sh
- takes gene symbols or entrez ids as input and returns RefSeq information (gene name, chromosome, start and stop position, gene synonyms, gene biotype, gene description, hgnc ids, and entrez ids)
- RefSeq files in GFF3 format (hg19 release 105.20220307, hg38 release 110) will be downloaded automatically from NCBI
- biotype.sh expects three arguments: genome build [hg19/hg38], type of input [symbol/entrez], and path to input file
- the input file (mygenes.txt) is expected to contain one gene symbol / entrez id per line
- the output file will have the same path/name as the input file with ending .refseq.[hg19/hg38].txt

Example for running biotype.sh:
```
printf '%s\n' 'MAPT' 'APOE' 'DRD2' > mygenes.txt
./biotype.sh hg19 symbol mygenes.txt
```
