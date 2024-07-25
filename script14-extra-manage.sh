#!/bin/bash

# EXTRA DATA MANAGEMENT #
# --------------------- #

## This code must be executed by lines or statements, do not execute it in one go.


#--- Search Habitat-generalists and habitat-specific from lists (used for the networks analysis):

# For Antarctic sponges:
while IFS= read -r line
do
	if grep -w "$line" asvs-spec.tsv
	then
	echo -e "$line\tHabitat-specific" >> tmp-spec-ant.txt
	else
	echo -e "$line\tOportunist/other" >> tmp-spec-ant.txt
	fi
done <  rownames-ant-netw.tsv

while IFS= read -r line
do
	if grep -w "$line" asvs-gen.tsv
	then
	echo -e "$line\tHabitat-generalist" >> tmp-gen-ant.txt
	else
	echo -e "$line\tOportunist/other" >> tmp-gen-ant.txt
	fi
done <  rownames-ant-netw.tsv

# For non-Antarctic sponges:
while IFS= read -r line
do
	if grep -w "$line" asvs-spec-noan.tsv
	then
	echo -e "$line\tHabitat-specific" >> tmp-spec-noa.txt
	else
	echo -e "$line\tOportunist/other" >> tmp-spec-noa.txt
	fi
done <  rownames-noa-netw.tsv

while IFS= read -r line
do
	if grep -w "$line" asvs-gen.tsv
	then
	echo -e "$line\tHabitat-generalist" >> tmp-gen-noa.txt
	else
	echo -e "$line\tOportunist/other" >> tmp-gen-noa.txt
	fi
done <  rownames-noa-netw.tsv

#--- Search for ASV sequences in "refseqs" file and put them in fasta files (Top 10 ASVs from SIMPER):
while IFS= read -r line
do
	a=$(grep -w "$line" refseqs.tsv | cut -f2)
	echo -e "> $line\n$a" >> drivers-asvs.fasta
done <  asvs_drivers.txt #list of top10 ASVs from SIMPER

sed -i "s/\"//g" drivers-asvs.fasta # TO UPLOAD IN BLAST NUCLEOTIDE (NCBI WEBSERVER) USING 16S rRNA FOR TAXONOMIC ANNOTATION

#--- Search for ASV sequences in "refseqs" file and put them in fasta files (habitat-specific and habitat-generalist ASVs):

for file in ant_specific noant_specific generalist
do
	sed -i "s/\"//g" asvs_$file.txt
	while IFS= read -r line
	do
		a=$(grep -w "$line" refseqs.tsv | cut -f2)
		echo -e "> $file | $line\n$a" >> habitat-asvs.fasta
	done <  asvs_$file.txt #list of habitat-specific and habitat-generalist ASVs
done
sed -i "s/\"//g" habitat-asvs.fasta # TO UPLOAD IN BLAST NUCLEOTIDE (NCBI WEBSERVER) USING 16S rRNA FOR TAXONOMIC ANNOTATION

#--- Search for ASV sequences in "refseqs" file and put them in fasta file (Top 10 high-degree ASVs in networks):

for file in ant noa
do
	sed -i "s/\"//g" tmp-code-$file.txt
	sed -i "s/-/_/g" tmp-code-$file.txt
	while IFS= read -r line
	do
		a=$(grep -w "$line" code-$file.tsv | cut -f1) # The translation of the numeric to alphabetic code
		b=$(grep -w "$a" refseqs.tsv | cut -f2)
		echo -e "> $file | $line | $a\n$b" >> hdegree-asvs.fasta #tmp-seqs-$file.txt
	done <  tmp-code-$file.txt # top10 high-degree ASVs
done
sed -i "s/\"//g" hdegree-asvs.fasta # TO UPLOAD IN BLAST NUCLEOTIDE (NCBI WEBSERVER) USING 16S rRNA FOR TAXONOMIC ANNOTATION


