#!/bin/bash

# EXPERIMENTAL DESIGN: SELECTING THE NON-ANTARCTIC SPONGE MICROBIOME DATASET FROM THE SPONGE MICROBIOME PROJECT (SMP) DATA AND STRUCTURING THE METADATA TABLE #
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #

## This code must be executed by lines or statements, do not execute it in one go.


## The metadata related to the SMP samples ("SraRunInfo.csv" and "SraRunTable.txt" files) and the SRA accesion list was downloaded from the SRA download options  and the Run Selector (NCBI).

#--- Isolate SRA and BioSample IDs:
cut -f1,26 tab-SraRunInfo.csv
cut -f1 tab-SraRunInfo.csv > biosample-ID.txt # save BioSample ID list

#---  Download BioSample full reports:
while IFS= read -r line
do
	wget https://www.ncbi.nlm.nih.gov/biosample/$line?report=full\&format=text -O $line.txt
done < biosample-ID.txt

#--- Create a customized metadata table from the BioSample reports:
echo -e "BioSample\tSample type\tLife stage\tHost health\tHost\tHost name\tCollection date\tGeographic location\tSampling site\tLatitude\tLongitude\tWater temperature\tDepth" >> tmp.out
for file in *.txt
do
	a=$( grep -o -m1 "\w*SAMEA\w*" $file )
	b=$( grep "sample type=\"" $file )
	c=$( grep "life stage=\"" $file )
	d=$( grep "host health state=\"" $file )
	e=$( grep "host=" $file )
	f=$( grep "host common name=" $file )
	g=$( grep "collection" $file )
	h=$( grep -m1 "geographic location" $file )
	i=$( grep "sampling_site=" $file )
	j=$( grep "latitude=" $file )
	k=$( grep "longitude=" $file )
	l=$( grep "water_temperature_degrees_c=" $file )
	m=$( grep "depth=" $file )
	echo -e "$a\t$b\t$c\t$d\t$e\t$f\t$g\t$h\t$i\t$j\t$k\t$l\t$m" >> tmp.out
done

#--- Concatenate the cutomized metadata table with desired columns of the "SraRunInfo.csv" file:
echo -e "SRA ID\tBases\tAverage length\tSize(mb)\tBiosample verif." >> tmp2.out
while IFS= read -r line
do
	grep $line SraRunInfo.csv | cut -d "," -f1,5,7,8,26 --output-delimiter=$'\t' >> tmp2.out
done < biosample-ID.out
paste -d "\t" tmp.out tmp2.out > metadata-SMP.tsv

## The cutomized metadata table was uploaded to Google Docs. There, the table was filtered by sample type as "sponge tissue", host health as "healthy", life stage as "adult" or "sissing_not provided", and depth as "<50m". Also, we discarded samples displaying the columns "host" and "host name" other than the host scientific name (e.g. "freshwater metagenome", "seawater", "sponges"). Ulva sp. was removed as it is classified as an algae according to NCBI Taxonomy. The "host y "host name" columns were merged as a new column identifiying the sponge species. The filtered table was saved as "metadata-SMP-filt-tmp.tsv".

## We realized many sponge orders reported did not coincide with the sponge species reported. Many sponge orders did not agree the World Porifera database (https://www.marinespecies.org/porifera/) or the NCBI taxonomy database (https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) for the given species. Therefore, we corrected the assignments of sponge order according to the sponge species reported (trusting the appropiate morphological sponge classification during the experimental procedures of the SMP preceding analyses at the host species level).

#--- Search the sponge order for each sponge species reported using the NCBI taxonomy (FTP download, /pub/taxonomy/taxdump_archive directory, file "rankedlineage.dmp", accessed on 01-08-2021):
echo -e "Order" > tmp3.out
while IFS= read -r line
do
	grep -P -m1 "\t$line\t" rankedlineage.dmp | cut -d "|" -f6 | tr -d "\t" >> tmp3.out # orders
	grep -P -m1 "\t$line\t" rankedlineage.dmp >> tmp4.out # complete assignments
	#echo -e "$line\t\t$(grep -P -m1 "\t$line\t" rankedlineage.dmp | cut -d "|" -f6 | tr -d "\t")" >> tmp4.out # To know the provenance of empty rows
done < list-species-SMP-filt.txt
# NOTE: Some uncorrect names of sponge species in the NCBI taxonomy assingments were corrected according to the World Porifera database.

## The new assigments and non-assigned orders from NCBI were verified manually using the World Porifera database. Verification was also applied for the taxonomy of the published Antarctic sponge samples.

#--- Paste new sponge order column with the 
paste -d "\t" metadata-SMP-filt-tmp.tsv tmp4.out > metadata-SMP-filt.tsv


#--- Filter the customized metadata table by sponge orders matching the taxonomy of the available Antarctic sponges:
echo -e "BioSample\tSample type\tLife stage\tHost health\tSpecies\tOrder\tCollection date\tGeographic location\tSampling site\tLatitude\tLongitude\tWater temperature\tDepth\tSRA ID\tBases\tAverage length\tSize(mb)\tBiosample verif." > desired-orders.tsv
while IFS= read -r line
do
	grep -w $line metadata-SMP-filt.tsv >> desired-orders.tsv
done < list-orders.txt

## The customized metadata table "metadata-SMP-filt.tsv" was then classified according to the marine environment, realm, province, and ecoregion using the coordinates tables provided by the Marine Ecoregions of the World (MEOW, Spalding et al., 2007, https://doi.org/10.1641/B570707) of the R package "meowR" (see Rscript "script2_metadata.R"). 

#--- Choose randomly 5 sponge species replicates per environment.
for envi in trop temp
do
	echo -e "BioSample\tSample type\tLife stage\tHost health\tSpecies\tOrder\tCollection date\tGeographic location\tSampling site\tLatitude\tLongitude\tWater temperature\tDepth\tSRA ID\tBases\tAverage length\tSize(mb)\tBiosample verif." > $envi.tsv
	while IFS= read -r line
	do
		grep "$line" $envi-metadata-SMP.tsv | rl -c5 >> $envi.tsv
	done < list-$envi.txt
done
# NOTE: For species with less than 5 replicates, the maximum number of replicates available was reptrieved.

#--- Verify sponge genera and sponge species matching the Antarctic sponges taxa in each environment-classified table:

# For Species:
while IFS= read -r line
do
	if grep -n "$line" *.tsv | cut -f1,5 >> verif-desired-species.csv
	then
		echo "$line found" >> tmp5-species.out
	else
		echo "$line not found" >> tmp5-species.out
	fi
done < desired-species.txt
grep "not found" tmp5-species.out
grep -c "not found" tmp5-species.out

# For genus:
while IFS= read -r line
do
	if grep -n -w "$line" *.tsv | cut -f1,5 >> verif-desired-genera.csv
	then
		echo "$line found" >> tmp6-genera.out
	else
		echo "$line not found" >> tmp6-genera.out
	fi
done < desired-genera.txt
grep "not found" tmp6-genera.out
grep -c "not found" tmp6-genera.out

#--- List the SRA IDs by environment-classified table:
cut -f14 trop.tsv > trop-sra-IDs.txt
cut -f14 temp.tsv > temp-sra-IDs.txt

#--- Download the datasets of the resulting customized metadata tables using the SRA IDs:
for file in trop temp 
do
	while IFS= read -r line
	do
		prefetch -s $line >> $file-size-SRA-datasets.txt
		prefetch $line -O $file-SRA-datasets
		cd $file-SRA-datasets/
		fasterq-dump $line.sra -x -O . >> $file-details-SRA-datasets.txt
		echo "----------------------------------------------------" >> $file-details-SRA-datasets.txt
		cd ..
	done < $file-sra-IDs.txt
done

#--- Create a file with the IDs, environemnt type, and Illumina sequencer (sequencing run):
for file in trop temp
do
	while IFS= read -r line
	do
		a=$(grep $line tab-SraRunInfo.csv | cut -f20)
		b=$(cat SraRunTable.tsv | cut -d "," -f37,38 | tr "," "\t")
		echo -e "$line\t$file\t$a\t$b" >> TS-by-seq-runs.txt
	done < $file-sra-IDs.txt
done

## IMPORTANT: The metadata related to Antarctic sponges was prepared manually and independently.

#--- Download SRA dataset from Sacristán-Soriano et al., 2020 (Antarctic sponge dataset):
while IFS= read -r line
do
	#prefetch -s $line >> size-SRA-datasets.txt
	#prefetch $line -O SRA-datasets
	cd SRA-datasets/
	fasterq-dump $line.sra -x -O . >> details-SRA-datasets.txt
	echo "----------------------------------------------------" >> details-SRA-datasets.txt
	cd ..
done < SraAccList.txt

## The datasets of the rest of Antarctic sponges were retrieved from our lab repository for practicity, though they were publicly available in SRA, with the exception of the newly generated dataset (1821 sequencing run).

## The Antarctic, tropical, and temperate sponge tables were merged in a single metadata one (metadata-sponges.tsv). The complete taxonomic assignments were added manually to the customized metadata table in LibreOffice. The sequencing run/lanes IDs were codified in the metadada from the "SraRunTable.tsv" as follows: "h2000" or "h2500" depending on the Illumina Hiseq sequencing platform used, followed by "l" from the word lane and the number of the lane used. The sequencing run IDs of Antarctic sponge datasets were codified according to the sampling year or publication year as the identifier.

## IMPORTANT: Only sponge samples matching the taxonomy of the available Antarctic sponges at least at the family level were kept at this point, the rest was manually removed from the selected dataset and metadata ("metadata-sponges.tsv" file).Then, the dataset (fastq files) was divided in different folders according to the sequencing run ID. 


#--- Copy fastq files from a list with the IDs classified by sequencing run to a folder (only for non-Antarctic sponges; the Antarctic sponge files were divided in different folders manually):
for num in h2000l7 h2000l8 h2500l3 h2500l4 h2500l5 h2500l7
do
	for file in $(cat ./TS-hiseq$num.txt)
	do 
		cp "$file" ./TS-hiseq$num/
	done
done

#--- Replace sequecing run IDs in r scripts:
for num in h2000l7 h2000l8 h2500l3 h2500l4 h2500l5 h2500l7 2013cgb 2014 2015 2020ss 1821 
do
	cp dada2-qc-2013.r dada2-qc-$num.r
	cp dada2-dn-2013.r dada2-dn-$num.r
	#sed -i "s/demux/sponges/g" dada2-qc-$num.r
	sed -i "s/2013/$num/g" dada2-qc-$num.r
	sed -i "s/2013/$num/g" dada2-dn-$num.r
done

#--- Tranform fastq files into required format of demultiplexed single-end reads to use QIIME2 software:
array=( sponges-h2000l7 sponges-h2000l8 sponges-h2500l3 sponges-h2500l4 sponges-h2500l5 sponges-h2500l6 sponges-h2500l7 )
a=1

for file in "${array[@]}"
do
	cd $file/
	for i in *.fastq.gz
	do
		rename "s/_L001/_${a}_L001/g" $i
		#echo "hola_${i}_$a"
		a=$((a+1))
	done
	cd ../
done

#--- Use QIIME2 software to visualize quality plots:
array=( h2000l7 h2000l8 h2500l3 h2500l4 h2500l5 h2500l7 ) #single-end

for folder in "${array[@]}"
do
	qiime tools import \
	--type 'SampleData[SequencesWithQuality]' \
	--input-path sponges-$folder \
	--input-format CasavaOneEightSingleLanePerSampleDirFmt \
	--output-path sponges-$folder.qza >> info-qplots.txt 2>> error-qplots.txt

	qiime demux summarize \
	--i-data sponges-$folder.qza \
	--o-visualization summarize-qc-$folder.qzv >> info-qplots.txt 2>> error-qplots.txt
done

#--- Use QIIME2 software to visualize quality plots:
array=( 2013 2013cgb 2014 2015 1821 2020ss ) #paired-end

for folder in "${array[@]}"
do 
	qiime tools import \
	--type EMPPairedEndSequences \
	--input-path sponges-$folder \
	--output-path imported-$folder.qza >> info-qplots.txt 2>> error-qplots.txt
	
	qiime demux emp-paired \
	--m-barcodes-file metadata-$folder.tsv \
	--m-barcodes-column barcode-sequence \
	--p-rev-comp-mapping-barcodes \
	--i-seqs imported-$folder.qza \
	--o-per-sample-sequences sponges-$folder.qza \
	--o-error-correction-details dm-details-$folder.qza >> info-qplots.txt 2>> error-qplots.txt

	qiime demux summarize \
	--i-data sponges-$folder.qza \
	--o-visualization summarize-qc-$folder.qzv >> info-qplots.txt 2>> error-qplots.txt
done

#--- Remove primers with Cutadapt (only for sequencing run IDs 1821 and 2020s):

for file in *_L001_R1_001.fastq.gz
do
	cutadapt -g GTGYCAGCMGCCGCGGTAA -o ./outs/out-$file $file --match-read-wildcards -O 19 2>> details-cutadapt.log
done

#--- Obtain sequence length distribution of fastq files:
array=( h2000l7 h2000l8 h2500l3 h2500l4 h2500l5 h2500l6 h2500l7 2013 2013cgb 2014 2015 1821 2020ss )

for folder in "${array[@]}"
do
	awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}; print "total reads: " counter}' ./sponges-$folder/*.fastq > ./read-lengths/read-length-distrib-$folder.csv
done


#--- Search host type in the metadata of Lurgi et al., 2019 from species name to join the metadata of the study:
while IFS= read -r line
do
	a=$(grep -m1 "$line" Lurgi2019-suppl-data.csv | cut -f15)
	echo -e "$line\t$a" >> host-type.txt
done < TS-by-sponge-species.txt


## From this part, the datasets were processed in R environment (Rscripts).


