# 1 / OBTAIN REFERENCE GENOME AND GENERATE A GENOME INDEX THROUGH STAR

mkdir /CULEX_TARSALIS_GENOME
cd /CULEX_TARSALIS_GENOME
# The genome should be here: https://osf.io/mdwqx/

# So I downloaded all files there:
# First I tried using wget: 
wget https://osf.io/mdwqx/

# But it would not work, it only downloads the html of the site.

# So I had to install a program specific for OSF:
pip install osfclient

# Then request to download files with the project ID, in this case: mdwqx
osf --project mdwqx clone

# Move the files downloaded and remove the download folder:
mv mdwqx/osfstorage/* .
rm -r mdwqx/

# Once we finally get the desired files:
gzip -dk Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3.gz
gzip -dk Culex-tarsalis_knwr_CONTIGS_CtarK1.fa.gz

# Now we have to generate an index file specific for STAR:

mkdir INDEXED_GENOME_CULEX_TARSALIS

STAR --runMode genomeGenerate --genomeDir ./INDEXED_GENOME_CULEX_TARSALIS --genomeFastaFiles Culex-tarsalis_knwr_CONTIGS_CtarK1.fa --sjdbGTFfile Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3 --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript ID --sjdbGTFtagExonParentGene Parent --runThreadN 12 --genomeChrBinNbits 12 --genomeSAindexNbases 13
# It took around 6 minutes: real 5m48.125s

# 2 / MAP RIBOSEQ READS TO THE GENOME THROUGH STAR

# The next step requires to align the reads from each file to the reference
# genome, using the index generated on the previous step
# The alignning program is STAR again and then we must specify several 
# parameters, such as the input files, the desired output folder and also 
# alignment-specific parameters such as how to deal with mismatches.

# The long command below performs three essential actions:

# 1/ stores on a variable the filename without extensions:
echo $filename | sed 's/\./ /g' | awk '{print $1}'
#    Please, be aware that when working with different datasets, the filtering
#    requirements may be different (sed/awk inputs)	

# 2/ generates a file-specific folder:
mkdir $output_folder/$current_folder

# 3/ Runs STAR with the desired parameters redirecting the output to the folder
#  generated on points 2 & 3
STAR --genomeDir $index_folder (...)

# The whole set of commands:
cd ..
input_folder=FASTQ_LINKS
index_folder=CULEX_TARSALIS_GENOME/INDEXED_GENOME_CULEX_TARSALIS/
mkdir STAR_ALIGNED_OUTPUT
output_folder=STAR_ALIGNED_OUTPUT
ls -1 FASTQ_LINKS/ | grep ribo | while read filename; do current_folder=$(echo $filename | sed 's/\./ /g' | awk '{print $1}'); mkdir $output_folder/$current_folder; STAR --genomeDir $index_folder --readFilesCommand zcat --quantMode GeneCounts --readFilesIn $input_folder/$filename --outFileNamePrefix $output_folder/$current_folder/$current_folder\_ --alignSJDBoverhangMin 1 --runThreadN 5 --seedSearchStartLmax 15 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 0 --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate; done

# With 'disown' at the end, so the process continue if we close the terminal:

ls -1 FASTQ_LINKS/ | grep ribo | while read filename; do current_folder=$(echo $filename | sed 's/\./ /g' | awk '{print $1}'); mkdir $output_folder/$current_folder; STAR --genomeDir $index_folder --readFilesCommand zcat --quantMode GeneCounts --readFilesIn $input_folder/$filename --outFileNamePrefix $output_folder/$current_folder/$current_folder\_ --alignSJDBoverhangMin 1 --runThreadN 5 --seedSearchStartLmax 15 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 0 --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate; done & disown 
