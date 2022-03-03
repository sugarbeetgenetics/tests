
#fastq_files dir
fastq_files=$1

#genome dir
genome=$2

#output dir
output=$3

###########################################################################################

#index genome
#bwa index $genome
#echo "echo "$(date) -------------------------- Indexing Done!"

###########################################################################################

#run trimmomatic
Trimmomatic="java -jar /software/7/apps/trimmomatic/0.38/trimmomatic-0.38.jar"

base_name=$(ls "$fastq_files"/*fastq.gz | awk '{split($1,a,"/");split(a[length(a)],b,"_"); print b[1]"_"b[2]"_"b[3]"_"b[4]}' | sort | uniq)
#echo ${base_name}

for file in $base_name
    do
        read1=${file}_1.fq.gz
        read2=${file}_2.fq.gz

        java -jar $Trimmomatic PE ${fastq_files}/${read1} ${fastq_files}/${read2}  \
        ${fastq_files}/output_forward_paired_${read1} ${fastq_files}/output_forward_unpaired_${read1} \
        ${fastq_files}/output_reverse_paired_${read2} ${fastq_files}/output_reverse_unpaired_${read2} \
        ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		echo "$(date) ------------- ${file} ------------- Trimming Done"
    done

############################################################################################

for file in $base_name
    do
		#mapping
		bwa-mem2 -t 40 -k 50 -M ${genome} ${fastq_files}/output_forward_paired_${read1} \
		${fastq_files}/output_reverse_paired_${read2} > ${output}/${file}.sam
		
		samtools view -@ 40 -S -b ${output}/${file}.sam -o ${output}/${file}.bam
		samtools sort -@ 40 ${output}/${file}.bam > ${output}/${file}.sorted.bam
		samtools index ${output}/${file}.sorted.bam
		
		#if sorted bam file exists, delete sam and unsorted bam
		sam=${output}/${file}.sam
		bam=${output}/${file}.bam
		sorted_bam=${output}/${file}.sorted.bam
		
		if test -f "$sorted_bam"; then
			echo "$sorted_bam exists. Deleting $sam and $bam"
			rm $sam $bam
		fi
		
		echo "$(date) ------------- ${file} ------------- Mapping Done"
    done
