
#bowtie2
bowtie="/home/ubuntu/softwares/bowtie/bowtie2-2.4.4-linux-x86_64/"

#fastq_files dir
fastq_files="/mnt/volume1/resistant_bulk_subset"

#genome dir
genome="/mnt/volume1/Genomes/EL10.2/Beta_vulgaris.fa"

#output dir
output="/mnt/volume1/test_output"

#trimmomatic

#index genome
#bwa index $genome
Trimmomatic="/home/ubuntu/softwares/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar"
base_name=$(ls "$fastq_files"/*fastq.gz | awk '{split($1,a,"/");split(a[length(a)],b,"_"); print b[1]}' | sort | uniq)
#echo ${base_name}

for file in $base_name
	do
		#echo "here1"
		read1=${file}_1_1000.fastq.gz
		read2=${file}_2_1000.fastq.gz
		
		java -jar $Trimmomatic PE ${fastq_files}/${read1} ${fastq_files}/${read2}  \
		${fastq_files}/output_forward_paired_${read1} ${fastq_files}/output_forward_unpaired_${read1} \
		${fastq_files}/output_reverse_paired_${read2} ${fastq_files}/output_reverse_unpaired_${read2} \
		ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	done

for file in $base_name
	do
		#mapping
		bwa mem -t 10 -k 50 -M ${genome} ${fastq_files}/output_forward_paired_${read1} \
		${fastq_files}/output_reverse_paired_${read2} > ${output}/${file}.sam
		samtools view -@ 10 -S -b ${output}/${file}.sam -o ${output}/${file}.bam
		samtools sort -@ 10 ${output}/${file}.bam > ${output}/${file}.sorted.bam
		samtools index ${output}/${file}.sorted.bam
	done

#variant calling
for file in $base_name
do
	bcftools mpileup -Ou -f ${genome} ${output}/${file}.sorted.bam \
	| bcftools call -Ou -mv | bcftools filter -s PASS -i \
	'%QUAL>20 && 15<DP<160' > ${output}/${file}.bcftools_minQ20_DP15_to_160_filter.vcf
done
