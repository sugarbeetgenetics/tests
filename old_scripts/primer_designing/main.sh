#!/bin/bash
#$1 fasta file
#$2 amplification size
#$3 reference genome



seq_id=$(grep '^>' "$1")

IFS='>' read -ra ADDR <<< "$seq_id"
sequence_id="${ADDR[1]}"

IFS=' ' read -ra ADDR <<< "$sequence_id"
sequence_id2="${ADDR[0]}"

unset IFS

mkdir $sequence_id2
cd $sequence_id2  

cwd=$pwd

python3.6 /media/pflanz/Hs1-2/Scripts/primer_designing/gap_finder.py $1 $2
if [ $? -eq 0 ]; then
    echo 'Defining regions to be used for primer designing finished successfully!'
else
    echo 'gap_finder.py failed'
    exit
fi

echo processing $sequence_id2

#generate sequence_target values for primer3
awk -F',' '{print $1","$2-$1}' $1"_region_info.csv" > $1"_input_intervals"
tr '\r\n' ' ' < $1"_input_intervals" > $1"_output_intervals"
target_intervals=$(<$1"_output_intervals")

awk -F',' '{print $1","$2-$1}' $1"_stretchedgaps_info.csv" > $1"_stretchedgaps_intervals"
tr '\r\n' ' ' < $1"_stretchedgaps_intervals" > $1"_excluded_intervals"
exclude_intervals=$(<$1"_excluded_intervals")

###get sequence from seq file
sequence=$(grep -v '^>' "$1")
complete_sequence=${sequence//$'\n'/}
####generate primer3_config file based on previous computation
echo "SEQUENCE_ID=$sequence_id2
SEQUENCE_TEMPLATE=$complete_sequence
SEQUENCE_TARGET=$target_intervals
SEQUENCE_INTERNAL_EXCLUDED_REGIONS=$exclude_intervals
PRIMER_TASK=pick_sequencing_primers
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
primer_num_return=10
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_PRODUCT_SIZE_RANGE=75-2000
PRIMER_EXPLAIN_FLAG=1
=" > $sequence_id2"_primer3_config"

/home/pflanz/softwares/primer3/src/primer3_core --format_output --output=$sequence_id2"_primer3_output" $sequence_id2"_primer3_config"

if [ $? -eq 0 ]; then
    echo "Done!"
else
    echo "Failed"
    exit
fi

primer_info=$sequence_id2"_PRIMER_INFO"

grep '[LEFT|RIGHT]_PRIMER' $sequence_id2"_primer3_output" > $primer_info
awk '{print ">"NR"\n"$10}' $sequence_id2"_PRIMER_INFO" > $sequence_id2"_PRIMER_SEQ"
echo 'Done'
#~/ncbi-blast-2.7.1+/bin/makeblastdb -in $3 -dbtype nucl

~/ncbi-blast-2.7.1+/bin/blastn -db $3 -query $primer_info -word_size 7 -ungapped -outfmt 6 -out $sequence_id2"_blastn_output" -perc_identity 100 -max_target_seqs 2 -num_threads 12

primers_info=$sequence_id2
primers=$sequence_id2"_PRIMER_SEQ"
blastn=$sequence_id2"_blastn_output"

echo 'Done78'
python3.6 /media/pflanz/Hs1-2/Scripts/primer_designing/best_primer.py $primers $blastn $sequence_id2 $primers_info
if [ $? -eq 0 ]; then
    echo "Done81!"
else
    echo "Failed"
    exit
fi