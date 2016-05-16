# Takes 3 arguments:
# $1 file containing full path to individual count files
# $2 allison's conversion file
# $3 transposed CDEs file with id, uuid, and gender

# Outputs 4 columns:
# Full path
# Filename
# File UUID
# Phen UUID

#i=files.txt
#j=htseq.txt
#k=cdes.txt

# First, pull out uuid from filename
if [ -f basename.txt ]; then 
	rm basename.txt;
fi

while read line; do
	a=$(basename $line) 
	echo $a >> basename.txt
done < $1

cut -f1 -d. basename.txt > file_id.txt

if [ -f conversion.txt ]; then 
	rm conversion.txt;
fi
# Next, match order of filename uuid with phenotype uuid
while read line; do
	grep $line $2 >> conversion.txt
done < file_id.txt

# Pull out resulting phenotype uuid column
cut -f3 conversion.txt > phen_uuid.txt

# Combine files to generate one conversion file
paste $1 basename.txt file_id.txt phen_uuid.txt > id_conversions.txt

rm basename.txt file_id.txt conversion.txt phen_uuid.txt
