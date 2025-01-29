#!/bin/bash

cat id.txt | 
while read LINE
do
	# Get sample ID and tumor type
	ID=$(echo "$LINE" | cut -f 1)
	TYPE=$(echo "$LINE" | cut -f 2)
	echo "$ID"
	echo "$TYPE"

	# Run cytospace
	if [[ $TYPE == "TNBC" ]]
	then
	  ST_data="${ID}_ST_data.txt"
          ST_coor="${ID}_Coordinates.txt"
	  cytospace --scRNA-path ~/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/cytospace/tnbc/latest/tnbc_scRNA_data.txt \
		        --cell-type-path ~/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/cytospace/tnbc/latest/tnbc_cell_type_labels.txt \
	        	--st-path ~/Desktop/CIIR/data/ST_data/cytospace/latest/$ST_data \
		    	--coordinates-path ~/Desktop/CIIR/data/ST_data/cytospace/latest/$ST_coor \
		    	-o ~/Desktop/CIIR/results/cytospace/${ID} -sm lap_CSPR
	elif [[ $TYPE == "HER2+" ]]
	then
	  ST_data="${ID}_ST_data.txt"
          ST_coor="${ID}_Coordinates.txt"
	  cytospace --scRNA-path ~/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/cytospace/her2/latest/her2_scRNA_data.txt \
                --cell-type-path ~/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/cytospace/her2/latest/her2_cell_type_labels.txt \
                --st-path ~/Desktop/CIIR/data/ST_data/cytospace/latest/$ST_data \
                --coordinates-path ~/Desktop/CIIR/data/ST_data/cytospace/latest/$ST_coor \
                -o ~/Desktop/CIIR/results/cytospace/${ID} -sm lap_CSPR
	fi
done
