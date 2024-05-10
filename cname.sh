#!/bin/bash

# Script to rename croco files to be compatible for ROMS OFFLINE


# Jonathan Rogerson
# August 2021

# SWITCH

# Set to 1 to name CROCO in line with ROFF by appending _N[value] to avg files
# Set to 0 to undo

SWITCH=1

# Set CROCO years

START=2002
END=2002

################################
# MAIN
################################

if [ ${SWITCH} -eq 1 ]

then

X=1;

for year in `seq $START $END`; do
 for month in `seq 1 12`; do
  for i in croco_avg_Y${year}M${month}.nc; do
   mv $i $(printf croco_N%02d.%s ${X%.*} ${i##*.})
   let X="$X+1"
   done
 done
done
 
   
elif [ ${SWITCH} -eq 0 ]; then

FILES=()
for f in `ls croco_N*.nc`; do
	FILES+=($(printf $f)) ;
done

result=()
for year in `seq $START $END`; do
  for month in `seq 1 12`; do
            result+=($(printf croco_avg_Y${year}M${month}.nc)) ;
 done
done

## get length of $FILES array
len=${#FILES[@]}

for (( i=0; i<$len; i++ )); do 
	mv ${FILES[$i]} ${result[$i]}  
done 

else

  echo "Not valid input"
 
fi 
