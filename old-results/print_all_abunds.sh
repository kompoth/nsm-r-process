#!/bin/bash
DIRS=$1

CMD="paste < ( awk 'printf (\"%4d\n\")'"
for d in $DIRS; do
  
done
paste <( awk '{print $1}' frdm/final_abundance.txt ) <( awk '{print $2}' frdm/final_abundance.txt ) <( awk '{print $2}' hfb/final_abundance.txt ) <( awk '{print $2}' ws4/final_abundance.txt ) <( awk '{print $2}' local2021jnt/final_abundance.txt )

{'printf ("%5s\t%s\n", $5, $1)'}
