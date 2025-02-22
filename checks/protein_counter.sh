#!/bin/bash

#grep ">" "$1" | sed -n 's/.*OS=\(.*\)OX.*/\1/p' | sort | uniq -c > protein_counts.csv

grep ">" "$1" | awk -F'OS=' '{print $2}' | awk -F'OX=' '{print $1}' | sort | uniq -c > protein_counts.csv

