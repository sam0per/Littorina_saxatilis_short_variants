#!/bin/bash

# request memory for job (default 6G, max 72G)
#$ -cwd
#$ -V

base_cov=$(zcat test/aln.bam.coverage.gz | awk 'NR>1 { x += $3; } END { print x }')
nchunks=1000
zcat aln.bam.coverage.gz | head -100000000 | awk 'BEGIN { bin='$base_cov' / '$nchunks' } NR==1 { next } NR==2 { chr=$1; pos=$2; last=$2; } ($1==chr && sum < bin) { sum += $3; last=$2 } ($1!=chr || sum > bin) { print chr":"pos"-"last; sum = $3; chr=$1; pos=$2; last=$2; } END { print chr":"pos"-"last; } ' > targets.regions
