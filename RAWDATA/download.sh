#!/bin/bash

#$ -q all.q
#$ -cwd

for url in $(grep "a href" data.html |sed -r 's/\s+<a href = //'|sed 's/>.*//'|grep fastq|sed "s/'//g")
do
 wget --user= --password= $url
done

