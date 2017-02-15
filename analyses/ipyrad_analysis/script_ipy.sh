#!/bin/bash
file_array=($(find . -name 'params*'))  # find all parameters files
for i in $file_array  # loop to run the demultiplexing step of ipyrad on each params file
	do
	/Users/srlab/miniconda2/bin/ipyrad -p $i -s 1
done