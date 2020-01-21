#!/bin/bash
for arg in "$@"
do 
	if [ "$arg" == "--inDir" ] || [ "$arg" == "-i"]
	then 
		inDir="$arg"
	fi
	if [ "$arg" == "--outDir"] || [ "$arg" == "-o"]
	then
		outDir="$arg"
done

if