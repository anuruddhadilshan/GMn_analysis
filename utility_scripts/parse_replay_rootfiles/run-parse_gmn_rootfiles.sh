#!/bin/bash


# List of arguments
passnum=$1
kinnum=$2
sbsfieldscale=$3
target=$4
output_dir=$5
output_rootfile_name=$6

root 'parse_gmn_rootfiles.C+('$passnum', '$kinnum', '$sbsfieldscale', '\"$target\"', '\"$output_dir\"', '\"$output_rootfile_name\"')'
