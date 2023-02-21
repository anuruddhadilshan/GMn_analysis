#!/bin/bash

kine_num=$1
configfilename=$2
outputfilename=$3 #Don't give .root/.pdf extensions. Just give the first part and the extentions will be added by the script.

root -l -b -q 'definecuts.C+('$kine_num','"$configfilename"','"$outputfilename"')'
