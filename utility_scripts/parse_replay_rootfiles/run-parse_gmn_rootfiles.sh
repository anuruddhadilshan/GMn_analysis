#!/bin/bash

# List of arguments
configfile=$1
scriptsdir=$2

root $scriptsdir'/parse_gmn_rootfiles.C('\"$configfile\"')'
