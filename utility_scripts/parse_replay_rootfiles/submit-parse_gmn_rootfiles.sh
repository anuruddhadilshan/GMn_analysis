#!/bin/bash

# Setting necessary envioronments / workflow names (ONLY user specific part)
# Simply set the abosolute path to the directory where all the shell scripts, parse script, and config files are located. These shell script assumes they are all in the same directory.
export SCRIPT_DIR=/work/halla/sbs/adr/GMn_analysis/utility_scripts/parse_replay_rootfiles
workflowname='parse-gmn-rootfiles'

#cd $SCRIPT_DIR

configfilename=$1
run_on_ifarm=$2 #input 1 to run on ifarm and any other number to run on batch mode.

if [[ "$#" -ne 2 ]]; then
	echo -e "\n--!--\n Invalid number of arguments!!"
	echo -e "This script expects 2 arguments: <configfilename> <run_on_ifarm>\n"
	exit;
fi

configfile=$SCRIPT_DIR'/'$1

script=$SCRIPT_DIR'/run-parse_gmn_rootfiles.sh'

if [[ $run_on_ifarm -ne 1 ]]; then
	swif2 create $workflowname
	echo "Adding requested parsing job to swif2 "
	swif2 add-job -workflow $workflowname -partition production -cores 1 -disk 100GB -ram 32GB $script $configfile $SCRIPT_DIR
	swif2 run $workflowname
	echo -e "\n Getting workflow status. May take a few minutes...\n"
	swif2 status $workflowname
else 
	echo 'Running the parsing locally on ifarm.'
	$script $configfile $SCRIPT_DIR
fi



 