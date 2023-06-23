**parse_replay_rootfiles** repository contains scripts to parse replayed ROOT files from the GMn/nTPE (should work for other SBS experiments such as GEn) to cut away obviously bad events for physics analysis from the main "T" event TTree. The motivation is to save time during the iterative physics analysis by not having to loop through all the "junk" data and loop through only a subset of events that have passed a set of relatively loose global cuts. Care should be take to choose these glbal cuts "loose enough", so that you will be not throwing away good/useful events for your analysis.

Credit goes to Sebastian Seeds for the initial idea of using parsed-down ROOT files for analysis and John Boyd for the first vesion of this script.

## 1.How it wokrs:
1. You will provide the program with a set of global cuts, along with the other useful information like run parameters and paths to input and output rootfiles and etc. within a configuration file.
2. The program will search and find the CODA runs for the specific run paramers and finds the replayed ROOT files either in the user provided locations or in the generic locations in /work or /volatile.
3. For each ROOT file, the program will make a copy of the "T" TTree parsed as per the input cut conditions and the other three TTres, "E", "TSLeft", and "TSsbs" are copied one to one, in a separate directory specified by the user. It is recommended to have this direcory in /volatile.

## 2. What each script does:
1. `parse_gmn_rootfiles.C`: This is the main script that does the parsing of ROOT files. User should not have to modify anything inside this file unless they are experimenting.
2. `<user given name>.cfg`: The configuration file where you will input all the necessary information like the global cuts, run parameters, and input and output ROOT file directories.
3. `run-parse_gmn_rootfiles.sh`: The shell script that runs the `parse_gmn_rootfiles.C`. Use should not have to moddify this file.
4. `submit-parse_gmn_rootfiles.sh`: The shell script that the user will have to execute to submit the parsing jobs either into the ifarm or into the batch farm.

## 3. Quick start guide:
1. Inside the directory where all the above mentioned scripts exists, create the `<user given name>.cfg` configuration file, following the provided `parse_gmn_rootfiles.cfg` file, and change the input fields as required.
2. Open the `submit-parse_gmn_rootfiles.sh` file and set the fields "SCRIPT_DIR" and "workflowname" with the user specific ones.
3. Run the command: `$ submit-parse_gmn_rootfiles.sh <configfilename> 1` if it needs to be run on ifarm or replace 1 with any other integer to run on batch farm.

## 4. Contact:
For any questions or concerns, contact,
>Anuruddha Rathnayake 
><adr4zs@virginia.edu>  