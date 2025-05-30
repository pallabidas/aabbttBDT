#!/bin/bash

#-----------------------------------------------
# This is the bash script called by job.sub
#-----------------------------------------------

export X509_USER_PROXY=$1
voms-proxy-info -all

export INPUT_FILE=$2
export OUTPUT_FILE=$3
export ISMC=$4
export ISEMBEDDED=$5

# Stage-in: eos cp the input files to the worker node
export EOS_MGM_URL=root://eosuser.cern.ch

COMMAND="bdteval inputFile=${INPUT_FILE} newOutputFile=1.0 newFile=${OUTPUT_FILE} isMC=${ISMC} isEmbedded=${ISEMBEDDED}"
echo "Attempting to execute: ${COMMAND}"
eval ${COMMAND}

