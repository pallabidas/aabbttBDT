#!/bin/bash

#--------------------------------------------------------
# Usage:
#   bash runCondorEval.sh
#
#--------------------------------------------------------
# path to .csv with list of samples to run on, formatted as: top-level directory to n-tuple (typically in EOS), sample name, config path, cross-section, and # of events in DAS
SAMPLES="hists-2017-benchmark.csv"
# optional description that only affects the job name in Condor
DESCRIPTION="May_skim-BDT-2017"
# skimpath: only used to find the number of jobs per sample, in postprocesssing some files might be missing (because of no trees) but we want the total number of files for the job indices: use the SVfit folder which has all files
SKIMPATH="/eos/cms/store/group/phys_susy/AN-24-166/pdas/HAA_svfit/2025-05-26-23h34m_2017/"
# submit: true (create .sub submit file and output directories, and submit with condor_submit) or false (only create .sub file and output directories, but do not submit)
SUBMIT=true

WORKING_DIR=${CMSSW_BASE}"/src/aabbttBDT/ROOT/bin"
EOS_BASE_DIR="/eos/cms/store/group/phys_susy/AN-24-166/${USER}/condorPostProcessing"

#--------------------------------------------------------
# Check if we used bash
#--------------------------------------------------------
if [[ "$0" != "$BASH_SOURCE" ]]; then
    echo ">>> ${BASH_SOURCE[0]}: Error: Script must be run with bash"
    return
fi

#--------------------------------------------------------
# Check if environments are set
#--------------------------------------------------------
if [[ -z ${CMSSW_BASE} ]]; then
    echo ">>> ${BASH_SOURCE[0]}: CMSSW environment is not set! Make sure to do cmsenv"
    exit 1
fi

#--------------------------------------------------------
# Get the voms-proxy-info certificate
#--------------------------------------------------------
export MYPROXYPATH="$(voms-proxy-info -path)"

# Use regex to extract the first letter of the username. ^ is the beginning of the string. \w means any alphanumeric character from the Latin alphabet. ( ) captures the variable
echo ${USER}
re="^(\w)"
# Bash's =~ operator matches the left hand side to the regex on the right hand side
[[ ${USER} =~ ${re} ]]
# The results of the regex can be extracted from the bash variable BASH_REMATCH. Index [0] returns all of the matches.
USER_FIRST_LETTER=${BASH_REMATCH[0]}


if [[ -f ${MYPROXYPATH} ]]; then
    echo ">>> runCondorEval.sh: Copying proxy from ${MYPROXYPATH} to /afs/cern.ch/user/${USER_FIRST_LETTER}/${USER}/private/x509up_file"
    cp ${MYPROXYPATH} /afs/cern.ch/user/${USER_FIRST_LETTER}/${USER}/private/x509up_file
else
    echo ">>> ${BASH_SOURCE[0]}: [ERROR]: x509 proxy not found on this machine, make sure voms-proxy-init was run, exiting"
    exit 1
fi

#--------------------------------------------------------
# First check if the input files exists
#--------------------------------------------------------
while IFS=, read -r SAMPLE_INFO PROCESS YEAR JOBSIZE ISMC ISEMBEDDED
do
    {
    # If stored in eos,
    if [[ ${SAMPLE_INFO} =~ ^/eos/ ]]; then
        if [[ ! $(ls -f ${SAMPLE_INFO}) ]]; then
            echo "${BASH_SOURCE[0]}: ERROR: File ${SAMPLE_INFO} is not there, aborting. (Make sure there are no extra newlines in ${SAMPLES}.)"
            exit 1
        fi
    else
    # Else use gfal-ls for remote files
        if [[ ! $(gfal-ls ${SAMPLE_INFO}) ]]; then
            echo "${BASH_SOURCE[0]}: ERROR: tried to gfal-ls file ${SAMPLE_INFO} but it is is not there, aborting. (Make sure there are no extra newlines in ${SAMPLES}.)"
            exit 1
         fi
    fi
    }
done < ${SAMPLES}

while IFS=, read -r SAMPLE_INFO SAMPLE YEAR JOBSIZE ISMC ISEMBEDDED
do
    {
    # Make a temporary .list of the input n-tuples
    INPUT_LIST="paths/${YEAR}/inputs_${SAMPLE}.list"
    echo ${INPUT_LIST}
    rm -r ${INPUT_LIST}
    mkdir -p "paths/${YEAR}/"
    # If the files are in /eos/, a simple "find" will work
    if [[ ${SAMPLE_INFO} =~ ^/eos/ ]]; then
        find $(dirname ${SAMPLE_INFO})/${SAMPLE}/*postprocessed*.root -type f > ${INPUT_LIST}
        cat ${INPUT_LIST}
    fi
}
done < ${SAMPLES}

#--------------------------------------------------------
# Then run the histogramming
#--------------------------------------------------------
start=$(date +%s)

DATETIME="$(date +"%Y-%m-%d-%Hh%Mm")"

# Create the Condor submit files
ABSOLUTE_PATH_TO_EXECUTABLE="${CMSSW_BASE}/bin/$SCRAM_ARCH/"

while IFS=, read -r SAMPLE_INFO SAMPLE YEAR ISMC ISEMBEDDED
do
{
    EOS_PREFIX="root://eosuser.cern.ch/"

    INPUT_LIST="paths/${YEAR}/inputs_${SAMPLE}.list"

    JOBLOG_DIR="logs/${DATETIME}_${DESCRIPTION}"

    EOS_DIR="${EOS_PREFIX}${EOS_BASE_DIR}/${DATETIME}_${DESCRIPTION}/${SAMPLE}/"

    # make the local submit directory (useful for checking the .sub file)
    # The EOS output directories are only created at the end of this .sh script if the execute flag is False.
    mkdir -p ${JOBLOG_DIR}/${SAMPLE}

    # Copy the template
    subfile="${JOBLOG_DIR}/h_${SAMPLE}.sub"
    cp jobTemplate.sub ${subfile}

    echo ">>> Writing to ${subfile}"

    JOB_FLAVOUR='"longlunch"'

    SAMPLE_FILE_NAME="postprocessed_ntuple_${SAMPLE}_\$(Process).root"
    INPUTDIR=$(dirname ${SAMPLE_INFO})
    INPUT_FILENAME="${EOS_PREFIX}${INPUTDIR}/${SAMPLE}/${SAMPLE_FILE_NAME}"
    OUTPUT_FILENAME="postprocessed_ntuple_${SAMPLE}_\$(Process).root"
    SUBLIST_NAME="paths/${YEAR}/inputs_${SAMPLE}.list"
    NJOBS=$(find ${SKIMPATH}${SAMPLE}/*.root | wc -l)
    echo "   >>> njobs: ${NJOBS}"

    # Edit parameters in the Condor .sub file
    sed -i "s|(job_flavour)|${JOB_FLAVOUR}|g" ${subfile}
    sed -i "s|(sample)|${SAMPLE}|g" ${subfile}
    sed -i "s|(input_file)|${INPUT_FILENAME}|g" ${subfile}
    sed -i "s|(timestamp)|${DATETIME}|" ${subfile}
    sed -i "s|(output_file)|${OUTPUT_FILENAME}|g" ${subfile}
    sed -i "s|(n_jobs)|${NJOBS}|g" ${subfile}
    sed -i "s|(isMC)|${ISMC}|g" ${subfile}
    sed -i "s|(isEmbedded)|${ISEMBEDDED}|g" ${subfile}

    # Paths to ship with the job
    sed -i "s|(absolute_path_to_executable)|${ABSOLUTE_PATH_TO_EXECUTABLE}|g" ${subfile}
    sed -i "s|(eos_dir)|${EOS_DIR}|g" ${subfile}

    mysubmitcommand="condor_submit -batch-name histo_${DESCRIPTION}_${DATETIME}_${SAMPLE} -file ${subfile} priority=100 "
    if [[ ${SUBMIT} == true ]]; then
        mkdir -p ${EOS_BASE_DIR}/${DATETIME}_${DESCRIPTION}/${SAMPLE}
        mkdir -p ${EOS_BASE_DIR}/${DATETIME}_${DESCRIPTION}/${SAMPLE}/logs/${DATETIME}_${DESCRIPTION}/${SAMPLE}
        eval ${mysubmitcommand}
    else
        echo ">>> Submit was false, did not submit"
    fi
}
done < ${SAMPLES}
wait
