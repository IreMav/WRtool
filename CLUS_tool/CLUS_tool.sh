#!/bin/bash

#Read inputs variables:
source inputsCLUS_tool.sh

#Build the name of the outputs with the given inputs
mkdir -p "${OUTPUT_PATH}Log/"
name_outputs="${varname}${level}_${freq}_${model}_${kind}_${res}_${numens}ens_${season}_${area}_${syr}-${eyr}"

# Redirect stdout ( > ) into a named pipe ( >() ) running "tee"
exec > >(tee -i "${OUTPUT_PATH}Log/${name_outputs}_logfile.txt")

echo "----------------------------------------------------------------------"
echo "CLUS_tool.sh is running using list of variables from inputCLUS_tool.sh"
echo "----------------------------------------------------------------------"

# Print I/O directories
echo "Input directory is ${INPUT_PATH}" 

#Compute absolute input file names
absfilenames=( "${filenames[@]/#/$INPUT_PATH}" )
echo
echo "Input ${#filenames[@]} files (absolute path):"
printf '%s\n' "${absfilenames[@]}"
echo

echo "Output directory will be in ${OUTPUT_PATH}" 
echo "Output file names will be: <varname>_${name_outputs}"
echo

#======================PRECOMPUTATION
python ${WRTOOL_PATH}precompute.py "$WRTOOL_PATH" "$OUTPUT_PATH" "$name_outputs" "$varname" "$level" "$filterwnd" "$numens" "$season" "$area" "${absfilenames[@]}" 
#rm -rf $INPUT_PATH0*.nc

source activate cdms2
echo '============ cdms2 environment activated ============'

#======================EOF COMPUTATION
python ${WRTOOL_PATH}compute.py "$WRTOOL_PATH" "$OUTPUT_PATH" "$name_outputs" "$numens" "$numpcs" "$perc" "$enstoselect"

#======================CLUSTER ANALYSIS COMPUTATION
python ${WRTOOL_PATH}clustering.py "$WRTOOL_PATH" "$OUTPUT_PATH" "$name_outputs" "$numens" "$numpcs" "$perc" "$enstoselect" "$numclus" 

#======================CLUSTER ANALYSIS COMPARISON
python ${WRTOOL_PATH}clusters_comparison.py "$WRTOOL_PATH" "$OUTPUT_PATH" "$name_outputs" "$numens" "$numpcs" "$perc" "$enstoselect" "$numclus" "$nameREF"

#======================SIGNIFICANCE OF CLUSTER PARTITION
#python ${WRTOOL_PATH}clusters_sig.py "$WRTOOL_PATH" "$OUTPUT_PATH" "$name_outputs" "$numens" "$numpcs" "$perc" "$enstoselect" "$season" "$freq"

#======================PLOTS
python ${WRTOOL_PATH}clusters_plots.py "$WRTOOL_PATH" "$OUTPUT_PATH" "$name_outputs" "$numens" "$numpcs" "$perc" "$enstoselect" "$numclus"


source deactivate
echo '============ cdms2 environment deactivated ============'


