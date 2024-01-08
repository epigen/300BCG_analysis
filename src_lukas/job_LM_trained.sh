#!/bin/bash

analysis="donor_fixed_effects"
coefs="SAMPLE:VISIT"

for Y_name in "CYTO" "CM" "WB_PER_ML"
  do
    if [[ ${Y_name} == "CYTO" ]]; then
        blood="BLOOD"
    elif [[ ${Y_name} == "CM" ]]; then
              blood="WHOLE_BLOOD"
    else
      blood="None"
    fi
    job="LM_${Y_name}_${analysis}.${coefs}"
    RUN="--${analysis} --Y_name ${Y_name} --coefs ${coefs} --model final.${Y_name}.${blood}.${coefs} --blood ${blood} --remove_evening --save_outputs --return_corrected_X"
#    echo "${job}"
#    echo "${RUN}"
#    python linear_models.py ${RUN} > logs/${job}.log 2>&1
#    echo ""
    sbatch --wrap "python linear_models.py ${RUN}" \
    --job-name=${job} --output=logs/${job}.log \
    --mem="16G" --cpus-per-task=16 --time=2:00:00 --partition=tinyq --qos=tinyq
done
