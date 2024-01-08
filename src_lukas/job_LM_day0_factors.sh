#!/bin/bash

analysis="V1"
extras="DONOR:AGE DONOR:SEX DONOR:BMI DONOR:oralContraceptivesIncludingMen SAMPLE:alcoholInLast24h SAMPLE:VISIT_TIME_REAL"

for Y_name in "CYTO" "CM" "WB_PER_ML"
do
  for factors in "include_host_factors" "include_blood_perc" "include_blood_per_ml"
  do
    if [[ ${factors} == "include_host_factors" ]]; then
      if [[ ${Y_name} == "CYTO" ]]; then
        blood="BLOOD"
      elif [[ ${Y_name} == "CM" ]]; then
        blood="WHOLE_BLOOD"
      elif [[ ${Y_name} == "WB_PER_ML" ]]; then
        blood="None"
      fi

    elif [[ ${factors} == "include_blood_perc" ]]; then
      if [[ ${Y_name} == "CYTO" ]]; then
        blood="None"
      else
        continue
      fi

    elif [[ ${factors} == "include_blood_per_ml" ]]; then
      if [[ ${Y_name} == "CM" ]]; then
        blood="None"
      else
        continue
      fi
    fi

    job="LM_${Y_name}_${analysis}${factors}"
    RUN="--${analysis} --${factors} --Y_name ${Y_name} --model final.${Y_name}.${blood} --blood ${blood} --extras ${extras}"
#    echo "${job}"
#    echo "${RUN}"
#    python linear_models.py ${RUN} > logs/${job}.log 2>&1
#    echo ""
    sbatch --wrap "python linear_models.py ${RUN}" \
    --job-name=${job} --output=logs/${job}.log \
    --mem="16G" --cpus-per-task=16 --time=2:00:00 --partition=tinyq --qos=tinyq
  done
done
