#!/bin/bash

analysis="visit_interaction"

for Y_name in "CYTO" "CM" "WB_PER_ML"
do
  for factors in "include_host_factors" "include_blood_per_ml" "include_CM"
  do
    if [[ ${factors} == "include_host_factors" ]]; then
      if [[ ${Y_name} == "CYTO" ]]; then
        blood="BLOOD"
      elif [[ ${Y_name} == "CM" ]]; then
        blood="WHOLE_BLOOD"
      elif [[ ${Y_name} == "WB_PER_ML" ]]; then
        blood="None"
      fi

    elif [[ ${factors} == "include_blood_per_ml" ]]; then
      if [[ ${Y_name} == "CYTO" || ${Y_name} == "CM" ]]; then
        blood="None"
      else
        continue
      fi

    elif [[ ${factors} == "include_CM" ]]; then
      if [[ ${Y_name} == "CYTO" ]]; then
        blood="None"
      else
        continue
      fi
    fi

    job="LM_${Y_name}_${analysis}${factors}"
#    --visit_interaction_and_correct_season corrects for seasonal interaction, can be used with --include_host_factors
    RUN="--${analysis} --${factors} --Y_name ${Y_name} --model final.${Y_name}.${blood} --blood ${blood} --LMM --remove_evening"
#    echo "${job}"
#    echo "${RUN}"
#    python linear_models.py ${RUN} > logs/${job}.log 2>&1
#    echo ""
    JOB_ID=$(sbatch --parsable --wrap "python linear_models.py ${RUN}" \
    --job-name=${job} --output=logs/${job}.log \
    --mem="16G" --cpus-per-task=16 --time=2:00:00 --partition=tinyq --qos=tinyq)
    printf "Submitted batch job ${JOB_ID}\n"

    if [[ ${factors} == "include_host_factors" ]]; then
      job="fix_interaxtion_${Y_name}"
      sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
      --wrap "python fix_interaxtion_linear_models.py ${Y_name}" \
      --job-name=${job} --output=logs/${job}.log \
      --mem="16G" --cpus-per-task=16 --time=2:00:00 --partition=tinyq --qos=tinyq
    fi

  done
done
