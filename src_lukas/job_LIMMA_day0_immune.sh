#!/bin/bash

snodes="s[001-004]"
bnodes="b[001-003],b[006-007],b[011-012],b[014-018]"
ncpus=8
mem="8G"
time="2:00:00"
partition="tinyq"

celltype="PBMC"
results="results_check"

for model in "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_GMCSF_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_MIP1b_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_24h_PBMC_IL.12_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_IL.10_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_IL.1b_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_IL.1ra_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_IL.2R_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_IL.6_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_IP10_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_MIG_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_MIP1a_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_TNF.a_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_24h_PBMC_lactate_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.C.albicans.yeast_7d_PBMC_IL.17_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.LPS.100ng_24h_PBMC_IL.1b_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.LPS.100ng_24h_PBMC_IL.6_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.LPS.100ng_24h_PBMC_TNF.a_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_24h_PBMC_IL.10_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_24h_PBMC_IL.1b_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_24h_PBMC_IL.1ra_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_24h_PBMC_IL.6_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_24h_PBMC_MIP1a_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_24h_PBMC_MIP1b_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_24h_PBMC_lactate_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_7d_PBMC_IFNg_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.MTB_7d_PBMC_IL.17_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.S.aureus_24h_PBMC_IL.1b_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.S.aureus_24h_PBMC_IL.6_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.S.aureus_24h_PBMC_TNF.a_good" "V1.batch.sex.age.blood.TSS_enr.visit_time.CYTO.S.aureus_7d_PBMC_IFNg_good" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.ADA_P00813" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.TNFSF14_O43557" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.TRAIL_P50591" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.TRANCE_O14788" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.TWEAK_O43508" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.VEGFA_P15692" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.XE.BP1_Q13541" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.uPA_P00749" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.ST1A1_P50225" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CXCL9_Q07325" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.AXIN1_O15169" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.Beta.NGF_P01138" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CASP.8_Q14790" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CCL11_P51671" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CCL19_Q99731" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CCL20_P78556" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CCL23_P55773" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CCL25_O15444" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CCL28_Q9NRJ3" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CCL3_P10147" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CCL4_P13236" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CD244_Q9BZW8" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CD40_P25942" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CD5_P06127" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CD6_Q8WWJ7" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CD8A_P01732" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CDCP1_Q9H5V8" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CSF.1_P09603" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CST5_P28325" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CX3CL1_P78423" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CXCL10_P02778" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CXCL11_O14625" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CXCL1_P09341" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CXCL5_P42830" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.CXCL6_P80162" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.DNER_Q8NFT8" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.EN.RAGE_P80511" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.FGF.19_O95750" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.FGF.21_Q9NSA1" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.FGF.23_Q9GZV9" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.Flt3L_P49771" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.GDNF_P39905" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.HGF_P14210" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL.10RA_Q13651" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL.10RB_Q08334" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL.12B_P29460" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL.17A_Q16552" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL.17C_Q9P0M4" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL.18R1_Q13478" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL10_P22301" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL18_Q14116" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL6_P05231" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL7_P13232" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.IL8_P10145" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.LAP.TGF.beta.1_P01137" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.LIF.R_P42702" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.MCP.1_P13500" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.MCP.2_P80075" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.MCP.3_P80098" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.MCP.4_Q99616" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.MMP.10_P09238" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.MMP.1_P03956" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.NT.3_P20783" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.OPG_O00300" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.OSM_P13725" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.PD.L1_Q9NZQ7" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.SCF_P21583" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.SIRT2_Q8IXJ6" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.SLAMF1_Q13291" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.STAMPB_O95630" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.TGF.alpha_P01135" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.TNFB_P01374" "V1.batch.sex.age.bmi.oralContra.blood.TSS_enr.visit_time.CM.TNFRSF9_Q07011"
do
   JOB_ID=1
   ls ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml

#  LIMMA
   job="limma_${celltype}.${model}"
   printf "${job}: "
   JOB_ID=$(sbatch --parsable --wrap "Rscript limma.R ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml $ncpus" \
   --job-name=${job} --output=logs/${job}.log \
   --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition --exclude=$bnodes)
   printf "Submitted batch job ${JOB_ID}\n"

#  VOLCANO
   job="volcano_${celltype}.${model}"
   printf "${job}: "
   JOB_ID=$(sbatch --parsable --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes --wrap "python DE_volcanos.py --config ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
   --job-name=${job} --output=logs/${job}.log \
   --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition)

#  ENRICHR
   job="Enrichr_${celltype}.${model}"
   printf "${job}: "
   sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes --wrap "python DE_enrichr.py --config ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml" \
   --job-name=${job} --output=logs/${job}.log \
   --mem=$mem --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition

#  LOLA
   lola_universe=resources/LOLA/lola_universe.bed  # resources/LOLA/lola_universe.TSS_PROXIMAL.bed
   region_filter=""  # "--regions TSS_PROXIMAL"
   job="LOLA_${celltype}.${model}"
   echo ${job}
   echo "Universe size:"
   wc -l ${lola_universe}
   sbatch --dependency=afterok:${JOB_ID} --kill-on-invalid-dep=yes \
   --wrap "python DE_lola.py --top_n 1000 --config ${results}/DE/${celltype}.${model}/de_config_${celltype}.${model}.yml --n_jobs ${ncpus} ${region_filter} --universe ${lola_universe}" \
   --job-name=${job} --output=logs/${job}.log \
   --mem=180G --cpus-per-task=$ncpus --time=$time --partition=$partition --qos=$partition --exclude=$bnodes

    echo ""
done

