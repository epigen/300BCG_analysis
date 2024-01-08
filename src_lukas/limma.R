library(ggplot2)
library(dplyr)
library(yaml)
library(edgeR)
library(splines)
library(OpenMPController)
omp_set_num_threads(4)

args = commandArgs(trailingOnly=TRUE)
yaml_fn = args[1]
print(yaml_fn)
config = yaml.load_file(yaml_fn)
list2env(config, environment())

if (exists("block_on_batch") && block_on_batch){
    block_on_batch = TRUE
} else {
    block_on_batch = FALSE
}

if ((exists("dream") && dream) || (exists("var_part") && var_part)){
    library(BiocParallel)
    library(variancePartition)

    if (exists("dream") && dream){
        stopifnot(length(contrasts) == 1)
    }
    ncpus = args[2]
    param = SnowParam(ncpus, "SOCK", progressbar=TRUE)
    register(param)
} else {
    library(data.table)
    library(limma)
}

file_suffix = sprintf("%s.%s", celltype, model)

annotations = read.csv(annot_fn, check.names=TRUE, na.strings="")
colnames(annotations) = make.names(colnames(annotations))
rownames(annotations) = annotations$SAMPLE.ID
print("Annotations")
print(dim(annotations))


if(exists("remove_samples") && !is.null(remove_samples)){
    print("Removing specified samples")
    annotations = annotations[!(annotations$SAMPLE.DONOR %in% remove_samples),]
    print(dim(annotations))
}

if(exists("keep_samples") && !is.null(keep_samples)){
    print("Keeping specified samples")
    annotations = annotations[annotations$SAMPLE.DONOR %in% keep_samples,]
    print(dim(annotations))
}

if(remove_wrong_batch){
    print("Removing wrong batch samples")
    annotations = annotations[annotations$LAB.WRONG_BATCH == "False",]
    print(dim(annotations))
}
if(remove_exclusions){
    print("Removing exclusions")
    annotations = annotations[is.na(annotations$SAMPLE.EXCLUSION),]
    print(dim(annotations))
}
if(exists("only_responders") && only_responders){
    print("Keeping only responders")
    annotations = annotations[annotations$thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.2_responder == "R",]
    print(dim(annotations))
}
if(exists("visits") && !is.null(visits)){
    print("Keeping only these visits")
    print(visits)
    annotations = annotations[annotations$SAMPLE.VISIT %in% visits,]
    print(dim(annotations))
}

print("Keeping only non-null annotations")
samples = annotations[, columns]
samples = na.omit(samples)
print(dim(samples))

if(remove_300BCG315){
    print("Removing 300BCG315")
    samples = samples[rownames(samples) != "300BCG315_V1_PBMC",]
    print(dim(samples))
}

if(exists("remove_300BCG085") && remove_300BCG085){
    print("Removing 300BCG085")
    samples = samples[rownames(samples) != "300BCG085T3m_RPMI",]
    samples = samples[rownames(samples) != "300BCG085T3m_LPS",]
    samples = samples[rownames(samples) != "300BCG085T0_LPS",]
    samples = samples[rownames(samples) != "300BCG085T0_RPMI",]
    print(dim(samples))
}

if(remove_evening){
    print("Removing evening cohort")
    samples = samples[samples$DONOR.IC_TIME_REAL <= 13,]
    print(dim(samples))
}

print("Loading counts")
counts = read.csv(counts_fn, check.names=FALSE)
rownames(counts) = counts$ID
counts["ID"] = NULL
stopifnot(isNumeric(counts))
print(dim(counts))

print("Removing without ATAC-seq")
samples = samples[rownames(samples) %in% colnames(counts),]
print(dim(samples))

if(complete_design){
    print("Complete design")
    visit_counts = count(samples, SAMPLE.DONOR)
    donors = levels(droplevels(visit_counts[visit_counts$n==3,]$SAMPLE.DONOR))
    samples = samples[samples$SAMPLE.DONOR %in% donors, ]
    print(dim(samples))
}

if(useful_samples){
    print("Useful design")
    visit_counts = count(samples, SAMPLE.DONOR)
    donors = levels(droplevels(visit_counts[(visit_counts$n==3 | visit_counts$n==2),]$SAMPLE.DONOR))
    samples = samples[samples$SAMPLE.DONOR %in% donors, ]
    print(dim(samples))
}

if("time" %in% columns){
    samples$time = droplevels(samples$time)
    samples$time = relevel(samples$time, "T0")
}

if("stim" %in% columns){
    samples$stim = droplevels(samples$stim)
    samples$stim = relevel(samples$stim, "RPMI")
}

if("ts" %in% columns){
    samples$ts = droplevels(samples$ts)
    samples$ts = relevel(samples$ts, "T0_RPMI")
}

if("SAMPLE.VISIT" %in% columns){
    samples$SAMPLE.VISIT = droplevels(samples$SAMPLE.VISIT)
    samples$SAMPLE.VISIT = relevel(samples$SAMPLE.VISIT, "V1")
}
if("DONOR.SEX" %in% columns){
    samples$DONOR.SEX = droplevels(samples$DONOR.SEX)
    samples$DONOR.SEX = relevel(samples$DONOR.SEX, "F")
}
if("SAMPLE.TISSUE" %in% columns){
    samples$SAMPLE.TISSUE = droplevels(samples$SAMPLE.TISSUE)
    if(grepl(celltype, design)){
        print(paste("Encoding", celltype, "vs. the rest"))
        samples[[celltype]] = as.integer(samples$SAMPLE.TISSUE == celltype)
    } else {
        print(paste("Selecting only", celltype))
        samples = samples[samples$SAMPLE.TISSUE == celltype, ]
        print(dim(samples))
    }
}
if("LAB.BATCH" %in% columns){
    samples$LAB.BATCH = droplevels(samples$LAB.BATCH)
}
if("SAMPLE.DONOR" %in% columns){
    samples$SAMPLE.DONOR = droplevels(samples$SAMPLE.DONOR)
}
if("donor" %in% columns){
    samples$donor = droplevels(samples$donor)
}
if("DONOR.CIRCAD_REPLIC" %in% columns){
    samples$DONOR.CIRCAD_REPLIC = relevel(samples$DONOR.CIRCAD_REPLIC, "MOR")
}
if("DONOR.IC_TIME_CAT" %in% columns){
    samples$DONOR.IC_TIME_CAT = droplevels(samples$DONOR.IC_TIME_CAT)
    samples$DONOR.IC_TIME_CAT = relevel(samples$DONOR.IC_TIME_CAT, "MOR8_9")
}
for (score in c(
    'thm.innate_nonspecific_24h_V3_FC1.2_responder',
    'thm.innate_nonspecific_24h_wo_LAC_V3_FC1.2_responder',
    'thm.adaptive_MTB_7d_V3_FC1.2_responder',
    'thm.heterologous_nonspecific_7d_V3_FC1.2_responder',
    'thm.innate_MTB_24h_V3_FC1.2_responder',
    'thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.2_responder',
    'thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.25_responder',
    'thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.3_responder',
    'thm.scaled.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.2_responder',
    'thm.scaled.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.25_responder',
    'thm.scaled.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.3_responder',
    'thm.innate_nonspecific_24h_wo_LAC_IL10_IL1ra_V3_FC1.6_responder',
    'thm.heterologous_nonspecific_7d_V3_FC2.7_responder',
    'thm.adaptive_MTB_7d_V3_FC3.1_responder',
    'thm.IFNg_MTB_7d_V3_FC1.2_responder'
)){
    if(score %in% columns){
        # print('Releveling')
        # print(score)
        samples[[score]] = relevel(samples[[score]], "N")
    }
}

# Unify counts and samples
counts = counts[,rownames(samples)]

if(debug){
    print('DEBUG mode')
    counts = counts[1:100,]
    print(dim(counts))
}

if (exists("bootstrap") && !is.null(bootstrap)){
    print('bootstraping')
    set.seed(bootstrap)
    idx = sample(1:nrow(samples), size=nrow(samples), replace=TRUE)
    samples = samples[idx,]
    counts = counts[,idx]
    print(rownames(samples))
    print(colnames(counts))
}

if("SAMPLE.VISIT" %in% columns){
    print(count(samples,SAMPLE.VISIT))
}

if(exists("combat_seq") && combat_seq){
    library(sva)
    counts = ComBat_seq(counts, batch=samples$LAB.BATCH)
    write.csv(counts, file=gzfile(sprintf("%s/DE/%s/de_batch_corr_counts_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=rownames(counts))
}

if(!is.null(splines)){
    print('Making splines')
    print(splines)
    SPLINES = ns(samples[,splines], df=splines_df)
}

if(!exists('normal_data') || !normal_data){
    if(batch_corrected && ((exists("dream") && dream) || (exists("var_part") && var_part))){
        print("Loading raw counts")
        raw_counts = read.csv('../data/DE/quantification_filtered_PBMC.csv.gz', check.names=FALSE)
        rownames(raw_counts) = raw_counts$ID
        raw_counts["ID"] = NULL
        stopifnot(isNumeric(raw_counts))
        print(dim(raw_counts))
        raw_counts = raw_counts[,rownames(samples)]
        if(debug){
            raw_counts = raw_counts[1:100,]
        }
        stopifnot(identical(rownames(counts), rownames(raw_counts)))
        stopifnot(identical(colnames(counts), colnames(raw_counts)))
    }
    if(!batch_corrected || ((exists("dream") && dream) || (exists("var_part") && var_part))){
        dge = DGEList(if(batch_corrected && ((exists("dream") && dream) || (exists("var_part") && var_part))) raw_counts else counts,
                      samples=samples,
                      genes=rownames(if(batch_corrected && ((exists("dream") && dream) || (exists("var_part") && var_part))) raw_counts else counts))

        # print("BEFORE")
        # print(dge$samples$lib.size)
        # print(dge$samples$norm.factors)
        # libsizes = read.csv("../data/DE/libsizes.csv.gz", check.names=FALSE, row.names=1)
        # dge$samples$lib.size = libsizes[rownames(dge$samples), "libsize"]

        print("Calculating Normalization Factors")
        dge = calcNormFactors(dge, method="TMM")

        # print("AFTER")
        # print(dge$samples$lib.size)
        # print(dge$samples$norm.factors)

        if(save_non_log_CPM){
            non_log_cpm = t(t(dge$counts + 0.5) / (dge$samples$lib.size * dge$samples$norm.factors + 1) * 1e6)
            stopifnot(max(non_log_cpm) > 50)
            write.csv(non_log_cpm, file=gzfile(sprintf("%s/DE/%s/de_non_log_cpm_counts_%s.csv.gz", results_dir, file_suffix, file_suffix)))
        }
    }
}
if((exists("dream") && dream) || (exists("var_part") && var_part)){
    # DREAM and variancePartition
    if (!exists('normal_data') || !normal_data){
        print("Vooming while dreaming")

        if(block_on_donor || block_on_batch){
            block=if (block_on_donor) "+ (1|SAMPLE.DONOR)" else "+ (1|LAB.BATCH)"
            design = paste(design, block)
        }
        print(design)

        v = voomWithDreamWeights(dge, formula(design), dge$samples)
        save(v, dge, design, results_dir, file_suffix, file=sprintf("%s/DE/%s/de_voom_%s.RData", results_dir, file_suffix, file_suffix))
        print("Saved the voom object")
    }

    if (exists("dream") && dream){
        print("Fitting my dreams")
        fitmm = dream(v, formula(design), dge$samples, useWeights=TRUE, REML=TRUE, computeResiduals=FALSE)
        save(fitmm, v, file=sprintf("%s/DE/%s/de_lmfit_%s.RData", results_dir, file_suffix, file_suffix))
        print("Saved the fit")

        model_matrix = fitmm$design
        head(model_matrix)
        write.csv(model_matrix, file=gzfile(sprintf("%s/DE/%s/de_design_%s.csv.gz", results_dir, file_suffix, file_suffix)))

        for (coef in contrasts[[1]]){
            print(coef)
            results = topTable(fitmm, coef=coef, number=Inf, sort.by="none")
            write.csv(results, file=gzfile(sprintf("%s/DE/%s/de_results_%s_%s.csv.gz", results_dir, file_suffix, coef, file_suffix)))
        }
        save(fitmm, v, model_matrix, file=sprintf("%s/DE/%s/de_lmfit_%s.RData", results_dir, file_suffix, file_suffix))

    } else if (exists("var_part") && var_part){
        if (exists('normal_data') && normal_data){
            vp = fitExtractVarPartModel(counts, formula(design), samples, useWeights=FALSE, REML=FALSE, showWarnings=FALSE)
        } else {
            if(!batch_corrected){
                modelFit = fitVarPartModel(v, ~ RUN.TSS_ENRICHMENT + (1|LAB.BATCH), dge$samples)
                res = residuals(modelFit)
                vp = fitExtractVarPartModel(res, formula(design), dge$samples, useWeights=TRUE, REML=FALSE, showWarnings=FALSE)
            } else {
                vp = fitExtractVarPartModel(counts, formula(design), dge$samples, useWeights=TRUE, REML=FALSE,
                                            weightsMatrix=v$weights, showWarnings=FALSE)
            }
        }
        write.csv(vp, file=gzfile(sprintf("%s/DE/%s/de_variance_%s.csv.gz", results_dir, file_suffix, file_suffix)))
        fig = plotVarPart(sortCols(vp))
        ggsave(sprintf("%s/DE/%s/de_plotVarPart_%s.pdf", results_dir, file_suffix, file_suffix), fig)
        fig = plotPercentBars(sortCols(vp)[1:10,])
        ggsave(sprintf("%s/DE/%s/de_plotPercentBars_%s.pdf", results_dir, file_suffix, file_suffix), fig)
    }

} else {
    # LIMMA
    model_matrix = model.matrix(formula(design), samples)
    colnames(model_matrix) = make.names(
    gsub("\\(", "__",
    gsub("\\)", "",
    gsub("SAMPLE.VISITV1", "V1",
    gsub("SAMPLE.VISITV2", "V2",
    gsub("SAMPLE.VISITV3", "V3",
    gsub("SAMPLE.TISSUEPBMC", "PBMC",
    gsub("SAMPLE.TISSUEnkcell", "nkcell",
    gsub("SAMPLE.TISSUEmonocyte", "monocyte",
    gsub("SAMPLE.TISSUEcd8t", "cd8t",
    gsub("DONOR.IC_TIME_CAT", "",
    gsub("DONOR.CIRCAD_REPLIC", "",
    gsub("_responderR", "_R",
    gsub("_responderN", "_N",
    colnames(model_matrix)))))))))))))))
    head(model_matrix)
    if(exists("drop_cols")){
        model_matrix = model_matrix[,!colnames(model_matrix) %in% drop_cols]
        head(model_matrix)
    }
    write.csv(model_matrix, file=gzfile(sprintf("%s/DE/%s/de_design_%s.csv.gz", results_dir, file_suffix, file_suffix)))
    stopifnot(all(round(svd(model_matrix)$d, 6) != 0))

    if(!batch_corrected){
        if(!exists("cons_correlation")){
            print("vooming")
            pdf(sprintf("%s/DE/%s/de_mean_var_trend_%s.pdf", results_dir, file_suffix, file_suffix))
            v = voom(dge, model_matrix, plot=TRUE)
            x = dev.off()
        }
    } else{
        v = counts
    }

    if(block_on_donor || block_on_batch){
        block=if (block_on_donor) samples$SAMPLE.DONOR else samples$LAB.BATCH
        if(!exists("cons_correlation")){
            print("duplicateCorrelation")
            corr_fit = duplicateCorrelation(v, model_matrix, block=block)
            print(corr_fit$consensus)
            cons_correlation = corr_fit$consensus
        }
        if(!batch_corrected){
            print("vooming")
            pdf(sprintf("%s/DE/%s/de_mean_var_trend2_%s.pdf", results_dir, file_suffix, file_suffix))
            v = voom(dge, model_matrix, block=block, correlation=cons_correlation, plot=TRUE)
            x = dev.off()
            # there used to be if(save_all_results) but this makes no sense, probably a bug, removed it on 18 August 2022
            print("duplicateCorrelation")
            corr_fit = duplicateCorrelation(v, model_matrix, block=block)
            print(corr_fit$consensus)
            cons_correlation = corr_fit$consensus
        }
        print("fitting")
        lmfit = lmFit(v, model_matrix, block=block, correlation=cons_correlation)
    } else {
        print("fitting")
        lmfit = lmFit(v, model_matrix)
    }
    pdf(sprintf("%s/DE/%s/de_fitted_mean_var_trend_%s.pdf", results_dir, file_suffix, file_suffix))
    plotSA(lmfit)
    x = dev.off()

    # save(lmfit, v, model_matrix, file=sprintf("%s/DE/%s/de_lmfit_%s.RData", results_dir, file_suffix, file_suffix))
    # save(lmfit, file=sprintf("%s/DE/%s/de_lmfit_%s.RData", results_dir, file_suffix, file_suffix))

    make_contrasts = function(...){makeContrasts(..., levels=colnames(model_matrix))}
    for (c in contrasts){
        contr.matrix = do.call(make_contrasts, args=as.list(c))
        if(batch_corrected){
            print('Running eBayes with trend=TRUE')
        }
        contr_fit = eBayes(contrasts.fit(lmfit, contrasts=contr.matrix), trend=batch_corrected, robust=TRUE)

        if (exists("treat_FC") && !is.null(treat_FC)){
            print(as.list(c))
            for (trFC in treat_FC){
                print(log2(trFC))
                tfit <- treat(contr_fit, lfc=log2(trFC))
                treatResults = topTreat(tfit, coef=c, number=Inf, sort.by="none")
                write.csv(treatResults, file=gzfile(sprintf("%s/DE/%s/de_treat_%s_results_%s.csv.gz", results_dir, file_suffix, trFC, file_suffix)))
            }
        }
        dt = decideTests(contr_fit, p=0.1)
        print(summary(dt))
        write.fit(contr_fit, file=gzfile(
            sprintf("%s/DE/%s/de_results_p5_%s%s.csv.gz", results_dir, file_suffix, file_suffix,
                    if (length(contrasts) > 1) sprintf("_%s", paste(c, collapse="_")) else "")), sep =",", results=dt)
    }

    if(!batch_corrected){
        write.csv(v$design, file=gzfile(sprintf("%s/DE/%s/de_design_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=rownames(v$targets))
        #colnames(v$weights) = colnames(v$E)
        #write.csv(v$weights, file=gzfile(sprintf("%s/DE/%s/de_weights_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=v$genes[,1])
        #write.csv(dge$samples[c('lib.size','norm.factors')], file=gzfile(sprintf("%s/DE/%s/de_tmm_factors_%s.csv.gz", results_dir, file_suffix, file_suffix)))
        #write.csv(v$E, file=gzfile(sprintf("%s/DE/%s/de_counts_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=v$genes[,1])
        #write.csv(residuals(lmfit, v$E), file=gzfile(sprintf("%s/DE/%s/de_residuals_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=v$genes[,1])
        # write.csv(lmfit$coefficients, file=gzfile(sprintf("%s/DE/%s/de_coefficients_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=v$genes[,1])
        #write.csv(lmfit$df.residual, file=gzfile(sprintf("%s/DE/%s/de_df_residual_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=v$genes[,1])
    } else {
        write.csv(model_matrix, file=gzfile(sprintf("%s/DE/%s/de_design_%s.csv.gz", results_dir, file_suffix, file_suffix)))
        #write.csv(residuals(lmfit, counts), file=gzfile(sprintf("%s/DE/%s/de_residuals_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=rownames(counts))
        # write.csv(lmfit$coefficients, file=gzfile(sprintf("%s/DE/%s/de_coefficients_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=rownames(counts))
        #write.csv(lmfit$df.residual, file=gzfile(sprintf("%s/DE/%s/de_df_residual_%s.csv.gz", results_dir, file_suffix, file_suffix)), row.names=rownames(counts))
    }
}
print('Done.')
