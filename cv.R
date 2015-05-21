library(signalHsmm)
library(seqinr)
library(cvTools)
library(hmeasure)
library(pbapply)

source("aa_groups.R")


# working directory ------------------------------

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"


pos_seqs <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)
neg_seqs <- read.fasta(paste0(pathway, "sept_neg.fasta"), seqtype = "AA")
#remove sequences with atypical aminoacids
atyp_aa <- which(sapply(neg_seqs, function(i) any(i %in% c("X", "J", "Z", "B", "U"))))
too_short <- which(sapply(neg_seqs, length) < 80)
neg_seqs <- neg_seqs[-unique(c(atyp_aa, too_short))]

too_short <- which(sapply(pos_seqs, length) < 80)
pos_seqs <- pos_seqs[-c(too_short)]

library(ROCR)

fold_res <- pblapply(1L:10, function(dummy) {
  pos_ids <- cvFolds(length(pos_seqs), K = 5)
  cv_neg <- neg_seqs[sample(1L:length(neg_seqs), length(pos_seqs))]
  lapply(all_groups, function(agg_group) {
    lapply(1L:5, function(fold) {
      model_cv <- train_hsmm(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] != fold]], agg_group)
      test_dat <- c(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] == fold]],
                    cv_neg[pos_ids[[4]][,][pos_ids[[5]] == fold]])
      preds <- cbind(t(sapply(predict(model_cv, test_dat), function(single_pred)
        c(prob = single_pred[["sp_probability"]], cs_pred = single_pred[["sp_end"]]))),
        cs_real = sapply(test_dat, function(i) 
          ifelse(is.null(attr(i, "sig")[2]), NA, attr(i, "sig")[2])))
      preds
    })
  })
})

save(fold_res, all_groups, file = paste0(pathway, "fold_res_df.RData"))


fold_res[[1]][[1]][[1]]

single_cv <- fold_res[[1]][[1]][[1]]

#ROCR prediction object
pred_obj <- prediction(single_cv[, "prob"], !is.na(single_cv[, "cs_real"]))

predicted_cl <- single_cv[, "prob"] > 0.2
real_cl <- !is.na(single_cv[, "cs_real"])

table(predicted_cl, real_cl)


HMeasure(!is.na(single_cv[, "cs_real"]), single_cv[, "prob"])[["metrics"]]
slot(prec, "x.values")[[1]][which.max(slot(prec, "y.values")[[1]] > 0.5)]
HMeasure(true.class, scores)