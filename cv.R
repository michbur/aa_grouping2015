library(signalHsmm)
library(seqinr)
library(cvTools)
library(hmeasure)
library(pbapply)

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

pos_ids <- cvFolds(length(pos_seqs), K = 5)
neg_ids <- cvFolds(length(neg_seqs), K = 5)

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



fold_res_df <- t(sapply(fold_res, function(i) rowMeans(matrix(unlist(i), ncol = 5))))
colnames(fold_res_df) <- c("AUC", "H", "Gini", "Sens", "Spec", "FPR", "mean_cs", "med_cs")

save(fold_res_df, all_groups, file = "fold_res_df.RData")


