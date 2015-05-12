#creation of aggregation groups
load("aa_nprop.RData")

traits <- list(size = c(30, 36, 54, 201),
               polarity = 202,
               pi = 203,
               hydroph = c(1, 26, 33, 57),
               alpha = c(12, 23, 155, 232))

grouping_properties <- t(aa_nprop[unlist(traits), ])
save(grouping_properties, file = "grouping_properties.RData")

#traits with polarity
all_traits_combn_polarity <- cbind(expand.grid(traits[["size"]], traits[["hydroph"]], traits[["alpha"]]),
                          traits[["polarity"]], traits[["pi"]])

colnames(all_traits_combn_polarity) <- c("size", "hydroph", "alpha", "polarity", "pi")

all_groups_polarity <- lapply(1L:nrow(all_traits_combn_polarity), function(single_trait_combn) {
  cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn_polarity[single_trait_combn, ]), ])))
  gr <- cutree(cl, k = 4)
  names(gr) <- tolower(names(gr))
  agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
  names(agg_gr) <- 1L:length(agg_gr)
  agg_gr
})


#traits without polarity
all_traits_combn_nonpol <- cbind(expand.grid(traits[["size"]], traits[["hydroph"]], 
                                             traits[["alpha"]]), traits[["pi"]])

colnames(all_traits_combn_nonpol) <- c("size", "hydroph", "alpha", "pi")


all_groups_nonpol <- lapply(1L:nrow(all_traits_combn_nonpol), function(single_trait_combn) {
  cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn_nonpol[single_trait_combn, ]), ])))
  gr <- cutree(cl, k = 4)
  names(gr) <- tolower(names(gr))
  agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
  names(agg_gr) <- 1L:length(agg_gr)
  agg_gr
})


aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("f", "w", "y", "s", "t", 
                                                                                                     "c", "n", "q"))

aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("s", "t", "c", "n", "q", "y", "w"))

library(signalHsmm)

all_groups <- c(list(aaaggregation), list(aa2), list(aa1), all_groups_polarity,
                all_groups_nonpol)
