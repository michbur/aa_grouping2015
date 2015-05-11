#creation of aggregation groups
load("aa_nprop.RData")

traits <- list(size = c(30, 36, 54, 201),
               polarity = 202,
               pi = 203,
               hydroph = c(1, 26, 33, 57),
               alpha = c(12, 23, 155, 232))

all_traits_combn <- cbind(expand.grid(traits[["size"]], traits[["hydroph"]], traits[["alpha"]]),
                          traits[["polarity"]], traits[["pi"]])

colnames(all_traits_combn) <- c("size", "hydroph", "alpha", "polarity", "pi")

grouping_properties <- t(aa_nprop[unlist(traits), ])
save(grouping_properties, file = "grouping_properties.RData")

all_groups <- lapply(1L:nrow(all_traits_combn), function(single_trait_combn) {
  cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn[single_trait_combn, ]), ])))
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

all_groups <- c(all_groups, list(aa2), list(aa1), list(aaaggregation))
