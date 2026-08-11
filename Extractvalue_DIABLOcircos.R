# 以下代码必须接着regenerate_DIABLO_circos_ManualColor.R运行才不会报错
# ============================================================
# Extract DIABLO Circos correlation:
# Plasma T.b.MCA <-> Ileum T.b.MCA
# ============================================================
my_color <- c("#FFB194","#BEA6A0","#E69A5A","#845EC2","#C197FF","#00C9A7","#005B44")
pal <- colorRampPalette(my_color)
image(x = 1:7, y =1, z=as.matrix(1:7),col = pal(7))
my_color_MCD <- my_color[1:3]


feature_name <- "TbetaMCA"
block1 <- "Plasma"
block2 <- "Ileum"

# SAME components used in your circosPlot()
comp_use <- seq_len(min(2L, fit$ncomp))
# 1. Extract the exact value used by the Circos plot
# ------------------------------------------------------------
# 1. Check that T.b.MCA exists in both blocks
# ------------------------------------------------------------

if (!feature_name %in% colnames(fit$X[[block1]])) {
  stop(feature_name, " not found in ", block1)
}

if (!feature_name %in% colnames(fit$X[[block2]])) {
  stop(feature_name, " not found in ", block2)
}


# ------------------------------------------------------------
# 2. Find variables selected by DIABLO on components 1-2
# ------------------------------------------------------------

selected_features <- lapply(names(fit$X), function(block) {
  
  loading_mat <- fit$loadings[[block]][, comp_use, drop = FALSE]
  
  rownames(loading_mat)[
    apply(abs(loading_mat), 1, max, na.rm = TRUE) > 0
  ]
})

names(selected_features) <- names(fit$X)


cat("\nT.b.MCA selected in Plasma: ",
    feature_name %in% selected_features[[block1]], "\n")

cat("T.b.MCA selected in Ileum: ",
    feature_name %in% selected_features[[block2]], "\n")


# ------------------------------------------------------------
# 3. Calculate correlation of every selected feature
#    with its DIABLO latent components
#
#    This follows the calculation used by circosPlot()
# ------------------------------------------------------------

feature_component_cor <- lapply(names(fit$X), function(block) {
  
  feats <- selected_features[[block]]
  
  if (length(feats) == 0) {
    return(NULL)
  }
  
  tmp <- cor(
    fit$X[[block]][, feats, drop = FALSE],
    fit$variates[[block]][, comp_use, drop = FALSE],
    use = "pairwise.complete.obs"
  )
  
  # Give every variable a unique Block::Feature name
  rownames(tmp) <- paste(block, rownames(tmp), sep = "::")
  
  tmp
})

names(feature_component_cor) <- names(fit$X)


# ------------------------------------------------------------
# 4. Combine correlations and construct the SAME
#    similarity matrix used by circosPlot()
# ------------------------------------------------------------

C <- do.call(rbind, feature_component_cor)

circos_similarity_matrix <- C %*% t(C)


plasma_feature <- paste(block1, feature_name, sep = "::")
ileum_feature  <- paste(block2, feature_name, sep = "::")


# ------------------------------------------------------------
# 5. Extract Plasma T.b.MCA <-> Ileum T.b.MCA
# ------------------------------------------------------------

if (!plasma_feature %in% rownames(circos_similarity_matrix)) {
  
  cat("\n", plasma_feature,
      " was NOT selected by DIABLO on components ",
      paste(comp_use, collapse = ", "), ".\n",
      sep = "")
  
} else if (!ileum_feature %in% rownames(circos_similarity_matrix)) {
  
  cat("\n", ileum_feature,
      " was NOT selected by DIABLO on components ",
      paste(comp_use, collapse = ", "), ".\n",
      sep = "")
  
} else {
  
  circos_r <- circos_similarity_matrix[
    plasma_feature,
    ileum_feature
  ]
  
  cat("\n========================================\n")
  cat("DIABLO CIRCOS RESULT\n")
  cat("========================================\n")
  
  cat("Plasma feature :", feature_name, "\n")
  cat("Ileum feature  :", feature_name, "\n")
  cat("Components     :", paste(comp_use, collapse = ", "), "\n")
  cat("Circos similarity =", circos_r, "\n")
  cat("|similarity|      =", abs(circos_r), "\n")
  
  cat(
    "Pass cutoff 0.7   =",
    abs(circos_r) >= 0.7,
    "\n"
  )
  
  cat(
    "Direction         =",
    ifelse(circos_r > 0, "Positive",
           ifelse(circos_r < 0, "Negative", "Zero")),
    "\n"
  )
}
# 2. See exactly how that Circos value was constructed
# ============================================================
# Show contribution from each DIABLO component
# ============================================================

if (plasma_feature %in% rownames(C) &&
    ileum_feature %in% rownames(C)) {
  
  plasma_comp_cor <- C[plasma_feature, ]
  ileum_comp_cor  <- C[ileum_feature, ]
  
  component_result <- data.frame(
    Component = paste0("comp", comp_use),
    Plasma_TbMCA_cor = as.numeric(plasma_comp_cor),
    Ileum_TbMCA_cor  = as.numeric(ileum_comp_cor),
    Product = as.numeric(plasma_comp_cor * ileum_comp_cor)
  )
  
  print(component_result)
  
  cat(
    "\nCircos similarity = sum(Product) =",
    sum(component_result$Product),
    "\n"
  )
}
# 3. If you actually want the ordinary Plasma-vs-Ileum Pearson correlation
raw_r <- cor(
  fit$X[["Plasma"]][, "TbetaMCA"],
  fit$X[["Ileum"]][, "TbetaMCA"],
  use = "pairwise.complete.obs",
  method = "pearson"
)

cat(
  "\nRaw Pearson correlation:",
  raw_r,
  "\n"
)





# ============================================================
# Direct correlation:
# Plasma T.b.MCA vs Ileum T.b.MCA
# ============================================================

feature <- "TbetaMCA"

# Extract values
plasma_tb <- fit$X[["Plasma"]][, feature]
ileum_tb  <- fit$X[["Ileum"]][, feature]

# ------------------------------------------------------------
# Make sure samples are matched by sample ID
# ------------------------------------------------------------

common_samples <- intersect(
  names(plasma_tb),
  names(ileum_tb)
)

# If vectors do not carry names, use matrix rownames instead
if (length(common_samples) == 0) {
  
  common_samples <- intersect(
    rownames(fit$X[["Plasma"]]),
    rownames(fit$X[["Ileum"]])
  )
  
  dat_cor <- data.frame(
    Sample = common_samples,
    Plasma_TbMCA = fit$X[["Plasma"]][common_samples, feature],
    Ileum_TbMCA  = fit$X[["Ileum"]][common_samples, feature]
  )
  
} else {
  
  dat_cor <- data.frame(
    Sample = common_samples,
    Plasma_TbMCA = plasma_tb[common_samples],
    Ileum_TbMCA  = ileum_tb[common_samples]
  )
}

# Remove samples with missing values
dat_cor <- dat_cor[
  complete.cases(dat_cor[, c("Plasma_TbMCA", "Ileum_TbMCA")]),
]

cat("Number of matched samples =", nrow(dat_cor), "\n\n")


# ============================================================
# Pearson correlation
# ============================================================

pearson_res <- cor.test(
  dat_cor$Plasma_TbMCA,
  dat_cor$Ileum_TbMCA,
  method = "pearson"
)


# ============================================================
# Spearman correlation
# ============================================================

spearman_res <- cor.test(
  dat_cor$Plasma_TbMCA,
  dat_cor$Ileum_TbMCA,
  method = "spearman",
  exact = FALSE
)


# ============================================================
# Print results
# ============================================================

cat("=====================================\n")
cat("PLASMA vs ILEUM T.b.MCA\n")
cat("=====================================\n\n")

cat("Pearson correlation\n")
cat("-------------------\n")
cat("r =", round(unname(pearson_res$estimate), 4), "\n")
cat("P =", format.pval(pearson_res$p.value, digits = 4), "\n")
cat(
  "95% CI =",
  round(pearson_res$conf.int[1], 4),
  "to",
  round(pearson_res$conf.int[2], 4),
  "\n\n"
)

cat("Spearman correlation\n")
cat("--------------------\n")
cat("rho =", round(unname(spearman_res$estimate), 4), "\n")
cat("P =", format.pval(spearman_res$p.value, digits = 4), "\n")



# Add experimental group
Y_char <- as.character(Y)
names(Y_char) <- rownames(fit$X[["Plasma"]])

dat_cor$Group <- Y_char[dat_cor$Sample]

table(dat_cor$Group)


# ============================================================
# Correlation separately within each experimental group
# ============================================================

group_results <- do.call(
  rbind,
  lapply(unique(dat_cor$Group), function(g) {
    
    tmp <- dat_cor[dat_cor$Group == g, ]
    
    # Need enough observations
    if (nrow(tmp) < 3) {
      return(NULL)
    }
    
    p_res <- cor.test(
      tmp$Plasma_TbMCA,
      tmp$Ileum_TbMCA,
      method = "pearson"
    )
    
    s_res <- cor.test(
      tmp$Plasma_TbMCA,
      tmp$Ileum_TbMCA,
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Group = g,
      N = nrow(tmp),
      
      Pearson_r =
        unname(p_res$estimate),
      
      Pearson_P =
        p_res$p.value,
      
      Spearman_rho =
        unname(s_res$estimate),
      
      Spearman_P =
        s_res$p.value
    )
  })
)

group_results$Pearson_r <-
  round(group_results$Pearson_r, 4)

group_results$Pearson_P <-
  signif(group_results$Pearson_P, 4)

group_results$Spearman_rho <-
  round(group_results$Spearman_rho, 4)

group_results$Spearman_P <-
  signif(group_results$Spearman_P, 4)

print(group_results)



library(ggplot2)

# ============================================================
# Linear regression:
# Plasma T.b.MCA ~ Ileum T.b.MCA
# ============================================================

fit_lm <- lm(
  Plasma_TbMCA ~ Ileum_TbMCA,
  data = dat_cor
)

lm_summary <- summary(fit_lm)

# Regression coefficients
intercept <- coef(fit_lm)[1]
slope     <- coef(fit_lm)[2]

# R-squared
R2 <- lm_summary$r.squared

# P value for regression slope
p_value <- coef(lm_summary)[2, 4]


# ============================================================
# Pearson correlation
# ============================================================

pearson_res <- cor.test(
  dat_cor$Plasma_TbMCA,
  dat_cor$Ileum_TbMCA,
  method = "pearson"
)

pearson_r <- unname(pearson_res$estimate)


# ============================================================
# Print numerical results
# ============================================================

cat("=====================================\n")
cat("Plasma vs Ileum T.b.MCA\n")
cat("=====================================\n")

cat("Pearson r =", round(pearson_r, 4), "\n")
cat("R² =", round(R2, 4), "\n")
cat("Slope =", round(slope, 4), "\n")
cat("Intercept =", round(intercept, 4), "\n")
cat(
  "P =",
  format.pval(
    p_value,
    digits = 4,
    eps = 0.0001
  ),
  "\n"
)


# ============================================================
# Create regression equation text
# ============================================================

if (intercept >= 0) {
  
  equation_text <- paste0(
    "y = ",
    round(slope, 3),
    "x + ",
    round(intercept, 3)
  )
  
} else {
  
  equation_text <- paste0(
    "y = ",
    round(slope, 3),
    "x - ",
    round(abs(intercept), 3)
  )
}


# ============================================================
# Format P value
# ============================================================

if (p_value < 0.0001) {
  
  p_text <- "P < 0.0001"
  
} else if (p_value < 0.001) {
  
  p_text <- "P < 0.001"
  
} else {
  
  p_text <- paste0(
    "P = ",
    signif(p_value, 3)
  )
}


# ============================================================
# Combine all statistics
# ============================================================

stat_text <- paste0(
  equation_text,
  "\n",
  "Pearson r = ", round(pearson_r, 3),
  "\n",
  "R² = ", round(R2, 3),
  "\n",
  p_text
)


# ============================================================
# Scatter plot
# ============================================================

p <- ggplot(
  dat_cor,
  aes(
    x = Ileum_TbMCA,
    y = Plasma_TbMCA,
    color = Group
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  # Overall regression line
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    linewidth = 0.9
  ) +
  
  # Statistics
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = stat_text,
    hjust = 1.1,
    vjust = 1.2,
    size = 4.5,
    lineheight = 1.15
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  labs(
    x = "Ileum T.b.MCA",
    y = "Plasma T.b.MCA",
    color = "Group"
  ) +
  
  theme(
    plot.margin = margin(
      10, 15, 10, 10
    )
  ) + scale_color_manual(values = my_color_MCD)

print(p)


# ============================================================
# Save figure
# ============================================================

ggsave(
  "Plasma_vs_Ileum_TbMCA_correlation.pdf",
  p,
  width = 7,
  height = 5
)

ggsave(
  "Plasma_vs_Ileum_TbMCA_correlation.png",
  p,
  width = 6,
  height = 5,
  dpi = 300
)


