library(ggplot2)
library(reshape2)
library(dplyr)
df_8 <- df_7
ct_names <- c("Nodule_CT","GGO_CT","Consolidation_CT","Halo_CT","Pleuraleffusion_CT","Cavity_CT")
df_8_ct <- df_8[,ct_names]
df_8_num <- as.data.frame(lapply(df_8_ct, function(x) {
    case_when(
        x == "Yes" ~ 1,
        x == "No" ~ 0,
        TRUE     ~ NA_real_
    )
}))
# -----------------------------
# 1. Load libraries
# -----------------------------
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

# -----------------------------
# 2. Define your CT columns
# -----------------------------
ct_cols <- c(
  "GGO_CT", 
  "Consolidation_CT", 
  "Pleuraleffusion_CT",
  "Nodule_CT", 
  "Halo_CT"
)


# -----------------------------
# 4. Create NA-safe co-occurrence matrix
# -----------------------------
n <- length(ct_cols)
co_counts <- matrix(0, nrow = n, ncol = n)
colnames(co_counts) <- rownames(co_counts) <- ct_cols

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    x <- df_8_num[[ct_cols[i]]]
    y <- df_8_num[[ct_cols[j]]]

    valid <- complete.cases(x, y) # ignore NA rows
    co_counts[i, j] <- sum(x[valid] == 1 & y[valid] == 1)
  }
}

# -----------------------------
# 5. Pairwise correlation matrix (also NA-safe)
# -----------------------------
co_corr <- cor(df_8_num, use = "pairwise.complete.obs")

# -----------------------------
# 6. Melt both matrices
# -----------------------------
counts_long <- melt(co_counts) %>%
  rename(CT1 = Var1, CT2 = Var2, Count = value)

corr_long <- melt(co_corr) %>%
  rename(CT1 = Var1, CT2 = Var2, Corr = value)

# Merge
heatmap_df <- left_join(counts_long, corr_long, by = c("CT1", "CT2"))

# Keep only lower triangle
heatmap_df <- heatmap_df %>%
  mutate(
    i = as.numeric(factor(CT1)),
    j = as.numeric(factor(CT2))
  ) %>%
  filter(i >= j)

# -----------------------------
# 7. Final triangular heatmap
# -----------------------------
ggplot(heatmap_df, aes(CT1, CT2, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%d\n(%.2f)", Count, Corr)),
            size = 4) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(
    title = "Co-occurrence of CT Findings",
    x = "", 
    y = "",
    fill = "Co-occurrence\nCount"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
