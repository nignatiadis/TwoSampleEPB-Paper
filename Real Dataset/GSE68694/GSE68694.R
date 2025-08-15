library(dplyr)
library(ggplot2)
library(GEOquery)
library(REBayes)
library(AnnotationDbi)
library(hgu133plus2.db)
library(matrixStats)
library(Biobase)

source('../../EPB main/EPB.R')
source('../../helper functions/Rejection Region helper.R')
source('../../helper functions/Distribution Visualization helper.R')

# Data Preprocessing

#gset = getGEO("GSE68694", GSEMatrix = TRUE)
#eset = gset[[1]]
#group_assignment = pData(eset)$title
#X = exprs(eset)   # probe expression matrix (probes x samples)

## 1) drop Affy control probes
#keep_ctrl = !grepl("^AFFX", rownames(X))
#X = X[keep_ctrl, , drop=FALSE]

## 2) Remove low-expression probes
## keep probes expressed above a global cutoff in at least half the samples
#cutoff = quantile(X, 0.05)                         # 5th percentile of GSE68694 intensities
#min_samples = ceiling(0.5 * ncol(X))               # expressed in >= half of samples
#keep_expr = rowSums(X > cutoff) >= min_samples
#expr_filtered = X[keep_expr, , drop=FALSE]         # probe-level filtered matrix

## 3) Merge probes into genes (GPL570) ---------------------------------------
## Map probes -> gene symbols
#ann = AnnotationDbi::select(hgu133plus2.db, keys     = rownames(expr_filtered), columns  = c("SYMBOL"), keytype  = "PROBEID")

## Keep one row per probe ID, aligned to expr_filtered
#ann = ann[!duplicated(ann$PROBEID), ]
#expr_filtered = expr_filtered[ann$PROBEID, , drop=FALSE]
#ann = ann[match(rownames(expr_filtered), ann$PROBEID), ]

## Keep only probes with a valid SYMBOL
#has_symbol = !is.na(ann$SYMBOL) & ann$SYMBOL != ""
#E = expr_filtered[has_symbol, , drop=FALSE]
sym = ann$SYMBOL[has_symbol]

## average probes per gene
#gene_sum    = rowsum(E, group = sym)
#gene_counts = as.integer(table(sym)[rownames(gene_sum)])
#gene_expr_pre   = sweep(gene_sum, 1, gene_counts, "/")   # rows = genes, cols = samples

#X1_pre = gene_expr_pre[, 1:3] # GSE68694 FATPAD
#X2_pre = gene_expr_pre[, 4:6] # GSE68694 MIND
#info_pre = information_extractor(X1_pre, X2_pre)
#filter_out_idx = c(which(info_pre$S1_list == 0), which(info_pre$S2_list == 0)) # = c(8988)
#expr = gene_expr_pre[-filter_out_idx, ]
#write.csv(expr, "Filtered GSE68694/GSE68694_gene_expression_matrix.csv", row.names = TRUE)

# DE analysis

expr = read.csv("Filtered GSE68694/GSE68694_gene_expression_matrix.csv")
m = nrow(expr) # 20915



K1 = 3
K2 = 3
X1 = expr[, 2:4] # filtered FATPAD
X2 = expr[, 5:7] # filtered MIND

# Run VREPB, DVEPB, Welch, B_F and VE_test

info = information_extractor(X1, X2)

alpha = 0.1
NPMLE_1D_parameter = c(1000, 0, 1.0)
NPMLE_2D_parameter = c(80, 80, 0.01, 1.0)
algorithm_list = c(1,2,3,4,5)

#result = solver(X1, X2, NPMLE_1D_parameter, NPMLE_2D_parameter, algorithm_list)

# Number of discoveries for each method

#length(my_BH(result$VREPB, alpha)) # 6805
#length(my_BH(result$DVEPB, alpha)) # 8303
#length(my_BH(result$Welch, alpha)) # 4970
#length(my_BH(result$B_F, alpha)) # 0
#length(my_BH(result$EV_test, alpha)) # 6869

#df_2D_GSE68694 = data.frame(x = result$'2D_grid'[,1], y = result$'2D_grid'[,2], prob = result$'2D_mass')
#write.csv(df_2D_GSE68694, './Result/2D_GSE68694.csv')

# Generate the H plot for GSE68694

df_2D_GSE68694 = read.csv('./Result/2D_GSE68694.csv')

Total_mass = 0.999999
plot_H_GSE68694 = plotter_2D_ALL(Total_mass, df_2D_GSE68694, '')
ggsave(filename = './Result/H(var1, var2) of GSE68694.jpg', plot = plot_H_GSE68694, width = 12, height = 10, dpi = 500)

# Generate the G plot for FSE68694

df_1D_GSE68694 = data.frame(x = df_2D_GSE68694$x/df_2D_GSE68694$y, prob = df_2D_GSE68694$prob)

plot_G_GSE68694 = plotter_1D(df_1D_GSE68694, '')
ggsave(filename = './Result/G(lambda) of GSE68694.jpg', plot = plot_G_GSE68694, width = 8, height = 6, dpi = 500)

# Generate the Rejection Region plot

grid = df_2D_GSE68694[, 2:3]
mass = df_2D_GSE68694[, 4]
grid_1D = grid[, 1] / grid[, 2]

tau_list = info$S1_list/info$S2_list
log_tau_list = log(tau_list)
log_tau_list = log_tau_list[log_tau_list >= -7.5 & log_tau_list <= 7.5]
tau_grid = exp(seq(min(log_tau_list), max(log_tau_list), length.out = 1000))

t_BF_list = c()
for (i in 1:m) {
  X1_i = as.numeric(X1[i,])
  X2_i = as.numeric(X2[i,])
  m1_i = mean(X1_i)
  m2_i = mean(X2_i)
  s1_i = var(X1_i)
  s2_i = var(X2_i)
  t_BF_i = (m1_i - m2_i) / sqrt(s1_i/K1 + s2_i/K2)
  t_BF_list = c(t_BF_list, t_BF_i)
}
t_BF_list = t_BF_list[log(tau_list) >= -7.5 & log(tau_list) <= 7.5] # Range of sample t_BF to determine range of root solver

## Rejection Region for VREPB

t_VREPB = c()

for (tau in tau_grid) {
  t_VREPB_i = BF_pvalue_VREPB_solver_positive(K1, K2, grid_1D, mass, tau, alpha)
  t_VREPB = append(t_VREPB, t_VREPB_i)
}

## Rejection Region for B_F
t_BF = c()

for (tau in tau_grid) {
  t_BF_i = BF_pvalue_BF_solver_positive(tau, K1, K2, alpha)
  t_BF = append(t_BF, t_BF_i)
}

## Rejection Region for Welch

t_Welch = c()

for (tau in tau_grid) {
  t_Welch_i = BF_pvalue_Welch_solver_positive(tau, K1, K2, alpha)
  t_Welch = append(t_Welch, t_Welch_i)
}

## Rejection Region for VE-test

t_EV = c()

for (tau in tau_grid) {
  t_EV_i = BF_pvalue_EV_solver_positive(tau, K1, K2, alpha)
  t_EV = append(t_EV, t_EV_i)
}

## Plot rejection region for VREPB, B_F, Welch, EV_test

size = 1

curve_colors = c(
  "VREPB"   = "#0072B2",   # blue
  "B-F"     = "#CC79A7",   # magenta
  "Welch"   = "#E69F00",   # gold
  "EV-test"  = "#D55E00"     # orange
)

plot_2D = function(u1, u2) {
  ggplot(data.frame(u1 = u1, u2 = u2), aes(x = u1, y = u2)) +
    stat_bin2d(bins = 80) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(x = expression(log(hat(lambda)[i])), y = expression(T[i]^{BF}), fill = "Count") +
    #scale_y_continuous(limits = c(-4, 4)) +
    scale_y_continuous(limits = c(-12, 12)) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major = element_line(colour = "gray80", size=0.5),
      panel.grid.minor = element_line(colour = "gray90", size=0.25),
      legend.position = "top",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      plot.margin = margin(12, 16, 12, 12)
    )
}

p = plot_2D(log_tau_list, t_BF_list) +
  geom_line(data = data.frame(x = log(tau_grid), y = t_VREPB), aes(x = x, y = y, color = "VREPB"), size = size) +
  geom_line(data = data.frame(x = log(tau_grid), y = -t_VREPB), aes(x = x, y = y, color = "VREPB"), size = size) +
  
  geom_line(data = data.frame(x = log(tau_grid), y = t_BF), aes(x = x, y = y, color = "B-F"), size = size) +
  geom_line(data = data.frame(x = log(tau_grid), y = -t_BF), aes(x = x, y = y, color = "B-F"), size = size) +
  
  geom_line(data = data.frame(x = log(tau_grid), y = t_EV), aes(x = x, y = y, color = "EV-test"), size = size) +
  geom_line(data = data.frame(x = log(tau_grid), y = -t_EV), aes(x = x, y = y, color = "EV-test"), size = size) +
  
  geom_line(data = data.frame(x = log(tau_grid), y = t_Welch), aes(x = x, y = y, color = "Welch"), size = size) +
  geom_line(data = data.frame(x = log(tau_grid), y = -t_Welch), aes(x = x, y = y, color = "Welch"), size = size) +
  
  scale_color_manual(
    name = NULL, 
    values = curve_colors,
    breaks = c("VREPB", "B-F", "Welch", "EV-test"),
    labels = c("VREPB", "B-F", "Welch", "EV-test")
  ) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.7), 
         color = guide_legend(override.aes = list(size = 2)))

# Save as JPG

ggsave("./Result/Rejection_Region_GSE68694.jpg", p, width = 7, height = 5.5, dpi = 400)