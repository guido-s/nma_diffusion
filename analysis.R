#### Hat matrix and geometric series ####

# Rücker G, Davies AL, Schwarzer G (2026):
# Network meta-analysis and diffusion.
# *Research Synthesis Methods*, accepted for publication

library("netmeta")

settings.netmeta(number.of.studies = FALSE)

source("functions.R")
source("data.R")

options(width = 100)

# Figure 1 ####

pdf("Figure1.pdf", width = 15.5)
par(mfrow = c(1, 3))
netgraph(nma2, cex = 2)
netgraph(nma_dong, seq = "optimal", rotate = 1 / n * 360, cex = 2)
labs <- c("Ante", "Hand", "Lido-prop", "Lido-pre", "Keta-pre",
  "NSAIDS-pre", "Opioid-pre")
netgraph(nma_jalota, start = "prcomp", seq = "optimal", rotate = 103,
  labels = labs, scale = 1.2, cex = 2)
dev.off()


# Example 1 ####

## Figure 2 ####

pdf("Figure2.pdf", width = 9, height = 6)
par(mfrow = c(2, 3))
diffusion(nma2, N = 5, adj = 0.5, col = "darkgray")
dev.off()

## Figure 4 ####

pdf("Figure4.pdf")
walk(nma2, c("A", "B", "C", "D", "E"))
dev.off()


# Example 2 ####

## Figure 5 ####

pdf("Figure5.pdf")
draw.TE(nma_dong, N = 10, text = TRUE)
dev.off()

## Figure 6 ####

pdf("Figure6.pdf", width = 9)
rbc <- walk(nma_dong, nma_dong$trts, 50)
dev.off()


# Example 3 ####

## Figure 7 ####
pdf("Figure7.pdf", width = 10)
dif <- diffusion(nma_jalota, N = 3000, equal = TRUE, proportions = TRUE)
dif
#
xpos <- 2950
ypos <- c(0.02, 0.30, 0.575, 0.725, 0.835, 0.9, 0.95)
#
props <- paste0(rownames(dif), " (", sprintf("%.1f", dif * 100), " %)")
for (i in seq_along(props))
  text(xpos, ypos[i], props[i], adj = 1)
dev.off()

## Figure 8 ####

pdf("Figure8.pdf", width = 13, height = 10)
par(mfrow = c(2, 2))
draw.TE(nma_jalota, N = 1000, diffusion.type = "simple")
draw.TE(nma_jalota, N = 10, text = TRUE) # lazy
draw.TE(nma_jalota, N = 10, diffusion.type = "absorbing",
  ref = "Lidocaine-propofol admixture", text = TRUE)
draw.TE(nma_jalota, N = 10, diffusion.type = "absorbing", text = TRUE)
dev.off()

## Figure 9 (Appendix F) ####

# Like Figure 8, but based on the random effects model 

pdf("Figure9.pdf", width = 13, height = 10)
par(mfrow = c(2, 2))
draw.TE(nma_jalota, N = 1000, diffusion.type = "simple",
  random = TRUE)
draw.TE(nma_jalota, N = 10, text = TRUE, random = TRUE)
draw.TE(nma_jalota, N = 10, diffusion.type = "absorbing",
  ref = "Lidocaine-propofol admixture", text = TRUE, random = TRUE)
draw.TE(nma_jalota, N = 10, diffusion.type = "absorbing", text = TRUE,
  random = TRUE)
dev.off()
