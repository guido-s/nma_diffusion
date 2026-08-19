#### Hat matrix and geometric series ####

# Rücker G, Davies AL, Schwarzer G (2026):
# Network meta-analysis and diffusion.
# *Research Synthesis Methods*, accepted for publication

library("netmeta")
library("MASS")

settings.netmeta(number.of.studies = FALSE)

source("functions.R")
source("data.R")

options(width = 100)

# Figure 1 ####

pdf("Figure1.pdf", width = 15)
par(mfrow = c(1, 3))
netgraph(nma2, cex = 2)
netgraph(nma.dong, seq = "optimal", rotate = 45, cex = 2)
labs <- c("Ante", "Hand", "Lido-prop", "Lido-pre", "Keta-pre",
  "NSAIDS-pre", "Opioid-pre")
netgraph(nma.jalota, start = "prcomp", seq = "optimal", rotate = 120, 
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
draw.TE(nma.dong, N = 10, text = TRUE)
dev.off()

## Figure 6 ####

pdf("Figure6.pdf", width = 9)
rbc <- walk(nma.dong, nma.dong$trts, 50)
dev.off()


# Example 3 ####

## Figure 7 ####
pdf("Figure7.pdf", width = 10)
diffusion(nma.jalota, N = 3000, equal = TRUE, proportions = TRUE)
text(2250, 0.02, "Antecubital vein (0.2 %)")
text(2250, 0.30, "Hand vein (49.9 %)")
text(2250, 0.575, "Lidocaine-propofol admixture (16.3 %)")
text(2250, 0.725, "Lidocaine pretreatment (11.8 %)")
text(2250, 0.835, "Pretreatment with ketamine (9.5 %)")
text(2250, 0.900, "Pretreatment with NSAIDS (4.2 %)")
text(2250, 0.950, "Pretreatment with opioids (8.0 %)")
dev.off()

## Figure 8 ####

pdf("Figure8.pdf", width = 13, height = 10)
par(mfrow = c(2, 2))
draw.TE(nma.jalota, N = 1000, diffusion.type = "simple", verbose = FALSE)
draw.TE(nma.jalota, N = 10, text = TRUE) # lazy
draw.TE(nma.jalota, N = 10, diffusion.type = "absorbing", ref = 3, text = TRUE)
draw.TE(nma.jalota, N = 10, diffusion.type = "absorbing", text = TRUE)
dev.off()
