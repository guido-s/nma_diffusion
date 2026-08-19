# All data sets ####

# Note that not all data sets given here are used in the article

# Triangle ####
studlab <- c("A:B", "A:C", "B:C")
TE <- rep(1, 3)
seTE <- rep(1, 3)
treat1 <- c("A", "A", "B")
treat2 <- c("B", "C", "C")
nma0 <- netmeta(TE, seTE, treat1, treat2, studlab)


# Bicycle ####
study <- c("A:B", "A:D", "B:C", "B:D", "C:D")
effect <- rep(1, 5)
se.effect <- rep(1, 5)
trt1 <- c("A", "A", "B", "B", "C")
trt2 <- c("B", "D", "C", "D", "D")
nma1 <- netmeta(effect, se.effect, trt1, trt2, study)


# Bicycle extended ####
trt1 <- c("A", "A", "B", "B", "B", "C", "D")
trt2 <- c("B", "E", "C", "D", "E", "D", "E")
effect <- rep(1, 7)
se.effect <- rep(1, 7)
study <- paste0(trt1, trt2)
nma2 <- netmeta(effect, se.effect, trt1, trt2, study)


# Bicycle extended with tail ####
trt1 <- c("A", "A", "B", "B", "B", "C", "D", "E")
trt2 <- c("B", "E", "C", "D", "E", "D", "E", "F")
effect <- rep(1, 8)
se.effect <- rep(1, 8)
study <- paste0(trt1, trt2)
nma3 <- netmeta(effect, se.effect, trt1, trt2, study)


# Example (Dependent paths) ####
study <- 1:9
t1 <- LETTERS[c(1, 1, 2, 3, 3, 4, 5, 5, 6)]
t2 <- LETTERS[c(2, 3, 3, 4, 5, 5, 6, 7, 7)]
effect <- rep(1, 9)
se.effect <- rep(1, 9)
nma4 <- netmeta(effect, se.effect, t1, t2, study, random = FALSE)


# Chain (bipartite) ####
trt1 <- c("A", "B", "C", "D", "E", "F")
trt2 <- c("B", "C", "D", "E", "F", "G")
effect <- rep(0, 6)
se.effect <- rep(1, 6)
study <- paste0(trt1, trt2)
chain <- netmeta(effect, se.effect, trt1, trt2, study)


# Weighted circle (odd) ####
trt1 <- c("A", "B", "C", "D", "E", "F", "G")
trt2 <- c("B", "C", "D", "E", "F", "G", "A")
effect <- rnorm(7)
se.effect <- 1:7
study <- paste0(trt1, trt2)
nma7 <- netmeta(effect, se.effect, trt1, trt2, study)


# Weighted circle (even) ####
trt1 <- c("A", "B", "C", "D", "E", "F", "G", "H")
trt2 <- c("B", "C", "D", "E", "F", "G", "H", "A")
effect <- rnorm(8)
se.effect <- 1:8
study <- paste0(trt1, trt2)
nma8 <- netmeta(effect, se.effect, trt1, trt2, study)


# Star (bipartite) ####
n <- 7
trt1 <- rep(1, n - 1)
trt2 <- 2:n
effect <- rep(0, n - 1)
se.effect <- rep(1, n - 1)
study <- 1:(n - 1)
star <- netmeta(effect, se.effect, trt1, trt2, study)


# Senn2013 ####
pw_senn <- pairwise(treatment, n = n, mean = mean, sd = sd,
  data = Senn2013, studlab = study, sm = "MD")
nma_senn <- netmeta(pw_senn, nchar.trts = 4)


# Dong 2013 ####
pw_dong <- pairwise(treatment, death, randomized, studlab = id, 
  data = Dong2013, sm = "OR", allstudies = TRUE)
nma_dong <- netmeta(pw_dong)


# Stowe 2010 ####
pw_stowe <- pairwise(treat = list(t1, t2, t3), n = list(n1, n2, n3), 
  mean = list(y1, y2, y3), sd = list(sd1, sd2, sd3), 
  studlab = study, data = Stowe2010, sm = "MD")
nma_stowe <- netmeta(pw_stowe)


# Jalota 2011 ####
pw_jalota <- pairwise(trt, pain, n, studlab = id, sm = "RR",
  data = Jalota2011)
nma_jalota <- netmeta(pw_jalota, reference.group = "Hand vein")
