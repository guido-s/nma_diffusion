# (1) Matrix powers and power series ####

# Matrix powers and their partial sums

# M is a square matrix
# k is a (large) number of iterations

Power <- function(M, k) {
  S <- P <- diag(dim(M)[1])
  if (k > 0) {
    for (i in seq_len(k)) {
      P.new <- M %*% P # P.new = next power of M
      S <- S + P.new   # next partial sum of power series
      P <- P.new
    }
  }
  list(P = P, S = S)
}


# (2) Visualize the exchange between the nodes ####

# x           netmeta object
# t0          starting vector of length n (number of nodes of x)
#             default is starting in node 1
# equal       if TRUE, starts with equal proportions at each node (1 / n)
# N           number of iterations
# size        maximum point size for network graphs
# cex         size of treatment labels
# only.result if TRUE, show only the result of the last iteration
# random      if TRUE, a random effects model is used
# proportions if TRUE, show proportional degrees over time
# diffusion.type "simple" (default), "lazy", or "absorbing"; see paper.
# p           if diffusion.type is "lazy", the proportion of lazy mass (default is 0.5)
# ref         if not NULL, specifies the reference node (numeric)
# title       the step number, if proportions = FALSE
# flow        if TRUE, distribution of mass over edges is shown by their thickness
# verbose     if TRUE, distribution of mass over nodes is printed for all iterations 

diffusion <- function(x, t0 = NULL, equal = FALSE, N = 20, size = 200, cex = 1,
                      random = FALSE, only.result = FALSE, proportions = FALSE,
                      diffusion.type = "simple", p = 0.5, ref = NULL,
                      title = TRUE, flow = FALSE, verbose = FALSE,
                      ...) {
  if (is.null(t0)) {
    t0 <- c(1, rep(0, x$n - 1))
    if (equal)
      t0 <- rep(1 / x$n, x$n)
  }
  #
  if (!title)
    main <- ""
  #
  if (random)
    L <- x$L.matrix.random
  else
    L <- x$L.matrix.common # Laplacian matrix
  #
  D <- diag(diag(L))       # Degree matrix
  A <- D - L               # Adjacency matrix
  I <- diag(x$n)
  Tr <- A %*% diag(1 / diag(D)) # Transition/Diffusion matrix (if diffusion.type = "simple")
  if (diffusion.type == "lazy")
    Tr <- p * Tr + (1 - p) * I
  if (diffusion.type == "absorbing" & is.null(ref))
    Tr[, x$n] <- c(rep(0, x$n - 1), 1)
  if (diffusion.type == "absorbing" & !is.null(ref))
    Tr[, ref] <- c(rep(0, ref - 1), 1, rep(0, x$n - ref))
  #
  if (!proportions) {
    rescale <- TRUE
    for (i in 0:N) {
      t <- Power(Tr, i)$P %*% t0
      if (verbose)
        print(t)
      if (flow) {
        thickness <- Tr %*% diag(as.vector(t)) + diag(as.vector(t)) %*% t(Tr)
        colnames(thickness) <- rownames(thickness) <- colnames(A)
        mx <- max(thickness)
      }
      else {
        thickness <- "number.of.studies"
        mx <- 0.05
      }
      if (!only.result) {
        if (title)
          main <- paste("Step", i)
        #
        netgraph(x, seq = "opt", cex = cex, cex.points = sqrt(size * t), 
                 points.max = size, rescale.pointsize = rescale, 
                 thickness = thickness, rescale.thickness = FALSE,
                 main = main, lwd.min = 20 * mx, lwd.max = 30 * mx,
                 ...)
      }
    }
    #
    if (only.result) {
      if (title) main <- paste(N, "Steps")
      netgraph(x, seq = "opt", cex = cex, cex.points = sqrt(size * t), 
               points.max = size, rescale.pointsize = rescale, 
               thickness = thickness, rescale.thickness = FALSE,
               main = main, lwd.min = 20 * mx, lwd.max = 30 * mx,
               ...)
    }
  }
  else {
    plot(0, 0, xlim = c(0, N), ylim = c(0, 1), cex = 0, las = 1,
         xlab = "", ylab = "")
    tsum <- t <- t0
    for (j in 2:x$n)
      tsum[j] <- tsum[j - 1] + t[j]
    polygon(c(0, 1, 1, 0), c(0, 0, tsum[1], tsum[1]), col = 1, border = NA)
    for (k in 2:x$n) {
      polygon(c(0, 1, 1, 0), c(tsum[k - 1], tsum[k - 1], tsum[k], tsum[k]), 
              col = k, border = NA)
    }
    for (i in seq_len(N)) {
      t <- Tr %*% t
      print(t)
      tsum <- t
      for (j in 2:x$n)
        tsum[j] <- tsum[j - 1] + t[j]
      polygon(c(i, i + 1, i + 1, i), c(0, 0, tsum[1], tsum[1]),
              col = 1, border = NA)
      for (k in 2:x$n) {
        polygon(c(i, i + 1, i + 1, i),
                c(tsum[k - 1], tsum[k - 1], tsum[k], tsum[k]),
                col = k, border = NA)
      }
    }
  }
}


# (3) Calculate the hat matrix using a matrix power series ####

# x is a netmeta object
# N is a (large) number of iterations
# ref is the number of a treatment that is to be used as a reference
# random if true, a random effects model is used

Hat <- function(x, N = 1000, ref = NULL, random = FALSE) {
  B <- x$B.matrix # Edge-vertex incidence matrix (design matrix)
  if (x$m == 1)
    B <- t(B)
  I <- diag(x$n)                # Identity matrix
  if (random)
    L <- x$L.matrix.random
  else
    L <- x$L.matrix.common      # Laplacian 
  D <- diag(diag(L))            # Diagonal matrix of Laplacian (degree matrix)
  A <- D - L                    # Generalized (weighted) adjacency matrix
  Tr <- A %*% diag(1 / diag(D)) # diffusion/transition matrix
  Tr1 <- (Tr + I) / 2           # Lazy walk diffusion matrix
  if (random) {
    W <- as.matrix(x$W.matrix.random)                 # Diagonal weight matrix
    H <- as.matrix(hatmatrix(x, method = "R")$random) # Hat matrix 
  }
  else {
    W <- as.matrix(x$W.matrix.common)                 # Diagonal weight matrix
    H <- as.matrix(hatmatrix(x, method = "R")$common) # Hat matrix
  }
  #
  H1 <- B %*% diag(1 / diag(D)) %*% Power(Tr, N)$S %*% t(B) %*% W    # Hat matrix, special
  H2 <- B %*% diag(1 / diag(D)) %*% Power(Tr1, N)$S %*% t(B) %*% W / 2 # Hat matrix, general
  H3 <- check3 <- Tr2 <- NA # Hat matrix, with reference node
  #
  if (!is.null(ref)) {
    Br <- as.matrix(B[, -ref])
    Tr2 <- as.matrix(Tr[-ref, -ref])
    Dr <-  as.matrix(D[-ref, -ref])
    H3 <- Br %*% solve(Dr) %*% Power(Tr2, N)$S %*% t(Br) %*% W 
    check3 <- all.equal(H, H3)
  }
  check1 <- all.equal(H, H1)
  check2 <- all.equal(H, H2)
  res <- list(x = x, Tr = Tr, Tr1 = Tr1, Tr2 = Tr2,
              H = H, H1 = H1, H2 = H2, H3 = H3,
              check = list(check1, check2, check3))
  res
}


# (4) Visualize the iterative course from TE to TE.nma.common ####

# x              netmeta object
# N              number of iterations
# diffusion.type "simple", "lazy" (default), or "absorbing"; see paper
# ref            if not NULL, specifies the reference node (numeric)
# text           if TRUE, show Q values in the graph
# random         if TRUE, a random effects model is used
# verbose        if TRUE, give longer output

draw.TE <- function(x, N = 10, diffusion.type = "lazy", ref = NULL,
                    random = FALSE, text = FALSE, verbose = FALSE) {
  col1 <- seq_len(x$m)
  y.new <- matrix(0, nrow = x$m, ncol = N + 2)
  y.new[, 1] <- x$TE
  Q.n <- vector("numeric", N + 2)
  Q.n[1] <- x$Q
  #
  if (random) {
    ylim <- c(min(x$TE.nma.random) - sd(x$TE.nma.random), 
              max(x$TE.nma.random) + sd(x$TE.nma.random))
    d <- diag(x$L.matrix.random)
  }
  else {
    ylim <- c(min(x$TE.nma.common) - sd(x$TE.nma.common), 
              max(x$TE.nma.common) + sd(x$TE.nma.common))
    d <- diag(x$L.matrix.common)
  }
  #
  plot(rep(-1, x$m), x$TE, xlim = c(-1, N), ylim = ylim,
       xlab = "Iteration step", ylab = "Treatment effect", 
       las = 1, pch = 16, cex = 0.5, col = col1)
  r <- which(d == max(d))[1]
  #
  for (n in 0:N) {
    if (diffusion.type == "simple")
      y.new[, n + 2] <- Hat(x, n, random = random)$H1 %*% x$TE 
    else if (diffusion.type == "lazy")
      y.new[, n + 2] <- Hat(x, n, random = random)$H2 %*% x$TE
    else if (diffusion.type == "absorbing" & is.null(ref))
      y.new[, n + 2] <- Hat(x, n, r, random)$H3 %*% x$TE
    else if (diffusion.type == "absorbing" & !is.null(ref))
      y.new[, n + 2] <- Hat(x, n, ref, random)$H3 %*% x$TE
    if (verbose)
      print(y.new[, n + 2])
    #
    if (random) {
      Q.n[n + 2] <-
        t(y.new[, n + 2] - x$TE.nma.random) %*%
        x$W.matrix.random %*% (y.new[, n + 2] - x$TE.nma.random)
    }
    else {
      Q.n[n + 2] <-
        t(y.new[, n + 2] - x$TE.nma.common) %*%
        x$W.matrix.common %*% (y.new[, n + 2] - x$TE.nma.common)
    }
    #
    points(rep(n, x$m), y.new[, n + 2], pch = 16, cex = 0.5, col = col1)
    #
    for (edge in seq_len(x$m)) {
      lines(c(n - 1, n), c(y.new[edge, n + 1], y.new[edge, n + 2]), 
            col = col1[edge], lwd = 2)
    }
    #
    if (text)
      text(n, ylim[1], round(Q.n[n + 2], 2), cex = 0.8)
  }
  if (text)
    text(-1, ylim[1], round(x$Q, 2), cex = 0.8)
  #
  res <- list(y = y.new, Q = Q.n, ref = ref)
  #
  invisible(res)
}


# (5) Visualize the iterative course of hat matrix diagonal or variances ####

# x              netmeta object
# N              number of iterations
# hat            if TRUE, show the course of the hat matrix diagonal,
#                if FALSE, show the course of the NMA variance estimates
# diffusion.type "simple", "lazy" (default), or "absorbing"; see paper
# ref            if not NULL, specifies the reference node (numeric)
# random         if TRUE, a random effects model is used
# verbose        if TRUE, give longer output

draw.hat <- function(x, N = 10, hat = TRUE, diffusion.type = "lazy",
                     ref = NULL, random = FALSE, verbose = FALSE) {
  col1 <- seq_len(x$m)
  y.new <- matrix(0, nrow = x$m, ncol = N + 2)
  if (random)
    d <- diag(x$L.matrix.random)
  else
    d <- diag(x$L.matrix.common)
  r <- (which(d == max(d)))[1]
  if (hat) {
    y.new[, 1] <- rep(1, x$m)
    ylim = c(0, 1)
    plot(rep(-1, x$m), rep(1, x$m), xlim = c(-1, N), ylim = ylim,
         xlab = "Iteration step", ylab = "Hat diagonal", 
         las = 1, pch = 16, cex = 0.5, col = col1)
    for (n in 0:N) {
      if (diffusion.type == "simple")
        y.new[, n + 2] <- diag(Hat(x, n, random = random)$H1) 
      else if (diffusion.type == "lazy")
        y.new[, n + 2] <- diag(Hat(x, n, random = random)$H2)
      else if (diffusion.type == "absorbing" & is.null(ref))
        y.new[, n + 2] <- diag(Hat(x, n, r, random)$H3)
      else if (diffusion.type == "absorbing" & !is.null(ref))
        y.new[, n + 2] <- diag(Hat(x, n, ref, random)$H3)
      if (verbose)
        print(y.new[, n + 2])
      points(rep(n, x$m), y.new[, n + 2], pch = 16, cex = 0.5, col = col1)
      for (edge in seq_len(x$m)) {
        lines(c(n - 1, n), c(y.new[edge, n + 1], y.new[edge, n + 2]), 
              col = col1[edge], lwd = 2)
      }
    }
  }
  else {
    if (random)
      w <- diag(x$W.matrix.random)
    else
      w <- diag(x$W.matrix.common)
    y.new[, 1] <- x$seTE^2
    ylim = c(0, max(x$seTE^2))
    plot(rep(-1, x$m), y.new[, 1], xlim = c(-1, N), ylim = ylim,
         xlab = "Iteration step", ylab = "Variance", 
         las = 1, pch = 16, cex = 0.5, col = col1)
    for (n in 0:N) {
      if (diffusion.type == "simple")
        y.new[, n + 2] <- diag(Hat(x, n, random = random)$H1) / w 
      else if (diffusion.type == "lazy")
        y.new[, n + 2] <- diag(Hat(x, n, random = random)$H2) / w
      else if (diffusion.type == "absorbing" & is.null(ref))
        y.new[, n + 2] <- diag(Hat(x, n, r, random)$H3) / w
      else if (diffusion.type == "absorbing" & !is.null(ref))
        y.new[, n + 2] <- diag(Hat(x, n, ref, random)$H3) / w
      if (verbose)
        print(y.new[, n + 2])
      points(rep(n, x$m), y.new[, n + 2], pch = 16, cex = 0.5, col = col1)
      for (edge in seq_len(x$m)) {
        lines(c(n - 1, n), c(y.new[edge, n + 1], y.new[edge, n + 2]), 
              col = col1[edge], lwd = 2)
      }
    }
  }
}


# (6) Function for weighted networks, eight colors ####

# x        netmeta object
# nodes    Up to 8 nodes that are to be presented
# N        Number of iterations
# random   if TRUE, a random effects model is used

walk <- function(x, nodes, N = 50, random = FALSE) {
  cols <- c("#FFFF0002", "#FFAA0002", "#FF000002", "#00FF0002",
            "#FF55FF02", "#0000FF02", "#00FFFF02", "#00000002")
  cols.full <- c("#FFFF00", "#FFAA00", "#FF0000", "#00FF00",
                 "#FF55FF", "#0000FF", "#00FFFF", "#000000")
  #
  nodes <- nodes[order(nodes)]
  #
  B <- x$B.matrix
  #
  if (random)
    L <- x$L.matrix.random
  else
    L <- x$L.matrix.common
  #
  d <- diag(L)
  D <- diag(d)        
  A <- D - L
  I <- diag(x$n)
  Tr <- (A %*% diag(1 / diag(D)) + I) / 2
  #
  W <- diag(-netmeta:::uppertri(L)) 
  if (dim(L)[1] == 2)
    W <- -as.matrix(netmeta:::uppertri(L))
  colnames(W) <- rownames(W) <- names(x$prop.direct.common)
  #
  nn <- length(nodes)
  cols <- cols[seq_len(nn)]
  cols.full <- cols.full[seq_len(nn)]
  take <- volume <- matrix(0, ncol = nn, nrow = x$n)
  for (i in seq_len(nn)) {
    take[, i] <- as.numeric(x$trts==nodes[i])
    volume[, i] <- rep(N, x$n) * nn * d / sum(d)
  }
  colnames(take) <- colnames(volume) <- nodes
  rownames(take) <- rownames(volume) <- x$trts
  # Step 1
  bars <- t(volume) %*% diag(1 / diag(D)) / 2
  colnames(bars) <- x$trts
  ylim <- c(0, max(bars))
  barplot(bars, beside = TRUE, col = cols, ylim = ylim, las = 1,
          ylab = "Volume")
  volume <- volume - take
  # Step 2 to N
  for (i in seq_len(N)) {
    take.new <- Tr %*% take
    volume <- volume - take.new
    bars <- t(volume) %*% diag(1 / diag(D)) / 2
    colnames(bars) <- x$trts
    barplot(bars, beside = TRUE, col = cols, border = NA, las = 1,
            xaxt = "n", yaxt = "n", add = TRUE)
    take <- take.new
  }
  #
  barplot(bars, beside = TRUE, col = cols.full, las = 1, legend.text = nodes,
          args.legend = list(x = "bottomright", bg = "white"),
          xaxt = "n", yaxt = "n", add = TRUE)
                                        #
  volume <- t(volume) %*% diag(1 / diag(D)) / 2
  colnames(volume) <- x$trts
                                        #
  B0 <- netmeta:::createB(ncol = nn)
  #
  comps <- rep(NA, nn * (nn - 1) / 2)
  k <- 1
  for (i in seq_len(nn - 1)) {
    for (j in ((i + 1):nn)) {
      comps[k] <- paste(nodes[i], nodes[j], sep = ":")
      k <- k + 1
    }
  }
  colnames(B0) <- nodes
  rownames(B0) <- comps
  #
  diff <- -t(volume) %*% t(B0)
  C <- B0 %*% diff[nodes,]
  H <- C %*% W[comps, comps]
  #
  res <- list(C = C, H = H)
  invisible(res)
}
