# nolint start
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
if (!require(mgcv)) {install.packages("mgcv"); library(mgcv)}
if (!require(geometry)) {install.packages("geometry"); library(geometry)}
# if (!require(uniformly)) {install.packages("uniformly"); library(uniformly)}
if (!require(furrr)) {install.packages("furrr"); library(furrr)}
if (!require(pracma)) {install.packages("pracma"); library(pracma)}
if (!require(combinat)) {install.packages("combinat"); library(combinat)}
# if (!require(igraph)) {install.packages("igraph"); library(igraph)}
if (!require(readxl)) {install.packages("readxl"); library(readxl)}
library(markovchain)

# ------ Parameter space geometry ------
# function that generates a table of all possible communities for a species pool
# Arguments: N = number of species, nv = initial notation for absence
# Return: the table of all possible communities (as a matrix)
get_compo <- function(num, nv = 0) {
  record <- matrix(nv, nrow = 2^num, ncol = num)
  k <- 2
  for (s in 1:num){
    for (i in 1:choose(num, s)){
      record[k, utils::combn(num, s)[, i]] <- 1
      k <- k + 1
    }
  }
  return(record)
}

# function that partitions the parameter space from a given interaction matrix
# Arguments: in_mat = interaction matrix (representing LV dynamics)
# Return: a list of all possible regions (as matrices)
get_region <- function(in_mat) {
  num <- ncol(in_mat)
  l <- 0
  A <- in_mat
  B <- eye(nrow(in_mat), ncol(in_mat))
  record0 <- get_compo(num, 0)
  region <- list()
  for (l in 1 : 2^num){
    region[[l]] <- A %*% diag(record0[l, ]) + B %*% diag(1 - (record0[l, ]))
  }
  return(region)
}

# function that computes the omega value of a region
# Arguments: vertex = vertex vectors of a region, nsamples = number of samples for Monte Carlo calculation
# Return: the omega value, in [0,1]
calculate_omega <- function(vertex, nsamples = 10000) {
  if (is.null(ncol(vertex))) {
    return(0.5)
  } else {
    num <- nrow(vertex)
    vertex <- apply(vertex, 2, norm2)
    set.seed(1010)
    vertex <- cbind(
      vertex, vertex %*% t(abs(runif_on_sphere(n = nsamples, d = ncol(vertex), r = 1)))
    )
    vertex <- apply(vertex, 2, norm2)
    vertex <- cbind(vertex, rep(0, num))
    vol_ori <- convhulln(t(vertex), output.options = TRUE)$vol
    vol_ball <- pi^(num / 2) / gamma(num / 2 + 1)
    return(vol_ori / vol_ball)
  }
}

# function that computes the omega values of each regions in a partition
# Arguments: region = a list of all possible regions (as matrices)
# Return: a vector of omega values as measurements for partitions
evaluate_region <- function(region) {
  omega_vec <- map_dbl(region, calculate_omega)
  return(omega_vec)
}

# function that calculates omega vectors from interaction matrices
# Arguments: in_mats = a list of interaction matrices
# Return: a list of omega vectors as measurements for partitions
ensem_geometry <- function(in_mats) {
  regionss <- future_map(in_mats, ~get_region(.))
  omegass <- future_map(regionss, evaluate_region, .options = furrr_options(seed = TRUE))
  omegass <- map(omegass, norm1)
  return(omegass)
}

# function that gets the centroids of each region in a partition
# Arguments: region = a list of all possible regions (as matrices)
# Return: a list of list of location vectors of the centroids in parameter space
get_centroids <- function(region) {
  if (is.list(region)) {
    map(region, get_centroids)
  } else {
    center <- region %*% c(rep(1, ncol(region)))
    centrd <- as.vector(center / sqrt(sum(center^2)))
    return(centrd)
  }
}

# function that generates feasibile partitions with uniform distrbuted omega values
# Arguments: num = number of species, steps = number of replications
# Return: a list of omega vectors
unif_omegass <- function(num, steps) {
  omegass <- list(); NUM <- 2^num
  for (i in seq_len(steps)) {
    omega <- runif(NUM, 0, 1)
    omegass[[i]] <- norm1(omega)
  }
  return(omegass)
}

# function that generates feasibile partitions with lognormal distrbuted omega values
# Arguments: num = number of species, steps = number of replications
# logmean = mean value of lognormal distribution, logsd = standard deviation of lognormal distribution
# Return: a list of omega vectors
lnorm_omegass <- function(num, steps, logmean = 0, logsd = 1) {
  omegass <- list(); NUM <- 2^num
  for (i in seq_len(steps)) {
    omega <- rlnorm(NUM, logmean, logsd) # standard logNormal; for mimick Unif, use logmean= -0.8369882, logsd = 0.53636
    omegass[[i]] <- norm1(omega)
  }
  return(omegass)
}

# function that samples centroids from uniformly distributed points on the unit sphere
# Arguments: num = number of species, steps = number of replications, NUM = number of centroids
# Return: a list of location vectors of the centroids in parameter space
random_sphere <- function(num, steps, NUM = 2^num) {
  cent_list <- list()
  for (i in seq_len(steps)) {
    cent_list[[i]] <- runif_on_sphere(NUM, num) %>% t() %>% data.frame() %>% unname() %>% as.list()
  }
  return(cent_list)
}

# function that generates a random interaction matrix from a distribution
# Arguments: num = number of species, stren = standard deviation of interaction
# conne = connectivity of the matrix, dist = distribution (normal or lognormal), mean = mean value of interaction
# Return: a random interaction matrix
inte_matrix_random <- function(num, stren, conne = 1, dist = "norm", mean = 0) {
  if (dist == "lnorm") {
    Inte <- -rlnorm(num * num, mean = mean, sd = stren)
  } else if (dist == "norm") {
    Inte <- rnorm(num * num, mean = mean, sd = stren)
  }
  zeroes <- sample(
    c(rep.int(1, floor(num * num * conne)), rep.int(0, (num * num - floor(num * num * conne))))
  )
  Inte[which(zeroes == 0)] <- 0
  Inte <- matrix(Inte, ncol = num, nrow = num)
  diag(Inte) <- -1
  return(Inte)
}

# function that checks if a matrix is Lotka-dissipiative and thus gloablly stable
# Arguments: A = interaction matrix, thre = threshold of detection
# Return: TRUE if the matrix is Lotka-dissipiative, FALSE otherwise
check_diss <- function(A, thre = 10^(-6)) {
  all(eigen(A + t(A))$values < thre)
}

# function that generates an ensemble of globally stable interaction matrices
# Arguments: dist = ensemble distribution (normal or lognormal); num = number of species;
# steps = number of matrices; stren = standard deviation of interaction; conne = connectivity of each matrix
ensem_inmat_rand <- function(dist, num, steps, stren, conne) {
  i <- 1; j <- 1
  in_mats <- list()
  while (i <= 10000 * steps) {
    in_m <- inte_matrix_random(num, stren, conne, dist)
    if (check_diss(in_m)) {
      in_mats[[j]] <- in_m
      j <- j + 1
    }
    if (j == steps + 1) break
    i <-  i + 1
  }
  return(in_mats)
}

# ------ Information theory ------
# function that computes the weighted entropy of a probability vector
# Arguments: x = probability vector, w = weight vector (the same size as x)
# Return: the weighted entropy of the probability vector
logsum <- function(x, w = c(rep(1, length(x)))) {
  sum <- 0
  for (i in seq_along(x)) {
    if (x[i] > 0) {
      sum <- sum + w[i] * x[i] * log(x[i], base = 2)
    }
  }
  return(sum)
}

# function that computes the entropy of a probability vector weighted by itself
# used for computing the switching capacity of an unstructured transition matrix
# Arguments: pii = probability vector
# Return: the weighted entropy of the probability vector
unst_ent_adap <- function(pii) {
  return(-logsum(pii, pii) - logsum(1 - pii, pii))
}

# function that analyzes a transition matrix with markov theory
# Arguments: t_mat = transition matrix, p_vec = weight vector (can be the stationary distribution s_vec)
# sd_calc = whether to calculate s_vec from t_mat, or use the given p_vec
# Return: a list of markovian properties, such as the stationary distribution and the switching capacity
markov_analysis <- function(t_mat, p_vec, sd_calc) {
  mkv <- list()
  if (!sd_calc) {
    mkv$sd <- norm1(p_vec)
  } else {
    sd <- (t_mat - diag(rep(1, ncol(t_mat)))) %>% t() %>% nullspace()
    mkv$sd <- norm1(sd) %>% as.vector()
  }
  mkv$ent_sd <- -logsum(mkv$sd) %>% as.numeric()
  entropy_adap <- function(trans_mat) {
    pii <- diag(trans_mat)
    return(-logsum(pii, mkv$sd) - logsum(1 - pii, mkv$sd))
  }
  mkv$ent_adap <- entropy_adap(t_mat) %>% as.numeric()
  return(mkv)
}

# function that analyzes a transition matrix with information theory
# Arguments: t_mat = transition matrix, p_vec = weight vector (can be the stationary distribution s_vec)
# sd_calc = whether to calculate s_vec from t_mat, or use the given p_vec
# Return: a vector of stationary distribution entropy and switching capacity
tmat_analysis <- function(t_mat, p_vec, sd_calc) {
  mkv_t <- markov_analysis(t_mat, p_vec, sd_calc)
  ent_sd <- mkv_t$ent_sd
  ent_adap <- mkv_t$ent_adap
  t_ii <- mean(diag(t_mat))
  # wrapping results
  results <- as.matrix(cbind(ent_sd, ent_adap, t_ii))
  colnames(results) <- c("ent_sd", "ent_adap", "<t_ii>")
  return(results)
}

# function that computes the entropy values of matrices
# Arguments: t_mats = a list of transition matrices, p_vecc = a list of weight vectors
# sd_calc = whether to calculate s_vec from t_mat, or use the given p_vec
# Return: a matrix of stationary distribution entropy and switching capacity, each row represents a transition matrix
ensem_trans_analysis <- function(t_mats, p_vecs, sd_calc) {
  if (is.matrix(t_mats) && is.vector(p_vecs)) {
    t_mats <- list(t_mats); p_vecs <- list(p_vecs)
  }
  meas_list <- map2(t_mats, p_vecs, ~tmat_analysis(.x, .y, sd_calc))
  meas_mat <- do.call(rbind, meas_list)
  return(meas_mat)
}

# ------ Gravity model ------
# function that generates transition matrix based on gravity model
# Arguments: omegas = omega vector for a partition, centrds = centroids of regions for a partiton
# beta = structural parameter
# Return: the transition matrix
gravity_tmat <- function(omegas, centrds, beta) {
  NUM <- length(omegas)
  gravi_mat <- matrix(0, NUM, NUM)
  for (li in 1:NUM) {
    for (lj in 1:NUM) {
      if (li == lj) {
        dij <- 0
      } else {
        dij <- acos(centrds[[li]] %*% centrds[[lj]])
      }
      gravi_mat[li, lj] <- omegas[li] * omegas[lj] * exp(-beta * dij)
    }
    gravi_mat[li, ] <- norm1(gravi_mat[li, ])
  }
  return(gravi_mat)
}

# function that generates gravity-model transition matrices from model-driven (LV) approach
# Arguments: in_mats = a list of interaction matrices, omegass = a list of the corrsponding omega values
# beta = structural parameter
# Return: a list of transition matrices
ensem_tmat_LV <- function(in_mats, omegass, beta) {
  regionss <- map(in_mats, ~get_region(.))
  centrdss <- get_centroids(regionss)
  map2(omegass, centrdss, ~gravity_tmat(.x, .y, beta = beta))
}

# function that generates gravity-model transition matrices from model-free (LV) approach
# Arguments: in_mats = a list of interaction matrices, omegass = a list of the corrsponding omega values
# beta = structural parameter
# Return: a list of transition matrices
ensem_tmat_rand <- function(rcentrds, omegass, beta) {
  map2(omegass, rcentrds, ~gravity_tmat(.x, .y, beta = beta))
}

# ------ Empirical time series ------
# function that empirically infers both transition matrix
# and stationary distribution from a presence-absence sequence
# Arguments: pa_seq = a vector of presence-absence sequence, each element is a string of "0" and "1".
# For example, pa_seq[2] = "0010" means a 4-sp community has only the 3rd sp at time stamp 2
# Return: a list of t_mat = emprical transition matrix, p_vec = empirical stationary distribution


seq_analysis <- function(pa_seq) {
  ids <- unique(pa_seq)
  l_NUM <- length(ids); l_time <- length(pa_seq)
  pa_seq <- c(pa_seq, pa_seq[1]) # periodic boundary condition

  freq_mat <- matrix(0, l_NUM, l_NUM)
  colnames(freq_mat) <- rownames(freq_mat) <- ids
  for (t in 1:l_time) {
    freq_mat[pa_seq[t], pa_seq[t + 1]] <- freq_mat[pa_seq[t], pa_seq[t + 1]] + 1
  }
  # colnames(freq_mat) <- rownames(freq_mat) <- NULL
  results <- list()
  results$t_mat <- t(apply(freq_mat, 1, norm1))
  results$p_vec <- norm1(rowSums(freq_mat))
  return(results)
}

# function that coverts the local tibble of abundances-time to the presence-absence sequence
# Arguments: data_local = abundances time series; a row for a taxa and a column for a time stamp,
# pa_thre = threshold of determining presence/absence
# Return: the presence-absence sequence vector
collect2seq <- function(data_local, pa_thre) {
  abun2comp <- function(abun, th = pa_thre) {
    paste(as.numeric(abun >= th), collapse = "")
  }
  return(apply(data_local[, -1], 2, abun2comp) %>% as.vector())
}

# function that coverts the total tibble of abundances-time to a list of transition measurements
# Arguments: row_ids = which row (taxa) to be analyzed together; length(row_ids[i,]) = size of taxon pool
# data_donor = total tibble of abundances-time; a row for a taxa and a column for a time stamp,
# pa_thre = threshold of determining presence/absence
# Return: a list, each element is a list of t_mat = emprical transition matrix and p_vec = empirical stationary distribution
ensem_data2trans <- function(row_ids, data_donor, pa_thre = 10^(-4)) {
  rows_list <- lapply(c(seq_len(nrow(row_ids))), function(i) row_ids[i, ])
  emp_tmats <- future_map(rows_list, ~data_donor[. , ] %>% collect2seq(. , pa_thre) %>% seq_analysis(), .options = furrr_options(seed = TRUE))
  return(emp_tmats)
}

# function that shuffules the local tibble of abundances-time by column and converts it to the presence-absence sequence
# Arguments: data_local = abundances time series; a row for a taxa and a column for a time stamp,
# pa_thre = threshold of determining presence/absence
# Return: the presence-absence sequence vector
rand2seq <- function(data_local, pa_thre) {
  abun2comp <- function(abun, th = pa_thre) {
    paste(sample(as.numeric(abun >= th)), collapse = "")
  }
  return(apply(data_local[, -1], 2, abun2comp) %>% as.vector())
}

# function that shuffules the total tibble of abundances-time by column and coverts it to a list of transition measurements
# Arguments: row_ids = which row (taxa) to be analyzed together; length(row_ids[i,]) = size of taxon pool
# data_donor = total tibble of abundances-time; a row for a taxa and a column for a time stamp,
# pa_thre = threshold of determining presence/absence
# Return: a list, each element is a list of t_mat = emprical transition matrix and p_vec = empirical stationary distribution
ensem_shulffle2trans <- function(row_ids, data_donor, pa_thre = 10^(-4)) {
  rows_list <- lapply(c(seq_len(nrow(row_ids))), function(i) row_ids[i, ])
  emp_tmats <- future_map(rows_list, ~data_donor[. , ] %>% rand2seq(. , pa_thre) %>% seq_analysis(), .options = furrr_options(seed = TRUE))
  return(emp_tmats)
}

telep_mat <- function(t_mat, alpha = 10^-4) {
  if(is.null(nrow(t_mat))) return(t_mat)
  else {
    np <- nrow(t_mat)
    return((1 - alpha) * t_mat + alpha * (1/np) * ones(np, np))
  }
}


# function that coverts the local tibble of abundances-time to the presence-absence sequence
# Arguments: data_local = abundances time series; a row for a taxa and a column for a time stamp,
# pa_thre = threshold of determining presence/absence
# Return: the presence-absence sequence vector
collect2prob <- function(data_local, pa_thre) {
  abun_prob_comp <- function(abun, th = pa_thre) {
    if(sum(abun) < th) {
      return(paste(as.numeric(abun >= th), collapse = ""))
    } else {
      abun <- norm1(abun)
      return(paste(sapply(abun, function(x) rbinom(1, 1, x)), collapse = ""))
    }
  }
  return(apply(data_local[, -1], 2, function(x) abun_prob_comp(x, pa_thre)) %>% as.vector())
}

collect2prob_max <- function(data_local, pa_thre) {
  data_loc <- as.matrix(data_local[, -1]); dimnames(data_loc) <- NULL

  normm <- function(vec) {
    if (max(abs(vec)) == 0) return(vec)
    else return(vec / max(vec)) # maybe can be normalized with mean(vec) or so?
  }
  data_loc <- t(apply(data_loc, 1, normm))

  abun_prob_comp <- function(abun, th = pa_thre) {
    if(sum(abun) < th) {
      return(paste(as.numeric(abun >= th), collapse = ""))
    } else {
      return(paste(sapply(abun, function(x) rbinom(1, 1, x)), collapse = ""))
    }
  }

  return(apply(data_loc, 2, function(x) abun_prob_comp(x, pa_thre)) %>% as.vector())
}

ensem_sto2ts <- function(row_ids, data_donor, pa_thre = 10^(-4)) {
  rows_list <- lapply(c(seq_len(nrow(row_ids))), function(i) row_ids[i, ])
  pa_tss <- future_map(rows_list, ~data_donor[. , ] %>% collect2prob(. , pa_thre), .options = furrr_options(seed = TRUE))
  return(pa_tss)
}

ensem_sto2rand_ts <- function(pa_tss) {
  shuffle_pa <- function(pa_ts) {
    if (length(pa_ts) > 1){
      sapply(pa_ts, shuffle_pa, USE.NAMES = FALSE)
    } else {
      return(paste(sample(strsplit(pa_ts, "")[[1]]), collapse = ""))
    }
  }
  sf_pa_tss <- future_map(pa_tss, ~shuffle_pa(.x), .options = furrr_options(seed = TRUE))
  return(sf_pa_tss)
}

ensem_pa2trans <- function(pa_tss) {
  emp_tmats <- future_map(pa_tss, ~seq_analysis(.x), .options = furrr_options(seed = TRUE))
  return(emp_tmats)
}

# function that implements the empirical analysis of time series data
# Arguments: sizes = a range vector of size of taxon pool; steps = number of replications
# data_ts = total time series of abundances-time; rep_num = num of replications for shuffling
# rep_num = num. of replications for shuffling
# can_thre = ratio of the maximum abundance to determine the candidate taxa
# Return: a tibble of empirical measurements
analyze_empirical_timeseries <- function(sizes, steps, data_ts, rep_num = 1, can_thre = 1e-6, pa_thre = 1e-4, tele_alpha = 1e-4) {
  emp_total <- tibble(size = double(), type = character(), ent_adap = double(), mean_t_ii = double())
  rowMax <- function(mat) {
    apply(mat, 1, max)
  }
  presence <- data_ts %>% transmute(total = rowMax(select(., -1)))
  candidates <- which(presence > can_thre)

  for (num in sizes) {
    row_ids <- c()
    set.seed(180)
    for (k in seq_len(steps)) {
      row_ids <- rbind(row_ids, sample(candidates, num))
    }

    # Regular empirical analysis
    emp <- ensem_data2trans(row_ids, data_ts, pa_thre)
    t_mats <- map(emp, "t_mat")
    f_vecs <- map(emp, "p_vec")
    results <- ensem_trans_analysis(t_mats, f_vecs, sd_calc = FALSE)
    emp_total <- rbind(
      emp_total,
      tibble(size = num, type = "empirical", ent_adap = results[, 2], mean_t_ii = results[, 3])
    )

    # # Teleportation
    # tt_mats <- map(t_mats, ~telep_mat(. , alpha = tele_alpha))
    # f_vecs <- map(emp, "p_vec")
    # results <- ensem_trans_analysis(tt_mats, f_vecs, sd_calc = FALSE)
    # emp_total <- rbind(
    #   emp_total,
    #   tibble(size = num, type = "telep", ent_adap = results[, 2], mean_t_ii = results[, 3])
    # )

    # Regular randomized analysis
    results <- list()
    for (r in 1:rep_num) {
      emppd <- ensem_shulffle2trans(row_ids, data_ts, pa_thre)
      t_mats <- map(emppd, "t_mat")
      f_vecs <- map(emppd, "p_vec")
      results[[r]] <- ensem_trans_analysis(t_mats, f_vecs, sd_calc = FALSE)
    }
    ave_results <- Reduce("+", results) / length(results)
    emp_total <- rbind(
      emp_total,
      tibble(size = num, type = "randomized", ent_adap = ave_results[, 2], mean_t_ii = ave_results[, 3])
    )

  #   # Stochastic empirical analysis
  #   results_sto <- list()
  #   results_sf_sto <- list()

  #   for (r in 1:rep_num) {
  #     sto_emp_ts <- ensem_sto2ts(row_ids, data_ts, pa_thre)
  #     emp_sto <- ensem_pa2trans(sto_emp_ts)
  #     t_mats <- map(emp_sto, "t_mat")
  #     f_vecs <- map(emp_sto, "p_vec")
  #     results_sto[[r]] <- ensem_trans_analysis(t_mats, f_vecs, sd_calc = FALSE)

  #     sf_sto_emp_ts <- ensem_sto2rand_ts(sto_emp_ts)
  #     emp_sf_sto <- ensem_pa2trans(sf_sto_emp_ts)
  #     t_mats <- map(emp_sf_sto, "t_mat")
  #     f_vecs <- map(emp_sf_sto, "p_vec")
  #     results_sf_sto[[r]] <- ensem_trans_analysis(t_mats, f_vecs, sd_calc = FALSE)
  #   }

  #   ave_results_1 <- Reduce("+", results_sto) / length(results_sto)
  #   ave_results_2 <- Reduce("+", results_sf_sto) / length(results_sf_sto)

  #   emp_total <- rbind(
  #     emp_total,
  #     tibble(size = num, type = "emp_sto", ent_adap = ave_results_1[, 2])
  #   )

  #   emp_total <- rbind(
  #     emp_total,
  #     tibble(size = num, type = "rand_sto", ent_adap = ave_results_2[, 2])
  #   )

  }
  return(emp_total)
}


# function that binarizes the total time series tibble, for PCoA analysis
# Arguments: data = time series of abundances-time
# Return: a matrix in the same shape, but with 0/1 values
bin_pa <- function(data) {
  bin <- function(value) {
    if (value > 10^-4) 1
    else 0
  }
  binpa <- matrix(0, nrow(data), ncol(data))
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(ncol(data))) {
      binpa[i, j] <- bin(data[i, j])
    }
  }
  return(binpa)
}

# ------ Plot and Snippets ------
# function that customizes the gplot theme
# Arguments ... =  other theme options
theme_my <- function(...) {
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    ...
  )
}

# sinppet for 1-norm and 2-norm normalization of a vector
norm1 <- function(vec) {
  vec / sum(vec)
}
norm2 <- function(vec) {
  vec / sqrt(sum(vec^2))
}
dist2 <- function(vec) {
  sum(vec^2)
}
runif_on_sphere <- function(n, d, r = 1) {
  sims <- matrix(rnorm(n * d), nrow = n, ncol = d)
  r * sims / sqrt(apply(sims, 1L, crossprod))
}
rowMax <- function(mat) {
  apply(mat, 1, max)
}

simulate_mc <- function(t_mat, nt, state0) {
  if (nrow(t_mat) == 1) {
    state1 <- colnames(t_mat)[1]
    return(replicate(nt, state1))
  }
  mc <- new("markovchain", states = rownames(t_mat), byrow = TRUE, transitionMatrix = t_mat, name = "test")
  set.seed(123)
  mcSim <- markovchain::rmarkovchain(n = nt, object = mc, t0 = state0)
  return(mcSim)
}

pa_ts_dist <- function(seq1, seq2) {
  if (length(seq1) != length(seq2)) stop("Lengths of sequences not equal")
  seq1_l <- strsplit(seq1, "") %>% map(. , as.integer)
  seq2_l <- strsplit(seq2, "") %>% map(. , as.integer)
  dist <- map2_dbl(seq1_l, seq2_l, ~dist2(.x - .y))
  sqrt(mean(dist)) # max: sqrt(num)
}

generate_null_pa <- function(num, prob = 0.5) {paste(rbinom(num, 1, prob),collapse = "")}

generate_rand_pa <- function(num, omegas) {
  sample(names(omegas), 1, prob = omegas)
}

generate_uns_pa <- function(num, fvec) {sample(names(fvec), 1, prob = fvec)}
generate_zero_pa <- function(num) {paste(c(rep(0, num)), collapse = "")}

rep_ts <- function(ts, nrep) {
  rep_ts <- c(rep(NA, length(ts) * nrep))
  for (k in (seq_len(nrep) - 1)) {
    rep_ts[nrep * seq_along(ts) - k] <- ts
  }
  return(rep_ts)
}
shuffle_ts <- function(ts) {
  strsplit(ts, "") %>% map(sample) %>% map_chr(paste, collapse = "")
}
copy_ts <- function(ts, ncopy) {
  c(rep(ts, ncopy))
}

# nolint end