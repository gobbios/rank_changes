# these are the helper functions to run the simulations
# partly they are written fairly awkwardly (= inefficient and slow), sorry...

# just returns the femdom index given an ordered(!) matrix and the index of 
#   females' positions in the matrix (i.e. their rank)
get_fem_dom <- function(m, fem_index) {
  dmat <- m > t(m)
  dmat[m == t(m)] <- NA
  # for females: how many males is each female dominant over
  fem_dom <- rowSums(dmat[fem_index, -fem_index, drop = FALSE], na.rm = TRUE)
  # get proportion and mean
  mean(fem_dom / (ncol(m) - length(fem_index)))
}


# creates matrix with fixed (within tolerance) femdom value, given total
#   number of ids and number of females
# this is an iterative algorithm, so it might not produce a result at all
#   and return 'NULL'
# females in the result have 3 letter codes starting with 'F'
# males have 4-letter codes starting with 'm'
mat_fem_dom <- function(n_ids, n_fem, fem_dom_val, tol = 0.05, max_tries = 10000) {
  m <- matrix(nrow = n_ids, ncol = n_ids, 0)
  m[upper.tri(m)] <- 1
  still_looking <- TRUE
  break_counter <- 0
  res <- NULL
  while (still_looking) {
    fem_ids <- sample(seq_len(n_ids), n_fem)
    xres <- get_fem_dom(m, fem_index = fem_ids)
    if (xres < fem_dom_val + tol & xres > fem_dom_val - tol) {
      still_looking <- FALSE
      res <- fem_ids
    } else {
      break_counter <- break_counter + 1
    }
    if (break_counter > max_tries) break
  }
  
  if (!is.null(res)) {
    # create ids and apply to column and row names
    #   and make sure there are no duplicates...
    mids <- unique(sapply(seq_len(1000), function(x) paste0("m", paste(sample(letters, 2, TRUE), collapse = ""), collapse = "")))
    fids <- unique(sapply(seq_len(1000), function(x) paste0("F", paste(sample(LETTERS, 3, TRUE), collapse = ""), collapse = "")))
    ids <- mids[seq_len(n_ids)]
    ids[res] <- fids[seq_len(n_fem)]
    
    colnames(m) <- ids
    rownames(m) <- ids
    return(m)
  } else {
    return(NULL)
  }
}

# returns 3 values for unknown relationships:
#   total (across all dyads)
#   mixed sex (female/male dyads)
#   within sex (female/female AND male/male combined)
split_prunks <- function(m) {
  m <- m[sort(colnames(m)), sort(colnames(m))]
  m <- m + t(m)
  
  f_index <- which(substr(colnames(m), 1, 1) == "F")
  m_index <- which(substr(colnames(m), 1, 1) == "m")
  
  res <- c(mixed_dyads = mean(m[f_index, m_index] < 1),
           within_sex_dyads = mean(c(m[f_index, f_index][upper.tri(m[f_index, f_index])], 
                                     m[m_index, m_index][upper.tri(m[m_index, m_index])]) < 1),
           total = mean(m[upper.tri(m)] < 1)
  )
  res
}

# remove interactions to increase unknown relationships to specified value
#   (within tolerance)
thin_matrix <- function(m, total_prunk, tol = 0.05, max_tries = 10000) {
  n_dyads <- ncol(m) * (ncol(m) - 1)  / 2
  total_prunks <- round(n_dyads * total_prunk)
  
  still_looking <- TRUE
  break_counter <- 0
  res <- NULL
  
  while (still_looking) {
    test_mat <- m
    test_mat[sample(which(upper.tri(m)), total_prunks)] <- 0 
    test_res <- split_prunks(test_mat)
    if (test_res["total"] < total_prunk + tol & test_res["total"] > total_prunk - tol) {
      res <- test_mat
      still_looking <- FALSE
      return(test_mat)
    } else {
      break_counter <- break_counter + 1
    }
    if (break_counter > max_tries) break
  }
  
  
}


# gets ranks and rank changes, first based on ground truth
# returns a list with 3 data frames: full (mixed sex), females, and males
# each list contains columns for true rank, DS rank, and whether
#   rank differs between truth and DS
analyze_matrix <- function(m) {
  f_index <- which(substr(colnames(m), 1, 1) == "F")
  m_index <- which(substr(colnames(m), 1, 1) == "m")
  
  # ground truth: taken from order matrix is supplied with
  full <- data.frame(id = colnames(m))
  full$truth <- seq_len(nrow(full))
  f_ranks <- data.frame(id = colnames(m)[f_index])
  f_ranks$truth <- seq_len(nrow(f_ranks))
  m_ranks <- data.frame(id = colnames(m)[m_index])
  m_ranks$truth <- seq_len(nrow(m_ranks))
  
  # results for 'empirical' after reordering according to DS
  full$ds <- NA
  f_ranks$ds <- NA
  m_ranks$ds <- NA
  
  
  # reorder according to DS
  x <- EloRating::DS(m, prop = "Dij")
  x <- x[order(x$normDS, decreasing = TRUE), ]  
  new_order <- x$ID
  for (i in seq_len(nrow(full))) {
    full$ds[i] <- which(new_order == full$id[i])
  }
  
  new_order_fem <- new_order[substr(new_order, 1, 1) == "F"]
  for (i in seq_len(nrow(f_ranks))) {
    f_ranks$ds[i] <- which(new_order_fem == f_ranks$id[i])
  }
  
  new_order_males <- new_order[substr(new_order, 1, 1) == "m"]
  for (i in seq_len(nrow(m_ranks))) {
    m_ranks$ds[i] <- which(new_order_males == m_ranks$id[i])
  }
  
  full$changed <- as.numeric(full$truth != full$ds)
  f_ranks$changed <- as.numeric(f_ranks$truth != f_ranks$ds)
  m_ranks$changed <- as.numeric(m_ranks$truth != m_ranks$ds)
  
  
  list(full = full, m_ranks = m_ranks, f_ranks = f_ranks)
}
