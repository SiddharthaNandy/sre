cholesky_embedding <- function(U, big_dim, TKR_Index, OKnot_Index, j0) {
  
  if (!is.matrix(U)) stop("U must be a matrix")
  p <- ncol(U)
  if (nrow(U) != p) stop("U must be square (nrow == ncol)")
  G <- length(TKR_Index)
  if (length(OKnot_Index) != G) stop("TKR_Index and OKnot_Index must have same number of groups")
  
  big_group_lens <- vapply(OKnot_Index, length, integer(1))
  big_group_starts <- cumsum(c(1, head(big_group_lens, -1)))
  
  small_group_lens <- vapply(TKR_Index, length, integer(1))
  if (sum(small_group_lens) != p) stop("Total columns implied by TKR_Index must equal ncol(U)")
  small_group_starts <- cumsum(c(1, head(small_group_lens, -1)))
  
  B <- matrix(0, big_dim, big_dim)
  
  for (g in seq_len(G)) {
    
    small_ids <- TKR_Index[[g]]
    big_ids   <- OKnot_Index[[g]]
    ns <- length(small_ids)
    nb <- length(big_ids)
    big_start <- big_group_starts[g]
    small_start <- small_group_starts[g]
    
    pos_in_group <- pmatch(small_ids, big_ids)
    if (any(is.na(pos_in_group))) stop(sprintf("Some TKR ids in group %d not found in OKnot_Index[[%d]]", g, g))
    
    for (k in seq_len(ns)) {
      
      sc <- small_start + (k - 1)
      pos <- pos_in_group[k]
      col_idx <- big_start + (pos - 1)
      
      nz_indices <- which(U[, sc] != 0)
      if (length(nz_indices) == 0) next
      nz_len <- max(nz_indices)
      
      if (nz_len <= j0) {
        rows <- seq_len(nz_len)
      } else {
        r <- nz_len - j0
        leading <- seq_len(j0)
        # Tail positions relative to the column in the big group
        tail_rows <- (big_start + pos - r):(big_start + pos - 1)
        rows <- c(leading, tail_rows)
      }
      
      vals <- U[seq_len(nz_len), sc]
      if (length(rows) != length(vals)) stop("Internal length mismatch!")
      B[rows, col_idx] <- vals
      
    }
    
  }
  
  B[lower.tri(B)] <- 0
  return(B)
  
}