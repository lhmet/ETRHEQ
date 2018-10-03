#https://sourceforge.net/p/octave/data-smoothing/ci/584d6abf26660549ff38c6a8d5fec86d2c1cef8b/tree/inst/smooth.m#l202
smooth_matlab <- function(x, k) {
  stopifnot(k > 0, k %% 2 != 0)
  
  sapply(
    1:length(x),
    function(i) {
      if (i <= (k - 1) / 2) {
        idx1 <- 1
        idx2 <- 2 * i - 1
      } else {
        if (i <= length(x) - (k - 1) / 2) {
          idx1 <- i - (k - 1) / 2
          idx2 <- i + (k - 1) / 2
        } else {
          idx1 <- i - (length(x) - i)
          idx2 <- i + (length(x) - i)
        }
      }
      res <- mean(x[idx1:idx2], na.rm = TRUE)
      res[is.na(res)] <- NA
      return(res)
    }
  )
}