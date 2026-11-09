# 三项配对率差的固定D排序精确区间，base R helper。
# 由paired_intervals.R调用，与PairedCI的归纳构造分别报告。
paired_exact <- function(nplus,nzero,nminus,conf_level=.95,grid_points=201) {
  counts <- c(nplus,nzero,nminus)
  if (length(counts)!=3 || any(!is.finite(counts)) || any(counts<0) || any(counts!=floor(counts)) || sum(counts)==0) stop("三类人数须为非负整数，总和须大于0。")
  if (length(conf_level)!=1 || !is.finite(conf_level) || conf_level<=0 || conf_level>=1) stop("置信水平须在0与1之间。")
  if (length(grid_points)!=1 || !is.finite(grid_points) || grid_points<3 || grid_points!=floor(grid_points)) stop("网格数须为至少3的整数。")
  n <- sum(counts)
  observed <- (nplus-nminus)/n
  alpha <- (1-conf_level)/2
  cells <- expand.grid(plus=0:n,minus=0:n)
  cells <- cells[cells$plus+cells$minus<=n,]
  cells$zero <- n-cells$plus-cells$minus
  cell_counts <- as.matrix(cells[,c("plus","zero","minus")])
  log_coefficient <- lgamma(n+1)-rowSums(lgamma(cell_counts+1))
  differences <- cells$plus-cells$minus
  upper_mask <- as.numeric(differences>=nplus-nminus)
  lower_mask <- as.numeric(differences<=nplus-nminus)
  maximum_tail <- function(theta,upper) {
    left <- 0
    right <- max(0,1-abs(theta))
    mask <- if (upper) upper_mask else lower_mask
    probability <- function(s) {
      probabilities <- cbind(pmin(1,pmax(0,(1-s+theta)/2)),s,pmin(1,pmax(0,(1-s-theta)/2)))
      log_mass <- matrix(log_coefficient,nrow=length(s),ncol=nrow(cell_counts),byrow=TRUE)
      for (category in 1:3) {
        terms <- outer(log(probabilities[,category]),cell_counts[,category])
        # 0乘log(0)在这里表示0次出现，该概率幂等于1。
        terms[is.nan(terms)] <- 0
        log_mass <- log_mass+terms
      }
      pmin(1,pmax(0,as.numeric(exp(log_mass)%*%mask)))
    }
    if (right<=left) return(probability(left))
    grid <- seq(left,right,length.out=grid_points)
    values <- probability(grid)
    best <- max(values)
    if (best>=1-1e-14) return(1)
    center <- values[2:(grid_points-1)]
    previous <- values[1:(grid_points-2)]
    following <- values[3:grid_points]
    peaks <- which(center>=previous & center>=following & (center>previous | center>following))+1
    peaks <- unique(c(peaks,which.max(values)))
    for (index in peaks) {
      if (index>1 && index<grid_points) {
        fitted <- optimize(function(s) -probability(s),c(grid[index-1],grid[index+1]),tol=1e-12)
        best <- max(best,-fitted$objective)
      }
    }
    best
  }
  endpoint <- function(left,right,upper) {
    for (iteration in 1:60) {
      middle <- (left+right)/2
      probability <- maximum_tail(middle,upper)
      if ((probability<alpha)==upper) left <- middle else right <- middle
      if (right-left<=1e-7) break
    }
    if (upper) left else right
  }
  lower <- if (observed==-1) -1 else endpoint(-1,observed,TRUE)
  upper <- if (observed==1) 1 else endpoint(observed,1,FALSE)
  p_value <- min(1,2*min(maximum_tail(0,TRUE),maximum_tail(0,FALSE)))
  data.frame(method="paired_fixed_D_central",estimate=observed,lower=lower,upper=upper,p_value=p_value,grid_points=grid_points)
}
