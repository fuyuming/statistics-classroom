options(width=42)
paired <- read.csv("../data/paired_vitaminA.csv")
X <- paired$X
Y <- paired$Y
d <- round(X - Y, 2)
one_sample <- wilcox.test(
  d,
  mu=0,
  alternative="two.sided",
  exact=TRUE,
  conf.int=TRUE
)
one_sample
