options(width=42)
paired <- read.csv("../data/paired_vitaminA.csv")
X <- paired$X
Y <- paired$Y
paired_result <- wilcox.test(
  X,
  Y,
  paired=TRUE,
  alternative="two.sided",
  exact=TRUE,
  conf.int=TRUE
)
paired_result
