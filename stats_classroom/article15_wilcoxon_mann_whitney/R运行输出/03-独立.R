options(width=42)
textbook <- read.csv("../data/textbook_example.csv")
A <- textbook$value[textbook$group == "A"]
B <- textbook$value[textbook$group == "B"]
independent <- wilcox.test(
  A,
  B,
  alternative="greater",
  exact=TRUE,
  conf.int=TRUE
)
independent
