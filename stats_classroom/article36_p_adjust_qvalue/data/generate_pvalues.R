# 统计课堂36：共享教学P值；这是模拟模型，不是患者原始资料。
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261121)
m <- 1000
p <- c(runif(800), rbeta(200, shape1 = 0.2, shape2 = 1))
data <- data.frame(id = seq_len(m), p = p, is_null = c(rep(1, 800), rep(0, 200)))
write.csv(data, "data/screening_pvalues.csv", row.names = FALSE)
