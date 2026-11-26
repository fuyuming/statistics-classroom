# 手动运行本文件可按固定种子重建同一份CSV；主入口只读取CSV，不会自动重建。
project_root <- normalizePath(getwd(), mustWork = FALSE)
if (basename(project_root) == "代码") project_root <- dirname(project_root)
set.seed(20261124)
z <- matrix(rnorm(15), nrow=5, ncol=3)
mu <- c(10,20,30); groups <- c("A1","A2","A3"); rows <- list(); j <- 1
for (s in "A") for (i in 1:3) for (r in 1:5) {
  scale <- 2
  rows[[j]] <- data.frame(scenario=s,group=groups[i],replicate=r,value=round(mu[i]+scale*z[r,i],2),mu=mu[i],z=z[r,i]); j <- j+1
}
out <- file.path(project_root,"数据","residual_cases.csv"); dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
write.table(do.call(rbind,rows),out,row.names=FALSE,quote=FALSE,sep=",",fileEncoding="UTF-8")
