# 统计课堂09：五个置信区间与覆盖率模拟（R总入口）
# 打开根目录的Rproj，在RStudio点Source完整运行。只用R自带包。
# 六个编号脚本按问题组织；本入口逐个调用，结果保存在运行结果/R。
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if (length(file_arg)>0) {
  root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else if (!is.null(source_file)) {
  root <- dirname(dirname(normalizePath(source_file)))
} else {
  root <- getwd()
}
# 先认识共享CSV中的字段及数据。
for (f in list.files(file.path(root,"数据"), pattern="csv$", full.names=TRUE)) {
  cat("\n数据：", basename(f), "\n")
  print(read.csv(f))
}
for (f in list.files(file.path(root,"代码"), pattern="^0[1-6].*\\.R$", full.names=TRUE)) {
  cat("\n正在运行：", basename(f), "\n")
  source(f, encoding="UTF-8")
}
# 集中保存五题结果，便于与Python、MATLAB核对。
result_files <- list.files(file.path(root,"运行结果","R"), pattern="^0[1-5].*csv$", full.names=TRUE)
all_results <- lapply(result_files, read.csv)
all_results <- Filter(function(x) all(c("estimate","lower","upper") %in% names(x)), all_results)
all_results <- do.call(rbind, lapply(all_results, function(x) x[1,c("estimate","lower","upper")]))
write.csv(all_results, file.path(root,"运行结果","R","all_intervals.csv"), row.names=FALSE)
