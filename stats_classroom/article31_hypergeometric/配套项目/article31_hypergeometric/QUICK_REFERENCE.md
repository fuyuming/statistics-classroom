# 三语函数速查

|任务|R|Python|MATLAB|
|---|---|---|---|
|组合数|`choose(n,k)`|`math.comb(n,k)`|`nchoosek(n,k)`|
|概率质量|`dhyper(x,m,n,k)`|`hypergeom.pmf(x,M,n,N)`|`hygepdf(x,M,K,N)`|
|左尾|`phyper(4,m,n,k)`|`hypergeom.cdf(4,M,n,N)`|`sum(hygepdf(0:4,M,K,N))`|
|右尾|`phyper(3,m,n,k,lower.tail=FALSE)`|`hypergeom.sf(3,M,n,N)`|`sum(hygepdf(4:9,M,K,N))`|
|Fisher 精确检验|`fisher.test(tab)`|`fisher_exact(tab)`|`[h,p]=fishertest(tab,'Tail','both')`|

这里 R 的 `m,n,k` 分别是总体成功数、失败数、抽取数；Python 的 `M,n,N` 和 MATLAB 的 `M,K,N` 只是参数命名不同，依次均对应 `33,9,22`；R的三个参数则为 `9,24,22`。
