#大小n,社区数K，实验重复次数rept
n <- 1200
K <- 2
rept <- 25

#概率矩阵P
a <- 0.1*c(3,4,5,6,7)
# # 情形1
# pdata=c(1,a[1],a[1],1)
# #情形2
# pdata=c(1,a[2],a[2],1)
# #情形3
# pdata=c(1,a[3],a[3],1)
# #情形4
# pdata=c(1,a[4],a[4],1)
#情形5
pdata=c(1,a[5],a[5],1)
P<-matrix(pdata,K,K)
P

#度异质参量theta
i <- 1:n
c0 <- 0.7
d0 <- 0.025
theta <- c0+(c0-d0)*(i/n)^2

#标签向量l
l<-rbinom(n,1,0.5)+1
l

#生成主要信息矩阵omega
omega<-matrix(0,n,n)
for (i in 1:n)
{
  for (j in 1:n)
  {  
    omega[i,j] <- theta[i]*P[l[i],l[j]]*theta[j]
  }
}
omega
#生成噪音矩阵W
W<-matrix(0,n,n)
for (i in 1:n)
{
  for (j in 1:n)
  { 
    if (j>i){
      W[i,j] <- rbinom(1,1,omega[i,j])-omega[i,j]
    }
    
  }
}
W<-W+t(W)
isSymmetric(W)
#生成N（V,E)的邻接矩阵AA,并找到其最大连通分量N0(V0,E0)的邻接矩阵A0

# PAA<-omega-diag(diag(omega))
# PAA
# AA<-matrix(0,n,n)
# for (i in 1:n)
# {
#   for (j in 1:n)
#   {
#     AA[i,j] <- rbinom(1,1,PAA[i,j])
#   }
# }
# AA<-AA+W

AA<-omega-diag(diag(omega))+W
library(igraph)
g1 <- graph_from_adjacency_matrix(AA,mode = "undirected", weighted = NULL,
                                  diag = TRUE, add.colnames = NULL, add.rownames = NA)
V(g1)$name <- as.character(1:n)
# plot(g1)

#找最大连通分量
if (!is.connected(g1))
{
  cg1<-decompose.graph(g1)#cg1为连通分量列表
  cgm1<-cg1[[which.max(sapply(cg1,vcount))]]#lst[[]]提取元素；lst[]提取子列表
}else
{
  cgm1<-g1
}
A0<-get.adjacency(cgm1)#此时A0为dgCMatrix类下的矩阵
A0<-as.matrix(A0)#转换为矩阵格式
n0<-dim(A0)[1]
n0

#取出对应连通分量的标签
if (n0 != n)
{
  vg1 <- as.vector(V(g1)$name)
  vcgm1 <- as.vector(V(cgm1)$name)
  uc <- as.integer(setdiff(vg1,vcgm1))
  l <- l[-uc]
}


source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/SCORE_change.R')
ex1_out <- SCORE(A0,K)
getError(ex1_out,l)

source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/PCA.R')
ex1_out <- SCORE(A0,K)
getError(ex1_out,l)
