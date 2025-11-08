library(igraph)
#参数
#大小n,社区数K，实验重复次数rept
n <- 1200
K <- 2
rept <- 100

#实验a
a <- 0.8
c0 <- 0.5
d0 <- 0.2
# #实验b
# a <- 0.6
# c0 <- 0.7
# d0 <- 0.025

#概率矩阵P
pdata=c(1,a,a,1)
P<-matrix(pdata,K,K)
P

#度异质参量theta
i <- 1:n
gamma <- c(0.5,1,1.5,2,2.5)
#情形1
theta1 <- c0+(c0-d0)*(i/n)^gamma[1]
#情形2
theta2 <- c0+(c0-d0)*(i/n)^gamma[2]
#情形3
theta3 <- c0+(c0-d0)*(i/n)^gamma[3]
#情形4
theta4 <- c0+(c0-d0)*(i/n)^gamma[4]
#情形5
theta5 <- c0+(c0-d0)*(i/n)^gamma[5]
theta <- list(theta1,theta2,theta3,theta4,theta5)

#标签向量l
l0<-rbinom(n,1,0.5)+1
l0


#生成主要信息矩阵omega
fomega <- function(t = 1)
{
  theta <- theta[[t]]
  omega<-matrix(0,n,n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {  
      omega[i,j] <- theta[i]*P[l0[i],l0[j]]*theta[j]
    }
  }
  return(omega)
}

#生成噪音矩阵W
fW <- function(omega)
{
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
  return(W)
}

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
fA <- function(omega,W)
{
  AA<-omega-diag(diag(omega))+W
  g1 <- graph_from_adjacency_matrix(AA,mode = "undirected", weighted = NULL,
                                    diag = TRUE, add.colnames = NULL, add.rownames = NA)
  V(g1)$name <- as.character(1:n)
  # plot(g1)
  
  #找最大连通分量
  if (!is.connected(g1))
  {
    cg1<-decompose.graph(g1)#cg1为连通分量列表
    cgm1<-cg1[[which.max(sapply(cg1,vcount))]]#lst[[]]提取元素；lst[]提取子列表
  }else{
    cgm1<-g1
  }
  A0<-get.adjacency(cgm1)#此时A0为dgCMatrix类下的矩阵
  A0<-as.matrix(A0)#转换为矩阵格式
  return(list(A0,g1,cgm1))
}

#取出对应连通分量的标签
fl <- function(A0,l0,g1,cgm1)
{
  n0<-dim(A0)[1]
  if (n0 != n)
  {
    vg1 <- as.vector(V(g1)$name)
    vcgm1 <- as.vector(V(cgm1)$name)
    uc <- as.integer(setdiff(vg1,vcgm1))
    l <- l0[-uc]
  }
  else{
    l <- l0
  }
  return(l)
}

ScoreError <- c()
OPCAError <- c()
for (t in 1:5)
{
  omega <- fomega(t)
  W <- fW(omega)
  lst <- fA(omega,W)
  A0 <- lst[[1]]
  g1 <- lst[[2]]
  cgm1 <- lst[[3]]
  l <- fl(A0,l0,g1,cgm1) 
  source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/SCORE_change.R')
  ex1_out <- SCORE(A0,K)
  ScoreErrorLst <- getError(ex1_out,l)
  ScoreError <- append(ScoreError,ScoreErrorLst$error.rate)
  source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/PCA.R')
  ex1_out <- SCORE(A0,K)
  OPCAErrorLst <- getError(ex1_out,l)
  OPCAError <- append(OPCAError,OPCAErrorLst$error.rate)
}
ScoreError
OPCAError
plot(gamma,ScoreError, type='o', lwd=2, col="red",pch = 19,ylab = 'ErrorRate')
lines(gamma,OPCAError,lty = 4, lwd=2, col="green")
legend(x='topleft',legend = c('SCORE','oPCA'),col = c('red','green'),lty = 1)