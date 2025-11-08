start_time <- Sys.time() # 记录初始时间
set.seed(1)
library(igraph)
library(EnvStats)
library(ScorePlus)
#参数
#大小n,社区数K，实验重复次数rept
n <- 1000
K <- 4
rept <- 2
#概率矩阵P
#情形1
pdata1=c(1,1/3,1/3,1/3,1/3,1,1/3,1/3,1/3,1/3,1,1/3,1/3,1/3,1/3,1)
P1<-matrix(pdata1,K,K)
#情形2
pdata2=c(1,2/3,0.1,0.1,2/3,1,0.5,0.5,0.1,0.5,1,0.5,0.1,0.5,0.5,1)
P2<-matrix(pdata2,K,K)
P <- list(P1,P2)
#度异质参量theta
alpha <- 5
beta <- 4/5
cn <- 3*log(n)/n
theta <- cn*rpareto(n,beta,alpha)
# i <- 1:n
# theta<-0.015+0.785*(i/n)^2

#标签向量l0
l0 <- 1:n
l0[1:(n/4)] <- 1
l0[(n/4+1):(n/2)] <- 2
l0[(n/2+1):(3*n/4)] <- 3
l0[(3*n/4+1):n] <- 4

#生成主要信息矩阵omega
fomega <- function(t = 1)
{
  theta <- theta
  omega<-matrix(0,n,n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {  
      omega[i,j] <- theta[i]*P[[t]][l0[i],l0[j]]*theta[j]
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

ScoreError <- list(c(),c())
OPCAError <- list(c(),c())
nPCAError <- list(c(),c())
RSCError <- list(c(),c())
SCORE_2_Error <- list(c(),c())
SCORE_1_Error <- list(c(),c())
SCOREPlusError <- list(c(),c())
meanScore <- c()
meanoPCA <- c()
meannPCA <- c()
meanSCORE_2 <- c()
meanSCORE_1 <- c()
meanSCOREPlus <- c()
meanRSC <- c()
varScore <- c()
varoPCA <- c()
varnPCA <- c()
varSCORE_2 <- c()
varSCORE_1 <- c()
varSCOREPlus <- c()
varRSC<- c()
for (t in 1:2)
{
  omega <- fomega(t)
  for (k in 1:rept)
  {
    W <- fW(omega)
    lst <- fA(omega,W)
    A0 <- lst[[1]]
    g1 <- lst[[2]]
    cgm1 <- lst[[3]]
    l <- fl(A0,l0,g1,cgm1) 
    source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/SCORE_change.R')
    ex1_out <- SCORE(A0,K)
    ScoreErrorLst <- getError(ex1_out,l)
    ScoreError[[t]] <- append(ScoreError[[t]],ScoreErrorLst$error.rate)
    source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/oPCA.R')
    ex2_out <- oPCA(A0,K)
    OPCAErrorLst <- getError(ex2_out,l)
    OPCAError[[t]] <- append(OPCAError[[t]],OPCAErrorLst$error.rate)
    source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/SCORE_q.R')
    ex3_out <- SCORE_q(A0,K,2)
    SCORE_2_ErrorLst <- getError(ex3_out,l)
    SCORE_2_Error[[t]] <- append(SCORE_2_Error[[t]],SCORE_2_ErrorLst$error.rate)
    source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/SCORE_q.R')
    ex4_out <- SCORE_q(A0,K,1)
    SCORE_1_ErrorLst <- getError(ex4_out,l)
    SCORE_1_Error[[t]] <- append(SCORE_1_Error[[t]],SCORE_1_ErrorLst$error.rate)
    source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/nPCA.R')
    ex5_out <- nPCA(A0,K)
    nPCAErrorLst <- getError(ex5_out,l)
    nPCAError[[t]] <- append(nPCAError[[t]],nPCAErrorLst$error.rate)
    source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/RSC.R')
    ex6_out <- RSC(A0,K)
    RSCErrorLst <- getError(ex6_out,l)
    RSCError[[t]] <- append(RSCError[[t]],RSCErrorLst$error.rate)
    source('D:/Document backup/OneDrive/常用文档/本研/社区发现/数据与程序/算法/代码/SCOREplus.R')
    # ex7_out <- SCOREplus(A0,K)
    ex7_out = ScorePlus::SCOREplus(A0,K)$labels

    SCOREPlusErrorLst <- getError(ex7_out,l)
    SCOREPlusError[[t]] <- append(SCOREPlusError[[t]],SCOREPlusErrorLst$error.rate)
  }
  meanScore <- c(meanScore,mean(ScoreError[[t]]))
  meanoPCA <- c(meanoPCA,mean(OPCAError[[t]]))
  meanSCORE_2 <- c(meanSCORE_2,mean(SCORE_2_Error[[t]]))
  meanSCORE_1 <- c(meanSCORE_1,mean(SCORE_1_Error[[t]]))
  meanSCOREPlus <- c(meanSCOREPlus,mean(SCOREPlusError[[t]]))
  meannPCA <- c(meannPCA,mean(nPCAError[[t]]))
  meanRSC <- c(meanRSC,mean(RSCError[[t]]))
  varScore <- c(varScore,var(ScoreError[[t]]))
  varoPCA <- c(varoPCA,var(OPCAError[[t]]))
  varSCORE_2 <- c(varSCORE_2,var(SCORE_2_Error[[t]]))
  varSCORE_1 <- c(varSCORE_1,var(SCORE_1_Error[[t]]))
  varSCOREPlus <- c(varSCOREPlus,var(SCOREPlusError[[t]]))
  varnPCA <- c(varnPCA,var(nPCAError[[t]]))
  varRSC <- c(varRSC,var(RSCError[[t]]))
}
ScoreError 
OPCAError 
nPCAError
RSCError 
SCORE_2_Error 
SCORE_1_Error 
SCOREPlusError 
meanScore 
meanoPCA 
meannPCA 
meanSCORE_2 
meanSCORE_1 
meanSCOREPlus 
meanRSC 
varScore 
varoPCA 
varnPCA 
varSCORE_2 
varSCORE_1 
varSCOREPlus 
varRSC

end_time <- Sys.time() # 记录终止时间
end_time - start_time # 计算时间差