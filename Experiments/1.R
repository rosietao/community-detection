source('D:/Extracurricular/community_detection/2_Code/algorithm/PCA.R')

#大小n,社区数K，实验重复次数rept
n <- 1000
K <- 2
rept <- 50
#概率矩阵P
pdata=c(1,0.5,0.5,1)
P<-matrix(pdata,K,K)
P
#度异质参量theta
theta <- rep(0.2,n)
theta

result_opca = c()
result_npca = c()
result_score = c()
result_scoreplus = c()

for(m in 1:rept){
  #标签向量l
  l<-rbinom(n,1,0.5)+1
  #生成主要信息矩阵omega
  omega<-matrix(0,n,n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {  
      omega[i,j] <- theta[i]*P[l[i],l[j]]*theta[j]
    }
  }
  #生成噪音矩阵W
  W<-matrix(0,n,n)
  for (p in 1:n)
  {
    for (q in 1:n)
    { 
      if (q>p){
        W[p,q] <- rbinom(1,1,omega[p,q])-omega[p,q]
      }
      
    }
  }
  W<-W+t(W)
  AA<-omega-diag(diag(omega))+W
  
  library(igraph)
  g1 <- graph_from_adjacency_matrix(AA,mode = "undirected", weighted = NULL,
                                    diag = TRUE, add.colnames = NULL, add.rownames = NA)
  V(g1)$name <- as.character(1:n)
  #plot(g1)
  
  #components(g1)
  #groups(components(g1))
  
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
  #ex1_out <- SCORE(A0,K)
  #getError(ex1_out,l)
  
  ans = answer(A0,K,l)$error.rate
  result_opca[m]  = ans[1]
  result_npca[m]  = ans[2]
  result_score[m]  = ans[3]
  result_scoreplus[m] = ans[4]
}

c(mean(result_opca),mean(result_npca),mean(result_score),mean(result_scoreplus))
