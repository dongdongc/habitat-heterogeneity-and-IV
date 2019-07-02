


{
  {
    
    # packages required#############
    
    
    # install.packages(c("cubature","pracma","lattice","mvtnorm","ggplot2","parallel",
    #                   "sp","foreach","doParallel","gstat","basicTrendline","optimx","numDeriv",
    #               "gplots","vioplot","magrittr","RColorBrewer","fields"))
    
    library(ExtDist) # required for four-parameter beta distribution
    library(cubature) # required for multivariate integration
    library(pracma) # alternative integration algorithms
    library(lattice)
    library(mvtnorm) 
    library(ggplot2)    
    library(parallel)
    library(foreach)
    library(doParallel)
    library(gplots)
    library(vioplot)
    require(magrittr)
    library(RColorBrewer)
    library(fields)
  }
  
  
  par(mar=c(4.5,4.5,4.5,4.5))
  par(mfrow=c(2,2))
  
  
  a<-rep(1:10,10)
  b<-sample(sample(sample(a)))
  random<-matrix(b,10,10)
  
  
  
  
  # dynamics with no variation
  
  comp.sens.no.var<-function(n1,n2,lambda1,lambda2,r1,r2,d1,d2,time.steps){
    
    library(ExtDist) # required for four-parameter beta distribution
    library(cubature) # required for multivariate integration
    library(pracma) # a
    
    generations=time.steps
    N1<-array(0,c(12,12,generations+1))
    N2<-array(0,c(12,12,generations+1))
    L1<-array(0,c(12,12,generations+1))
    L2<-array(0,c(12,12,generations+1))
    
    N1[2:11,2:11,1]<-n1
    N2[2:11,2:11,1]<-n2
    for (t in 1:time.steps){
      
      N1[1,,t]<-N1[11,,t]
      N1[12,,t]<-N1[2,,t]
      N1[,1,t]<-N1[,11,t]
      N1[,12,t]<-N1[,2,t]
      N2[1,,t]<-N2[11,,t]
      N2[12,,t]<-N2[2,,t]
      N2[,1,t]<-N2[,11,t]
      N2[,12,t]<-N2[,2,t]
      
      #######seeds production number
      for (i in 2:11){
        for (j in 2:11) {
          lambda01<-lambda1[i-1,j-1]
          lambda02<-lambda2[i-1,j-1]
          
          L1[i,j,t]=N1[i,j,t]*lambda01/(1+r1*(N1[i,j,t]+N2[i,j,t]))
          L2[i,j,t]=N2[i,j,t]*lambda02/(1+r2*(N2[i,j,t]+N1[i,j,t]))
        }
      }
      
      for (i in 2:11) {
        for (j in 2:11) {
          
          N1[i,j,(t+1)]=((1-4*d1)*L1[i,j,t]+d1*(L1[(i-1),j,t]+L1[i,(j-1),t]+L1[i,(j+1),t]+L1[(i+1),j,t]))
          N2[i,j,(t+1)]=((1-4*d2)*L2[i,j,t]+d2*(L2[(i-1),j,t]+L2[i,(j-1),t]+L2[i,(j+1),t]+L2[(i+1),j,t]))
        }
      }
    }
    N<-array(0,c(10,20,time.steps+1))
    N[1:10,1:10,]<-N1[2:11,2:11,]
    N[1:10,11:20,]<-N2[2:11,2:11,]
    
    
    nn1<-N[1:10,1:10,]
    nn2<-N[1:10,11:20,]
    t<-dim(nn1)[3]
    nn1_t<-c(1:t)
    for (i  in 1:(t)) {
      nn1_t[i]<-mean(nn1[,,i])
    }
    nn2_t<-c(1:t)
    for (i in 1:(t)) {
      nn2_t[i]<-mean(nn2[,,i])
    }
    
    return(cbind(nn1_t,nn2_t))
  }
  
  
  
  
  
  
  
  # functions required for generating dynamics when individuals vary 
  # function describing dynamics when variation between individuals in competitive sensitivity is beta distributed 
  
  comp.sens.beta<-function(n1,n2,lambda1,lambda2,r1,r2,r.var1,r.var2,d1,d2,time.steps){
    
    library(ExtDist) # required for four-parameter beta distribution
    library(cubature) # required for multivariate integration
    library(pracma) 
    
    
    generations=time.steps
    N1<-array(0,c(12,12,generations+1))
    N2<-array(0,c(12,12,generations+1))
    L1<-array(0,c(12,12,generations+1))
    L2<-array(0,c(12,12,generations+1))
    
    N1[2:11,2:11,1]<-n1
    N2[2:11,2:11,1]<-n2
    
    
    r1.min=r1-(min(c(r2,r1))-0.0001)
    r1.max=r1+(min(c(r2,r1))-0.0001)
    r2.min=r2-(min(c(r2,r1))-0.0001)
    r2.max=r2+(min(c(r2,r1))-0.0001)
    
    
    # function to estimate alpha and beta parameters of a four-parameter beta distribution with given mean and variance
    beta.4.par <- function(mean,variance,min,max){
      alpha <- (-(mean-min) * (variance + mean^2 - max*mean - min*mean + min*max))/((max-min)*variance)
      beta <- ((alpha*(max-min))/(mean-min)) - alpha
      if(alpha<=0|beta<=0){
        stop("Distribution not defined because shape parameters are <= 0. Try reducing the required variance, and/or make sure the mean does not equal the maximum or minimum values")
      } else {
        return(c(alpha = alpha, beta = beta))
      }
    }
    
    # shape parameters for beta distributions describing intraspecific variation in both species
    
    sp1.par<-beta.4.par(mean=r1,variance=r.var1,min=r1.min,max=r1.max)
    sp2.par<-beta.4.par(mean=r2,variance=r.var2,min=r2.min,max=r2.max)
    
    
    
    for (t in 1:time.steps){
      #######roll the landscape
      N1[1,,t]<-N1[11,,t]
      N1[12,,t]<-N1[2,,t]
      N1[,1,t]<-N1[,11,t]
      N1[,12,t]<-N1[,2,t]
      
      N2[1,,t]<-N2[11,,t]
      N2[12,,t]<-N2[2,,t]
      N2[,1,t]<-N2[,11,t]
      N2[,12,t]<-N2[,2,t]
      
      for (i in 2:11){
        for (j in 2:11) {
          lambda01<-lambda1[i-1,j-1]
          lambda02<-lambda2[i-1,j-1]
          
          
          sp1.growth.fun <- function(x){
            (1*lambda01/(1+x*(N1[i,j,t]+N2[i,j,t]))) * # the annual plant model
              dBeta_ab(x,shape1 = sp1.par[1],shape2 = sp1.par[2],a = r1.min,b = r1.max) # distribution of competitive sensitivity β分布的概率密度函数
            
          }
          # species 2 growth function (= integrand of equation 3 in the manuscript)
          
          sp2.growth.fun <- function(x){
            (1*lambda02/(1+x*(N1[i,j,t]+N2[i,j,t]))) *
              dBeta_ab(x,shape1 = sp2.par[1],shape2 = sp2.par[2],a = r2.min,b = r2.max)}
          
          # integration gives annual growth rate#
          (sp1.growth<-integrate(sp1.growth.fun,lower=r1.min,upper=r1.max)$value)
          (sp2.growth<-integrate(sp2.growth.fun,lower=r2.min,upper=r2.max)$value)
          
          
          L1[i,j,t]<-N1[i,j,t]*sp1.growth
          L2[i,j,t]<-N2[i,j,t]*sp2.growth
          
        }
      }
      
      for (i in 2:11) {
        for (j in 2:11) {
          
          N1[i,j,(t+1)]=((1-4*d1)*L1[i,j,t]+d1*(L1[(i-1),j,t]+L1[i,(j-1),t]+L1[i,(j+1),t]+L1[(i+1),j,t]))
          N2[i,j,(t+1)]=((1-4*d2)*L2[i,j,t]+d2*(L2[(i-1),j,t]+L2[i,(j-1),t]+L2[i,(j+1),t]+L2[(i+1),j,t]))
        }
      }
    }
    N<-array(0,c(10,20,time.steps+1))
    N[1:10,1:10,]<-N1[2:11,2:11,]
    N[1:10,11:20,]<-N2[2:11,2:11,]
    nn1<-N[1:10,1:10,]
    nn2<-N[1:10,11:20,]
    t<-dim(nn1)[3]
    nn1_t<-c(1:t)
    for (i  in 1:(t)) {
      nn1_t[i]<-mean(nn1[,,i])
    }
    nn2_t<-c(1:t)
    for (i in 1:(t)) {
      nn2_t[i]<-mean(nn2[,,i])
    }
    
    return(cbind(nn1_t,nn2_t))###
    # return(N)
  }
  
  
  noplot<-function(no){
    plot(no[,1],type="l",pch=20,xlab="Time",ylab = "Population density",ylim = c(0,400),lwd=1.5,cex.lab=1)
    #lines(no[,1],lwd=4,xlab=" ",ylab = "",ylim = c(0,400))
    par(new=T)
    plot(no[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,400),col="red",lwd=1.5,cex.lab=1)
    #lines(no[,2],type="l",lwd=4,xlab = "",ylab = "",ylim = c(0,400),col="red")
  }
  
  
  
  co_heatmap<-function(v_d) {
    
    time<-nrow(v_d[[1]])
    len1<-length(v_d)
    len<-sqrt(len1)
    b<-c(1:len1)
    for (i in 1:len1) {
      a<-v_d[[i]][time,]
      {
        if(a[1]<=2 && a[2]<=2)
          b[i]=0
        else if((a[1]<=2 && a[2]>2))
          b[i]=1   
        else if((a[1]>2 && a[2]>2))  
          b[i]=2 
        else if((a[1]>2 && a[2]<=2))  
          b[i]=3   
      }
    }
    
    dim(b)<-c(len,len)
    library("gplots")
    library(RColorBrewer)
    
    cols<-brewer.pal(n=3,name="Paired")
    
    
    b1<-b
    for (i in 1:len) {
      b1[i,]<-b[len+1-i,]
    } 
    colnames(b1)<-seq(0.5,10,by=0.5) 
    rownames(b1)<-seq(10,0.5,by=-0.5)
    heatmap.2(b1,Rowv = FALSE,Colv = FALSE,col=cols,scale = "none",key = TRUE,
              symkey = FALSE,density.info = "none",trace="none",cexRow = 0.8,xlab ="r.var1",ylab = "r.var2" )
    
  }
  
  
  
  
  
  
  
  
  
  
  #parallel parameters########################################################################
  mm<-rep(1:50,50)
  jjj<-function(x){
    b<-0
    for (i in 1:x) {
      a<-rep(i,x)
      b<-c(b,a)
    }
    c<-b[-1]
    return(c)
  }
  nn<-jjj(50)##########
  
}

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################


#generating autocorrelated landscapes###########################################################################################################
{
  
  library(sp)
  library(ExtDist) # required for four-parameter beta distribution
  library(cubature) # required for multivariate integration
  library(pracma) # alternative integration algorithms
  
  library(lattice)
  library(mvtnorm) 
  library(ggplot2)    
  library(parallel)
  library(foreach)
  library(doParallel)
  library("gplots")
  library(vioplot)
  require(magrittr)
  
  library(RColorBrewer)
  library(gstat)
  # create structure
  xy <- expand.grid(1:10, 1:10)
  names(xy) <- c("x","y")
  
  par(mfrow=c(3,4))
  par(mar=c(4,4,4,4))
  
  landscape<-array(NA,c(10,10,10))
  
  dim(landscape)
  
  for (j in 1:10) {
    # define the gstat object (spatial model)
    g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=5.5, model=vgm(psill=0.5,model="Gau",range=3), nmax=4)
    
    # make four simulations based on the stat object
    yy <- predict(g.dummy, newdata=xy, nsim=1)
    
    # show one realization
    gridded(yy) = ~x+y
    #spplot(yy[1])
    
    # show all four simulations:
    #spplot(yy)
    
    
    aa<-yy$sim1
    
    as.vector(aa)
    
    
    rowname<-c(1:100)
    
    rank(aa)
    
    ab<-cbind(rowname,rank(aa),aa)
    ac<-c(1:100)
    
    for (i in 1:100) {
      if (ab[i,2]<=10)
        ac[i]<-1
      else if(ab[i,2]<=20)
        ac[i]<-2
      else if(ab[i,2]<=30)
        ac[i]<-3
      else if(ab[i,2]<=40)
        ac[i]<-4
      else if(ab[i,2]<=50)
        ac[i]<-5
      else if(ab[i,2]<=60)
        ac[i]<-6
      else if(ab[i,2]<=70)
        ac[i]<-7
      else if(ab[i,2]<=80)
        ac[i]<-8
      else if(ab[i,2]<=90)
        ac[i]<-9
      else if(ab[i,2]<=100)
        ac[i]<-10
    }
    
    ad<-matrix(ac,10,10)
    
    landscape[,,j]<-ad
    quilt.plot(Y,X,ad,nrow=10,ncol=10,add=F,col = cols)
    
  }
  
  landscape
}


#initial #################################
initial.beta=175

n1<-matrix(c(rep(initial.beta,100)),10,10)
n2<-matrix(c(rep(initial.beta,100)),10,10)

a<-rep(1:10,10)

b<-sample(a)

random<-matrix(b,10,10)






################################################################################## ############################################################### 
################################################################################## ############################################################### 
################################################################################## ############################################################### 
################################################################################## ############################################################### 
################################################################################## ############################################################### 




#competitive dynamics########################################################################################################## ############################################################### 



{
  #homo,no iv######################################################################################
  l1<-matrix(rep(5.5,100),10,10)
  l2=l1/1.3
  #x11()
  
  no1<-comp.sens.no.var(n1,n2,lambda1 = l1,lambda2 = l2,r1=0.012,r2=0.011,d1=0.01,d2=0.01,time.steps = 1000)
  
  
  #homo, IV######################################################################################
  l1<-matrix(rep(5.5,100),10,10)
  l2=l1/1.3
  
  no2<-comp.sens.beta(n1,n2,lambda1 =l1,lambda2 =l2,r1=0.012,r2=0.011,r.var1 =5e-05,r.var2 = 5e-05,d1=0.01,d2=0.01,time.steps =1000)
  
  
  #hetero, no IV############################################################################################
  l1<-as.matrix(random)
  l2<-as.matrix(random/1.3)  
  
  no3<-comp.sens.no.var(n1,n2,lambda1 = l1,lambda2 = l2,r1=0.012,r2=0.011,d1=0.01,d2=0.01,time.steps = 1000)
  
  
  #hetero, IV####################################################################################
  
  
  no4<-comp.sens.beta(n1,n2,lambda1 =l1,lambda2 =l2,r1=0.012,r2=0.011,r.var1 =5e-05,r.var2 = 5e-05,d1=0.01,d2=0.01,time.steps =1000)
  
}


#example of competitive dynamics in heterogeneous landscape##################################
{
  
  auto<-landscape[,,1]
  class(auto)
  #homo, no IV######################################################################################
  l1<-matrix(rep(5.5,100),10,10)
  l2=l1/1.3
  
  
  no5<-comp.sens.no.var(n1,n2,lambda1 = l1,lambda2 = l2,r1=0.012,r2=0.011,d1=0.01,d2=0.01,time.steps = 1000)
  
  
  #homo, IV#####################################################################################
  l1<-matrix(rep(5.5,100),10,10)
  l2=l1/1.3
  
  no6<-comp.sens.beta(n1,n2,lambda1 =l1,lambda2 =l2,r1=0.012,r2=0.011,r.var1 =5e-05,r.var2 = 5e-05,d1=0.01,d2=0.01,time.steps =1000)
  
  
  #hetero, no IV############################################################################################
  l1<-as.matrix(auto)
  l2<-as.matrix(auto/1.3)  
  
  no7<-comp.sens.no.var(n1,n2,lambda1 = l1,lambda2 = l2,r1=0.012,r2=0.011,d1=0.01,d2=0.01,time.steps = 1000)
  
  
  #hetero, IV###################################################################################
  
  l1<-as.matrix(auto)
  l2<-as.matrix(auto/1.3)
  
  no8<-comp.sens.beta(n1,n2,lambda1 =l1,lambda2 =l2,r1=0.012,r2=0.011,r.var1 =5e-05,r.var2 = 5e-05,d1=0.01,d2=0.01,time.steps =1000)
  
}
{
  
  
  # pdf(file = "figure1bar.pdf",width =7.5 ,height = 7.5)
  
  #  tiff(filename = "fig1.tiff",width =108 ,height = 108,units ="mm",res = 600,compression = 'lzw')
  par(mar=c(2.5,2.5,1,1))
  par(mfrow=c(2,2)) 
  par(mgp=c(1.5,0.6,0))
  
  library(fields)
  auto<-landscape[,,1]
  y <- rep(1:10, each=10)
  x <- rep(1:10, 10)
  cols<-brewer.pal(n=10,name="RdYlBu")
  quilt.plot(x,y,random,nrow=10,ncol=10,xaxt="n",yaxt="n",add=F,main="Random",col = cols,add.legend=FALSE,cex.lab=0.6,cex.axis=0.6,cex.main=0.7)
  
  
  mtext("(a)",side = 3,at=-1,font = 2,cex=0.7)
  
  quilt.plot(x,y,auto,nrow=10,ncol=10,add=F,xaxt="n",yaxt="n",main="Autocorrelated",col = cols,add.legend=FALSE,cex.lab=0.6,cex.axis=0.6,cex.main=0.7)
  
  mtext("(b)",side = 3,at=-1,font = 2,cex=0.7)
  
  
  
  {
    #Homogeneous landscape, no IV
    plot(no1[,1],type="l",pch=20,xlab="Time",ylab = "Population density",ylim = c(0,900),main="Random",lwd=1.5,lty=1,cex.lab=0.6,cex.axis=0.6,cex.main=0.7)
    par(new=T)
    plot(no1[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,900),col="red",lwd=1.5,cex.lab=1,xaxt="n",yaxt="n",lty=1)
    #Homogeneous landscape, IV
    par(new=T)
    plot(no2[,1],type="l",pch=20,xlab="",ylab = "",ylim = c(0,900),main="",lwd=1.5,cex.lab=1,xaxt="n",yaxt="n",lty=2)
    par(new=T)
    plot(no2[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,900),col="red",lwd=1.5,cex.lab=1,xaxt="n",yaxt="n",lty=2)
    #Heterogeneous landscape, no IV
    par(new=T)
    plot(no3[,1],type="l",pch=20,xlab="",ylab = " ",ylim = c(0,900),main="",cex.lab=1,xaxt="n",yaxt="n",lty=3,lwd=3)
    par(new=T)
    plot(no3[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,900),col="red",cex.lab=1,xaxt="n",yaxt="n",lty=3,lwd=3)
    #Heterogeneous landscape, IV
    par(new=T)
    plot(no4[,1],type="l",pch=20,xlab="",ylab = "",ylim = c(0,900),main="",lwd=1.5,cex.lab=1,cex.axis=1,xaxt="n",yaxt="n",lty=4)
    par(new=T)
    plot(no4[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,900),col="red",lwd=1.5,cex.lab=1,cex.axis=1,xaxt="n",yaxt="n",lty=4)
    
    legend("topleft",legend = c(expression(paste(italic(sigma[h]^2==0),",   ",italic(sigma[r]^2==0))),
                                expression(paste(italic(sigma[h]^2==0),",   ",sigma[italic(r)]^2==5%*%10^-5)),
                                expression(paste(italic(sigma[h]^2==8.7),", ",italic(sigma[r]^2==0))),
                                expression(paste(italic(sigma[h]^2==8.7),", ",sigma[italic(r)]^2==5%*%10^-5))),
           lty = c(1,2,3,4),lwd = c(1.5,1.5,3,1.5),bty = "n",col = c("black","black","black","black"),cex = 0.6)
    mtext("(c)",side = 3,at=-200,font = 2,cex=0.7)
  }
  
  {
    #Homogeneous landscape, no IV
    plot(no5[,1],type="l",pch=20,xlab="Time",ylab = "Population density",ylim = c(0,900),main="Autocorrelated",lwd=1.5,lty=1,cex.lab=0.6,cex.axis=0.6,cex.main=0.7)
    par(new=T)
    plot(no5[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,900),col="red",lwd=1.5,cex.lab=1,xaxt="n",yaxt="n",lty=1)
    #Homogeneous landscape, IV
    par(new=T)
    plot(no6[,1],type="l",pch=20,xlab="",ylab = "",ylim = c(0,900),main="",lwd=1.5,cex.lab=1,xaxt="n",yaxt="n",lty=2)
    par(new=T)
    plot(no6[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,900),col="red",lwd=1.5,cex.lab=1,xaxt="n",yaxt="n",lty=2)
    #Heterogeneous landscape, no IV
    par(new=T)
    plot(no7[,1],type="l",pch=20,xlab="",ylab = "",ylim = c(0,900),main="",cex.lab=1,lty=3,xaxt="n",yaxt="n",lwd=3)
    par(new=T)
    plot(no7[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,900),col="red",cex.lab=1,xaxt="n",yaxt="n",lty=3,lwd=3)
    #Heterogeneous landscape, IV
    par(new=T)
    plot(no8[,1],type="l",pch=20,xlab="",ylab = "",ylim = c(0,900),main="",lwd=1.5,cex.lab=1,cex.axis=1,xaxt="n",yaxt="n",lty=4)
    par(new=T)
    plot(no8[,2],type="l",pch=20,xlab = "",ylab = "",ylim = c(0,900),col="red",lwd=1.5,cex.lab=1,cex.axis=1,xaxt="n",yaxt="n",lty=4)
    
    #legend("topleft",legend = c("Homo + no IV","Homo + IV","Heter + no IV", "Heter + IV"),lty = c(1,2,3,4),lwd = c(2,2,4,2),bty = "n",col = c("black","black","black","black"))
    #legend("topleft",legend = c("Homo + no IV","Homo + IV","Heter + no IV", "Heter + IV"),lty = c(1,2,3,4),lwd = c(2,2,4,2),bty = "n",col = c("black","black","black","black"))
    mtext("(d)",side = 3,at=-200,font = 2,cex=0.7)
  }
  # dev.off()
}
}
  
  
  #dispersal rate and final population density########################################
  dispersal_density<-matrix(NA,10,82)
  
  
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  for (i in 1:10) {
    print(i)
    
    
    mvd02<-foreach(ii=0:40,.combine = cbind)  %dopar% comp.sens.beta(n1,n2,lambda1 =l1,lambda2 =l2,r1=0.012,r2=0.011,r.var1 =5e-05,r.var2 = 5e-05,d1=0.005*ii,d2=0.005*ii,time.steps =1000)
    #write.csv(mvd2,file = "b2.csv")
    
    dispersal_density[i,]<-mvd02[1000,]
    
    
  }
  stopCluster(cl)
  
  
  
  
  
  #effects of IV between spcecies##################################################################################################
  
  co_heatmap<-function(v_d) {
    
    time<-nrow(v_d[[1]])
    len1<-length(v_d)
    len<-sqrt(len1)
    b<-c(1:len1)
    for (i in 1:len1) {
      a<-v_d[[i]][time,]
      {
        if(a[1]<=2 && a[2]<=2)
          b[i]=0
        else if((a[1]<=2 && a[2]>2))
          b[i]=1   ####
        else if((a[1]>2 && a[2]>2))  
          b[i]=2   ##########
        else if((a[1]>2 && a[2]<=2))  
          b[i]=3   ###
      }
    }
    
    dim(b)<-c(len,len)
    library("gplots")
    library(RColorBrewer)
    
    cols<-brewer.pal(n=3,name="Paired")
    
    
    b1<-b
    for (i in 1:len) {
      b1[i,]<-b[len+1-i,]
    } 
    colnames(b1)<-seq(0.5,10,by=0.5) 
    rownames(b1)<-seq(10,0.5,by=-0.5)
    heatmap.2(b1,Rowv = FALSE,Colv = FALSE,col=cols,scale = "none",key = FALSE,
              symkey = FALSE,density.info = "none",trace="none",cexRow = 1,xlab =expression("Species 1 variance in r  " (x10^-5)),ylab = expression("Species 2 variance in r  "(x10^-5)) )
    return(b1)
  }
  
  
  
  
  
  
  {#random landscape###################################################################################################
    randomIV<-array(NA,c(50,50,50))
    
    for (i in 1:50) {
      
      l1<-as.matrix(random)
      l2<-as.matrix(random/1.3)  
      
      no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      registerDoParallel(cl)
      
      mvd22<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% comp.sens.beta(n1,n2,lambda1 =l1,lambda2 =l2,r1=0.012,r2=0.011,r.var1 =(0.2e-05)*ii,r.var2 = (0.2e-05)*jj,d1=0.01,d2=0.01,time.steps =1000)
      #write.csv(mvd2,file = "b2.csv")
      
      stopCluster(cl)
      
      #png(filename = "variance random.png",width = 150,height = 150,units ="mm",res = 400)
      res<-co_heatmap(mvd22)
      randomIV[,,i]<-res
      
      
    }
    #dev.off()
    write.csv(randomIV,file = "1randomIV.csv")
  }
  
  
  
  {##autocorrelated######################################################################################################
    autoIV<-array(NA,c(50,50,50))
    
    for (i in 1:50) {
      auto<-landscape[,,i]
      l1<-as.matrix(auto)
      l2<-as.matrix(auto/1.3) 
      
      no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      registerDoParallel(cl)
      
      mvd2<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% comp.sens.beta(n1,n2,lambda1 =l1,lambda2 =l2,r1=0.012,r2=0.011,r.var1 =(0.2e-05)*ii,r.var2 = (0.2e-05)*jj,d1=0.01,d2=0.01,time.steps =1000)
      #write.csv(mvd2,file = "b2.csv")
      
      stopCluster(cl)
      
      #png(filename = "variance auto.png",width = 150,height = 150,units ="mm",res = 400)
      res1<-co_heatmap(mvd2)
      autoIV[,,i]<-res1
      
    }
    
    write.csv(autoIV,file = "1autoIV.csv")
    #dev.off()
  }
  
  
  #dispersal rate varies between species######################################################
  
  
  co_heatmap<-function(v_d) {
    
    time<-nrow(v_d[[1]])
    len1<-length(v_d)
    len<-sqrt(len1)
    b<-c(1:len1)
    for (i in 1:len1) {
      a<-v_d[[i]][time,]
      {
        if(a[1]<=2 && a[2]<=2)
          b[i]=0
        else if((a[1]<=2 && a[2]>2))
          b[i]=1  
        else if((a[1]>2 && a[2]>2))  
          b[i]=2  
        else if((a[1]>2 && a[2]<=2))  
          b[i]=3  
      }
    }
    
    dim(b)<-c(len,len)
    library("gplots")
    library(RColorBrewer)
    
    cols<-brewer.pal(n=3,name="Paired")
    
    
    b1<-b
    for (i in 1:len) {
      b1[i,]<-b[len+1-i,]
    } 
    colnames(b1)<-seq(0.01,0.2,by=0.01) 
    rownames(b1)<-seq(0.2,0.01,by=-0.01) 
    heatmap.2(b1,Rowv = FALSE,Colv = FALSE,col=cols,scale = "none",key = FALSE,
              symkey = FALSE,density.info = "none",trace="none",cexRow = 1,xlab ="d1",ylab = "d2" )
    return(b1)
  }
  
  
  {#random landscape############################################################
    randomdispersal<-array(NA,c(50,50,10))
    
    for (i in 1:10) {
      
      
      
      random1<-as.matrix(random)
      random2<-as.matrix(random1/1.3) 
      
      no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      registerDoParallel(cl)
      
      mvd31<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% comp.sens.beta(n1,n2,lambda1 =random1,lambda2 =random2,r1=0.012,r2=0.011,r.var1 =5e-05,r.var2 = 5e-05,d1=0.004*ii,d2=0.004*jj,time.steps =1000)
      #write.csv(mvd31,file = "mvd31.csv")
      
      stopCluster(cl)
      
      
      res<-co_heatmap(mvd31)
      
      randomdispersal[,,i]<-res
      
    }
    write.csv(randomdispersal,file = "1randomdispersal.csv")
  }
  
  {
    #autocorrelated landscape######################################################################################
    autodispersal<-array(NA,c(50,50,10))
    
    for (i in 1:10) {
      auto<-landscape[,,i]
      auto1<-as.matrix(auto)
      auto2<-as.matrix(auto/1.3) 
      
      no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      registerDoParallel(cl)
      
      mvd32<-foreach(ii=nn,jj=mm,.multicombine = TRUE)  %dopar% comp.sens.beta(n1,n2,lambda1 =auto1,lambda2 =auto2,r1=0.012,r2=0.011,r.var1 =5e-05,r.var2 = 5e-05,d1=0.004*ii,d2=0.004*jj,time.steps =1000)
      write.csv(mvd31,file = "mvd32.csv")
      
      stopCluster(cl)
      
      
      res<-co_heatmap(mvd31)
      
      autodispersal[,,i]<-res
      
    }
    write.csv(autodispersal,file = "1autodispersal.csv")
    
  }
  
  
  
  
  #################################################################
  #effect of variaiton in patch quality on competitive dynamics##########################
  
  
  
  #gamma distribution #########################################
  lambdagamma<-function(mean.lambda,var.lambda){
    # function for deriving shape and rate parameters of a gamma distribution given known mean and variance
    gamma.par<-function(mean,variance){
      shape<-(mean^2)/variance
      rate<-shape/mean
      return(c(shape,rate))
    }
    if(var.lambda>0){
      (shape1<-gamma.par(mean=mean.lambda,variance=var.lambda)[1]) 
      (shape2<-gamma.par(mean=mean.lambda,variance=var.lambda)[2])
    } else {shape1<-0;shape2<-0}
    if (shape1==0){individual.lambdas <- rep(mean.lambda,100)}
    else{individual.lambdas <- rgamma(100,shape = shape1,rate = shape2)}
    
    return(individual.lambdas)
  }
  lambdagamma(5.5,5)
  
  #################################################################################################
  co_heter<-function(var){
    a<-lambdagamma(5.5,var)
    b<-sample(sample(sample(a)))
    random<-matrix(b,10,10)
    l1<-as.matrix(random)
    l2<-as.matrix(random/1.3)  
    initial.beta=100
    n1<-matrix(c(rep(initial.beta,100)),10,10)
    n2<-n1
    ##异质生境、个体有差异
    no<-comp.sens.beta(n1,n2,lambda1 =l1,lambda2 =l2,r1=0.012,r2=0.011,r.var1 =5e-05,r.var2 = 5e-05,d1=0.01,d2=0.01,time.steps =1000)
    return(no[1000,])
  }
  
  co_heter(5)
  
  
  {
    # reptites#####################################################################
    res<-matrix(NA,200,50)
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    for (i in 1:50) {
      co_heter1<-foreach(i1=1:100,.combine = rbind)  %dopar% co_heter(var = (i1/10))
      res[1:100,i]<-co_heter1[,1]
      res[101:200,i]<-co_heter1[,2]
    }
    
    vargamma<-seq(0.1,10,0.1)
    stopCluster(cl)
    res12<-apply(res, 1, mean)
    
    
    #png(filename = "heterogeneity_density.png",width = 120,height = 120,units ="mm",res = 400)
    par(mar=c(5,5,5,5))
    
    par(mfrow=c(1,1))
    plot(vargamma,res12[1:100],xlab = "Heterogeneity of patch quality",ylab = "Final population density",pch=20,ylim = c(0,600),xlim = c(0,10))
    par(new=T)
    plot(vargamma,res12[101:200],xlab = " ",ylab = " ",pch=20,ylim = c(0,600),col="red",xaxt="n",yaxt="n",xlim = c(0,10))
    
    legend("topright",legend = c("Species 1","Species 2"),lty = c(0,0),pch=c(20,20),bty = "n",col = c("black","red"))
  }
  
  heterogeneity<-cbind(res12[1:100],res12[101:200])
  write.csv(heterogeneity,file = "heterogeneity.csv")
  
  
  
  
  
  