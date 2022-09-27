#We start with a bunch of facial photographs, 5 photos of 10 persons each. 
#All the images are 200x180 in size. So they are 36000 dimensional 
#vectors. Our aim is to reduce the dimension as much as possible without 
#losing our ability to recognise the faces.

#We start with loading the images.

library(jpeg)
#------
loadImages = function() {
  name = c("male.pacole","male.pspliu","male.sjbeck","male.skumar","male.rsanti",
           "female.anpage","female.asamma","female.klclar","female.ekavaz","female.drbost")
  #------  
  x = matrix(0,nrow=10*5,ncol=200*180)
  #------
  k = 0
  
  for(i in 1:10) {
    for(j in 1:5) {
      k = k + 1
      #The first argument below should be the location of the face folder.
      filename = paste("E:/ISI Classes/Semester II/Arnab Sir_Vector and Matrices II/Eigen face project/grayfaces/",name[i],".",j,".jpg",sep="")
      x[k,] = as.vector(readJPEG(filename))
    }
  }
  #------
  return(x)
}
#---------------
#Next, we perform PCA. The following function does precisely this, but without using R's built-in tools like
#prcomp or princomp. We perform the operations explicitly so that you can appreciate what goes on
#"behind the scene".

process = function(x) {
  meanx = apply(x,2,mean)
  
  y = scale(x,scale=F) #Rows of y are the cases
  A = y %*% t(y)
  eig = eigen(A)    
  
  P = t(y) %*% eig$vec[,-50] #Columns of P are e vectors of Y'Y
  
  Q = apply(P,2,function(x) x/sqrt(sum(x*x))) #Columns of Q form onb for rowspace of Y
  
  
  scores = y %*% Q
  #scores[i,] is for i-th image
  return(list(centre=meanx,onb=Q,scores=scores,values=eig$values))
}
#------------------
#Each principal component is again a 200*180 dimensional vector, and hence may
#be considered as an image. It might be instructive to take a look at these.
#The following function helps you to just that.

showFace = function(newCoord, i) {
  plot(1:2,ty='n',main="0")
  y = abs(newCoord$onb[,i])
  extreme = range(y)
  y = (y-extreme[1])/(extreme[2]-extreme[1])
  dim(y) = c(200,180)
  rasterImage(as.raster(y),1,1,2,2)
}
#----------------
#The next function reconstucts the faces from the principal components.
#It starts with the "mean face" and then gradually adds the details.
#You'll need to hit "enter" to step through the process.

showSteps = function(newCoord,i) {
  meanx = newCoord$centre
  Q = newCoord$onb
  scores = newCoord$scores
  values = newCoord$values
  expl = 100*cumsum(newCoord$values)/sum(newCoord$values)
  
  coeff = as.vector(scores[i,])
  
  plot(1:2,ty='n',main="0")
  y = meanx
  dim(y) = c(200,180)
  plot(as.raster(y))
  readline()
  for(k in 1:49) {
    if(k==1)
      temp = Q[,1]*coeff[1]
    else
      temp=Q[,1:k] %*% as.vector(coeff[1:k])
    recons = meanx + temp
    recons[recons<0]=0
    recons[recons>1]=1
    dim(recons) = c(200,180)
    plot(as.raster(recons),main=paste(k,":",values[k],", ",expl[k]))
    readline()
  }
}

filename = "male.pacole.1.jpg" #Use appropriate location of file in device with filename here
tmp=readJPEG(filename)

#rasterImage(tmp[90:130,70:115],-3,0,3,25)

x = loadImages()
newcoord = process(x)
#showSteps(newcoord,26)

mat_score=newcoord$scores[,-11:-49]

#Within matrix formation
W_list <- list()
for (i in 1:10)
{
  W_i=cov(mat_score[((5*(i-1))+1):((5*(i-1))+5),])#cov of the i-th 5 points (cluster i)
  W_list[[i]]=W_i
}
W=(5-1)*Reduce('+',W_list)/(50-1)#combined within matrix

#Between matrix formation
B1=apply(mat_score[1:5,],2,mean)#Mean of the first cluster
for (i in 2:10)
{
  B_i=apply(mat_score[((5*(i-1))+1):((5*(i-1))+5),],2,mean)#Mean of the i-th cluster
  B1=rbind(B1,B_i)
}
B=cov(B1)#Covariance matrix when all points in a cluster equals the cluster mean

C=solve(W)%*%B
main_vec=eigen(C)$vectors[,1:2]
comp_score=mat_score%*%main_vec

plot(comp_score[,1])
plot(comp_score[,2])
plot(comp_score[,1], comp_score[,2])

range_list=list()
for (i in 1:10)
{
  range1=range(comp_score[((5*(i-1))+1):((5*(i-1))+5),1])
  range2=range(comp_score[((5*(i-1))+1):((5*(i-1))+5),2])
  range=list(range1,range2)
  range_list[[i]]=range
}

plot(comp_score[,1])
plot(comp_score[,2])
plot(comp_score[,1], comp_score[,2]) #plotting different clusters