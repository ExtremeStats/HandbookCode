

#all data must be in row form, that is each row corresponds to one covariate. 

frechet2<-function(x){
  initdat<-matrix(0, nrow=length((x[,1])), ncol=length(x[1,]))
  l<-length(x[1,])
  for (i in 1:length(x[,1])){            
    distf<-ecdf(x[i,])
    initdat[i,]<-(-log(pmin(l/(l+1)*(distf)(x[i,]))))^(-1/2)
  }
  return(initdat)
}

### angles will compute the observations with the largest euclidead norm, or radial components.

angles<-function(x, thresh_k){
  
  #y<-frechet2(x) ### if you haven’t already standardized to Frechet(2) rv.
  y<-(x)
  csums<-apply(rbind(y,0),2,function(x)sum(x^2))
  ind<-which(csums >=sort(csums,T)[thresh_k])
  z<-(apply(rbind(y,0)[,ind],2,function(x)x^2/sum(x^2)))
  return(z)
}

###computes scaling

scaling<-function(x,ind, thresh_k){
  
  #z<-frechet2(x) ### if you haven’t already standardized to Frechet(2) rv.
  z<-x
  y<-angles(z, thresh_k)[ind,]
  m<-ifelse(min(ind)<0, length(y[,1])-1, length(ind))
  p<-vector()
  p<-length(x[,1])*mean(apply(rbind(y,-1), 2, function(x)max(x)))
  return(p)
}


### scaling step towards the goal of obtaining the initial nodes “at once”, via scaling by the scalar a, according to Algo. 4.

initialnodesa<-function(z, a, thresh_k){
  
  v<-matrix(nrow=length(z[,2]), ncol=length(z[,2]))
  q<-matrix(nrow=length(z[,2]), ncol=length(z[,2]))
  
  #z<-frechet2(z) ##### if you haven’t already standardized to Frechet(2) rv.
  
  for (i in 1:length(z[,2])){
    for (j in 1:length(z[,2])){
      v[j,i]<-scaling(rbind((z)[c(i,j),]),c(1,2), thresh_k)
      q[j,i]<-(1+a^2)/2*scaling(c(a,1)*rbind((z)[c(i,j),]), c(1,2), thresh_k)-v[j,i]-(a^2-1)
    }
  }
  return(q)
}


### scaling step towards finding one descendant at a time according to Algo. 5.
###found- all found indices taken from the initial ordering of x !!
### if X is d=10 dimensional, and you have found the first 2 nodes, say  5 and 6 (from the initial ordering of X), then tofind must be set to 8; that is (d-number of found nodes).

descnodes<-function(x, found, tofind, a, thresh_k){
  
  v<-vector()
  len<-length(x[,1])
  lefound<-length(found)
  
  for (i in c(1:length(x[,1]))[-found]){
    v[i]<-(a^2*(lefound+1)+tofind-1)/(tofind+lefound)*scaling(rbind((rbind(a*x[(i),],a*x[found,],x[-c(i,found),]))),c(1:length((rbind(x[(i),],x[found,],x[-c(i,found),]))[,1])), thresh_k)-scaling(rbind((x)),c(1:len), thresh_k)-(a^2-1)*scaling(rbind(x[(i),],x[found,]),c(1:length((rbind(x[(i),],x[found,]))[,1])), thresh_k)
  }
  return(v)
}


### Below both algorithms above are iteratively combined to return a causal order. epsilon_1, and epsilon_2 are the same as in Algo. 5


causalorder<-function(x, thresh_k, epsilon_1, epsilon_2, a){
  x<-rbind(x)
  #x<-frechet2(x) ##### if you haven’t already standardized to Frechet(2) rv.
  
  d<-length(x[,1])
  nv<-rep(0,d)
  
  #initial nodes found according to Algo. 4.
  
  
  pairs<-initialnodesa(x, a,thresh_k)
  sc<-apply(pairs,2,function(x)min(x))
  
  indd<-which(sc== max(sc));
  
  #nv[(d-length(indd)+1):d]<- -scalingsq[1,indd]
  nv[(d-length(indd)+1):d]<- indd
  
  #descendants, found iteratively by applying the descnodes function, corresponding to Algo. 5.
  
  while(length(which(nv==0))>=1){
    scdif<-(descnodes(x, nv[which(nv>0)], (d-length(which(nv>0))), a, thresh_k))
    
    scdif[which(is.na(scdif))]<- -10^10  ####just a very small value, which will not be reached by the differences.
    
    #### this will pick the argument with the biggest positive difference… it may happen, very rarely, that you might get an error, as this may return two indices, but as long as there is an ordering in the end, it should be fine.
    
    nv[(d-length(which(nv>0)))]<-which(scdif==max(scdif)) 
    
  }
  return(nv)
}


### the x argument below (standardized to Frechet2) must be reordered according to its causal order. B is computed according to Algo. 1.

Bmatrix<-function(x, thresh_k){
  
  k<-length(x[,1])
  b=matrix(0, nrow= k, ncol= k)
  for (i in k:2){
    b[i,i]<-ifelse(i<=(k-1), scaling(x,c(i:k), thresh_k),scaling(x, k, thresh_k))-ifelse(i<=(k-1), ifelse(i<=(k-2), scaling(x, c((i+1):k), thresh_k), scaling(x, k, thresh_k)),0)
  }
  b[1,1]<-(scaling(x,c(1:k), thresh_k)-sum(diag(b)))
  for (i in 1:(k-1)){
    for (j in (i+1):(k)){
      if(i==1)
        b[i,j]<-((scaling(x,-c((i+1):j), thresh_k)-sum(b[i,c(i:(j-1))])-ifelse(j<k,scaling(x,-c(i:(j)), thresh_k),0)))
      else  
        b[i,j]<-((scaling(x,-c((1:(i-1)),((i+1):j)), thresh_k)-sum(b[i,c(i:(j-1))])-ifelse(j<k,scaling(x,-c(1:j), thresh_k),0)))
    } 
  }
  b<- pmax(b,0)^.5
  
  return(b)
}


#### working example with j2dn with 21 variables #####

###  !!!   compute ordering for DATA  (say columns correspond to variables, since this is typical !!)
# (You might want to check different values of the eps’ terms, for instance, say, if eps_2=.00005 doesn’t, you might select a larger value, i.e.  eps_2=.000075, and so on. )

# orderj2dn<-causalorder(frechet2(t(DATA)), 139, .00001, .00007, 1.0001)
                       
#compute ML coefficient matrix, and standardize rows
                       
# MLmatData<-t(apply(Bmatrix(frechet2(t(DATA)[orderData,]),139),1,function(x)x/sum(x^2)^.5))
                       
                       
                       
                       
                       