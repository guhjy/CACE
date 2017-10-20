

get_probs<-function(z,zRange,mat){
  out = c(rep(NA,length(z)))
  for(i in 1:length(z)){
    r = rank(append(z[i],zRange),ties.method = 'last')[1]
    if(r>length(zRange)){out[i]=mat[i,(r-1)]}
    else if(r==1){out[i]=mat[i,1]}
    else{
      fl_val = mat[i,(r-1)]
      cl_val = mat[i,r]
      out[i] = mean(c(fl_val,cl_val))
    }
  }
  return(out)
}




