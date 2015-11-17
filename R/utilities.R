
listadd = function(a,b){
   b = list(b)
   if (is.null(a)) return(b)
   return(c(a,b))
}

list2matrix = function(vlist,bNAasZero=TRUE) {
   if(class(vlist)=="list"){
      maxlen = max(sapply(vlist,length))
      vlist = lapply(vlist,function(x){extend(x,maxlen,bNAasZero)})
      res = matrix(0,maxlen,length(vlist))
      for(i in 1:length(vlist)) {
         res[,i] = vlist[[i]]
      }
      colnames(res) = names(vlist)
      rownames(res) = 1:nrow(res)
      return(res)
   }else{
      return(vlist)
   }
}

extend = function(vec,len,bNAasZero=TRUE) {
   if(bNAasZero){
      vec[is.na(vec)]=0
   }
   vlen = length(vec)
   if(vlen < len) {
      toadd = rep(0,len-vlen)
      vec = c(vec,toadd)
   }
   return(vec)
}

QCminimize = function(obj) {
   obj@annotation = list(version=QCannotation(obj))
   obj@DBA$binding=NULL
   obj@DBA$merged=NULL
   obj@DBA$called=NULL
   
   gc()
   res = obj
   return(res)
}


