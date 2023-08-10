generate.valueMatrix<-function(ds,matrix,pars,time){
	
	# ds=data_col.ts
	# matrix=SAM
	# year=2010
	
	results=matrix(NA,nrow=nrow(matrix),ncol=ncol(matrix),dimnames = dimnames(matrix))
	period=which(ds$time==time)
	attach(ds[period,])
	attach(pars)
	for(r in 1:nrow(matrix)){
		for(c in 1:ncol(matrix)){
			cellVal<-as.character(matrix[r,c])
			if(nchar(cellVal)>0){
				#If there's an entry in the cell
				eq<-eval(parse(text=cellVal))
				results[r,c]=eq
			}
		}
	}
	detach(pars)
	detach(ds[period,])
	return(results)
}

mymatplot<-function(dataset,varnames,location){
	dstemp=sapply(varnames,function(x) return(eval(parse(text=x),envir=dataset)))
	matplot(1:nrow(dstemp),dstemp, main="" ,type="l", ylab="", xlab="", lwd=2,lty=1)
	legend(location,legend=varnames,lty=1,lwd=2,col=1:length(varnames),bty='n')
}

mymatplotcompare<-function(datasets,varnames,location){
	dstemp=sapply(varnames,function(x) return(eval(parse(text=x),envir=datasets[[1]])))
	ltys=rep(1,length(varnames))
	namesScen<-names(datasets)
	if(length(datasets)>1){
		for(i in 2:length(datasets)){
			dstemp=cbind(dstemp,sapply(varnames,function(x) return(eval(parse(text=x),envir=datasets[[i]]))))
			ltys=c(ltys,rep(i,length(varnames)))
		}
	}
	matplot(1:nrow(dstemp),dstemp, main="" ,type="l", ylab="", xlab="", lwd=2,col=1:length(varnames),lty=ltys)
	legend(location,legend=c(varnames,namesScen),lty=c(rep(1,length(varnames)),seq(1,length(namesScen))),lwd=2,col=c(1:length(varnames),rep(1,length(namesScen))),bty='n')
}

growth<-function(var){
	return(100*(var[-1]/var[-length(var)]-1))
}