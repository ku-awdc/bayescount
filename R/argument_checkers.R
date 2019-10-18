checksingleprob <- function(x, na.ok=FALSE){
	if(length(x)!=1){
		stop(paste(substitute(x), 'must be of length 1'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(is.na(x)){
		if(!na.ok){
			stop(paste(substitute(x), 'must be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.numeric(x)){
			stop(paste(substitute(x), 'must be numeric'), call.=FALSE)
		}
		if(! x>=0){
			stop(paste(substitute(x), 'must be >=0'), call.=FALSE)
		}
		if(! x<=1){
			stop(paste(substitute(x), 'must be <=1'), call.=FALSE)
		}
	}
}

checksinglelogical <- function(x, na.ok=FALSE){
	if(length(x)!=1){
		stop(paste(substitute(x), 'must be of length 1'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(is.na(x)){
		if(!na.ok){
			stop(paste(substitute(x), 'must be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.logical(x)){
			stop(paste(substitute(x), 'must be logical'), call.=FALSE)
		}
		if(! x>=0){
			stop(paste(substitute(x), 'must be >=0'), call.=FALSE)
		}
		if(! x<=1){
			stop(paste(substitute(x), 'must be <=1'), call.=FALSE)
		}
	}
}

checksingleposint <- function(x, na.ok=FALSE){
	if(length(x)!=1){
		stop(paste(substitute(x), 'must be of length 1'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(is.na(x)){
		if(!na.ok){
			stop(paste(substitute(x), 'must be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.numeric(x)){
			stop(paste(substitute(x), 'must be numeric'), call.=FALSE)
		}
		if(! x>0){
			stop(paste(substitute(x), 'must be >0'), call.=FALSE)
		}
		if(x%%1 !=0){
			stop(paste(substitute(x), 'must be an integer'), call.=FALSE)
		}
	}
}

checksingleposdouble <- function(x, na.ok=FALSE){
	if(length(x)!=1){
		stop(paste(substitute(x), 'must be of length 1'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(is.na(x)){
		if(!na.ok){
			stop(paste(substitute(x), 'must be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.numeric(x)){
			stop(paste(substitute(x), 'must be numeric'), call.=FALSE)
		}
		if(! x>0){
			stop(paste(substitute(x), 'must be >0'), call.=FALSE)
		}
	}
}

checklen2posdouble <- function(x, na.ok=FALSE){
	if(length(x)!=2){
		stop(paste(substitute(x), 'must be of length 2'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(any(is.na(x))){
		if(!na.ok){
			stop(paste(substitute(x), 'must all be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.numeric(x)){
			stop(paste(substitute(x), 'must be numeric'), call.=FALSE)
		}
		if(! all(x>0)){
			stop(paste(substitute(x), 'must all be >0'), call.=FALSE)
		}
	}
}

checkintlimit <- function(x){
	if(x > .Machine$integer.max || x < -.Machine$integer.max){
		stop(paste(substitute(x), 'is beyond the integer limits of your machine so cannot be passed to the underlying C++ code'), call.=FALSE)
	}
	return(as.integer(x))
}
