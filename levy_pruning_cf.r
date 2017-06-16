##########some characteristic exponents#########

get_cf = function(m) {
    if (m=="BM")    return(ebm)
    if (m=="OU")    return(ebm)
    if (m=="EB")    return(ebm)
    if (m=="JN")    return(ejn)
    if (m=="VG")    return(evg)
    if (m=="NIG")   return(enig)
    if (m=="BMJN")  return(ebmjn)
    if (m=="BMVG")  return(ebmvg)
    if (m=="BMNIG") return(ebmnig)
}


####################
# GAUSSIAN PROCESS #
####################

# BM: Brownian motion
ebm = function(k,t,sigma_bm,singularity=FALSE,variance=FALSE) {
	if (singularity) {
		return(c(-Inf,Inf))
	}
	if (variance) {
		return(sigma_bm^2)
	}
	if (sigma_bm < 0) {
		-Inf	
	} else {
		-sigma_bm^2*k^2/2*t
	}
}


####################
# COMPOUND POISSON #
####################

# JN: compound Poisson process w/ normal jumps
ejn = function(k,t,lambda_jn,sigma_jn,singularity=FALSE,variance=FALSE) {
	if (singularity) {
		return(c(-Inf,Inf))
	}
	if (variance) {
		return(lambda_jn*sigma_jn^2)
	}
	if (lambda_jn < 0 || sigma_jn < 0){
		-Inf	
	} else {
		lambda_jn*t*(exp(-sigma_jn^2*k^2/2)-1)
	}
}

# BM+JN: Brownian motion + compound Poisson process w/ normal jumps
ebmjn = function(k,t,sigma_bm,lambda_jn,sigma_jn,singularity=FALSE,variance=FALSE) {
	if (singularity) {
		return(c(-Inf,Inf))
	}
	if (variance) {
		return(sigma_bm^2+lambda_jn*sigma_jn^2)
	}
	if (sigma_bm < 0 || lambda_jn < 0 || sigma_jn < 0){
		-Inf	
	} else {
		-sigma_bm^2*k^2/2*t + lambda_jn*t*(exp(-sigma_jn^2*k^2/2)-1)
	}
}


##################
# VARIANCE GAMMA #
##################

# VG: variance Gamma process
evg = function(k,t,sigma_vg,nu_vg,mu_vg=0,singularity=FALSE,variance=FALSE) {
	#parameterized as in Madan, Carr and Chang (1998)
    # sigma_vg : \sigma, volatility of BM
    # nu_vg    : \nu, variance rate of gamma time change
    # mu_vg    : \theta, directional trend
    # k        : u, the cf variable
	if (variance) {
		return(sigma_vg^2)
	}
	if (singularity) {
		v = sqrt(as.complex(-mu_vg^2-2*sigma_vg^2/nu_vg))
		return(c(1i*mu_vg+v,1i*mu_vg-v)/sigma_vg^2)
	}
	if (nu_vg < 0 || sigma_vg < 0){
		-Inf
	} else {
		-t/nu_vg*log(1-(1i)*mu_vg*nu_vg*k + k^2*sigma_vg^2*nu_vg/2 )
	}
}

# BM+VG: Brownian motion + variance Gamma process
ebmvg = function(k,t,sigma_bm,sigma_vg,nu_vg,mu_vg=0,singularity=FALSE,variance=FALSE) {
	#parameterized as in Madan, Carr and Chang (1998)
    # sigma_vg : \sigma, volatility of BM
    # nu_vg    : \nu, variance rate of gamma time change
    # mu_vg    : \theta, directional trend
    # k        : u, the cf variable
	if (variance) {
		return(sigma_bm^2+sigma_vg^2)
	}
	if (singularity) {
		v = sqrt(as.complex(-mu_vg^2-2*sigma_vg^2/nu_vg))
		return(c(1i*mu_vg+v,1i*mu_vg-v)/sigma_vg^2)
	}
	if (sigma_bm < 0 || nu_vg < 0 || sigma_vg < 0){
		-Inf
	} else {
		 -sigma_bm^2*k^2/2*t - t/nu_vg*log(1-(1i)*mu_vg*nu_vg*k + k^2*sigma_vg^2*nu_vg/2 )
	}
}


###########################
# NORMAL INVERSE GAUSSIAN #
###########################

# NIG: normal inverse Gaussian process
enig = function(k,t,alpha_nig,delta_nig,beta_nig=0,mu_nig=0,singularity=FALSE,variance=FALSE) {
	if (variance) {
		return(t*alpha_nig^2*delta_nig/(alpha_nig^2-beta_nig^2)^(3/2))
	}
	if (singularity) {
		return(c(-Inf,Inf))
	}
	
    if (alpha_nig < abs(beta_nig)) {
		-Inf
	} else if (delta_nig < 0) {
		-Inf
	} else { 
		gamma_nig = sqrt(alpha_nig^2-beta_nig^2)
		t*((1i)*mu_nig*k+delta_nig*(gamma_nig-sqrt(alpha_nig^2-(beta_nig+(1i)*k)^2)))
	}
}

# BM+NIG: Brownian motion + normal inverse Gaussian process
ebmnig = function(k,t,sigma_bm,alpha_nig,delta_nig,beta_nig=0,mu_nig=0,singularity=FALSE,variance=FALSE) {
	if (variance) {
		return(sigma_bm^2 + t*alpha_nig^2*delta_nig/(alpha_nig^2-beta_nig^2)^(3/2))
	}
	if (singularity) {
		return(c(-Inf,Inf))
	}
	if (sigma_bm < 0) {
		-Inf
	} else if (alpha_nig < abs(beta_nig)) {
		-Inf
	} else if (delta_nig < 0) {
		-Inf
	} else { 
		gamma_nig = sqrt(alpha_nig^2-beta_nig^2)
		-sigma_bm^2*k^2/2*t + t*((1i)*mu_nig*k+delta_nig*(gamma_nig-sqrt(alpha_nig^2-(beta_nig+(1i)*k)^2)))
	}
}



