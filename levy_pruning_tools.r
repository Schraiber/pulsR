library(ape)
library(geiger)

#####################
# solve for moments #
#####################

get_params_for_var = function(process_var=1,process_kurt=0,frac_of_var=0,halflife=0,decay=0,tip=0) {
    ret = list()
    ret$bm = list(sigma.bm=sqrt(process_var))
    ret$ou = get_params_OU(process_var,halflife)
    ret$eb = get_params_EB(process_var,decay)
    ret$jn = get_params_JN(process_var,process_kurt,1)
    ret$vg = get_params_VG(process_var,process_kurt,1)
    ret$nig = get_params_NIG(process_var,process_kurt,1)
    ret$bmjn = get_params_JN(process_var,process_kurt,frac_of_var)
    ret$bmvg = get_params_VG(process_var,process_kurt,frac_of_var)
    ret$bmnig = get_params_NIG(process_var,process_kurt,frac_of_var)
    ret$tip = tip
    return(ret)
}

get_params_JN = function(process_var,process_kurt,frac_of_var) {
    lambda = 3*frac_of_var^2/process_kurt
    delta2 = (process_kurt*process_var)/(3*frac_of_var)
    sigma2 = process_var-frac_of_var*process_var
    return(list(lambda.jn=lambda, delta.jn=sqrt(delta2), sigma.bm=sqrt(sigma2)))
}

get_params_VG = function(process_var,process_kurt,frac_of_var) {
    nu = process_kurt/(3*frac_of_var^2)
    tau2 = frac_of_var*process_var
    sigma2 = process_var-frac_of_var*process_var
    return(list(nu.vg=nu, sigma.vg=sqrt(tau2), sigma.bm=sqrt(sigma2)))
}

get_params_OU = function(process_var,halflife) {
    return(list(sigma.bm=sqrt(2*process_var*log(2)/halflife),alpha.ou=log(2)/halflife))
}

get_params_EB = function(process_var,decay) {
    # stationary variance of process is zero as t->Inf
    return(list(sigma.bm=sqrt(process_var),
                decay.eb=decay))
}

get_params_NIG = function(process_var,process_kurt,frac_of_var) {
	alpha = sqrt(3/(frac_of_var*process_kurt*process_var))
	delta = sqrt(3*frac_of_var*process_var/process_kurt)
	sigma2 = process_var - frac_of_var*process_var
	return(list(sigma.bm = sqrt(sigma2), alpha.nig = alpha, delta.nig = delta))
}

###################
# compute moments #
###################

get_moments_JN = function(lambda.jn, delta.jn, sigma.bm) {
	v = lambda.jn*delta.jn^2 + sigma.bm^2
	k = 3*lambda.jn*delta.jn^4/(sigma.bm^2+lambda.jn*delta.jn^2)^2
	return(list(var=v,kurt=k))
}

get_moments_VG = function(nu.vg, sigma.vg, sigma.bm) {
	v = sigma.bm^2+sigma.vg^2
	k = 3*nu.vg
	return(list(var = v, kurt = k))
}

get_moments_NIG = function(alpha.nig, delta.nig, sigma.bm) {
	v = sigma.bm^2 + delta.nig/alpha.nig
	k = 3*delta.nig/(alpha.nig*(sigma.bm^2*alpha.nig+delta.nig)^2)
	return(list(var = v, kurt = k))
}

#################
# Levy measures #
#################

JN_measure = function(x, lambda, delta, log = FALSE) {
	ret = log(lambda) + dnorm(x,sd=delta,log=TRUE)
	if (log) {
		return(ret)
	} else {
		return(exp(ret))
	}
}

VG_measure = function(x, nu, sigma, mu = 0 ,log=FALSE) {
	ret = -log(nu*abs(x)) + mu*x/sigma^2 - sqrt(2/nu-mu^2/sigma^2)/sigma*abs(x)
	if (log) {
		return(ret)
	} else {
		return(exp(ret))
	}
} 
	
NIG_measure = function(x, alpha, delta, beta = 0, log=FALSE) {
	ret = log(alpha) + log(delta) - log(pi) - log(abs(x)) + beta*x + log(besselK(alpha*abs(x),1))
	if (log) {
		return(ret)
	} else {
		return(exp(ret))
	}
}

#################
# rate of jumps #
#################

jump_rate = function(x, levy_measure, ...) {
    if (x == 0) {
	    total = integrate(levy_measure, -Inf, Inf, ..., stop.on.error=FALSE)$value
    } else {

        negative = integrate(levy_measure, -Inf, -x, ..., stop.on.error=FALSE)$value
        positive = integrate(levy_measure, x, Inf, ..., stop.on.error=FALSE)$value
        total = negative + positive
    }
	return(total)
}

###########################
# moment from Levy measure#
###########################

moment_from_measure = function(k, levy_measure, ...) {
	integrate(function(y){y^k*levy_measure(y,...)},-Inf,Inf)
}


#########################
# brlen transformations #
#########################

OU.brlen = function(phy,theta=1e-6) {
  phy = reorder(phy,"postorder")
  n_tip = phy$Nnode + 1
  
  # get raw ages/times
  a = branching.times(phy)
  T = max(a)
  t = T - a
  
  # get OU-scaled ages/times
  t.ou = 1/(2*theta) * exp(-2*theta*(T-t)) * (1-exp(-2*theta*t))
  h.ou = 1/(2*theta) * exp(-2*theta*(T-T)) * (1-exp(-2*theta*T))
  a.ou = h.ou - t.ou
  a.ou = c( rep(0, n_tip), a.ou)
  
  # assign OU-scaled times to tree
  for (i in 1:nrow(phy$edge))
  {
      phy$edge.length[i] = a.ou[phy$edge[i,1]] - a.ou[phy$edge[i,2]]
  }
  return(phy)
}

EB.brlen = function(phy,r=1e-6) {
    phy = reorder(phy,"postorder")
    n_tip = length(phy$tip.label)

    # get raw ages/times
    a = branching.times(phy)
    T = max(a)
    t = T - a
    t = c(rep(T, n_tip), t)

    # assign EB-scaled times to tree
    for (i in 1:nrow(phy$edge))
    {
        t_pa = t[phy$edge[i,1]]
        t_ch = t[phy$edge[i,2]]
        dx = exp(r*t_ch) - exp(r*t_pa)
        phy$edge.length[i] = dx/r
    }
    
    return(phy)
}

####################
# parameter labels #
####################

# parameter labeling
format_params = function(p, m) {
    param_names = c("sigma_bm",
                    "lambda_jn",
                    "delta_jn",
                    "sigma_vg",
                    "nu_vg",
                    "alpha_nig",
                    "delta_nig",
                    "alpha_ou",
                    "decay_eb",
                    "sigma_tip")
    n_param = length(param_names)
    x = rep(0, n_param)
    names(x) = param_names

    x[n_param] = p[length(p)]
    if (m=="BM") {
        x[1]=p[1]
    } else if (m=="BMJN") {
        x[1]=p[1]
        x[2]=p[2]
        x[3]=p[3]
    } else if (m=="BMVG") {
        x[1]=p[1]
        x[4]=p[2]
        x[5]=p[3]
    } else if (m=="BMNIG") {
        x[1]=p[1]
        x[6]=p[2]
        x[7]=p[3]
    } else if (m=="JN") {
        x[2]=p[1]
        x[3]=p[2]
    } else if (m=="VG") {
        x[4]=p[1]
        x[5]=p[2]
    } else if (m=="NIG") {
        x[6]=p[1]
        x[7]=p[2]
    } else if (m=="OU") {
        x[1]=p[1]
        x[8]=p[2]
    } else if (m=="EB") {
        x[1]=p[1]
        x[9]=p[2]
    }
    return(x)
}


#################
# data cleaning #
#################

drop.outlier = function(phy,dat,n=1,drop_zero=T,verbose=F)
{

    # expects dat to be a named vector
    if (!is.vector(dat))
    {
        stop("ERROR: dat is not a vector")
    }
    if (is.null(names(dat)))
    {
        stop("ERROR: dat does not have names")
    }

    to.drop = c()
    to.drop.zero = c()
    to.drop.outlier = c()

    if (drop_zero)
    {
        td = drop.zero(phy,dat)
        phy = td$phy
        dat = td$dat
        if (length(td$dropped) > 0)
        {
            to.drop.zero = td$dropped
        }
    }
    
    phy = reorder(phy,'postorder')
    dat = dat[phy$tip.label]
   
    # only drop outliers if n > 0 
    if (n > 0) {
        contrast = pic(dat,phy)
        nodes = order(abs(contrast),decreasing=TRUE)[1:n]
        nodes = length(phy$tip.label)+nodes
        for (node in nodes)
        {
            clade = extract.clade(phy,node)
            to.drop.outlier = clade$tip.label
        }
    }
    to.drop = unique(c(to.drop.outlier, to.drop.zero))

    #to.drop = unique(to.drop)
    #cat("to.drop\n")
    #print(to.drop)
    #cat("phy$tip.labels\n")
    #print(phy$tip.label)
    #print(c(length(to.drop), length(phy$tip.label)))
    
    if (length(to.drop) != length(phy$tip.label))
    {
	newPhy = drop.tip(phy, to.drop)
	newDat = dat[newPhy$tip.label]
    } else {
        newPhy = phy
        newDat = dat[newPhy$tip.label]
    }
    if (verbose) {
        cat("Dropped taxa\n")
        print(to.drop)
        cat("Dropped contrasts\n")
        print(contrast[order(abs(contrast),decreasing=TRUE)[1:n]])
        cat("Dropped taxa (drop_zero)    =",length(to.drop.zero),"\n")
        cat("Dropped taxa (drop_outlier) =",length(to.drop.outlier),"\n")
        cat("\n")
    }
    return(list(phy=newPhy,dat=newDat))
}

drop.zero = function(phy,dat,eps=0) {
	#expects dat to be a named vector
	if (!is.vector(dat)) {
		stop("ERROR: dat is not a vector")
	}
	if (is.null(names(dat))) {
		stop("ERROR: dat does not have names")
	} 
	phy = reorder(phy,"postorder")
	dat = dat[phy$tip.label]
	bad_nodes = c()
	pseudo_obs = rep(0,Nnode(phy,internal.only=FALSE))
	for (i in 1:length(phy$tip.label)) {
		pseudo_obs[i] = dat[i]
	}
	for (i in seq(1,nrow(phy$edge),2)) {
		cur_node = phy$edge[i,1]
		cur_left = phy$edge[i,2]
		cur_right = phy$edge[i+1,2]
		t_l = phy$edge.length[i]
		t_r = phy$edge.length[i+1]
		contrast = pseudo_obs[cur_left]-pseudo_obs[cur_right]
		if (!is.nan(contrast) && abs(contrast) <= eps) {
			bad_nodes = c(bad_nodes,cur_node)
		}
		pseudo_obs[cur_node] = t_r*pseudo_obs[cur_left]+t_l*pseudo_obs[cur_right]
		pseudo_obs[cur_node] = pseudo_obs[cur_node]/(t_l+t_r)

	}
	to.drop = c()
	for (node in bad_nodes) {
		clade = extract.clade(phy,node)
		to.drop = c(to.drop,clade$tip.label)
	}
	to.drop = unique(to.drop)
	newPhy = drop.tip(phy,to.drop)
	newDat = dat[newPhy$tip.label]
	return(list(phy=newPhy,dat=newDat,bad_nodes=bad_nodes,dropped=to.drop))
}


### plot PIC outliers for BM component of BM+LP

pic_outliers = function(phy,dat,sigma_bm) {
	cur_pic = pic(dat,phy)
	z_vals = cur_pic
	p_vals = pnorm(abs(z_vals),sd=sigma_bm,lower.tail=FALSE)
	return(-log10(p_vals))
}

get_p_cols = function(p_vals,phy) {
	normed = p_vals/max(p_vals)
	cols_rgb = colorRamp(c("white","red"),bias=100)(normed)
	cols_hex = apply(cols_rgb,1,function(x){rgb(x[1],x[2],x[3],maxColorValue=255)})
	names(cols_hex) = names(p_vals)
	return(cols_hex)
}

plot_jumps = function(x,cutoff=-log10(0.05),cex=.5,adj=.5,main="") {

    if (is.null(x$phy)) stop("x does not contain phy object!")
    if (is.null(x$dat)) stop("x does not contain dat object!")
    if (is.null(x$params)) stop("x does not contain params object!")
    if (!("sigma_bm" %in% names(x$params))) stop("x does not contain sigma_bm parameter!")
    
	sigma = x$params["sigma_bm"]
	dat_p_vals = pic_outliers(x$phy, x$dat, sigma)
    max_p_val = max(dat_p_vals[dat_p_vals != Inf])
    if (any(dat_p_vals==Inf)) {
        inf_idx = dat_p_vals==Inf
        dat_p_vals[inf_idx] = 1.5 * max_p_val
    }
	filt_p_vals = dat_p_vals[dat_p_vals>=cutoff]

	p_cols = get_p_cols(filt_p_vals)
	plot(x$phy,adj=adj,cex=cex,main=main)
	nodelabels(pch=16,node = as.numeric(names(p_cols)), col=p_cols,frame="circle")
	invisible(dat_p_vals)
}