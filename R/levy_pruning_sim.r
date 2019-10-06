require(fBasics) #to simulate alpha stable guys, ?dstable
require(ape)
require(gsl) #for better bessel functions
require(stabledist)
require(statmod)

#source("levy_pruning_tools.r")

rlevy = function(phy, model, par, n=1) {

    # format data
    par = check_args(model, par)
    print(par)
    phy = reorder(phy, "cladewise")
    
    # simulate
    if (model=="BM") {
        x = sim_bm(phy=phy,sigma_bm=par$sigma_bm,sigma_tip=par$sigma_tip)
    } else if (model=="OU") {
        phy2 = OU.brlen(phy=phy,theta=par$theta_ou)
        phy2 = reorder(phy2, "cladewise")
        x = sim_bm(phy=phy2,sigma_bm=par$sigma_bm,sigma_tip=par$sigma_tip)
    } else if (model=="EB") {
        phy2 = EB.brlen(phy=phy, r=par$decay_eb)
        phy2 = reorder(phy2, "cladewise")
        x = sim_bm(phy=phy2,sigma_bm=par$sigma_bm,sigma_tip=par$sigma_tip)
    } else if (model=="JN") {
        x = sim_cpp(phy=phy,
                    sigma=0,
                    lambda=par$lambda_jn,
                    kernel=rnorm,
                    mean=0,
                    sd=par$delta_jn,
                    sigma_tip=par$sigma_tip)
    } else if (model=="VG") {
        x = sim_vg(phy=phy,
                   sigma_bm=0,
                   sigma_vg=par$sigma_vg,
                   nu_vg=par$nu_vg,
                   mu_vg=0,
                   sigma_tip=par$sigma_tip)
    } else if (model=="NIG") {
       x = sim_nig(phy=phy,
                    sigma_bm=0,
                    alpha_nig=par$alpha_nig,
                    beta_nig=0,
                    delta_nig=par$delta_nig,
                    mu_nig=0,
                    sigma_tip=par$sigma_tip) 
    } else if (model=="BMJN") {
        x = sim_cpp(phy=phy,
                    sigma=par$sigma_bm,
                    lambda=par$lambda_jn,
                    kernel=rnorm,
                    mean=0,
                    sd=par$delta_jn,
                    sigma_tip=par$sigma_tip)
    } else if (model=="BMVG") {
        x = sim_vg(phy=phy,
                   sigma_bm=par$sigma_bm,
                   sigma_vg=par$sigma_vg,
                   nu_vg=par$nu_vg,
                   mu_vg=0,
                   sigma_tip=par$sigma_tip)
    } else if (model=="BMNIG") {
        x = sim_nig(phy=phy,
                    sigma_bm=par$sigma_bm,
                    alpha_nig=par$alpha_nig,
                    beta_nig=0,
                    delta_nig=par$delta_nig,
                    mu_nig=0,
                    sigma_tip=par$sigma_tip)
    } else if (model=="AS") {
        x = sim_stable(phy=phy,
                       sigma=0,
                       alpha=par$alpha_as,
                       cee=par$c_as,
                       sigma_tip=par$sigma_tip)
    } else if (model=="BMAS") {
         x = sim_stable(phy=phy,
                       sigma=par$sigma_bm,
                       alpha=par$alpha_as,
                       cee=par$c_as,
                       sigma_tip=par$sigma_tip)
    }
    
    return(x)
}


jump_dif = function(sigma,lambda,kernel,...) {
	#sigma = rate*t of brownian motion
	#lambda = rate*t of poisson process
	#kernel = function to draw from jump kernel
	#... = parameters for kernel (MUST BE IN CORRECT ORDER!)
	numJumps = rpois(1,lambda)
	#Jbranch <<- append(Jbranch,numJumps)
	curState = 0
	curState = rnorm(1,0,sd=sigma)
	if (numJumps > 0) {
		curState = curState + sum(kernel(numJumps,...))
	}
	return(curState)	
}

sim_bm = function(phy, sigma_bm, sigma_tip = 0) {
	#phy is an ape-format tree
	#sigma is rate of brownian motion
	#nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
    nodes = rep(NA, nrow(phy$edge)+1)
	nodes[phy$edge[1,1]] = 0 # root value
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
        dx = rnorm(n=1, mean=0, sd=sigma_bm*sqrt(curLen)) #,lambda*curLen,kernel,...)
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + dx
        names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
	}
    nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
    return( nodes + rnorm(length(phy$tip.label),mean=0,sd=sigma_tip) )
}

sim_stable = function(phy, sigma, alpha, cee, sigma_tip=0) {
	#phy is ape-format tree
	#sigma is rate of brownian motion
	#cee is the scale parameter
	#alpha is the stability parameter
	nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
	nodes[phy$edge[1,1]] = 0
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, mean=0, sd = sigma*sqrt(curLen)) + rstable(1,alpha=alpha,beta=0,gamma=(curLen)^(1/alpha)*cee)
        names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
	}
	#print(nodes)
    nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
    return( nodes + rnorm(length(phy$tip.label),0,sigma_tip) )
	
}

sim_vg = function(phy,sigma_bm,sigma_vg,nu_vg,mu_vg,sigma_tip=0) {
	#parameterized as in Madan, Carr, Chang (1998)
	#phy is ape-format tree
	#sigma_bm is brownian motion rate
	#nu_vg controls kurtosis of VG
	#sigma_vg controls rate of VG
	#mu_vg controls skewness of VG
	
	#implements as a time-changed brownian motion.
	nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
	nodes[phy$edge[1,1]] = 0
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
		curTimeChange = rgamma(1,curLen/nu_vg,scale=nu_vg)
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, sd = sigma_bm*sqrt(curLen)) + rnorm(1, sd = sigma_vg*sqrt(curTimeChange), mean = mu_vg*curTimeChange)
        	names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
	}
	#return(rnorm(length(phy$tip.label),nodes[1:length(phy$tip.label)],sigma_tip))
    nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
    return( nodes + rnorm(length(phy$tip.label),0,sigma_tip) )
}

#implements as a time-changed brownian motion.
rVGBM = function(n,t,kap,sig_vg,sig_bm) {
    return(rVG(n,t,kap,sig_vg)+rnorm(n,sd=sig_bm*sqrt(t)))
}
rVG = function(n,t,kap,sig) {
	gammas = rgamma(n,t/kap,scale=kap)
	to_return = rnorm(n,sd=sig*sqrt(gammas))
	to_return[to_return==0] = .Machine$double.xmin
	return(to_return)
}

sim_nig = function(phy,sigma_bm,alpha_nig,beta_nig,delta_nig,mu_nig,sigma_tip=0) {
	#parametrized as in PhD seminar for modeling normal inverse gaussian processses
	#www2.math.uni-wuppertal.de/~ruediger/pages/vortraege/ws1112/nig_5.pdf
	#phy is ape-format tree
	#sigma_bm is brownian motion rate
	#alpha_nig controls kurtosis of NIG
	#beta_nig controls skewness of NIG
	#delta_nig controls the scale of NIG
	#mu_nig is location parameter of NIG

	#implements as time-changed brownian motion
	gamma_nig = sqrt(alpha_nig^2 - beta_nig^2)
	nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
	nodes[phy$edge[1,1]] = 0
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
		curTimeChange = rinvgauss(1,delta_nig/gamma_nig*curLen,delta_nig^2*curLen^2)
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, sd = sigma_bm*sqrt(curLen)) + rnorm(1, sd = sqrt(curTimeChange), mean = mu_nig*curLen + beta_nig*curTimeChange)
		names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
	}
	nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
	return(nodes + rnorm(length(phy$tip.label),0,sigma_tip) )
}

rNIG = function(n=1,sigma_bm,alpha_nig,beta_nig,delta_nig,mu_nig,t=1) {
	#parametrized as in PhD seminar for modeling normal inverse gaussian processses
	#www2.math.uni-wuppertal.de/~ruediger/pages/vortraege/ws1112/nig_5.pdf
	#phy is ape-format tree
	#sigma_bm is brownian motion rate
	#alpha_nig controls kurtosis of NIG
	#beta_nig controls skewness of NIG
	#delta_nig controls the scale of NIG
	#mu_nig is location parameter of NIG

	#implements as time-changed brownian motion
	gamma_nig = sqrt(alpha_nig^2 - beta_nig^2)
    curLen = t
    curTimeChange = rinvgauss(n,delta_nig/gamma_nig*curLen,delta_nig^2*curLen^2)
    x = rnorm(n, sd=sigma_bm*sqrt(curLen)) + rnorm(n, sd=sqrt(curTimeChange), mean=mu_nig*curLen + beta_nig*curTimeChange)
    return(x)
}

#this simulates a compound poisson model with jump kernel given by kernel and parameters given by ...
sim_cpp = function(phy, sigma, lambda, kernel,...,sigma_tip = 0) {
	#phy is an ape-format tree
	#sigma is rate of brownian motion
	#lambda is rate of poisson process
	#... are kernel parameters
	#nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
    nodes = rep(NA, nrow(phy$edge)+1)
	nodes[phy$edge[1,1]] = 0 #runif(1,-5,5) #sets the root
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
        dx = jump_dif(sigma*sqrt(curLen),lambda*curLen,kernel,...)
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + dx
        names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
		
	}
    nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
    return( nodes + rnorm(length(phy$tip.label),0,sigma_tip) )
}

check_args = function(model, par) {
    model = toupper(model)
    if ( !(model %in% c("BM","OU","EB","JN","VG","AS","NIG","BMJN","BMVG","BMNIG","BMAS")) ) {
        stop("Provided model type invalid. Valid models: BM, OU, EB, JN, VG, NIG, AS, BMJN, BMVG, BMNIG, BMAS")
    }
    
    if (!("sigma_tip" %in% names(par)))
        par$sigma_tip = 0
     
    # things with a Brownian component
    if (model=="BM" || model=="OU" || model=="EB" || 
        model=="BMJN" || model=="BMVG" || model=="BMNIG" || model=="BMAS") {
        
        if (!("sigma_bm" %in% names(par)))
            stop("Value for par$sigma_bm missing!")
    }
    
    # Ornstein-Uhlenbeck model
    if (model=="OU") {
        if (!("theta_ou" %in% names(par)))
            stop("Value for par$theta_ou missing!")
    }
    
    # Early Burst model
    if (model=="EB") {
        if (!("decay_eb" %in% names(par)))
            stop("Value for par$decay_eb missing!")
    }
    
    # jump normal models
    if (model=="JN" || model=="BMJN") {
        if (!("lambda_jn" %in% names(par)))
            stop("Value for par$lambda_jn missing!")
        if (!("delta_jn" %in% names(par)))
            stop("Value for par$delta_jn missing!")
    }
    
    # variance gamma models
    if (model=="VG" || model=="BMVG") {
        if (!("nu_vg" %in% names(par)))
            stop("Value for par$nu_vg missing!")
        if (!("sigma_vg" %in% names(par)))
            stop("Value for par$sigma_vg missing!")
    }
    
    # normal inverse gaussian models
    if (model=="NIG" || model=="BMNIG") {
        if (!("alpha_nig" %in% names(par)))
            stop("Value for par$alpha_nig missing!")
        if (!("delta_nig" %in% names(par)))
            stop("Value for par$delta_nig missing!")
    }
    
    # alpha stable models
    if (model=="AS" || model=="BMAS") {
        if (!("c_as" %in% names(par)))
            stop("Value for par$c_as missing!")
        if (!("alpha_as" %in% names(par)))
            stop("Value for par$alpha_as missing!")
    }
    
    # return parameters
    return(par)
}


.simulate_test = function(phy) {
    models = c("BM","OU","EB","JN","VG","NIG","AS","BMJN","BMVG","BMNIG","BMAS")
    par = list(sigma_bm=1,
               theta_ou=0.1,
               decay_eb=-0.1,
               lambda_jn=0.5,
               delta_jn=2,
               nu_vg=5,
               sigma_vg=2,
               alpha_nig=0.8,
               delta_nig=5,
               c_as=2,
               alpha_as=1.5,
               sigma_tip=0.1)
    
    phy = reorder(phy, "cladewise")
    x = list()
    for (m in models) {
        x[[m]] = rlevy(phy=phy, model=m, par=par)
        
        if (any(is.na(x[[m]]))) {
            stop(paste("simulations using",m,"generated NA values!"))
        }
    }
    
    return(x)
}
