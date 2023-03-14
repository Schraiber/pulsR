require(neldermead)
require(ape)
require(geiger)

#source("levy_pruning_prob.r")
#source("levy_pruning_cf.r")

################
# OPTIMIZATION #
################

# ensures initial values are within tolerance limits
tol_safe = function(x, tol) {
    # some redundancy, but modify if needed
    if (tol < 0 && x > tol) {
        x = tol - rexp(1,100)
    } else if (tol > 0 && x < tol) {
        x = tol + rexp(1,100)
    }
    return(x)
}

# Arguments
# phy               : an ape tree object
# dat               : a vector of univariate trait data with names
#                     matching phy$tip.label
# model             : a string indicating the model (BM, OU, EB, JN, VG, NIG, BMJN, BMVG, BMNIG)
# par               : initial parameters (must match model dimensions/bounds)
# sigma_tip         : estimate sigma_tip value? (default==TRUE)
# tol               : tolerance to terminate optimization (when delta-lnL or delta-param < tol)
# maxfun            : max number of function calls during optimization
# maxiter           : max number of iterations during optimization
# weight_scaling    : weight scaling to improve integration
# silent            : do not print text during optimization
# development       : use untested features
#
# Returns list object with the following elements
# params            : MLE parameters
# n_params          : number of parameters
# lnL               : MLE log likelihood
# AIC               : Akiake Information Criterion score
# optim             : neldermead optimization object

fit_reml_levy = function(phy,dat,model,par=NA,sigma_tip=T,tol=1e-4,maxfun=2000,maxiter=1000,weight_scaling=1e5,silent=T,development=F) {

    # validate model choice
    model = toupper(model)
    if (development) {
        if ( !(model %in% c("BM","OU","EB","JN","VG","NIG","BMJN","BMVG","BMNIG","EBJN","EBVG","EBNIG")) ) {
            stop("Provided model type invalid. Valid models: BM, OU, EB, JN, VG, NIG, BMJN, BMVG, BMNIG, EBJN, EBVG, EBNIG")
        }
    } else {
        if ( model %in% c("EBJN","EBVG","EBNIG") ) {
            stop("Set development=true to use these untested models: EBJN, EBVG, EBNIG")
        }
        if ( !(model %in% c("BM","OU","EB","JN","VG","NIG","BMJN","BMVG","BMNIG")) ) {
            stop("Provided model type invalid. Valid models: BM, OU, EB, JN, VG, NIG, BMJN, BMVG, BMNIG")
        }
    }

    # model tip noise?
    tip_noise = NA # estimate tipnoise
    tip_noise_lower = 1e-9
    if (!sigma_tip) {
        tip_noise = 0
        tip_noise_lower = 0
    }

    # initial parameter guess
    s = exp( runif(3, -0.5, 0.5) )
    alpha_ou = 1e-4
    decay_eb = -1e-4
    if (model=="OU" && runif(1) < 0.5) {
        # OU
        phy2 = chronos(phy, quiet=T)
        par_init = fitContinuous(phy2,dat,model="OU",SE=tip_noise)$opt
        alpha_ou = tol_safe( par_init$alpha * rbeta(1,1,5), tol )
    } else if (model%in%c("EB","EBJN","EBVG","EBNIG") && runif(1) < 0.5) {
        # EB
        phy2 = chronos(phy, quiet=T)
        par_init = fitContinuous(phy2,dat,model="EB",SE=tip_noise)$opt
        decay_eb = tol_safe( par_init$a * rbeta(1,1,5), -tol )
    } else {
        # all other models
        par_init = fitContinuous(phy,dat,model="BM",SE=tip_noise)$opt
    }
    sig_bm = tol_safe( sqrt(par_init$sigsq) * s[1], tol )
    if (sigma_tip) {
        sig_tip = tol_safe( par_init$SE * s[2], tol )
    } else {
        sig_tip = 0
    }
    proc_var = sig_bm^2
    if (runif(1) < 0.5) {
        proc_kurt = rgamma(1, shape=2, rate=0.2)
    } else {
        proc_kurt = runif(1,0.5,5)
    }
    frac_of_var = rbeta(1,2,2)


    # BROWNIAN MOTION
    if (model=="BM") {
        # set initial params
        if (any(is.na(par))) {
            par = c(sig_bm, sig_tip)
        }
        lower = rep(1e-9,2)
        upper = rep(1e+3,2)
        lower[2] = tip_noise_lower

        fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
            if (!sigma_tip) x[2] = 0
            if (any(x<0)) {
                ret = -Inf
            } else {
                ret = llik_reml_levy(phy,dat,"BM",sigma_bm=x[1],theta_ou=0,decay_eb=0,sigma_tip=x[2],weight_scaling=weight_scaling,silent=silent)
            }
            if (any(ret == -Inf)) ret = -(0.1/tol)
            return(list(f=-sum(ret),
            		g=c(),
            		c=c(),
            		gc=c(),
            		index=index,
            		this=list(costfargument=fmsfundata)))
        }
    }
    # ORNSTEIN-UHLENBECK
    #else if (all.equal(phi,ebm)==T && is.ou && !is.eb) {
    else if (model=="OU") {
      # set initial params
      if (any(is.na(par))) {
          par = c(sig_bm, alpha_ou, sig_tip)
      }
      lower = rep(1e-9,3)
      upper = rep(1e+3,3)
      lower[3] = tip_noise_lower
      fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
        if (!sigma_tip) x[3] = 0
        if (any(x<0)) {
            ret = -Inf
        } else {
            ret = llik_reml_levy(phy,dat,"OU",sigma_bm=x[1],theta_ou=x[2],decay_eb=0,sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
        }
        if (any(ret == -Inf)) ret = -(0.1/tol)
        return(list(f=-sum(ret),
                    g=c(),
                    c=c(),
                    gc=c(),
                    index=index,
                    this=list(costfargument=fmsfundata)))
      }
    }
    # EARLY BURST
    else if (model=="EB") {
      # set initial params
      if (any(is.na(par))) {
          par = c(sig_bm, decay_eb, sig_tip)
      }
      lower = c(1e-9,-1e+3,1e-9)
      upper = c(1e+3,-1e-9,1e+3)
      lower[3] = tip_noise_lower
      fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
        if (!sigma_tip) x[3] = 0
        if (x[1] < 0 || x[2] > 0 || x[3] < 0) {
            ret = -Inf
        } else {
            ret = llik_reml_levy(phy,dat,"EB",sigma_bm=x[1],theta_ou=0,decay_eb=x[2],sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
        }
        if (any(ret == -Inf)) ret = -(0.1/tol)
        return(list(f=-sum(ret),
                    g=c(),
                    c=c(),
                    gc=c(),
                    index=index,
                    this=list(costfargument=fmsfundata)))
      }
    }
    # BROWNIAN MOTION + COMPOUND POISSON WITH NORMAL JUMPS
    else if (model=="BMJN") {
        # set initial params
        if (any(is.na(par))) {
            par_tmp = get_params_JN(proc_var, proc_kurt, frac_of_var)
            sigma_bm  = par_tmp$sigma.bm
            lambda_jn = par_tmp$lambda.jn
            delta_jn  = par_tmp$delta.jn
            par = c(sigma_bm, lambda_jn, delta_jn, sig_tip)
        }
        lower = rep(1e-9,4)
        upper = rep(1e+6,4)
        lower[4] = tip_noise_lower
        fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
            if (!sigma_tip) x[4] = 0
            if (any(x<0)) {
                ret = -Inf
            } else {
                ret = llik_reml_levy(phy,dat,"BMJN",sigma_bm=x[1],lambda_jn=x[2],sigma_jn=x[3],
                                theta_ou=0,decay_eb=0,sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
            }
            if (any(ret == -Inf)) ret = -(0.1/tol)
            return(list(f=-sum(ret),
            		g=c(),
            		c=c(),
            		gc=c(),
            		index=index,
            		this=list(costfargument=fmsfundata)))
        }
    }
    # PURE COMPOUND POISSON WITH NORMAL JUMPS
    else if (model=="JN") {
        # set initial params
        if (any(is.na(par))) {
            par_tmp = get_params_JN(proc_var, proc_kurt, 1)
            lambda_jn = par_tmp$lambda.jn
            delta_jn  = par_tmp$delta.jn
            par = c(lambda_jn, delta_jn, sig_tip)
        }
        lower = rep(1e-9,3)
        upper = rep(1e+3,3)
        lower[3] = tip_noise_lower
        fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
            if (!sigma_tip) x[3] = 0
            if (any(x<0)) {
                ret = -Inf
            } else {
                ret = llik_reml_levy(phy,dat,"JN",lambda_jn=x[1],sigma_jn=x[2],
                                theta_ou=0,decay_eb=0,sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
            }
            if (any(ret == -Inf)) ret = -(0.1/tol)
            return(list(f=-sum(ret),
            		g=c(),
            		c=c(),
            		gc=c(),
            		index=index,
            		this=list(costfargument=fmsfundata)))
        }
    }
    # BROWNIAN MOTION + VARIANCE GAMMA
    else if (model=="BMVG") {
        # set initial params
        if (any(is.na(par))) {
            par_tmp = get_params_VG(proc_var, proc_kurt, frac_of_var)
            par = c(par_tmp$sigma.bm, par_tmp$nu.vg, par_tmp$sigma.vg, sig_tip)
        }
        lower = rep(1e-9,4)
        upper = rep(1e+3,4)
        lower[4] = tip_noise_lower
        fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
            if (!sigma_tip) x[4] = 0
            if (any(x<0)) {
                ret = -Inf
            } else {
                ret = llik_reml_levy(phy,dat,"BMVG",sigma_bm=x[1],sigma_vg=x[2],nu_vg=x[3],
                                theta_ou=0,decay_eb=0,sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
            }
            if (any(ret == -Inf)) ret = -(0.1/tol)
            return(list(f=-sum(ret),
            		g=c(),
            		c=c(),
            		gc=c(),
            		index=index,
            		this=list(costfargument=fmsfundata)))
        }
    }
    # PURE VARIANCE GAMMA
    else if (model=="VG") {
        # set initial params
        if (any(is.na(par))) {
            par_tmp = get_params_VG(proc_var, proc_kurt, 1)
            par = c(par_tmp$nu.vg, par_tmp$sigma.vg, sig_tip)
        }
        lower = rep(1e-9,3)
        upper = rep(1e+3,3)
        lower[3] = tip_noise_lower
        fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
            if (!sigma_tip) x[3] = 0
            if (any(x<0)) {
                ret = -Inf
            } else {
                ret = llik_reml_levy(phy,dat,"VG",sigma_vg=x[1],nu_vg=x[2],
                                theta_ou=0,decay_eb=0,sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
            }
            if (any(ret == -Inf)) ret = -(0.1/tol)
            return(list(f=-sum(ret),
            		g=c(),
            		c=c(),
            		gc=c(),
            		index=index,
            		this=list(costfargument=fmsfundata)))
        }
    }
    # BROWNIAN MOTION + NORMAL INVERSE GAUSSIAN
    else if (model=="BMNIG") {
        # set initial params
        if (any(is.na(par))) {
            par_tmp = get_params_NIG(proc_var, proc_kurt, frac_of_var)
            par = c(par_tmp$sigma.bm, par_tmp$alpha.nig, par_tmp$delta.nig, sig_tip)
        }
        lower = rep(1e-9,4)
        upper = rep(1e+3,4)
        lower[4] = tip_noise_lower
        fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
            if (!sigma_tip) x[4] = 0
            if (any(x<0)) {
                ret = -Inf
            } else {
                ret = llik_reml_levy(phy,dat,"BMNIG",sigma_bm=x[1],alpha_nig=x[2],delta_nig=x[3],
                                theta_ou=0,decay_eb=0,sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
            }
            if (any(ret == -Inf)) ret = -(0.1/tol)
            return(list(f=-sum(ret),
            		g=c(),
            		c=c(),
            		gc=c(),
            		index=index,
            		this=list(costfargument=fmsfundata)))
        }
    }
    # PURE NORMAL INVERSE GAUSSIAN
    else if (model=="NIG") {
        # set initial params
        if (any(is.na(par))) {
            par_tmp = get_params_NIG(proc_var, proc_kurt, frac_of_var)
            par = c(par_tmp$alpha.nig, par_tmp$delta.nig, sig_tip)
        }
        lower = rep(1e-9,3)
        upper = rep(1e+3,3)
        lower[3] = tip_noise_lower
        fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
            if (!sigma_tip) x[3] = 0
            if (any(x<0)) {
                ret = -Inf
            } else {
                ret = llik_reml_levy(phy,dat,"NIG",alpha_nig=x[1],delta_nig=x[2],
                                theta_ou=0,decay_eb=0,sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
            }
            if (any(ret == -Inf)) ret = -(0.1/tol)
            return(list(f=-sum(ret),
            		g=c(),
            		c=c(),
            		gc=c(),
            		index=index,
            		this=list(costfargument=fmsfundata)))
        }
    }
    # EARLY BURST w/ RATE DECAYING JN
    else if (model=="EBJN") {
      # set initial params
      if (any(is.na(par))) {
        par_tmp = get_params_JN(proc_var, proc_kurt, 1)
        lambda_jn = par_tmp$lambda.jn
        delta_jn  = par_tmp$delta.jn
        par = c(decay_eb, lambda_jn, delta_jn, sig_tip)
      }
      lower = c(-1e+3,1e-9,1e-9,1e-9)
      upper = c(-1e-9,1e+3,1e+3,1e+3)
      lower[3] = tip_noise_lower
      fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
        if (!sigma_tip) x[4] = 0
        if (x[1] > 0 || x[2] <  0 || x[3] < 0 || x[4] < 0) {
            ret = -Inf
        } else {
            ret = llik_reml_levy(phy,dat,"EBJN",theta_ou=0,decay_eb=x[1],lambda_jn=x[2],sigma_jn=x[3],sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
        }
        if (any(ret == -Inf)) ret = -(0.1/tol)
        return(list(f=-sum(ret),
                    g=c(),
                    c=c(),
                    gc=c(),
                    index=index,
                    this=list(costfargument=fmsfundata)))
      }
    }
    
    # EARLY BURST w/ TIME DECAYING VG
    else if (model=="EBVG") {
      # set initial params
      if (any(is.na(par))) {
        par_tmp = get_params_VG(proc_var, proc_kurt, 1)
        par = c(decay_eb, par_tmp$nu.vg, par_tmp$sigma.vg, 1)
      }
      lower = c(-1e+3,1e-9,1e-9,1e-9)
      upper = c(-1e-9,1e+3,1e+3,1e+3)
      lower[3] = tip_noise_lower
      fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
        if (!sigma_tip) x[4] = 0
        if (x[1] > 0 || x[2] <  0 || x[3] < 0 || x[4] < 0) {
            ret = -Inf
        } else {
            ret = llik_reml_levy(phy,dat,"EBVG",theta_ou=0,decay_eb=x[1],nu_vg=x[2],sigma_vg=x[3],sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
        }
        if (any(ret == -Inf)) ret = -(0.1/tol)
        return(list(f=-sum(ret),
                    g=c(),
                    c=c(),
                    gc=c(),
                    index=index,
                    this=list(costfargument=fmsfundata)))
      }
    }

    # EARLY BURST w/ TIME DECAYING NIG
    else if (model=="EBNIG") {
      # set initial params
      if (any(is.na(par))) {
        par_tmp = get_params_NIG(proc_var, proc_kurt, 1)
        par = c(decay_eb, par_tmp$alpha.nig, par_tmp$delta.nig, sig_tip)
      }
      lower = c(-1e+3,1e-9,1e-9,1e-9)
      upper = c(-1e-9,1e+3,1e+3,1e+3)
      lower[3] = tip_noise_lower
      fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
        if (!sigma_tip) x[4] = 0
        if (x[1] > 0 || x[2] <  0 || x[3] < 0 || x[4] < 0) {
            ret = -Inf
        } else {
            ret = llik_reml_levy(phy,dat,"EBNIG",theta_ou=0,decay_eb=x[1],alpha_nig=x[2],delta_nig=x[3],sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
        }
        if (any(ret == -Inf)) ret = -(0.1/tol)
        return(list(f=-sum(ret),
                    g=c(),
                    c=c(),
                    gc=c(),
                    index=index,
                    this=list(costfargument=fmsfundata)))
      }
    }
    else {
        return(NA)
    }

    x0 = transpose(par)
    nm = neldermead()
    nm = neldermead.set(nm,'numberofvariables',length(par))
    nm = neldermead.set(nm,'function',fn)
    nm = neldermead.set(nm,'x0',x0)
    nm = neldermead.set(nm,'verbose',FALSE)
    nm = neldermead.set(nm,'storehistory',TRUE)
    nm = neldermead.set(nm,'verbosetermination',FALSE)
    nm = neldermead.set(nm,'method','box')
    nm = neldermead.set(nm,'boundsmin',lower)
    nm = neldermead.set(nm,'boundsmax',upper)
    nm = neldermead.set(nm,'tolfunmethod', TRUE)
    nm = neldermead.set(nm,'tolfunrelative',tol)
    nm = neldermead.set(nm,'tolxmethod',TRUE)
    nm = neldermead.set(nm,'tolxrelative',tol)
    nm = neldermead.set(nm,'maxfunevals',maxfun)
    nm = neldermead.set(nm,'maxiter',maxiter)
    nm = neldermead.set(nm,'boxtermination',TRUE)
    nm = neldermead.set(nm,'boxtolf',tol*1e-2)
    nm = neldermead.search(nm)
    nm$gc = gc()
    nm$phy = phy

    if (!sigma_tip) {
        tip_noise_idx = length( nm$optbase$xopt[,1] )
        p = nm$optbase$xopt[,1]
        p[length(p)]=0
        nm$optbase$xopt[,1] = p
    }
    
    if (!silent) cat ("...done!\n")
    
    # return
    results = list()
    results$model = model
    results$dat = dat
    results$phy = phy
    results$params = format_params(nm$optbase$xopt[,1], model)
    results$n_params = sum( results$params != 0 )
    results$lnL = -nm$optbase$fopt
    results$AIC = 2*(results$n_params-results$lnL)
    results$optim = nm
    return(results)
}
 

get_best_idx = function(x) {
    n = length(x)
    best = which.min(unlist(sapply(1:n,
        function(i)
        {
            z=x[[i]]$optbase$fopt
            if (!is.numeric(z))
                z = NA
            return(z)
        }
    )))
    return(best)
}


