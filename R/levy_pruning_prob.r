#source("levy_pruning_tools.r")
 
rescale_error_cf = function(left_scale,right_scale,left_times,right_times,sigma_left,sigma_right,t_l,t_r,v_l,v_r,sig_p2) {
	#rescale left by multiplying by RIGHT branch rescaling factor
	new_left_scale = (sig_p2*t_r+v_r)*c(1,left_scale)/
		(sig_p2*(t_l+t_r)+v_l+v_r)
	#rescale right by multiplying by LEFT branch rescaling factor
	new_right_scale = (sig_p2*t_l+v_l)*c(1,right_scale)/
		(sig_p2*(t_l+t_r)+v_l+v_r)
	#add the new branch lengths
	new_left_times = c(t_l,left_times)
	new_right_times = c(t_r,right_times)
	#add the new sigmas
	new_sigma_left = c(0,sigma_left)
	new_sigma_right = c(0,sigma_right)
	return(list(scales=c(new_left_scale,new_right_scale),times=c(new_left_times,new_right_times),sigma=c(new_sigma_left,new_sigma_right)))
}

#NB: this function expects phi to be the characteristic exponent!
error_cf.vectorized = function(k,scale,times,phi,...) {
	if (length(times)!=length(scale)) {
		stop("length(times) = ",length(times)," != length(scale) = ",length(scale))
 	} 
 	if (length(times) == 0) {
 		return(rep(1,length(k)))
 	}
 	#vectorize the computation
 	scaleKmat = scale%*%t(k)
 	timeMat = matrix(nrow=length(times),ncol=length(k),times)
 	errMat = phi(scaleKmat,timeMat,...)
	cur_err = colSums(errMat)
	return(cur_err)
}

#NB: This function expects phi to be the characteristic exponent!
error_cf = function(k,sigma_tip,scale,times,phi,...) {

 	# here be the tipward pass to compute the cf
 	cur_err = 0
	#TODO: THIS FLAGGING ISN'T WORKING
	for (i in 1:length(scale)) {
		#TODO: for the tips, need to add a *different* phi
		if (sigma_tip[i] == 0) {
			cur_err = cur_err + phi(scale[i]*k,times[i],...)
		} else {
			#cur_err = cur_err - sigma_tip[i]^2*scale[i]^2*k^2/2
			cur_err = cur_err - sigma_tip[i]*sigma_tip[i]*scale[i]*scale[i]*k*k/2
		}
	}
	return(cur_err)
}

#NB: This function expects phi to be the characteristic exponent!
contrast_cf = function(k,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,...) {
	#compute the cf of the change to the left
	cur_cf = phi(k,t_l,...)
	#compute the cf of the change to the right (evaluated at negative)
	cur_cf = cur_cf+phi(-k,t_r,...)
	#compute the cf of the error to the left
	cur_cf = cur_cf+error_cf(k,sigma_left,left_scale,left_times,phi,...)
	#compute the cf of the error to the right (evaluated at negative)
	cur_cf = cur_cf+error_cf(-k,sigma_right,right_scale,right_times,phi,...)
	return(exp(cur_cf))
}

weight = function(k,p=8,q=4) {
    # returns 0.5*erfc(k/p-q)
    return(pnorm(-sqrt(2)*(k/p-q)))
}

fourier_real = function(x,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,..., weight_scaling=100000) {
	#gets the real part of the fourier integral
	return(function(k) {
		cur_cf = contrast_cf(k,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,...)
		return(weight(k,sqrt(weight_scaling/abs(x)),sqrt(weight_scaling*abs(x)/4))*(Re(cur_cf)*cos(-k*x)-Im(cur_cf)*sin(-k*x)))
	})
}

contrast_likelihood = function(x,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,...,weight_scaling=100000) {
	contrast_re = fourier_real(x,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,...,weight_scaling=weight_scaling)
	#compute integral
	real_part = 1/(pi)*integrate(contrast_re,0,Inf,subdivisions=1000000L,rel.tol=1e-50,abs.tol=1e-50,stop.on.error=F)$value
	return(real_part)
}

llik_reml_levy = function(phy,dat,model,...,return_inf=TRUE,theta_ou=0,decay_eb=0,sigma_tip=0,weight_scaling=100000,silent=T) {

    ### INITIALIZE DATA ###
	#sort the tree so that it can be pruned
	phy = reorder(phy,"postorder")
  
	#ensure that data is a named vector
	if (!is.vector(dat)) {
		stop("ERROR: dat is not a vector")
	}
	if (is.null(names(dat))) {
		stop("ERROR: dat does not have names")
	} 
    # label and order data
    old_length = length(dat)
    dat = dat[phy$tip.label]
    if (length(dat)!=old_length) {
    	warning("There are discrepancies between dat and phy")
    }

    ### INITIALIZE MODEL ###
    model = toupper(model)
    if ( !(model %in% c("BM","OU","EB","JN","VG","NIG","BMJN","BMVG","BMNIG")) ) {
        stop("Provided model type invalid. Valid models: BM, OU, EB, JN, VG, NIG, BMJN, BMVG, BMNIG")
    }
    phi = get_cf(model)
	
	# rescale by EB/OU parameters for use with ebm
    if (theta_ou > 0 && model=="OU") {
        phy = OU.brlen(phy=phy, theta=theta_ou)
 	}
    else if (decay_eb != 0 && model=="EB") {
        phy = EB.brlen(phy=phy, r=decay_eb)
    }

    ### INITIALIZE WORKSPACE ####
	scales = list()
	times = list()
	sigma_tip_list = list()
	pseudo_obs = rep(0,Nnode(phy,internal.only=FALSE))
	contrast_lnL = rep(0,Nnode(phy,internal.only=FALSE))
	contrast_val = rep(0,Nnode(phy,internal.only=FALSE))
	error_var = rep(0,Nnode(phy,internal.only=FALSE))
    #get the process variance
	sig_p2 = phi(1,1,...,variance=TRUE)


    ### COMPUTE LIKELIHOOD ###
	#first initialize the tips
	for (i in 1:length(phy$tip.label)) {
		scales[[i]] = 1 #initalize this at 1 so that the tip errors are included
		times[[i]] = 0
		sigma_tip_list[[i]] = sigma_tip
		pseudo_obs[i] = dat[i]
		contrast_lnL[i] = 0
		error_var[i] = sigma_tip^2
	}
	#then loop over the internal nodes
	for (i in seq(1,nrow(phy$edge),2)) {
		#get the nodes
		cur_node = phy$edge[i,1]
		cur_left = phy$edge[i,2]
		cur_right = phy$edge[i+1,2]
		#get the relevant scalings and times
		left_scale = scales[[cur_left]]
		right_scale = scales[[cur_right]]
		sigma_left = sigma_tip_list[[cur_left]]
		sigma_right = sigma_tip_list[[cur_right]]
		left_time = times[[cur_left]]
		right_time = times[[cur_right]]
		t_l = phy$edge.length[i]
		t_r = phy$edge.length[i+1]
		v_l = error_var[cur_left]
		v_r = error_var[cur_right]
		#compute the contrast likelihood
		contrast = pseudo_obs[cur_left]-pseudo_obs[cur_right]
		cur_like = contrast_likelihood(contrast,sigma_left,sigma_right,left_scale,right_scale,left_time,right_time,t_l,t_r,phi,...,weight_scaling=weight_scaling)
		contrast_val[cur_node] = contrast
		if (cur_like > 0) {
			contrast_lnL[cur_node] = log(cur_like)			
		} else if (!return_inf) {
			contrast_lnL[cur_node] = -300
		} else {
            return(-Inf)
		}
		#compute the new error var
		new_error_var = (sig_p2*t_l+v_l)*(sig_p2*t_r+v_r)/(sig_p2*(t_l+t_r)+v_l+v_r)
		error_var[cur_node] = new_error_var
		#compute the new error scaling
		new_error = rescale_error_cf(left_scale,right_scale,left_time,right_time,sigma_left,sigma_right,t_l,t_r,v_l,v_r,sig_p2)
		scales[[cur_node]] = new_error$scales
		times[[cur_node]] = new_error$times
		sigma_tip_list[[cur_node]] = new_error$sigma
		#compute the new pseudo-observation
		pseudo_obs[cur_node] = ((sig_p2*t_r+v_r)*pseudo_obs[cur_left]+
			(sig_p2*t_l+v_l)*pseudo_obs[cur_right])/
			(sig_p2*(t_l+t_r)+v_l+v_r)
		#pseudo_obs[cur_node] = pseudo_obs[cur_node]/(t_l+t_r)
	}
	
    results = c(sum(contrast_lnL),...)
    names(results)[1] = "lnL"
    if (decay_eb!=0) {
        names(decay_eb)="decay_eb"
        results = c(results, decay_eb)
    }
    if (theta_ou>0) {
        names(theta_ou)="theta_ou"
        results = c(results, theta_ou)
    }
    names(sigma_tip)="sigma_tip"
    results = c(results,sigma_tip)
    
    if (!silent) print(results)

	return(contrast_lnL)
}

