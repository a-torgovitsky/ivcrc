capture program drop ivcrc_estimate
capture mata: mata drop gridparse()

program ivcrc_estimate, eclass		

				
		syntax [anything] [if] [in]  				///
                [,NOConstant						///
				  GENerate(string)					///
				  INTegral(string)					///				  
				  QUantiles(integer 100)			///
				  BANDwidth(numlist)				///
				  BOOTstrap(string)					///
				  Kernel(string) *]

				  
		_iv_parse `0'
        
        local lhs `s(lhs)'
        local endog `s(endog)'
        local exog `s(exog)'
        local inst `s(inst)'
        local 0 `s(zero)'	
		
		tempvar rankstat weigh
		tempname beta tempbeta bsub bbw		  
	  
		marksample touse
        markout `touse' `lhs' `exog' `inst' `endog'
		qui count if `touse'
		ereturn scalar N = r(N)
		
		local k : word count `kernel'
		if `k' > 1 {
			display as err "Only one kernel may be specified"
			error 198
		}
		
		** Define integral set & subsets
		if trim("`integral'") == "" {
			local integral ".01(.01).99"
		}
		mata: gridparse("`integral'")
		if `gridflag'==1 {
			dis as err "Integral subsets must contain 1600 or fewer grid points, one or more grid steps automatically reduced"	
		}
		if `ascendflag'>0 {
			dis as err "Integral subsets must not overlap"
			error 198
		}
		forvalues s = 1(1)`nsubsets' {
			numlist "`int`s''"
			local intgrid`s' `r(numlist)'	  
			tempname bbw_s`s'
		}
		
		** Bandwidth specification
		if trim("`bandwidth'") == "" {
			local bandwidth .05
		}
		local numbw : word count `bandwidth'
		
		** Conditional rank statistic computation
		qui forval s = 1(1)`= `quantiles' - 1' {
			local ss = `s'/`quantiles'
			qui qreg `endog' `inst' `exog' if `touse', quantile(`ss') 
			tempvar q_`s' ind_`s'
			qui predict `q_`s'' if `touse'
			qui gen `ind_`s'' = `q_`s''<=`endog' if `touse'
			local ind_list `ind_list' `ind_`s''
		}
		qui egen `rankstat' = rowmean(`ind_list') 
		
		tokenize "`generate'", parse(",")
		if "`generate'" != "" & trim("`bootstrap'") == "" {
			if "`2'" == "" {
				qui gen `1' = `rankstat' if `touse'
				display "Conditional rank statistic has been storred to `1'"
			}
			else if "`2'" == "," & "`3'" == "" {
				qui gen `1' = `rankstat' if `touse'
				display "Conditional rank statistic has been storred to `1'"
			}
			else if "`2'" == "," & ("`3'" == "replace" | "`3'" == "rep" | "`3'" == "repl" | "`3'" == "re" | "`3'" == "r") {
				qui replace `1' = `rankstat' if `touse'
				display "Conditional rank statistic has been storred to `1'"
			} 
			else {
				display as err "Syntax error with generate()"
				error 198
			}			
		}
				
		** Beta(r) computation via weighted least squares
		capture mata: mata drop beta*
		qui gen `weigh' = . if `touse'		
		qui forvalues s = 1(1)`nsubsets' {	
			qui foreach r of local intgrid`s' {
				qui foreach h of local bandwidth {					
					if  "`kernel'" == "" |			/// 
						"`kernel'" == "uniform" |	///
						"`kernel'" == "unif" |	 	///
						"`kernel'" == "un" | 		///
						"`kernel'" == "rectangle" |	///
						"`kernel'" == "rec"	{
						qui replace `weigh' = (1/2)*(abs((`rankstat' - `r')/`h')<=1) if `touse'
					}
					else if "`kernel'" == "biweight" | 	///
							"`kernel'" == "biweigh" | 	///
							"`kernel'" == "biw" | 		///
							"`kernel'" == "bi" {
						qui replace `weigh' = (15/16)*(1 - ((`rankstat' - `r')/`h')^2)^2*(abs((`rankstat' - `r')/`h')<=1) if `touse'
					}
					else if "`kernel'" == "triweight" | 		///
							"`kernel'" == "triweigh" | 			///
							"`kernel'" == "triw" {
						qui replace `weigh' = (35/32)*(1 - ((`rankstat' - `r')/`h')^2)^3*(abs((`rankstat' - `r')/`h')<=1) if `touse'
					}
					else if "`kernel'" == "cosine" |		///
							"`kernel'" == "cos" {
						qui replace `weigh' = (_pi/4)*cos((_pi/2)*(`rankstat' - `r')/`h')*(abs((`rankstat' - `r')/`h')<=1) if `touse'
					}
					else if "`kernel'" == "epanechnikov" |		///
							"`kernel'" == "epanech" |			///
							"`kernel'" == "epan" |				///
							"`kernel'" == "ep" {
						qui replace `weigh' = (3/4)*(1 - ((`rankstat' - `r')/`h')^2)*(abs((`rankstat' - `r')/`h')<=1) if `touse'
					}
					else if "`kernel'" == "gaussian" | 		///
							"`kernel'" == "gauss" | 		///
							"`kernel'" == "gaus" | 			///
							"`kernel'" == "gau" {
						qui replace `weigh' = (1/sqrt(2*_pi))*exp(-(1/2)*((`rankstat' - `r')/`h')^2) if `touse'
					}
					else if "`kernel'" == "triangle" | 	///
							"`kernel'" == "triangular" | 	///
							"`kernel'" == "trian" | 	///
							"`kernel'" == "tri" {
						qui replace `weigh' = (1 - abs((`rankstat' - `r')/`h'))*(abs((`rankstat' - `r')/`h')<=1) if `touse'
					}
					else {
						dis as err "Error specifying kernel function"
						error 198
					}				
					
					qui _regress `lhs' `endog' `exog' [w = `weigh'] if `touse', `noconstant'
					
					local numbwremain : word count `ferest()'
					local hiter `= `numbw' - `numbwremain''
					capture mata: betar_s`s'_bw`hiter' = betar_s`s'_bw`hiter' \ st_matrix("e(b)")
					if _rc!=0 {
							mata: betar_s`s'_bw`hiter' = st_matrix("e(b)")
					}				
				}
			}
		}

		** Beta, subset Beta(R), and bandwidth Beta(h) computations
		local vnames: coln e(b)	
		qui mata: vnames = st_local("vnames")
		qui forvalues s = 1(1)`nsubsets' {
			qui forvalues h = 1(1)`numbw' {			
				if `numbw'==1 {
					capture mata: betasub = (betasub,mean(betar_s`s'_bw`h')')
					if _rc!=0 {
						mata: betasub = mean(betar_s`s'_bw`h')'
					}				
				}											
				else {				
					capture mata: beta_s`s' = (beta_s`s',mean(betar_s`s'_bw`h')')
					if _rc!=0 {
						mata: beta_s`s' = mean(betar_s`s'_bw`h')'
					}
					capture mata: beta_bw`h' = (beta_bw`h',mean(betar_s`s'_bw`h')')
					if _rc!=0 {
						mata: beta_bw`h' = mean(betar_s`s'_bw`h')'
					}
				}
			}			
		}	
		if `numbw'==1 {
			qui mata: st_matrix("`beta'",rowsum(betasub:*subsweight)') 
			qui mata: st_matrix("`bsub'",betasub')
			*mata: st_local("bsubcoln",invtokens(J(1,cols(intvect),"R")+strofreal(1..cols(intvect))))
						
			qui matrix colnames `beta' = `vnames'			
			qui matrix colnames `beta' = Estimates:
			qui matrix colnames `bsub' = `vnames'
			
			if `subsetcoef'==1 {
				qui forvalues i = 1(1)`nsubsets' {
					matrix `tempbeta' = `bsub'[`i',1...]
					matrix colnames `tempbeta' = subset`i':
					matrix `beta' = `beta' , `tempbeta'			
				}
			}

			ereturn post `beta', esample(`touse') obs(`e(N)') depname("`lhs'")	
			
			*qui matrix rownames `bsub' = `vnames'
			*qui matrix colnames `bsub' = `bsubcoln'
			*ereturn matrix bsub = `bsub'
			
		}	
		else {
			** Computes beta(h)
			qui forvalues h = 1(1)`numbw' {				
				capture mata: betabw = (betabw,rowsum(beta_bw`h':*subsweight))
				if _rc!=0 {
					mata betabw = rowsum(beta_bw`h':*subsweight)
				}						
			}
			qui mata: st_matrix("`bbw'",betabw')
			
			qui matrix colnames `bbw' = `vnames'
			forvalues i = 1(1)`numbw' {
				matrix `tempbeta' = `bbw'[`i',1...]
				matrix colnames `tempbeta' = bw`i':
				if `i'==1 {
					matrix `beta' = `tempbeta'
				}
				else {
					matrix `beta' = `beta' , `tempbeta'
				}			
			}			
									
			if `subsetcoef'==1 {
				** Compute beta(h,R)
				qui forvalues s = 1(1)`nsubsets' {					
					mata: st_matrix("`bbw_s`s''",beta_s`s'')					
					matrix colnames `bbw_s`s'' = `vnames'
					qui forvalues i = 1(1)`numbw' {
						matrix `tempbeta' = `bbw_s`s''[`i',1...]
						matrix colnames `tempbeta' = bw`i'_subset`s':
						matrix `beta' = `beta' , `tempbeta'					
					}					
				}				
			}			

			ereturn post `beta', esample(`touse') obs(`e(N)') depname("`lhs'")	
		}		
	


				
end
	

mata: 
	function gridparse(integralstr)
	{	
		external intvect
		intvect = select(tokens(strtrim(integralstr),","),strmatch(tokens(strtrim(integralstr),","),","):==0)
		if (intvect[1,cols(intvect)]=="report") {
			external report
			report = 1
			st_local("subsetcoef",strofreal(1))
			intvect = intvect[1,1::(cols(intvect)-1)]
		}
		else {
			external report
			report = 0
			st_local("subsetcoef",strofreal(0))
		}
		st_local("nsubsets",strofreal(cols(intvect)))		
		subsupper = substr(intvect,strpos(intvect,")")+J(1,cols(intvect),1),strlen(intvect)-strpos(intvect,")"))
		subslower = substr(intvect,J(1,cols(intvect),1),strpos(intvect,"(")-J(1,cols(intvect),1))
		external lebsubm 
		external lebm 
		external subsweight
		lebsubm = (strtoreal(subsupper) - strtoreal(subslower))
		lebm = sum(lebsubm)
		subsweight = lebsubm:/lebm
		gstep = strtoreal(substr(intvect,strpos(intvect,"(")+J(1,cols(intvect),1),strpos(intvect,")")-strpos(intvect,"(")-J(1,cols(intvect),1)))
		numlistsize = lebsubm:/gstep
		gridflag = 0
		ascendflag = J(1,cols(intvect)-1,0)
		if (any(numlistsize :>= 1600)) {
			gridflag = 1
			bigg = (numlistsize :>= 1600)
			minindex(bigg,-sum(bigg),biggindex=.,biggmult=.)
			gstep[sort(biggindex,1)'] = lebsubm[sort(biggindex,1)']:/1599
			gstep = strofreal(gstep)
			intvect = subslower + J(1,cols(intvect),"(") + gstep + J(1,cols(intvect),")") + subsupper
		}
		for(i=1;i<=cols(intvect);i++) {
			if (i<cols(intvect)) ascendflag[1,i] = (strtoreal(subsupper[1,i])>=strtoreal(subslower[1,i+1]))
			st_local("int" + strofreal(i),intvect[1,i])
		}
		st_local("gridflag",strofreal(gridflag))
		st_local("ascendflag",strofreal(sum(ascendflag)))
	}
end

