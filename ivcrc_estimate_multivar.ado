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
				  Dendog(varlist)					///
				  VARCOEF(string)					/// 
				  Kernel(string) *]

				  
		_iv_parse `0'
        
        local lhs `s(lhs)'
		local endog `s(endog)'
        local exog `s(exog)'
        local inst `s(inst)'
        local 0 `s(zero)'	
		
		tokenize `endog'
		local k : word count `endog'
		forvalues v = 1(1)`k' {
			tempvar endog`v'
			gen `endog`v'' = ``v''
			tempvar rankstat`v' 
		}
		tempvar weigh subsampr
		tempname beta betasub betar betasubr 		  

		marksample touse
        markout `touse' `lhs' `exog' `inst' `endog' `dendog'
		qui count if `touse'
		ereturn scalar N = r(N)
		
		local kercount : word count `kernel'
		if `kercount' > 1 {
			display as err "Only one kernel may be specified"
			error 198
		}
		
		
		
		
		if trim("`varcoef'")!="" {
		
		if trim("`integral'") == "" {
			dis as err "Option integral() must be specified with varcoef()"
		}
		if trim("`bandwidth'") == "" {
			dis as err "Option bandwidth() must be specified with varcoef()"
		}
		** Define integral set & subsets (create upper and lower cutoffs for R)
		qui if trim("`integral'") != "" {
			tokenize "`integral'", parse(",")
			numlist "`1'"
			local minr `r(numlist)'
			numlist "`3'"
			local maxr `r(numlist)'
		}
		
		local numbw : word count `bandwidth'
		tokenize `varcoef'
		local k : word count `varcoef'
		** Loop over covariates `v' to compute `rankstat`v'' variables for kernel weights 
		qui gen `subsampr' = `touse'
		qui forvalues v = 1(1)`k' {		
			qui egen `rankstat`v'' = ``v''
			if "`integral'"!="" {
				qui replace `subsampr' = ((`rankstat`v'' >= `minr') & (`rankstat`v'' <= `maxr')) if `touse' 
			}
		}
		
		}
		
		
		
		if trim("`varcoef'")=="" {	
		
		
		** Define integral set & subsets (create upper and lower cutoffs for R)
		qui if trim("`integral'") != "" {
			tokenize "`integral'", parse(",")
			numlist "`1'"
			local minr `r(numlist)'
			numlist "`3'"
			local maxr `r(numlist)'
		}
		
		** Bandwidth specification
		if trim("`bandwidth'") == "" {
			local bandwidth .05
		}
		local numbw : word count `bandwidth'
		
		** Loop over basic endogenous vars `v' to compute `rankstat`v'' rank statistics 
		qui gen `subsampr' = `touse'
		qui forvalues v = 1(1)`k' {		
			** Conditional rank statistic computation
			qui forval s = 1(1)`= `quantiles' - 1' {
				local ss = `s'/`quantiles'
				qui qreg `endog`v'' `inst' `exog' if `touse', quantile(`ss') 
				tempvar q_`v'_`s' ind_`v'_`s'
				qui predict `q_`v'_`s'' if `touse'
				qui gen `ind_`v'_`s'' = `q_`v'_`s''<=`endog`v'' if `touse'
				local ind_list_`v' `ind_list_`v'' `ind_`v'_`s''
			}
			qui egen `rankstat`v'' = rowmean(`ind_list_`v'') if `touse'
			if "`integral'"!="" {
				qui replace `subsampr' = ((`rankstat`v'' >= `minr') & (`rankstat`v'' <= `maxr')) if `touse' 
			}
		}	
				
		tokenize "`generate'", parse(",")
		if "`generate'" != "" & trim("`bootstrap'") == "" {
			if "`2'" == "" {
				forvalues v = 1(1)`k' {
					qui gen `1'_`v' = `rankstat`v'' if `touse'
				}
				display "Conditional rank statistic has been saved"
 			}
			else if "`2'" == "," & "`3'" == "" {
				forvalues v = 1(1)`k' {
					qui gen `1'_`v' = `rankstat`v'' if `touse'
				}
				display "Conditional rank statistic has been saved"
			}
			else if "`2'" == "," & ("`3'" == "replace" | "`3'" == "rep" | "`3'" == "repl" | "`3'" == "re" | "`3'" == "r") {
				forvalues v = 1(1)`k' {
					qui replace `1'_`v' = `rankstat`v'' if `touse'
				}
				display "Conditional rank statistic has been saved"
			} 
			else {
				display as err "Syntax error with generate()"
				error 198
			}			
		}
				
		}
		
		** Product kernel and beta(r) at each r[i] in the empirical support of `rankstat`v'' 
		qui gen `weigh' = 1 if `touse'
		qui count if `touse'
		local N = `r(N)'
		qui foreach h of local bandwidth { 
		qui forvalues i = 1(1)`N' {
			 
					**Loop over basic edogenous vars to construct product kernel;
					qui forvalues v = 1(1)`k' { 					
					if  "`kernel'" == "" |			/// 
						"`kernel'" == "uniform" |	///
						"`kernel'" == "unif" |	 	///
						"`kernel'" == "un" | 		///
						"`kernel'" == "rectangle" |	///
						"`kernel'" == "rec"	{
						qui replace `weigh' = `weigh'*(1/2)*(abs((`rankstat`v'' - `rankstat`v''[`i'])/`h')<=1) if `touse'
					}
					else if "`kernel'" == "biweight" | 	///
							"`kernel'" == "biweigh" | 	///
							"`kernel'" == "biw" | 		///
							"`kernel'" == "bi" {
						qui replace `weigh' = `weigh'*(15/16)*(1 - ((`rankstat`v'' - `rankstat`v''[`i'])/`h')^2)^2*(abs((`rankstat`v'' - `rankstat`v''[`i'])/`h')<=1) if `touse'
					}
					else if "`kernel'" == "triweight" | 		///
							"`kernel'" == "triweigh" | 			///
							"`kernel'" == "triw" {
						qui replace `weigh' = `weigh'*(35/32)*(1 - ((`rankstat`v'' - `rankstat`v''[`i'])/`h')^2)^3*(abs((`rankstat`v'' - `rankstat`v''[`i'])/`h')<=1) if `touse'
					}
					else if "`kernel'" == "cosine" |		///
							"`kernel'" == "cos" {
						qui replace `weigh' = `weigh'*(_pi/4)*cos((_pi/2)*(`rankstat`v'' - `rankstat`v''[`i'])/`h')*(abs((`rankstat`v'' - `rankstat`v''[`i'])/`h')<=1) if `touse'
					}
					else if "`kernel'" == "epanechnikov" |		///
							"`kernel'" == "epanech" |			///
							"`kernel'" == "epan" |				///
							"`kernel'" == "ep" {
						qui replace `weigh' = `weigh'*(3/4)*(1 - ((`rankstat`v'' - `rankstat`v''[`i'])/`h')^2)*(abs((`rankstat`v'' - `rankstat`v''[`i'])/`h')<=1) if `touse'
					}
					else if "`kernel'" == "gaussian" | 		///
							"`kernel'" == "gauss" | 		///
							"`kernel'" == "gaus" | 			///
							"`kernel'" == "gau" {
						qui replace `weigh' = `weigh'*(1/sqrt(2*_pi))*exp(-(1/2)*((`rankstat`v'' - `rankstat`v''[`i'])/`h')^2) if `touse'
					}
					else if "`kernel'" == "triangle" | 	///
							"`kernel'" == "triangular" | 	///
							"`kernel'" == "trian" | 	///
							"`kernel'" == "tri" {
						qui replace `weigh' = `weigh'*(1 - abs((`rankstat`v'' - `rankstat`v''[`i'])/`h'))*(abs((`rankstat`v'' - `rankstat`v''[`i'])/`h')<=1) if `touse'
					}
					else {
						dis as err "Error specifying kernel function"
						error 198
					}					
					}  
					
					qui _regress `lhs' `endog' `dendog' `exog' [w = `weigh'] if `touse', `noconstant' 

					capture mata: betar = betar \ st_matrix("e(b)") 
					if _rc!=0 { 
							mata: betar = st_matrix("e(b)") 
					}	
					if "`integral'"!="" & `subsampr'[`i']==1 {
						capture mata: betasubr = betasubr \ st_matrix("e(b)") 
						if _rc!=0 { 
							mata: betasubr = st_matrix("e(b)") 
						}
					}	
					replace `weigh' = 1 											
		}
		}
		

		** Beta, and subsample Beta(R) from integral option 
		local vnames: coln e(b)	
		qui mata: vnames = st_local("vnames")
		qui mata: st_matrix("`beta'",mean(betar))			
		qui matrix colnames `beta' = `vnames'
		qui matrix colnames `beta' = Estimates: 

		if "`integral'"!=""  {
			qui mata: st_matrix("`betasub'",mean(betasubr))	
			qui matrix colnames `betasub' = `vnames'
			qui matrix colnames `betasub' = Subset:
			qui matrix `beta' = `beta' , `betasub'
		}
		ereturn post `beta', esample(`touse') obs(`e(N)') depname("`lhs'")	
			
			
end
	
