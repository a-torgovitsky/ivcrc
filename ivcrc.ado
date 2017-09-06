capture program drop ivcrc 

program ivcrc, eclass
		version 12

		syntax [anything] [if] [in]  				///
                [,NOConstant						///
				  GENerate(string)					///
				  INTegral(string)					///				  
				  QUantiles(integer 100)			///
				  BANDwidth(numlist)				///
				  BOOTstrap(string)					///
				  NORMal							///		
  				  DENDOG(varlist)					///
				  VARCOEF(string)					/// 
				  Kernel(string) *]

			
		** Estimation

		if (regexm("`0'","boot")==0 & regexm("`0'","bootstrap")==0) {
			display in gr "Default setting is no standard errors"
			ivcrc_estimate `0'
			ereturn display			
		}		
		else {
			if "`generate'" != "" & "`bootstrap'" != "" {
				display in gr "Generate option ignored while bootstrapping..." 
			}
			
			qui bootstrap , `bootstrap' : ivcrc_estimate `0'
			estat bootstrap , p `normal' `bc' `bca'
					
		}

end		
			
		
