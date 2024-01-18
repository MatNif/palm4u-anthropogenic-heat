{phstat.spc
Former revisions
----------------
 $Id: salsa+phstat.spc 2459 2017-09-13 14:10:33Z forkel $
}
#include atoms

 #DEFVAR  
    O       = O ;      		{oxygen atomic ground state (3P)}
    O3		= 3O ;          {ozone}  
    NO		= N + O ;       {nitric oxide}  
    NO2		= N + 2O ;      {nitrogen dioxide} 
    HNO3    = H + N + 3O ;  { nitric acid }
    H2SO4   = 2H + S +4O ;  {sulfuric acid}
    NH3     = 3H + N ;      {ammonia}
    OCNV    = ignore ;      {non-volatile OC}
    OCSV    = ignore ;      {semi-volatile OC}


#DEFFIX
    H2O		= H + 2O ;      {water}
    H2		= 2H ;          {molecular hydrogen}
    O2      = 2O ;          {molecular oxygen}              
    N2      = 2N ;          {molecular nitrogen}              

