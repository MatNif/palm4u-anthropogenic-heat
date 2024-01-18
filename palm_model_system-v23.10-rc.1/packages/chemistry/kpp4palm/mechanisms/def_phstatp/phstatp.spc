{phstatp.spc
Former revisions
----------------
 $Id: smog.spc 2459 2017-09-13 14:10:33Z forkel $
}
#include atoms

 #DEFVAR  
    O       = O ;      		{oxygen atomic ground state (3P)}
    O3		= 3O ;          {ozone}  
    NO		= N + O ;       {nitric oxide}  
    NO2		= N + 2O ;      {nitrogen dioxide} 
    PM10 	= ignore ;      {PM10} 

#DEFFIX
    H2O		= H + 2O ;      {water}
    H2		= 2H ;          {molecular hydrogen}
    O2      = 2O ;          {molecular oxygen}              
    N2      = 2N ;          {molecular nitrogen}              

