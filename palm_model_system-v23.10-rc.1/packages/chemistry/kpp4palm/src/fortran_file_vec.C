
// ############################################################################
//
//     create_mz_kpp_module                       
//
//     create vectorcode from .90 sources created by KPP to be used in MESSy
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################
//
//Current revisions:
//------------------
//
//
//Former revisions:
//-----------------------
//$Id: fortran_file_vec.C 3327 2018-10-09 19:55:00Z forkel $
// Added vector switch for INTEGER         :: j,k                       (08.03.2019, forkel)
//
// Line 112: do k=is,ie bydo k=1,vl; line 156 ff: replaced index k by j (18.09.2018, ketelsen)
//
// added edit_Initialize, changed loop direction for update_rconst in edit_Update_RCONST (Sept. 2018, forkel)
//
// initial version       (Nov. 2016, ketelsen)
//

#include <fstream>
#include "fortran_file.h"

#include "utils.h"


void fortran_file::edit_inc_vec (vector<string> &gvl) {
  vector<program_line>::iterator     ip;
  vector<string>::iterator           ig;

  cout << "Handling include: " <<name <<endl;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    if(ip->get_token(0) == "REAL")  {
      int pos = ip->get_token_number_from_string("::");
      if(pos > 1) { 
        pos ++;
        for (ig=gvl.begin(); ig != gvl.end(); ig++) {
          string var_name = *ig;
          if(ip->get_token(pos) == var_name)  {
            ip->update_token(pos+1,"(VL_DIM,"); 
          }
        }
      }
    }
    if(ip->get_token(0) == "EQUIVALENCE" && ip->get_token(2) == "C" )  {
      ip->update_token(3,"(1,");
      ip->update_token(9,"(1,");
    }
  }

  return;
}

void fortran_file::global_variables2vector (vector<string> &gvl) {
  vector<program_line>::iterator     ip;
  vector<string>::iterator           ig;

  cout << "Handling subroutine: " <<name <<endl;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    for (ig=gvl.begin(); ig != gvl.end(); ig++) {
      string var_name = *ig;
      ip->change_variable_to_vector (var_name);
    }
  }

  return;
}

void fortran_file::edit_Update_RCONST (vector <Vvar> &var_list) { 
  vector<program_line>::iterator     ip;
  vector<Vvar>::iterator             iv;
  string                             lo_line; 

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    ip->global_substitute ("*"," *");
    ip->global_substitute ("* *","**");
    ip->global_substitute (","," , ");
    ip->global_substitute ("/"," / ");
    ip->global_substitute ("1:VL","j");
  }
  for (ip=pline.begin(); ip != pline.end(); ip++) {
    for (iv=var_list.begin(); iv != var_list.end(); iv++) {
      ip->change_variable_to_vector_g (*iv);
    }
  }

  for (ip=pline.begin(); ip != pline.end(); ip++) {
      if(kpp_switches.is_vector() ) {
    	  ip->global_substitute ("phot (","phot(j,");
      }
  }

  ip = pline.begin()+1 ;
  lo_line = ip->get_line() ;
  lo_line.erase();
  if(kpp_switches.is_vector() ) {
    lo_line = " INTEGER         :: j,k";
  } else {
    lo_line = " INTEGER         :: k";
  }
  ip->set_line(lo_line);
  
  if(kpp_switches.is_vector() ) {
    ip = pline.begin()+2 ;
    lo_line = ip->get_line() ;
    lo_line.erase();
    lo_line = " do k=1,vl";
    ip->set_line(lo_line);

    ip = pline.begin()+3 ;
    lo_line = ip->get_line() ;
    lo_line.erase();
    lo_line = "  j = k";
    ip->set_line(lo_line);

    ip = pline.end()-2 ;
    lo_line = ip->get_line() ;
    lo_line.erase();
    lo_line = " end do";
    ip->set_line(lo_line);
  } else {
    ip = pline.begin()+3 ;
    lo_line = ip->get_line() ;
    lo_line.erase();
    lo_line = "  k = is";
    ip->set_line(lo_line);
  }

  return;
}

void fortran_file::edit_Initialize (vector <Vvar> &var_list) {
  vector<program_line>::iterator     ip;
  vector<Vvar>::iterator             iv;
  string                             lo_line;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    ip->global_substitute ("*"," *");
    ip->global_substitute ("* *","**");
    ip->global_substitute (","," , ");
    ip->global_substitute ("/"," / ");
    ip->global_substitute ("1:VL","j");
  }
  for (ip=pline.begin(); ip != pline.end(); ip++) {
    for (iv=var_list.begin(); iv != var_list.end(); iv++) {
      ip->change_variable_to_vector_g (*iv);
    }
  }

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    if(kpp_switches.is_vector() ) {
          ip->global_substitute ("qvap","qvap(j)");
    }
  }
  for (ip=pline.begin(); ip != pline.end(); ip++) {
    if(kpp_switches.is_vector() ) {
          ip->global_substitute ("fakt","fakt(j)");
    }
  }

  ip = pline.begin()+3 ;
  lo_line = ip->get_line() ;
  lo_line.erase();
  if(kpp_switches.is_vector() ) {
    lo_line = " INTEGER         :: j,k";
  } else {
    lo_line = " INTEGER         :: k";
  }
  ip->set_line(lo_line);

  if(kpp_switches.is_vector() ) {
    ip = pline.begin()+9 ;
    lo_line = ip->get_line() ;
    lo_line.erase();
    lo_line = " do k = is,ie";
    ip->set_line(lo_line);

    ip = pline.begin()+10 ;
    lo_line = ip->get_line() ;
    lo_line.erase();
    lo_line = "  j = k - is +1";
    ip->set_line(lo_line);

    ip = pline.end()-2 ;
    lo_line = ip->get_line() ;
//  lo_line.erase();
    lo_line = " end do";
    ip->set_line(lo_line);
  } else {
    ip = pline.begin()+7 ;
    lo_line = ip->get_line() ;
    lo_line.erase();
    lo_line = "  k = is";
    ip->set_line(lo_line);
  }

  return;
}


void fortran_file::edit_KppDecomp () {

  vector<program_line>::iterator     ip;
  string                             lo_line;

  bool declaration = true;

  for (ip=pline.begin(); ip != pline.end(); ip++) {

    if(ip->get_token(0) == "IER") {
      declaration = false;
    }
    if ( declaration ) {
      if(ip->get_token(0) == "REAL") {
        lo_line = ip->get_line() ;
        lo_line.erase();
         if(kpp_switches.de_indexing () == 2) {
           lo_line = "      REAL(kind=dp) :: JVS(:,:), a(VL)";
         } else {
           lo_line = "      REAL(kind=dp) :: JVS(:,:), W(VL,NVAR), a(VL)";
         }
        ip->set_line(lo_line);
      }

    } else {
      ip->change_variable_to_vector ("W");
      ip->change_variable_to_vector ("JVS");

      if(ip->get_token(0) == "IF" || ip->get_token(1) == "IF") {
        lo_line = ip->get_line() ;
        lo_line.insert(0,"! Not in vector Mode ");
        ip->set_line(lo_line);
      }
      if(ip->get_token(0) == "IER" && ip->get_token(2) == "k") {
        lo_line = ip->get_line() ;
        lo_line.insert(0,"! Not in vector Mode ");
        ip->set_line(lo_line);
      }
      if(ip->get_token(0) == "RETURN") {
        lo_line = ip->get_line() ;
        lo_line.insert(0,"! Not in vector Mode ");
        ip->set_line(lo_line);
      }
    }
  }

  return;
}

void fortran_file::edit_KppSolve () {

  vector<program_line>::iterator     ip;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    if(ip->get_token(0) == "REAL") {
      ip->substitute("NVAR",":,:");
      ip->substitute("LU_NONZERO",":,:");
    } else {
      ip->change_variable_to_vector ("JVS");
      ip->change_variable_to_vector ("X");
    }
  }

  return;
}

void fortran_file::edit_Jac_SP () {

  vector<program_line>::iterator     ip;
  string                             lo_line;

  for (ip=pline.begin(); ip != pline.end(); ip++) {
    if(ip->get_token(0) == "REAL") {
      ip->substitute("NVAR",":,:");
      ip->substitute("NFIX",":,:");
      ip->substitute("NREACT",":,:");
      ip->substitute("LU_NONZERO",":,:");
      if(ip->get_token(5) == "B") {
        lo_line = ip->get_line() ;
        lo_line.erase();
        lo_line = " REAL(kind=dp) :: B(VL," +ip->get_token(7)  +")";
        ip->set_line(lo_line);
      }
    } else {
      ip->change_variable_to_vector ("V");
      ip->change_variable_to_vector ("F");
      ip->change_variable_to_vector ("B");
      ip->change_variable_to_vector ("RCT");
      ip->change_variable_to_vector ("JVS");
    }
  }

  return;
}

void fortran_file::edit_Fun () {

  vector<program_line>::iterator     ip;

  bool declaration = true;

  for (ip=pline.begin(); ip != pline.end(); ip++) {

    if(ip->get_token(1) == "Computation" || ip->get_token(1) == "Told" ) {
      declaration = false;
    }
    if ( declaration ) {
      if(ip->get_token(0) == "REAL") {
        ip->substitute("NVAR",":,:");
        ip->substitute("NFIX",":,:");
        ip->substitute("NREACT",":,:");
      }
    } else {
      ip->change_variable_to_vector ("V");
      ip->change_variable_to_vector ("F");
      ip->change_variable_to_vector ("RCT");
      ip->change_variable_to_vector ("Vdot");
    }
  }

  return;
}

void fortran_file::edit_WAXPY () {

  vector<program_line>::iterator     ip;

  cout << "Handling subroutine: WAXPY" <<endl;

  for (ip=pline.begin(); ip != pline.end(); ip++) {

    if(ip->get_token(0) == "REAL") {
       ip->substitute("N",":,:");
       ip->substitute("N",":,:");
       ip->substitute("Alpha","Alpha(:)");
    } else {
        ip->change_variable_to_vector ("Y");
        ip->change_variable_to_vector ("X");
        if(ip->get_token(0) != "SUBROUTINE") {
        	if(ip->get_token(0) == "IF")  {
        		ip->substitute("Alpha","SUM(alpha(1:VL))");
        	} else {
        		ip->substitute("Alpha","alpha(1:VL)");
        	}
        }
    }


  }

  return;
}

void fortran_file::edit_FunTemplate () {

  vector<program_line>::iterator     ip;

  cout << "Handling subroutine: FunTemplate" <<endl;

  for (ip=pline.begin(); ip != pline.end(); ip++) {

    if(ip->get_token(0) == "REAL") {
       ip->substitute("NVAR",":,:");
       ip->substitute("T,","T(:),");
       ip->substitute("Told","Told(size(T)),Time(size(T))");
    }


  }

  return;
}

void fortran_file::edit_JacTemplate () {

  vector<program_line>::iterator     ip;

  cout << "Handling subroutine: JacTemplate" <<endl;

  for (ip=pline.begin(); ip != pline.end(); ip++) {

    if(ip->get_token(0) == "REAL") {
        ip->substitute("NVAR",":,:");
        ip->substitute("T,","T(:),");
        ip->substitute("Told","Told(size(T)),Time(size(T))");
        ip->substitute("LU_NONZERO",":,:");
    }


  }

  return;
}

