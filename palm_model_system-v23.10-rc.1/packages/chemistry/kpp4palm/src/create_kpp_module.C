
// ############################################################################
//
//     create_kpp_module 
//
//     create scalar code from .f90 sources created by KPP
//
//     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
//
// ############################################################################
//
//
//Former revisions:
//-----------------
// Update of copyright statements in templates/module_header              (12.01.2021, forkel)
//
// alphabetic ordering of PUBLIC statments in templates/module_header     (28.03.2019, forkel)
//
// renamed get_mechanismname to get_mechanism_name                        (26.03.2019, forkel)
//
//
// Deleted $Id since document_changes does not work for C and C++         (15.03.2019, forkel)
//
// OpenMP version    (15.03.2019, ketelsen)
//
// Added vector switch Kacc,Krej,IERRV, Commented add_line for istatf,    (05.03.2019, forkel)
//      added ,pe after ierr_u,         
//
// Added create_set_cs and cs_mech and get_mechanismname in module_header (05.03.2019, forkel)
//
// exclude kco_compress from handling by global_variables2vector (30.10.2018, forkel)
//
// Added  automatic line with mechanism name (read mech_list) (25.09.2018, forkel)
//
// Added  vl_glo = size(tempi,1) (20.09.2018, forkel) 
//
// Removed creation of fill_ Subroutine and creation of calls thereof (18.09.2018, ketelsen)
//
// Fix in order not to loose the values of qvap and fakt (12.09.2018, forkel)
//
// Bug fixes: moved kppi.add_line("    CALL initialize after fakt = fakti(is)
// Deleted definition of qvap,fakt in create_kpp_integrate again (03.09.2018, forkel)
//
// Changes for vector mode (edit_WAXPY, edit_FunTemplate, edit_JacTemplate, 
// some cleanup of comments, various changes in create_kpp_integrate) (July 2018, ketelsen)
//
//
// Added qvap and fakt                            (June 2018, forkel)
// --> Change in module_header: qvap, fakt added  (June 2018, forkel)
//
// re-established original uppercase/lowercase     (June 2018, forkel)
// --> Change in module_header: reset case in  Initialize, Integrate, and 
//     Update_rconst                               (June 2018, forkel)
//
// Removed preprocessor directive __chem again (2017-09-14, forkel)
// 
// Added phot                                                 (2017-09-14, forkel)
// --> Change in module_header: Variables for photolyis added (2017-09-14, forkel)
//
// change of some output to lowercase with uppercase Fortran (2017, forkel)
//
// Intial version of KP4 adapted to PALM                      (Nov. 2016, ketelsen)
//
#include <stdio.h>
// stdlib is necessary to define getenv:
#include <stdlib.h>

#include "create_kpp_module.h"
#include "utils.h"

void create_kpp_module::do_work (string s) {
   vector<fortran_file>::iterator  it;
   vector<string>::iterator        ic;
   vector<Vvar>::iterator          iv;

   expand_decomp                   exp_de;

   prefix = s;
   module_name = prefix;

   cout << "Create " << module_name << " from kpp Fortran sources" <<endl;
   cout << "Vector mode " << kpp_switches.is_vector() <<endl;
   cout << "De_indexing " << kpp_switches.de_indexing() <<endl;

   create_fortran_files_and_read();

// Generate first module lines 

     string first_line="MODULE " + module_name;
   mz_kpp.add_line(first_line);
   mz_kpp.add_line(" ");

//    string e5_line = first_line +"_e5";
//    e5_kpp.add_line(e5_line);
//    e5_line = "  USE             " + module_name;
//    e5_kpp.add_line(e5_line);
//    e5_kpp.add_line(" ");

// edit include files

   for(it=kpp_includes.begin();it!=kpp_includes.end();it++) {
     it->edit_inc(header_variables);

//   Create variable Species list and vector variable list

     if(it->get_name() == module_name + "_Parameters") {
       it->create_species_list(species_list);
     }
     if(it->get_name() == module_name + "_Global") {
       it->vector_variable_list(Vvar_list);
     }
   }

// Prepare expansion of decomposition subroutine

   if(kpp_switches.de_indexing () > 0 ) {
     exp_de.create_sparse_info (kpp_includes, module_name);
   }

// edit FORTRAN files

   for(it=kpp_files.begin();it!=kpp_files.end();it++) {
     it->edit_fortran ();
   }

// Generate a list of single subroutines from kpp-files
// kpp files are modules containing several subroutines

   copy_files_to_subroutines ();

// All header_variables to include list

   kpp_includes.push_back(header_variables);

// Create decomposition subroutine
   if(kpp_switches.de_indexing () > 0 ) {
     exp_de.create_routine (kpp_subroutines);
   }

   if(kpp_switches.is_vector()) {

   cout << "##### Hier kpp_switches.is_vector          " <<endl;
//   Change header section
     for(it=kpp_includes.begin();it!=kpp_includes.end();it++) {
       it->edit_inc_vec(global_variable_list);
     }

//   Change global variables to vector (except for kp4_compress, which has already the right form)
     
     for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
       if(it->get_name() != "kco_compress" ) {
         it->global_variables2vector (global_variable_list);
       }
     }

//   Edit individual subroutines 

     for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
       if(it->get_name() == "KppDecomp") {
         it->edit_KppDecomp();
       }
       if(it->get_name() == "KppSolve") {
         it->edit_KppSolve();
       }
       if(it->get_name() == "Jac_SP" ) {
         it->edit_Jac_SP();
       }
       if(it->get_name() == "Fun" ) {
         it->edit_Fun();
       }
       if(it->get_name() == "WAXPY" ) {
         it->edit_WAXPY();
       }
       if(it->get_name() == "FunTemplate" ) {
         it->edit_FunTemplate();
       }
       if(it->get_name() == "JacTemplate" ) {
         it->edit_JacTemplate();
       }
     }
   }

// Update_RCONST has to be changed also in scalar mode

   for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
     if(it->get_name() == "Update_RCONST") {
       it->edit_Update_RCONST(Vvar_list);
     }

     if(it->get_name() == "Initialize") {
       it->edit_Initialize(Vvar_list);
     }

   }

// Add Solver template to subroutine list
   if(kpp_switches.is_vector()) {
     add_solver_to_subroutine_list ();
   }

// The module header will be taken from ../templates/module_header.
// Please edit if header has to be changed.

   generate_module_header();

// Create subroutine to communicate mechanism name to other modules
   create_set_cs();

// Create kpp_integrate subroutine (chem_gasphase_integrate) for skalar and vector mode

   create_kpp_integrate();
// Copy include files

   for(it=kpp_includes.begin();it!=kpp_includes.end();it++) {
     it->copy_to_MZ_KPP(mz_kpp);
   }

   mz_kpp.add_line(" ");
   mz_kpp.add_line("! Interface Block ");
   mz_kpp.add_line(" ");
   for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
     string          buf;

     string prefix = "  ";
     for(ic=interface_ignore.begin();ic!=interface_ignore.end();ic++) {
       if(it->get_name() == *ic) {
         prefix = "!interface not working  ";
         break;
       }
     }

     buf = prefix + "interface            " + it->get_name() ;
     mz_kpp.add_line(buf);
     buf = prefix + "  module procedure   " + it->get_name();
     mz_kpp.add_line(buf);
     buf = prefix + "end interface        " + it->get_name();
     mz_kpp.add_line(buf);
     mz_kpp.add_line(" ");
   }

   mz_kpp.add_line(" ");

// Declare variables THREADPRIVATE for OpenMP version

   mz_kpp.add_line("  ! OpenMP directives generated by kp4 ");
   mz_kpp.add_line(" ");
   mz_kpp.add_line("  !$OMP THREADPRIVATE (vl,vl_glo,is,ie,data_loaded)");
   mz_kpp.add_line("  !$OMP THREADPRIVATE (c,var,fix,rconst,time,temp,stepmin,cfactor)");
   mz_kpp.add_line("  !$OMP THREADPRIVATE (qvap,fakt,cs_mech,a,icntrl,rcntrl)");
   mz_kpp.add_line(" ");
   if(kpp_switches.is_vector()) {
      mz_kpp.add_line("  ! Vector mode Only ");
	  mz_kpp.add_line(" ");
	  mz_kpp.add_line("  !$OMP THREADPRIVATE (kacc,krej,ierrv)");
	  mz_kpp.add_line("  !$OMP THREADPRIVATE (kpoints,kpoints_SAVE,index_org,done_check,index_step,cell_done)");
	  mz_kpp.add_line("  !$OMP THREADPRIVATE (f_done,kacc_done,krej_done,ierr_done,compress_done)");
	  mz_kpp.add_line(" ");
   }


// Copy FORTRAN subroutines to mz_kpp

   mz_kpp.add_line(" CONTAINS");
   
   for(it=kpp_subroutines.begin();it!=kpp_subroutines.end();it++) {
     mz_kpp.add_line(" ");
     it->copy_to_MZ_KPP(mz_kpp);
   }

// Finish module

   string last_line="end module " + module_name;
   mz_kpp.add_line("");
   mz_kpp.add_line(last_line);

// Write the complete module to file: mz_kpp.f

   write_module_file();

   return;
}

void create_kpp_module::create_fortran_files_and_read() {

   string                          name;
   ifstream                        in,in_c,in_b,in_i;
   fortran_file                    f_file;
   vector<fortran_file>::iterator  it;

// Open file with list of FORTRAN routines

   in.open("file_list");
   if( !in ) {
      cout << "cannot open " << endl; my_abort("file_list");
   }
   
// Create kpp_fortran routines
   while ( 1 ) {
     in >> name;
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_1");
     f_file.set_name(name);
     kpp_files.push_back(f_file);
   }
   in.close();

// Read FORTRAN code

   for(it=kpp_files.begin();it!=kpp_files.end();it++) {
     it->read();
   }

// Open file with list of include files

   in_c.open("include_list");
   if( !in_c ) {
      cout << "cannot open " << endl; my_abort("include_list");
   }

// Create kpp_includes vector
   while ( 1 ) {
     in_c >> name;
     if( in_c.eof() ) break;
     if( in_c.bad() ) my_abort("ERROR_READ_3");
     f_file.set_name(name);
     kpp_includes.push_back(f_file);
   }
   in_c.close();

// Read include files

   for(it=kpp_includes.begin();it!=kpp_includes.end();it++) {
     it->read();
   }

// Read Ignore list

   in_i.open("interface_ignore_list");
   if( !in_i ) {
      cout << "cannot open " << endl; my_abort("include_list");
   }

// Create kpp_includes vector
   while ( 1 ) {
     in_i >> name;
     if( in_i.eof() ) break;
     if( in_i.bad() ) my_abort("ERROR_READ_4");
     interface_ignore.push_back(name);
   }
   in_c.close();

}

void create_kpp_module::copy_files_to_subroutines () {
   string                          name;
   ifstream                        in;
   fortran_file                    s_file;
   vector<fortran_file>::iterator  it;

// Open file with list of FORTRAN routines

   in.open("subroutine_list");
   if( !in ) {
      cout << "cannot open " << endl; my_abort("subroutine_list");
   }

// Create vector kpp_subroutines

   while ( 1 ) {
     in >> name;
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_S1");
     s_file.set_name(name);
     kpp_subroutines.push_back(s_file);
   }
   in.close();

   header_variables.add_line(" ");
   header_variables.add_line("!  variable definations from  individual module headers ");
   header_variables.add_line(" ");

//  Loop over all FORTRAN Files

   for(it=kpp_files.begin();it!=kpp_files.end();it++) {
     it->copy_to_subroutine_vector(kpp_subroutines, header_variables);
   }
}

void create_kpp_module::add_solver_to_subroutine_list () {
   fortran_file                    s_file;

   string solver_name = getenv("KPP_SOLVER");
   cout << "KPP_SOLVER " <<solver_name <<endl;
   
   s_file.set_name(solver_name);
   s_file.read();
   kpp_subroutines.push_back(s_file);

   return;
}

void create_kpp_module::generate_module_header() {

   string                          buf;
   ifstream                        in;
   ifstream                        in_e5;
   program_line                    line;
   vector<fortran_file>::iterator  it;
   char                            distr[2];
   string                          diline;

// Read mechanism from mech_list

   in.open("mech_list");
   if( !in ) {
      cout << "cannot open " << endl; my_abort("mech_list");
   }

   while ( 1 ) {
     getline (in, buf);
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_4");
     line.set_line(buf);
     mz_kpp.add_line(line);
   }
   in.close();


// Read Modul Header from file $MZ_KPP_HOME/templates/module_header

   in.open("module_header");
   if( !in ) {
      cout << "cannot open " << endl; my_abort("module_header");
   }

   while ( 1 ) {
     getline (in, buf);
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_4");
     line.set_line(buf);
     mz_kpp.add_line(line); 
   }
   mz_kpp.add_line("                                                                 "); 
   mz_kpp.add_line("! Variables used for vector mode                                 "); 
   mz_kpp.add_line("                                                                 "); 
   if(kpp_switches.is_vector()) {
       mz_kpp.add_line("  logical,parameter          :: L_VECTOR = .TRUE.             ");
   } else {
       mz_kpp.add_line("  logical,parameter          :: L_VECTOR = .FALSE.            ");
   }
//  mz_pj_20070531+
   sprintf(distr,"%i",kpp_switches.de_indexing());
   diline = distr ;
   mz_kpp.add_line("  integer,parameter          :: I_LU_DI = " + diline );
//  mz_pj_20070531-

   mz_kpp.add_line("  integer,parameter          :: VL_DIM = " 
                 + kpp_switches.get_vector_length() ); 
   mz_kpp.add_line("  integer                     :: vl                              "); 
   mz_kpp.add_line("                                                                 "); 
   mz_kpp.add_line("  integer                     :: VL_glo                          "); 
   mz_kpp.add_line("  integer                     :: is,ie                           "); 
   mz_kpp.add_line("                                                                 "); 
   mz_kpp.add_line("                                                                 "); 
   if(kpp_switches.is_vector()) {
      mz_kpp.add_line("  integer, dimension(VL_dim)   :: Kacc,Krej                       "); 
      mz_kpp.add_line("  integer, dimension(VL_dim)   :: IERRV                           "); 
   }
   mz_kpp.add_line("  logical                     :: data_loaded = .false.             "); 
   if(kpp_switches.is_vector()) {
	  mz_kpp.add_line("  REAL(dp),POINTER,DIMENSION(:,:),CONTIGUOUS    :: var           ");
   } else {
      mz_kpp.add_line("  REAL(dp),POINTER,DIMENSION(:),CONTIGUOUS    :: var             ");
   }
   in.close();

   return;
}

void create_kpp_module::write_module_file() {
   ofstream                    out;
   ofstream                    out_e5;

   string out_file  = "kk_kpp.f90";
   out.open(out_file.c_str(), ios::out);
   if( !out ) {
      cout << "cannot open " << endl; my_abort(out_file);
   }

   mz_kpp.write_file (out);

   out.close();
   

   return;
}

void create_kpp_module::create_set_cs() {
   fortran_file          kppi;          
   vector<Vvar>::iterator               iv;
   string                               xline;
     
   string                          buf;
   ifstream                        in;
   program_line                    line;

   kppi.set_name("get_mechanism_name");
   kppi.add_line("SUBROUTINE get_mechanism_name                                       ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  IMPLICIT NONE                                                     ");
// Read mechanism from set_cm
// Tis got an own own subroutine to aviod being called at each timestep

   in.open("set_cm");
   if( !in ) {
      cout << "cannot open " << endl; my_abort("set_cm");
   }

   while ( 1 ) {
     getline (in, buf);
     if( in.eof() ) break;
     if( in.bad() ) my_abort("ERROR_READ_4");
     line.set_line(buf);
     kppi.add_line(line);
   }
   in.close();

   kppi.add_line("                                                                    ");
   kppi.add_line("  return                                                            ");
   kppi.add_line("END SUBROUTINE get_mechanism_name                                   ");
   kppi.add_line("                                                                    ");
   kpp_subroutines.push_back(kppi);

   return;
}


void create_kpp_module::create_kpp_integrate() {
   fortran_file          kppi;
   vector<Vvar>::iterator               iv;
   string                               xline;


   kppi.set_name("chem_gasphase_integrate");

   kppi.add_line("SUBROUTINE chem_gasphase_integrate (time_step_len, conc, tempi, qvapi, fakti, photo, ierrf, xnacc, xnrej, istatus, l_debug, pe, icntrl_i, rcntrl_i )  ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  IMPLICIT NONE                                                     ");
   kppi.add_line("                                                                    ");

   kppi.add_line("  REAL(dp), INTENT(IN)                   :: time_step_len           ");
   kppi.add_line("  REAL(dp),  DIMENSION(:,:),  INTENT(INOUT) :: conc                    ");
   kppi.add_line("  REAL(dp),  DIMENSION(:,:),  INTENT(IN)    :: photo                   ");
   kppi.add_line("  REAL(dp),  DIMENSION(:),  INTENT(IN)      :: tempi                   ");
   kppi.add_line("  REAL(dp),  DIMENSION(:),  INTENT(IN)      :: qvapi                   ");
   kppi.add_line("  REAL(dp),  DIMENSION(:),  INTENT(IN)      :: fakti                   ");
   kppi.add_line("  INTEGER,  INTENT(OUT), OPTIONAL        :: ierrf(:)                ");
   kppi.add_line("  INTEGER,  INTENT(OUT), OPTIONAL        :: xNacc(:)                ");
   kppi.add_line("  INTEGER,  INTENT(OUT), OPTIONAL        :: xNrej(:)                ");
   kppi.add_line("  INTEGER,  INTENT(INOUT), OPTIONAL      :: istatus(:)              ");
   kppi.add_line("  INTEGER,  INTENT(IN), OPTIONAL         :: PE                      ");
   kppi.add_line("  LOGICAL,  INTENT(IN), OPTIONAL         :: l_debug                 ");
   kppi.add_line("  INTEGER,  DIMENSION(nkppctrl),INTENT(IN), OPTIONAL  :: icntrl_i         ");
   kppi.add_line("  REAL(dp), DIMENSION(nkppctrl),INTENT(IN), OPTIONAL  :: rcntrl_i         ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  INTEGER                                 :: k   ! loop variable     ");
   kppi.add_line("  REAL(dp)                                :: dt                      ");
   kppi.add_line("  integer, dimension(20)                 :: istatus_u               ");
   kppi.add_line("  integer                                :: ierr_u                  ");
// kppi.add_line("  integer                                :: istatf                  ");
   kppi.add_line("  integer                                :: vl_dim_lo               ");
   kppi.add_line("                                                                    ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  if (present (istatus) )   istatus = 0                             ");
   kppi.add_line("  if (present (icntrl_i) )  icntrl  = icntrl_i                      ");
   kppi.add_line("  if (present (rcntrl_i) )  rcntrl  = rcntrl_i                      ");
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
      kppi.add_line("  IF ( PRESENT(l_debug) .AND. PRESENT(pe) ) CONTINUE                ");
      kppi.add_line("                                                                    ");
      kppi.add_line("  var => c(:,1:nvar)                                                ");
   } else {
      kppi.add_line("  var => c(1:nvar)                                                  ");
   }
   kppi.add_line("                                                                    ");
   kppi.add_line("  vl_glo = size(tempi,1)                                            ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  vl_dim_lo = VL_DIM                                                ");
   kppi.add_line("  DO k=1,VL_glo,vl_dim_lo                                           ");
   kppi.add_line("    is = k                                                          ");
   kppi.add_line("    ie = min(k+vl_dim_lo-1,VL_glo)                                  ");
   kppi.add_line("    vl = ie-is+1                                                    ");

   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    C(1:vl,:) = Conc(is:ie,:)                                     ");
   } else {
     kppi.add_line("    C(:) = Conc(is,:)                                             ");
   }

   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    temp(1:vl) = tempi(is:ie)                                     ");
   } else {
     kppi.add_line("    temp = tempi(is)                                              ");
   }
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    qvap(1:vl) = qvapi(is:ie)                                     ");
   } else {
     kppi.add_line("    qvap = qvapi(is)                                              ");
   }
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    fakt(1:vl) = fakti(is:ie)                                     ");
   } else {
     kppi.add_line("    fakt = fakti(is)                                              ");
   }

   kppi.add_line("                                                                    ");
   kppi.add_line("    CALL initialize                                                 ");

   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    phot(1:vl,:) = photo(is:ie,:)                                     ");
   } else {
     kppi.add_line("    phot(:) = photo(is,:)                                             ");
   }
   kppi.add_line("                                                                    ");
   kppi.add_line("    CALL update_rconst                                              ");
   kppi.add_line("                                                                    ");
   kppi.add_line("    dt = time_step_len                                              ");
   kppi.add_line("                                                                    ");
   kppi.add_line("    ! integrate from t=0 to t=dt                                    ");
   kppi.add_line("    CALL integrate(0._dp, dt, icntrl, rcntrl, istatus_u = istatus_u, ierr_u=ierr_u)");
   kppi.add_line("                                                                    ");
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    Conc(is:ie,:) = C(1:VL,:)                                     ");

   } else {
     kppi.add_line("   IF ( PRESENT(l_debug) .AND. PRESENT(PE) ) THEN                       ");
     kppi.add_line("      IF ( l_debug ) CALL error_output(Conc(is,:),ierr_u,pe)            ");
     kppi.add_line("   ENDIF                                                              ");
     kppi.add_line("                                                                      ");
     kppi.add_line("    Conc(is,:) = C(:)                                                 ");
   }

   kppi.add_line("                                                                    ");
   kppi.add_line("    ! Return Diagnostic Information                                 ");
   kppi.add_line("                                                                    ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("    if ( Present(ierrf) )    ierrf(is:ie) = IERRV(1:VL)              ");
     kppi.add_line("    if ( Present(xNacc) )    xNacc(is:ie) = Kacc(1:VL)               ");
     kppi.add_line("    if ( Present(xNrej) )    xNrej(is:ie) = Krej(1:VL)               ");
   } else {
     kppi.add_line("    if ( Present(ierrf) )    ierrf(is) = IERR_U                      ");
     kppi.add_line("    if ( Present(xNacc) )    xNacc(is) = istatus_u(4)                ");
     kppi.add_line("    if ( Present(xNrej) )    xNrej(is) = istatus_u(5)                ");
   }
   kppi.add_line("                                                                    ");
   kppi.add_line("    if ( present(istatus) )  then                                   ");
   if(kpp_switches.is_vector()) {
     kppi.add_line("      istatus(4) =   istatus(4) + sum(Kacc(1:VL))                  ");
     kppi.add_line("      istatus(5) =   istatus(5) + sum(Krej(1:VL))                  ");
     kppi.add_line("      istatus(3) =   istatus(4) + istatus(5)                       ");
     kppi.add_line("      istatus(6) =   istatus(6) + istatus_u(6)                     ");
     kppi.add_line("      istatus(7) =   istatus(7) + istatus_u(7)                     ");
   } else {
     kppi.add_line("      istatus(1:8) = istatus(1:8) + istatus_u(1:8)                 ");
   }
   kppi.add_line("    end if                                                          ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  END DO                                                            ");
   kppi.add_line(" ");

   kppi.add_line("                                                                    ");
   kppi.add_line("! Deallocate input arrays                                           ");
   kppi.add_line("                                                                    ");
   for(iv=Vvar_list.begin();iv!=Vvar_list.end();iv++) {
//     kppi.add_line("  if ( allocated("+ iv->name +") )   deallocate("+ iv->name +" )    ");
   }

   kppi.add_line("                                                                    ");
   kppi.add_line("  data_loaded = .false.                                             ");
   kppi.add_line("                                                                    ");
   kppi.add_line("  return                                                            ");
   kppi.add_line("END SUBROUTINE chem_gasphase_integrate                              ");

//   e5_subroutines.push_back(kppi);
   kpp_subroutines.push_back(kppi);

   return;
}

