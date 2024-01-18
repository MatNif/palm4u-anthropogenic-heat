//> @file fastv8_server.c
//--------------------------------------------------------------------------------------------------!
// This file is part of the PALM model system.
//
// PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along with PALM. If not, see
// <http://www.gnu.org/licenses/>.
//
// Copyright 2017-2023 Carl von Ossietzky Universität Oldenburg
// Copyright 2022-2023 pecanode GmbH
//--------------------------------------------------------------------------------------------------!
//
// Description:
// ------------
//> C Codebase of FASTv8 Coupler
//--------------------------------------------------------------------------------------------------!

#include<stdio.h>
#include<string.h>
#include<pthread.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/types.h> 
#include<sys/socket.h>
#include<netinet/in.h>
#include<netdb.h> 

//function headers that bind to Fortan (fastv8_updata.f90)
void palm_bld_data(int *turbine_id, int *numbld, int *numbldelem, double *rotrad);
void palm_sim_env(int *turbine_id, double *dtfast, double *shaft_height_fast_in);
void palm_turb_par(int *turbine_id, double *simtimef, double *vel_rot, double *xs_x, double *xs_y, double *xs_z);
void palm_force_table_bld(int *turbine_id, int *max_blade_elem_id, int *blade_elem_id, double *force_x, double *force_y, double *force_z, int *err);
void palm_position_table(int *turbine_id, int *max_blade_elem_id, int *blade_elem_id, double *pos_x, double *pos_y, double *pos_z, int *err);
void palm_vel_value(int *turbine_id, int *bladeid, int *elemid, double *compu, double *compv, double *compw, int *err);

/*------------------------------------------------------------------------------
| debug message
------------------------------------------------------------------------------*/
#define debug_msg(...) {FILE *fp; \
                        fp=fopen("DEBUG_CCODE","a");\
                        fprintf(fp, __VA_ARGS__); fprintf(fp, "\n");\
                        fclose(fp); }

#define error_msg(...) {FILE *fp; \
                        fp=fopen("DEBUG_CCODE","a");\
                        fprintf(fp, __VA_ARGS__); fprintf(fp, "\n");\
                        fclose(fp); }

#define status_msg(...) {FILE *fp; \
                        fp=fopen("DEBUG_CCODE","a");\
                        fprintf(fp, __VA_ARGS__); fprintf(fp, "\n");\
                        fclose(fp); }

//ssh for HLRN
#define SSHTUNNEL 0

//show output messages
#define CMSGOUT 1

//definition of precision
#define PREC 8

//definition of string array length
#define SOB 131072
#define SD 200
#define MTU 1400

//client targets
#define PALM 1
#define FAST 0

//communication control variables
char **port_fast;
char **host_fast;
int *fastturbdata;
char host_palm[SD];
char filename[] = "./commdata.srv"; 

//simulation control variables
int inumturb;
int ibldelemp;
int   simstepp;
double simtimep;
double simendtimep;
double simdtp;
double simturbonp;
double simdurationp;
int msgcp = 1;
double *adtfast;
double *atfast;
double *rhoairfast;
double *shaft_height_fast;
double *vel_rot;
double *simtimef;
double *xs_x;
double *xs_y;
double *xs_z;

//blade data
int *anumbld;
int *anumbldelem;
double *aradturb;

//data control variables
char delimiter[] = "#";
char msgrsp[SOB];
int *bytein_count;
int *byteout_count;
//int *cgelem;
//int *cgelemmax;

//variable to check whether everything has been sent
int *totalsent;
int *loopcounter;
int *stilltosend;

//debug variables
char command[8192];

//thread control
int *retptr;
pthread_t *tid;
int *retptr;
int msgid;

//---------------------------------------------------------------------------------------------------------------------------------
// function for cleaning up after simulation
//----------------------------------------------------------

int clean_exit(void)
{

  int i;

  free(anumbld);
  free(anumbldelem);
  free(aradturb);
  free(adtfast);
  free(atfast);
  for (i = 0; i < inumturb; i++){
    free(port_fast[i]);
  }
  free(port_fast);
//  free(cgelemmax);
//  free(cgelem);
  for (i = 0; i < inumturb; i++){
    free(host_fast[i]);
  }
  free(host_fast);

  return 0;

}

//---------------------------------------------------------------------------------------------------------------------------------
//exchange of time related variables between Fortran and C
//----------------------------------------------------------

int time_data(double simulated_time_palm, double curdt, int curstep, double turbon, double endtime)
{

  //sprintf(command, "echo \"%s%f:%f:%d:%f:%f\" >> DEBUG_CCODE", "Time PALM = ", simulated_time_palm, curdt, curstep, turbon, endtime);
  //system(command);	

  simstepp    = curstep;
  simdtp      = curdt;
  simtimep    = simulated_time_palm;
  simturbonp  = turbon;
  simendtimep = endtime;
 
  simdurationp = simendtimep-simturbonp;
  
//  if(CMSGOUT > 0){
//    debug_msg("[Debug] Time PALM = ", simtimep, simdtp, simstepp, simturbonp, simendtimep); 
//  }

  return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------
// byte to double 
//----------------------------------------------------------
double bfconv(unsigned char *af2b, int cbyte){
  double *val = (double*) &af2b[cbyte];

  return *val;

}

//---------------------------------------------------------------------------------------------------------------------------------
// byte to integer 
//----------------------------------------------------------

int biconv(unsigned char *af2b, int cbyte){
  int *biconvert = (int*) &af2b[cbyte];

  return *biconvert; 
}

//---------------------------------------------------------------------------------------------------------------------------------
// double to byte 
//----------------------------------------------------------

int fbconv(unsigned char *af2b, double fval, int cbyte){

  unsigned char *bytes = (unsigned char*) &fval;

  af2b[cbyte]   = bytes[0];
  af2b[cbyte+1] = bytes[1];
  af2b[cbyte+2] = bytes[2];
  af2b[cbyte+3] = bytes[3];
  af2b[cbyte+4] = bytes[4];
  af2b[cbyte+5] = bytes[5];
  af2b[cbyte+6] = bytes[6];
  af2b[cbyte+7] = bytes[7];

  return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------
// integer to byte 
//----------------------------------------------------------
int ibconv(unsigned char *af2b, int ival, int cbyte){
  unsigned char *bytes = (unsigned char *) &ival;

  af2b[cbyte]   = bytes[0];
  af2b[cbyte+1] = bytes[1];
  af2b[cbyte+2] = bytes[2];
  af2b[cbyte+3] = bytes[3];

  return 0;

}

//---------------------------------------------------------------------------------------------------------------------------------
// analyse FAST message 
//----------------------------------------------------------

/*
msgtype
a
-1 error in FAST
 1 connection check response (sending of blade data)
 2 end of FAST simulation (sending of stop signal)
 3 sending of positions
 4 sending of positions and forces
 5 sending a response: 1

*/

int analyse_msg(unsigned char *velmsg, int ccid){

  //data analysis variables
  int cpack, err;
  double bepos[3];
  double befor[3];
  double tvel_rot;
  double tsimtimef;
  double txs_x;
  double txs_y;
  double txs_z;

  int turbine_id_f;

  int ccidpone = ccid+1;
 
  err = 0;

 // int tmpelemid;   

  int msgcf;

  //int start_signal;
  int msgtypef;
  int cbyte = 0;

  //at least a header (start_signal, turbine id, message counter, message id) should be present
  if(bytein_count[ccid] < 16){
    error_msg("[Error] FAST message is corrupted (%d).", ccidpone);
    return -1;
  }

  //start_signal
  //start_signal = biconv(velmsg, cbyte);
  cbyte = cbyte+4;
  
  //turbine id
  turbine_id_f = biconv(velmsg, cbyte);
  cbyte = cbyte+4;
  if(turbine_id_f != ccidpone){
    error_msg("[Error] ID of responding turbine (%d) does not match intended addressee (%d).", turbine_id_f, ccidpone);
    return -1;
  }  
  
  //message counter
  msgcf = biconv(velmsg, cbyte);
  cbyte = cbyte+4;
  if(msgcf != msgcp){
    error_msg("[Error] Synchronization problem between FAST and PALM! Message counter mismatch (%d:%d).", msgcf, msgcp);
    return -1;
  }  

  //message type  
  msgtypef = biconv(velmsg, cbyte);
  cbyte = cbyte+4;

  if(CMSGOUT > 0){
    debug_msg("[Debug] Timep: %d", simtimep); 
    debug_msg("[Debug] Timestepp: %d", simstepp); 
    debug_msg("[Debug] Header: %d:%d:%d", turbine_id_f, msgcf, msgtypef);
  }

  
  if((msgtypef == -1) || ((msgtypef > 0) && (msgtypef < 6))){

    if(msgtypef == -1){
      //error in FAST
      return 2;
  
    }else if(msgtypef == 1){
      //connection check response


      if(bytein_count[ccid] != 52+8){
        error_msg("[Error] Unexpected number of bytes received %d:%d.", bytein_count[ccid], (ibldelemp*3+5)*4+8);
        return -1;
      }

      //time step of FAST instance
      adtfast[ccid] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
      status_msg("[Status] [adtfast] %f", adtfast[ccid]);

      //number of blades
      anumbld[ccid] = biconv(velmsg, cbyte);
      cbyte = cbyte+4;
      status_msg("[Status] [anumbld] %d", anumbld[ccid]);

      //number of blade elements
      anumbldelem[ccid] = biconv(velmsg, cbyte);
      cbyte = cbyte+4;
      status_msg("[Status] [anumbldelem] %d", anumbldelem[ccid]);

      //rotor radius
      aradturb[ccid] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
      status_msg("[Status] [aradturb] %f", aradturb[ccid]);

      ibldelemp = anumbld[ccid]*anumbldelem[ccid]+1;
      palm_bld_data(&ccidpone, &anumbld[ccid], &anumbldelem[ccid], &aradturb[ccid]);
    
      //air density of FAST instance
      rhoairfast[ccid] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
      status_msg("[Status] [rhoairfast] %f", rhoairfast[ccid]);
      
      //Shaft height
      shaft_height_fast[ccid] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
      status_msg("[Status] [shaft_height_fast] %f", shaft_height_fast[ccid]);

      return 0;

    }else if(msgtypef == 2){
      //end of simulation in FAST
   
      return 1;


    }else if(msgtypef == 3){
      //positions
      //check length of message
      if(bytein_count[ccid] != ((ibldelemp*3)*8+20)){
        error_msg("[Error] Unexpected number of bytes received %d:%d.", bytein_count[ccid], ((ibldelemp*3)*8+20)); //(ibldelemp*3+5)*4);
        return -1;
      }

    }else if(msgtypef == 4){
      //position and forces
      //check length of message
      if(bytein_count[ccid] != ((((ibldelemp-1)*2+1)*3*8)+20+8+24+8)){
        error_msg("[Error] Unexpected number of bytes received %d:%d.", bytein_count[ccid], ((((ibldelemp-1)*2+1)*3*8)+20+8+24+8)); //(((ibldelemp-1)*2+1)*3+5)*4);
        return -1;
      }
    }else if(msgtypef == 5){
      //response of FAST to sending velocities
   
      return 0;
      
    }
  }else{
    //unknown message type
    return -1;
  }
  

  if((msgtypef == 3) || (msgtypef == 4)){
    //positions of blade elements + hub from FAST 
    cpack = 0;
    
    while(cpack < ibldelemp) {

      bepos[0] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
      bepos[1] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
      bepos[2] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
     
      cpack++;
      //update blade element position in PALM
      palm_position_table(&ccidpone, &ibldelemp, (int *)&cpack, &bepos[0], &bepos[1], &bepos[2], &err);
      if(CMSGOUT > 0){
        debug_msg("[Debug] positions: %d %d %d %f %f %f", ccidpone, ibldelemp, cpack, bepos[0], bepos[1], bepos[2]);
      }
      
      if(err != 0){
        error_msg("[Error] Unable to update blade element positions in PALM (id: %d).", ccidpone);
        return -1;
      }
    }
  }
 
  if(msgtypef == 4){   
    //forces at blade elements from FAST
    cpack = 0;
    debug_msg("[Debug] sizeof velmsg: %d", sizeof(velmsg)); 
    while(cpack < ibldelemp-1){

      befor[0] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
      befor[1] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;
      befor[2] = bfconv(velmsg, cbyte);
      cbyte = cbyte+8;

      cpack++;
      //update forces at blade elements in PALM
      palm_force_table_bld(&ccidpone, &ibldelemp, (int *)&cpack, &befor[0], &befor[1], &befor[2], &err);
      if(CMSGOUT > 0){
        debug_msg("[Debug] forces:%d %d %d %f %f %f ", ccidpone, ibldelemp, cpack, befor[0], befor[1], befor[2]);
      }
      

      if(err != 0){
        error_msg("[Error] Unable to update forces at blade elements in PALM (id: %d).", ccidpone);
        return -1;
      }
    }
    
    //rotational velocity of FAST instance
    vel_rot[ccid] = bfconv(velmsg, cbyte);
    cbyte = cbyte+8;
    tvel_rot = vel_rot[ccid]; 
    //simulation time of FAST instance
    simtimef[ccid] = bfconv(velmsg, cbyte);
    cbyte = cbyte+8;
    tsimtimef = simtimef[ccid]; 
    //shaft coordinate system - xs direction
    xs_x[ccid] = bfconv(velmsg, cbyte);
    cbyte = cbyte+8;
    txs_x = xs_x[ccid]; 
    xs_y[ccid] = bfconv(velmsg, cbyte);
    cbyte = cbyte+8;
    txs_y = xs_y[ccid]; 
    xs_z[ccid] = bfconv(velmsg, cbyte);
    cbyte = cbyte+8;
    txs_z = xs_z[ccid]; 
    status_msg("[Status] [xs] %f, %f, %f", xs_x, xs_y, xs_z);
    
    palm_turb_par(&ccidpone, &tsimtimef, &tvel_rot, &txs_x, &txs_y, &txs_z); // sending FAST rotation velocity to fortran part of PALM
  }

return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------
// prepare PALM message 
//----------------------------------------------------------

int build_msg(unsigned char *msg2f, int msgtype, int ccid){

/*
msgtype

 1 connection initialisation
 2 PALM is ready
 3 sending of velocities
 4 resume simulation
 5 sending of velocities without response
 6 asking for positions

*/

  int i, j, err; 

  double velu, velv, velw;
  int ccidpone = ccid+1;
  err = 0;

  int cbld, cbldelem;
  double fval;
  int ival;
  int cval = 0;

  if(CMSGOUT > 0){
    debug_msg("[Debug] msg type %d", msgtype);
  }

  //valid message type?
  if((msgtype > 0) || (msgtype < 7)){

    //message header
    ival = 69;
    ibconv(msg2f, ival, cval);
    cval = cval+4;
    
    //turbine id
    ival = ccidpone;
    ibconv(msg2f, ival, cval);
    cval = cval+4;

    //message id
    ival = msgcp;
    ibconv(msg2f, ival, cval);
    cval = cval+4;

    //message type
    ival = msgtype;
    ibconv(msg2f, ival, cval);
    cval = cval+4;

    //message type = connection initialisation
    if(msgtype == 1){
 
      //host
      ival = 999;
      fbconv(msg2f, ival, cval);
      cval = cval+4;
      ival = 999;
      ibconv(msg2f, ival, cval);
      cval = cval+4;
      ival = 999;
      fbconv(msg2f, ival, cval);
      cval = cval+4;
      ival = 999;
      fbconv(msg2f, ival, cval);
      cval = cval+4;

      //port
      ival = 65535;
      ibconv(msg2f, ival, cval);
      cval = cval+4;

      //duration of simulation
      fval = simdurationp;
      fbconv(msg2f, fval, cval);
      cval = cval+8;
      
    //message type = palm ready
    }else if(msgtype == 2){


    //message type = regular type (send velocities)
    }else if(msgtype == 3){

      if(CMSGOUT > 0){
        debug_msg("[Debug] Sending 3*%d velocities to turbine %d, timestep: .",simstepp,  anumbld[ccid]*anumbldelem[ccid]+1, ccidpone);
      }

      //loop blades
      for (j = 0; j < anumbld[ccid]; j++){    
        //loop blade elements
        for (i = 0; i < anumbldelem[ccid]; i++){
          cbld = j+1;
          cbldelem = i+1;
          palm_vel_value(&ccidpone, &cbld, &cbldelem, &velu, &velv, &velw, &err);
          if(err != 0){
            error_msg("[Error] Unable to retrieve PALM velocities (id: %d).", ccidpone);
            return -1;
          }
          debug_msg("[Debug] velocities: [%d][%d] %f:%f:%f", cbld, cbldelem, velu, velv, velw); 

          fval = velu;
          fbconv(msg2f, fval, cval);
          cval = cval+8;
          fval = velv;
          fbconv(msg2f, fval, cval);
          cval = cval+8;
          fval = velw;
          fbconv(msg2f, fval, cval);
          cval = cval+8;
        }
      }

      //add hub velocities
      cbld = 0;
      cbldelem = 0;
      palm_vel_value(&ccidpone, &cbld, &cbldelem, &velu, &velv, &velw, &err); 
      debug_msg("[Debug] hub velocity: [%d][%d] %f:%f:%f", cbld, cbldelem, velu, velv, velw);

      fval = velu;
      fbconv(msg2f, fval, cval);
      cval = cval+8;
      fval = velv;
      fbconv(msg2f, fval, cval);
      cval = cval+8;
      fval = velw;
      fbconv(msg2f, fval, cval);
      cval = cval+8;

    //message type = resume simulation
    }else if(msgtype == 4){

    //message type = regular type (sending velocities without response)
    }else if(msgtype == 5){

      if(CMSGOUT > 0){
        debug_msg("[Debug] Sending 3*%d velocities to turbine %d, timestep: .",simstepp,  anumbld[ccid]*anumbldelem[ccid]+1, ccidpone);
      }

      //loop blades
      for (j = 0; j < anumbld[ccid]; j++){    
        //loop blade elements
        for (i = 0; i < anumbldelem[ccid]; i++){
          cbld = j+1;
          cbldelem = i+1;
          palm_vel_value(&ccidpone, &cbld, &cbldelem, &velu, &velv, &velw, &err);
          if(err != 0){
            error_msg("[Error] Unable to retrieve PALM velocities (id: %d).", ccidpone);
            return -1;
          }
          debug_msg("[Debug] velocities: [%d][%d] %f:%f:%f", cbld, cbldelem, velu, velv, velw); 

          fval = velu;
          fbconv(msg2f, fval, cval);
          cval = cval+8;
          fval = velv;
          fbconv(msg2f, fval, cval);
          cval = cval+8;
          fval = velw;
          fbconv(msg2f, fval, cval);
          cval = cval+8;
        }
      }

      //add hub velocities
      cbld = 0;
      cbldelem = 0;
      palm_vel_value(&ccidpone, &cbld, &cbldelem, &velu, &velv, &velw, &err); 
      debug_msg("[Debug] hub velocity: [%d][%d] %f:%f:%f", cbld, cbldelem, velu, velv, velw);

      fval = velu;
      fbconv(msg2f, fval, cval);
      cval = cval+8;
      fval = velv;
      fbconv(msg2f, fval, cval);
      cval = cval+8;
      fval = velw;
      fbconv(msg2f, fval, cval);
      cval = cval+8;

    //message type = asking for positions
    }else if(msgtype == 6){
   
    }    
    //end of message 
    ival = 69;
    ibconv(msg2f, ival, cval);
    cval = cval+4;

    byteout_count[ccid] = cval;

  }else{
    error_msg("[Error] Unknown message type (id: %d).", ccidpone);
    return -1;
  }

return 0;

}

//---------------------------------------------------------------------------------------------------------------------------------
// PALM client communication function 
//----------------------------------------------------------
void *commfunc(void *param){

  //socket variables
  int sockfd, cflag, i;
  struct addrinfo hints;
  struct addrinfo *result, *rp;
  struct sockaddr_in serv_addr;
  struct hostent *server;

  //message strings (könnte vlt auch nur einer sein)
  unsigned char msgout[SOB];
  unsigned char msgrsp[SOB];
  unsigned char msgrsppack[MTU];

  //communication control
  unsigned char omsglen[4];
  int imsglen;

  //function parameter variables
  int ccid = *(int *) param;
  free(param);

  //return value
  //retptr[ccid] = -1; 

  int bytesread = 0;

  memset(&hints, 0, sizeof(hints));
  hints.ai_family = AF_INET;    /* Allow IPv4 */
  hints.ai_socktype = 0;
  hints.ai_flags = 0;
  hints.ai_protocol = 0;          /* Any protocol */

  cflag = getaddrinfo(host_fast[ccid], port_fast[ccid], &hints, &result);
  if (cflag != 0) {
    error_msg("[Error] getaddrinfo failed (id: %d).", ccid+1);
    retptr[ccid] = -1;
    pthread_exit((void*) &retptr[ccid]);
  }

  for (rp = result; rp != NULL; rp = rp->ai_next) {
    sockfd = socket(rp->ai_family, rp->ai_socktype, rp->ai_protocol);
    if (sockfd == -1)
        continue;

    if (connect(sockfd, rp->ai_addr, rp->ai_addrlen) != -1)
        break;  /* Success */

    close(sockfd);
  }

  freeaddrinfo(result);           /* No longer needed */

  if (rp == NULL) {               /* No address succeeded */
    error_msg("[Error] Unable to open socket (id: %d).", ccid+1);
    retptr[ccid] = -1;
    pthread_exit((void*) &retptr[ccid]);
  }

  //prepare message string
  if(CMSGOUT > 0){
    debug_msg("[Debug] [build message] %d", msgid);
  }


  cflag = build_msg(msgout, msgid, ccid);
  if (cflag != 0){ 
    error_msg("[Error] Unable to prepare message for FAST (id: %d).", ccid+1);
    retptr[ccid] = -1;
    pthread_exit((void*) &retptr[ccid]);
  }

  //send data string to FAST
  if(CMSGOUT > 0){
    debug_msg("[Debug] PALM client sending data to: %s:%s (mc = %d)", host_fast[ccid], port_fast[ccid], msgcp);
    debug_msg("[Debug] Bytes: %d", byteout_count[ccid]);
  }

  //sending total number of bytes that will be sent shortly
  ibconv(omsglen, byteout_count[ccid], 0);

  totalsent[ccid] = 0;
  stilltosend[ccid] = 4;
  loopcounter[ccid] = 0;
  while(totalsent[ccid] < 4){
      cflag = write(sockfd, omsglen, 4);
      if (cflag < 0){ 
         error_msg("[Error] Unable to write to socket (id: %d).", ccid+1);
         retptr[ccid] = -1;
         pthread_exit((void*) &retptr[ccid]);
      }
      totalsent[ccid] += cflag;
      stilltosend[ccid] -= cflag;
      loopcounter[ccid] = loopcounter[ccid]+1;
      debug_msg("[Debug] Zu senden PALM A: 4");
      debug_msg("[Debug] Gesendet PALM A: %d", totalsent[ccid]);
      debug_msg("[Debug] Loopcounter A: %d", loopcounter[ccid]);
      //printf("[Information] Zu senden PALM A 4:\n");
      //printf("[Information] Zu senden PALM A:\n", totalsent);
  }

  //sending of main message (data)
  totalsent[ccid] = 0;
  stilltosend[ccid] = byteout_count[ccid];
  loopcounter[ccid] = 0;
  while(totalsent[ccid] < byteout_count[ccid]){
      cflag = write(sockfd, msgout+totalsent[ccid], stilltosend[ccid]);
      if (cflag < 0){ 
         error_msg("[Error] Unable to write to socket (id: %d).", ccid+1);
         retptr[ccid] = -1;
         pthread_exit((void*) &retptr[ccid]);
      }
      totalsent[ccid] += cflag;
      stilltosend[ccid] -= cflag;
      loopcounter[ccid] = loopcounter[ccid]+1;
      debug_msg("[Debug] Zu senden PALM B: %d", byteout_count[ccid]);
      debug_msg("[Debug] Gesendet PALM B: %d", totalsent[ccid]);
      debug_msg("[Debug] Loopcounter B: %d", loopcounter[ccid]);
      //      printf("[Information] Zu senden PALM B 4:\n");
      //      printf("[Information] Zu senden PALM B:\n", totalsent);
  }

  if(CMSGOUT > 0){
    debug_msg("[Debug] Sending completed (id: %d).", ccid+1);

    debug_msg("[Debug] Waiting for response (id: %d)...", ccid+1);

    debug_msg("[Debug] Response of server (id: %d):", ccid+1);
  }

  memset(msgrsp, 0, SOB);
  memset(msgrsppack, 0, MTU);
  bytesread = 0;

 
  //read number of bytes that FAST is going to send
  cflag = read(sockfd,msgrsppack, 4);
  
  
    
  if (cflag != 4){ 
    error_msg("[Error] Unable to read msg length from socket (id: %d).", ccid+1);
    error_msg("Returned length: %d, Returned value: %d.", cflag, msgrsppack[0]);
    retptr[ccid] = -1;
    pthread_exit((void*) &retptr[ccid]);
  } 
  	
  

  imsglen = biconv(msgrsppack, 0);
  
  debug_msg("Message length to be read: %d.",  imsglen); 

  if(CMSGOUT > 0){
    debug_msg("[Debug] Message length: %d (id: %d):", imsglen, ccid+1);
  }

  //read data
  while(bytesread < imsglen){
  
    cflag = read(sockfd,msgrsppack, MTU);

    if (cflag < 0){ 
      error_msg("[Error] Unable to read from socket (id: %d).", ccid+1);
      retptr[ccid] = -1;
      pthread_exit((void*) &retptr[ccid]);
    }

    //attach current package to the array that holds the whole response
    memcpy((msgrsp+bytesread), msgrsppack, cflag); 
    bytesread = bytesread + cflag;

    memset(msgrsppack, 0, MTU);
  }

  bytein_count[ccid] = bytesread;

  if(CMSGOUT > 0){
    debug_msg("[Debug] Bytes(Total): %d (%d)", bytein_count[ccid], imsglen);
  }

  //output response string
  if(CMSGOUT > 0){

    char* bufout_str = (char*) malloc (3*bytein_count[ccid] + 1);
    if (bufout_str == NULL){
      error_msg("[Error] Out of memory [bufout_ptr].");
      retptr[ccid] = -1;
      pthread_exit((void*) &retptr[ccid]);
    }

    char* bufout_ptr = bufout_str;
  
    for (i = 0; i < bytein_count[ccid]; i++){
      bufout_ptr += sprintf(bufout_ptr, "%02X", msgrsp[i]);
      if(i != bytein_count[ccid]-1){
        bufout_ptr += sprintf(bufout_ptr, ":");
      }
    }

    *(bufout_ptr + 1) = '\0';

    debug_msg("[Debug] %s", bufout_str);
    free(bufout_str);
  }

  //analyse response string
  cflag = analyse_msg(msgrsp, ccid); 
  if(cflag < 0){
    error_msg("[Error] Unable to understand FAST response (id: %d).", ccid+1);
    retptr[ccid] = -1;
    pthread_exit((void*) &retptr[ccid]);
  }else if(cflag == 1){
    status_msg("[Status] At least one FAST instance finished the simulation (id: %d).",ccid+1);
    retptr[ccid] = 1; 
    pthread_exit((void*) &retptr[ccid]);
  }else if(cflag == 2){
    error_msg("[Error] An error in at least one FAST instance occurred (id: %d).",ccid+1);
    retptr[ccid] = -1;
    pthread_exit((void*) &retptr[ccid]);
  }

  if(CMSGOUT > 0){
    debug_msg("[Debug] Fertig");
  }

  close(sockfd); 
  return NULL;
}

//---------------------------------------------------------------------------------------------------------------------------------
// PALM client for communication with FAST/PALM server
//----------------------------------------------------------
int commclient(int target, int ccmsgid){

  int err, iturb, thread_join_res, errjoin;
  int *tnum;
  double tdtfast;
  double trhoairfast;
  double shaft_height_fast_tmp;
  int ccidpone;

  msgid = ccmsgid;
  err = 0;

  //valid message type?
  if((msgid < 1) || (msgid > 6)){
    /*
    1 connection initialisation
    2 PALM is ready
    3 sending of velocities
    4 resume simulation
    */
    error_msg("[Error] Unknown message id in call to communication client.");
    return -1;   
  }

  status_msg("------------------------------------------------------------------------------");

  //identify communication partner
  if(target == FAST){

    //loop through all FAST server
    for (iturb = 0; iturb < inumturb; iturb++){
      tnum = malloc(sizeof(int*));
      if (tnum == NULL){
        error_msg("[Error] Out of memory [tnum].");
        return -1;
      }
      *tnum = iturb;

      //create a separate communication thread for every FAST server
      err = pthread_create(&tid[iturb], NULL, &commfunc, (void *) tnum);
      if (err != 0){
        error_msg("[Error] Unable to create client thread [%s]!", strerror(err));
        clean_exit();
        return -1;
      }else{
        status_msg("[Status] Client thread %d for communication with %s:%s running (mc=%d)...", iturb, host_fast[iturb], port_fast[iturb], msgcp);
      }
    }

    //join all threads
    errjoin = 0;
    for (iturb = 0; iturb < inumturb; iturb++){
      thread_join_res = pthread_join(tid[iturb], NULL);
      if(thread_join_res != 0){
        errjoin = 1;
        error_msg("[Error] Unable to join thread %d.", iturb);
      }
    }
    if(errjoin != 0){
      clean_exit();
      return -1;
    }

    //check thread results
    for (iturb = 0; iturb < inumturb; iturb++){
      status_msg("[Status] Thread %d returned: %d", iturb, retptr[iturb]);
      if(retptr[iturb] < 0){
        error_msg("[Error] Thread %d returned: %d.", iturb, retptr[iturb]);
        clean_exit(); 
        return -1;
      }

      //at least one FAST instance finished simulation
      if(retptr[iturb] == 1){
        status_msg("[Status] At least one FAST instance finished the simulation!");
        clean_exit();
        return 1;
      }
    }

    status_msg("[Status] All threads returned successfully.");

    //verfify that the time steps of all FAST instances are the same    
    if(msgcp == 1){
      for (iturb = 0; iturb < inumturb; iturb++){
        if(adtfast[iturb] != adtfast[0]){
          error_msg("[Error] All FAST instances must use the same time step.");
          clean_exit();
          return 1;
        }

        if(rhoairfast[iturb] != rhoairfast[0]){
          error_msg("[Error] All FAST instances must run with the same air density.");
          clean_exit();
          return 1;
        }
      } 

      //set time step PALM
      tdtfast = adtfast[0];  
      status_msg("[Status] [FAST time step size] %f", tdtfast);

      //set air density PALM
      trhoairfast = rhoairfast[0];  
      status_msg("[Status] [FAST air density] %f", trhoairfast);
      
      for (iturb = 0; iturb < inumturb; iturb++){
        //set shaft height
        shaft_height_fast_tmp = shaft_height_fast[iturb];  
        status_msg("[Status] [FAST Shaft height] %f", shaft_height_fast_tmp);
        ccidpone = iturb+1;

        palm_sim_env(&ccidpone, &tdtfast, &shaft_height_fast_tmp); // sending FAST time step and Shaft height to fortran part of PALM
      } 
    }

    msgcp++;
  }else{
    error_msg("[Error] Unknown communication target.");
    clean_exit();
    return 1;
  }

return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------
// allocate memory
//----------------------------------------------------------
int alloc_mem(void)
{

  int i;

  //communication data
    //array of FAST ports
    port_fast = malloc(inumturb*sizeof(char*));
    if (port_fast == NULL){
      error_msg("[Error] Out of memory [port_fast].");
      return -1;
    }

    for (i = 0; i < inumturb; i++){
      port_fast[i] = malloc(SD*sizeof(char));
      if (port_fast[i] == NULL){
        error_msg("[Error] Out of memory [port_fast].");
        return -1;
      }
    }
   
    //array of FAST hosts
    host_fast = malloc(inumturb*sizeof(char*));
    if (host_fast == NULL){
      error_msg("[Error] Out of memory [host_fast].");
      return -1;
    }

    for (i = 0; i < inumturb; i++){
      host_fast[i] = malloc(SD*sizeof(char));
      if (host_fast[i] == NULL){
        error_msg("[Error] Out of memory [host_fast].");
        return -1;
      }
    }

  //time data 
    //array with time steps of FAST instances
    adtfast = malloc(inumturb*sizeof(double));
    if (adtfast == NULL){
      error_msg("[Error] Out of memory [adtfast].");
      return -1;
    }

  //model data 
    //array with air densities of FAST instances
    rhoairfast = malloc(inumturb*sizeof(double));
    if (rhoairfast == NULL){
      error_msg("[Error] Out of memory [rhoairfast].");
      return -1;
    }
    
    vel_rot = malloc(inumturb*sizeof(double));
    if (vel_rot == NULL){
      error_msg("[Error] Out of memory [vel_rot].");
      return -1;
    }
    
    simtimef = malloc(inumturb*sizeof(double));
    if (vel_rot == NULL){
      error_msg("[Error] Out of memory [simtimef].");
      return -1;
    }
    
    shaft_height_fast = malloc(inumturb*sizeof(double));
    if (shaft_height_fast == NULL){
      error_msg("[Error] Out of memory [shaft_height_fast].");
      return -1;
    }
    
    xs_x = malloc(inumturb*sizeof(double));
    if (xs_x == NULL){
      error_msg("[Error] Out of memory [xs_x].");
      return -1;
    }    
    xs_y = malloc(inumturb*sizeof(double));
    if (xs_y == NULL){
      error_msg("[Error] Out of memory [xs_y].");
      return -1;
    }
    xs_z = malloc(inumturb*sizeof(double));
    if (xs_z == NULL){
      error_msg("[Error] Out of memory [xs_z].");
      return -1;
    }


  //blade data 
    //array with number of blades of each turbine
    anumbld = malloc(inumturb*sizeof(int));
    if (anumbld == NULL){
      error_msg("[Error] Out of memory [anumbld].");
      return -1;
    }

    //array with number of blade elements of each turbine
    anumbldelem = malloc(inumturb*sizeof(int));
    if (anumbldelem == NULL){
      error_msg("[Error] Out of memory [anumbldelem].");
      return -1;
    }

    //array with radius of turbine rotors
    aradturb = malloc(inumturb*sizeof(double));
    if (aradturb == NULL){
      error_msg("[Error] Out of memory [aradturb].");
      return -1;
    }

  //thread data
    //thread id 
    tid = malloc(inumturb*sizeof(pthread_t));
    if (tid == NULL){
      error_msg("[Error] Out of memory [tid].");
      return -1;
    }

    //return value
    retptr = malloc(inumturb*sizeof(int));
    if (retptr == NULL){
      error_msg("[Error] Out of memory [retptr].");
      return -1;
    }
    for (i = 0; i < inumturb; i++){
      retptr[i] = 0;
    }

    //return value (ein rückgabewert eventuell unnötig)
 /*   retptr = malloc(inumturb*sizeof(int));
    if (retptr == NULL){
      error_msg("[Error] Out of memory.");
      return -1;
    }
*/

  //data handling
    //max elements of response string
  /*  cgelemmax = malloc(inumturb*sizeof(int));
    if (cgelemmax == NULL){
      error_msg("[Error] Out of memory [cgelemmax].");
      return -1;
    }

    //current element number of response string
    cgelem = malloc(inumturb*sizeof(int));
    if (cgelem == NULL){
      error_msg("[Error] Out of memory [cgelem].");
      return -1;
    }
*/
    //packages in message to FAST
    bytein_count = malloc(inumturb*sizeof(int));
    if (bytein_count == NULL){
      error_msg("[Error] Out of memory [bytein_count].");
      return -1;
    }

    //packages in message to FAST
    byteout_count = malloc(inumturb*sizeof(int));
    if (byteout_count == NULL){
      error_msg("[Error] Out of memory [byteout_count].");
      return -1;
    }

    totalsent = malloc(inumturb*sizeof(int));
    if (totalsent == NULL){
      error_msg("[Error] Out of memory [totalsent].");
      return -1;
    }

    stilltosend = malloc(inumturb*sizeof(int));
    if (stilltosend == NULL){
      error_msg("[Error] Out of memory [stilltosend].");
      return -1;
    }

    loopcounter = malloc(inumturb*sizeof(int));
    if (loopcounter == NULL){
      error_msg("[Error] Out of memory [loopcounter].");
      return -1;
    }


  return 0;
}


//---------------------------------------------------------------------------------------------------------------------------------
// init communication
//----------------------------------------------------------
int init_comm(int ifnumturb, char* fast_addr[], char* fast_port[])
{

  int cflag, i;
  inumturb = ifnumturb;

  status_msg("[Status] Initializing communication routine.");

  status_msg("[Status] Number of turbines in simulation: %d", inumturb);

  //allocate memory
  cflag = alloc_mem();
  if (cflag < 0){ 
    error_msg("[Error] Unable to allocate memory!");
    return -1;
  }

  status_msg("[Status] Importing communication data.");
  for (i = 0; i < inumturb; i++){
    strcpy(host_fast[i], fast_addr[i]);
    strcpy(port_fast[i], fast_port[i]);
    status_msg("[Status] FAST INSTANCE NO. %d: %s:%s\n", i+1, host_fast[i], port_fast[i]);
    printf("[Status] FAST DATA%d: %s:%s\n", i+1, host_fast[i], port_fast[i]);
  }

  //initialise blade data arrays
  for (i = 0; i < inumturb; i++){
    anumbld[i] = -1;
    anumbldelem[i] = -1;
    aradturb[i]    = -1;
    adtfast[i]     = -1;
  }

  return 0;
}
