!> @file pmc_general_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
!> @author Klaus Ketelsen (no affiliation)
!
! Description:
! ------------
!> Structure definition and utilities of Palm Model Coupler
!--------------------------------------------------------------------------------------------------!
 MODULE pmc_general

#if defined( __parallel )
    USE, INTRINSIC ::  ISO_C_BINDING

    USE kinds

    USE MPI

    IMPLICIT NONE

    INTEGER(iwp) ::  pmc_max_array  !< max # of arrays which can be coupled
                                    !< - will be determined dynamically in pmc_interface

    INTEGER(iwp), PARAMETER ::  da_desclen       =  8  !<
    INTEGER(iwp), PARAMETER ::  da_namelen       = 16  !<
    INTEGER(iwp), PARAMETER ::  pmc_da_name_err  = 10  !<
    INTEGER(iwp), PARAMETER ::  pmc_max_models   = 64  !<
    INTEGER(iwp), PARAMETER ::  pmc_status_ok    =  0  !<
    INTEGER(iwp), PARAMETER ::  pmc_status_error = -1  !<

    TYPE ::  xy_ind  !< pair of indices in horizontal plane
       INTEGER(iwp) ::  i  !<
       INTEGER(iwp) ::  j  !<
    END TYPE

    TYPE ::  arraydef
       CHARACTER(LEN=da_namelen) ::  Name  !< name of array

       INTEGER(iwp) ::  coupleindex  !<
       INTEGER(iwp) ::  dimkey       !< key for NR dimensions and array type
       INTEGER(iwp) ::  nrdims       !< number of dimensions
       INTEGER(iwp) ::  RecvSize     !< size in receive buffer
       INTEGER(iwp) ::  SendSize     !< size in send buffer
       INTEGER(idp) ::  RecvIndex    !< index in receive buffer
       INTEGER(idp) ::  SendIndex    !< index in send buffer
       INTEGER(iwp) ::  ks           !< Start index in z direction (3d arrays on parent only)
       INTEGER(iwp) ::  ke           !< end   index in z direction (3d arrays on parent only)

       INTEGER(iwp), DIMENSION(4) ::  a_dim  !< size of dimensions

       TYPE(C_PTR) ::  data     !< pointer of data in parent space
       TYPE(C_PTR) ::  SendBuf  !< data pointer in send buffer
       TYPE(C_PTR) ::  RecvBuf  !< data pointer in receive buffer

       TYPE(arraydef), POINTER ::  next  !<

       TYPE(C_PTR), DIMENSION(2) ::  po_data  !< base pointers, pmc_s_set_active_data_array
                                              !< sets active pointer
    END TYPE arraydef

    TYPE ::  pedef
       INTEGER(iwp) ::  nr_arrays = 0  !< number of arrays which will be transfered
       INTEGER(iwp) ::  nrele          !< number of elements, same for all arrays

       TYPE(arraydef), POINTER, DIMENSION(:) ::  array_list  !< list of data arrays to be transfered

       TYPE(xy_ind), POINTER, DIMENSION(:) ::  locInd  !< xy index local array for remote PE
    END TYPE pedef

    TYPE ::  childdef
       INTEGER(iwp) ::  inter_comm        !< inter communicator model and child
       INTEGER(iwp) ::  inter_npes        !< number of PEs child model
       INTEGER(iwp) ::  intra_comm        !< intra communicator model and child
       INTEGER(iwp) ::  intra_rank        !< rank within intra_comm
       INTEGER(iwp) ::  model_comm        !< communicator of this model
       INTEGER(iwp) ::  model_npes        !< number of PEs this model
       INTEGER(iwp) ::  model_rank        !< rank of this model
       INTEGER(idp) ::  totalbuffersize   !<
       INTEGER(iwp) ::  win_parent_child  !< MPI RMA for preparing data on parent AND child side

       TYPE(pedef), DIMENSION(:), POINTER ::  pes  !< list of all child PEs
    END TYPE childdef

    TYPE ::  da_namedef  !< data array name definition
       CHARACTER(LEN=da_desclen) ::  childdesc     !< child array description
       CHARACTER(LEN=da_namelen) ::  nameonchild   !< name of array within child
       CHARACTER(LEN=da_namelen) ::  nameonparent  !< name of array within parent
       CHARACTER(LEN=da_desclen) ::  parentdesc    !< parent array description

       INTEGER(iwp) ::  couple_index  !< unique number of array
    END TYPE da_namedef

    TYPE(arraydef), POINTER ::  next  !<

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC pmc_g_setname

!
!-- Public variables, constants and types
    PUBLIC arraydef,                                                                               &
           childdef,                                                                               &
           da_desclen,                                                                             &
           da_namedef,                                                                             &
           da_namelen,                                                                             &
           next,                                                                                   &
           pedef,                                                                                  &
           pmc_da_name_err,                                                                        &
           pmc_max_array,                                                                          &
           pmc_max_models,                                                                         &
           pmc_status_error,                                                                       &
           pmc_status_ok,                                                                          &
           xy_ind

    INTERFACE pmc_g_setname
       MODULE PROCEDURE pmc_g_setname
    END INTERFACE pmc_g_setname


 CONTAINS


!---------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Add array to list of arraydef structure. No the arra "name" is schedules for parent child transfer
!---------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_g_setname( mychild, couple_index, aname )

    CHARACTER(LEN=*), INTENT(IN) ::  aname  !<

    INTEGER(iwp) ::  i  !<

    INTEGER(iwp), INTENT(IN) ::  couple_index  !<

    TYPE(childdef), INTENT(INOUT) ::  mychild  !<

    TYPE(pedef), POINTER ::  ape  !<


!
!-- Assign array to next free index in array list.
!-- Set name of array in arraydef structure
    DO  i = 1, mychild%inter_npes
       ape => mychild%pes(i)
       ape%nr_arrays = ape%nr_arrays + 1
       ape%array_list(ape%nr_arrays)%name        = aname
       ape%array_list(ape%nr_arrays)%coupleindex = couple_index
    ENDDO

 END SUBROUTINE pmc_g_setname
#endif


 END MODULE pmc_general
