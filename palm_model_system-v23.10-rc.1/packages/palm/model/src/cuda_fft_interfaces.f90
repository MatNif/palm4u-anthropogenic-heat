!> @file cuda_fft_interfaces.f90
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
! Description:
! ------------
!> FORTRAN interfaces for the CUDA fft
!> Routines for the fft along x and y (forward/backward) using the CUDA fft
!--------------------------------------------------------------------------------------------------!
 MODULE cuda_fft_interfaces

#if defined ( __cuda_fft )

    USE kinds

    INTEGER(iwp) ::  CUFFT_C2C = Z'29'    !< Complex to Complex, interleaved
    INTEGER(iwp) ::  CUFFT_C2R = Z'2c'    !< Complex (interleaved) to Real
    INTEGER(iwp) ::  CUFFT_D2Z = Z'6a'    !< Double to Double-Complex
    INTEGER(iwp) ::  CUFFT_FORWARD = -1   !<
    INTEGER(iwp) ::  CUFFT_INVERSE =  1   !<
    INTEGER(iwp) ::  CUFFT_R2C = Z'2a'    !< Real to Complex (interleaved)
    INTEGER(iwp) ::  CUFFT_Z2D = Z'6c'    !< Double-Complex to Double
    INTEGER(iwp) ::  CUFFT_Z2Z = Z'69'    !< Double-Complex to Double-Complex

    PUBLIC


!
!-- cufftPlan1d( cufftHandle *plan, int nx, cufftType type, int batch )
    INTERFACE CUFFTPLAN1D

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
       SUBROUTINE CUFFTPLAN1D( plan, nx, type, batch ) BIND( C, name='cufftPlan1d' )

          USE ISO_C_BINDING

          INTEGER(C_INT)        ::  plan   !<
          INTEGER(C_INT), VALUE ::  batch  !<
          INTEGER(C_INT), VALUE ::  nx     !<
          INTEGER(C_INT), VALUE ::  type   !<
       END SUBROUTINE CUFFTPLAN1D

    END INTERFACE CUFFTPLAN1D

!
!-- cufftDestroy( cufftHandle plan )  !!! remove later if not really needed !!!
    INTERFACE CUFFTDESTROY

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
       SUBROUTINE CUFFTDESTROY( plan ) BIND( C, name='cufftDestroy' )

          USE ISO_C_BINDING

          INTEGER(C_INT), VALUE ::  plan

       END SUBROUTINE CUFFTDESTROY

    END INTERFACE CUFFTDESTROY


    INTERFACE CUFFTEXECZ2D

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
       SUBROUTINE CUFFTEXECZ2D( plan, idata, odata ) BIND( C, name='cufftExecZ2D' )

          USE ISO_C_BINDING
          USE kinds

          COMPLEX(dp), DEVICE   ::  idata(:,:,:)  !<

          INTEGER(C_INT), VALUE ::  plan          !<

          REAL(dp), DEVICE      ::  odata(:,:,:)  !<

       END SUBROUTINE CUFFTEXECZ2D

    END INTERFACE CUFFTEXECZ2D


    INTERFACE CUFFTEXECD2Z

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
       SUBROUTINE CUFFTEXECD2Z( plan, idata, odata ) bind( C, name='cufftExecD2Z' )

          USE ISO_C_BINDING

          USE kinds

          COMPLEX(dp), DEVICE   ::  odata(:,:,:)  !<

          INTEGER(C_INT), VALUE ::  plan          !<

          REAL(dp), DEVICE      ::  idata(:,:,:)  !<

       END SUBROUTINE CUFFTEXECD2Z

    END INTERFACE CUFFTEXECD2Z

#else

    INTERFACE CUFFTdummy

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Dummy interface to avoid compiler warnings in case of no bublic objects declared.
!--------------------------------------------------------------------------------------------------!
       SUBROUTINE CUFFTdummy( dummy )

          USE kinds

          REAL(wp) ::  dummy  !<

       END SUBROUTINE CUFFTdummy

    END INTERFACE CUFFTdummy

#endif

 END MODULE cuda_fft_interfaces
