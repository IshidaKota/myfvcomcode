










!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!==============================================================================|
!  Calculate Bottom Drag Coefficient based on Bottom Roughness                 !
!   note:                                                                      !
!   when the log function derived from the constant stress log-viscous         !
!   layer is applied to an estuary, if the value of z0 is close to             !
!   (zz(kbm1)-z(kb)*dt1, drag coefficient "cbc" could become a huge            !
!   number due to near-zero value of alog function. In our application         !
!   we simply cutoff at cbc=0.005. One could adjust this cutoff value          !
!   based on observations or his or her experiences.                           !   
!   CALCULATES:   WUBOT(N), WVBOT(N) : BOTTOM SHEAR STRESSES                   !
!==============================================================================|

SUBROUTINE BOTTOM_ROUGHNESS

  !==============================================================================!
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_WD
  USE MOD_PAR

  IMPLICIT NONE
  INTEGER :: I,II
  REAL(SP), PARAMETER  :: KAPPA = .40_SP   !!VON KARMAN LENGTH SCALE
  REAL(SP), PARAMETER  :: VK2   = .160_SP  !!KAPPA SQUARED
  REAL(SP)             :: ZTEMP,BTPS,RR,U_TAUB,Z0B_GOTM,CFF

! USED IN 2D MODEL ONLY
!   REAL(SP), PARAMETER  :: CONST_CD=.0015_SP      !! CD SET CONSTANT TO THIS VALUE
   REAL(SP), PARAMETER  :: ALFA =  .166667_SP, &   !! POWER OF WATER DEPTH
                           NN   = 0.02_SP          !! FACTOR TO DIVIDE
!   REAL(SP), PARAMETER  :: CFMIN   = .0025_SP, &  !! DEEP WATER VALUE
!                           H_BREAK = 1.0_SP,    & !! 
!                           THETA   = 10._SP,   &  !!
!                           LAMB    = 0.3333333333_SP
  !==============================================================================!

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: BOTTOM_ROUGHNESS"

  !
  !  SET CONSTANTS
  !

  SELECT CASE(BOTTOM_ROUGHNESS_TYPE) 
  !==============================================================================|
  CASE(BR_ORIG) !USE ORIGINAL FVCOM FORM FOR BOTTOM FRICTION  |
  !==============================================================================|

     ! SET A EFFECTIVE MAXIMUM FOR CBC USING THE DEPTH
!!     WHERE (DT1 > 3.0_SP)
!!        CBC = VK2/(LOG((ZZ1(:,KBM1)-Z1(:,KB))*DT1/CC_Z0B))**2
!!     ELSEWHERE
!!        CBC = VK2/(LOG((ZZ1(:,KBM1)-Z1(:,KB))*3.0/CC_Z0B))**2
!!     END WHERE
     WHERE (DT1(1:NT) > 3.0_SP)
        CBC(1:NT) = VK2/(LOG((ZZ1(1:NT,KBM1)-Z1(1:NT,KB))*DT1(1:NT)/CC_Z0B(1:NT)))**2
     ELSEWHERE
        CBC(1:NT) = VK2/(LOG((ZZ1(1:NT,KBM1)-Z1(1:NT,KB))*3.0/CC_Z0B(1:NT)))**2
     END WHERE
     
     ! SET A MINIMUM VALUE FOR CBC
     WHERE (CBC < CBCMIN)
        CBC=CBCMIN
     END WHERE


  !==============================================================================|
  CASE(BR_GOTM) !GOTM FORMULATION FOR BOTTOM FRICTION    |
  !==============================================================================|

     !----Convert Input Z0B to GOTMS H0B
     !     H0B = Z0B/.03  
     ! DAS fixed bug to match gotm's friction.f90
     DO I=1,N
        U_TAUB = 0.0_SP
        DO II=1,40       
           IF (UMOL <= 0.0_SP) THEN
              Z0B_GOTM=CC_Z0B(I)   !0.03*H0B 
           ELSE
              Z0B_GOTM=0.1_SP*UMOL/MAX(UMOL,U_TAUB)+CC_Z0B(I) !0.03*H0B
           END IF
           ztemp=(zz1(I,kbm1)-z1(I,kb))*dt1(i)
           RR=KAPPA/(LOG((Z0B_GOTM+ZTEMP)/Z0B_GOTM))
           U_TAUB = RR*SQRT( U(I,KBM1)*U(I,KBM1) + V(I,KBM1)*V(I,KBM1) )
        END DO
        CBC(I) =   RR*RR
     END DO

  ! Added by Siqi Li@2015/08/29
  !==============================================================================|
  CASE(BR_UDEF) !READ DRAG COEFFICENT FROM THE BOTTOM ROUGHNESS FILE    |
  !==============================================================================|

     CBC = CBC_UD !; DEALLOCATE(CBC_UD)
  ! Siqi Li

  CASE DEFAULT
     CALL FATAL_ERROR ("BROUGH: UNKNOWN BOTTOM_ROUGHNESS_TYPE:"&
          & ,TRIM(BOTTOM_ROUGHNESS_TYPE) )
  END SELECT


  !==============================================================================|
  !  CALCULATE SHEAR STRESS ON BOTTOM  --> WUBOT/WVBOT                           |
  !==============================================================================|
  DO  I = 1, N
     BTPS= CBC(I)*SQRT(U(I,KBM1)**2+V(I,KBM1)**2)
     WUBOT(I) = -BTPS * U(I,KBM1)
     WVBOT(I) = -BTPS * V(I,KBM1)
     !for plb case only
  
  END DO


  !==============================================================================|
  !  Calculate shear stress on nodes (x-component, y-component, magnitude) 
  !==============================================================================|
    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUBOT,WVBOT)
  TAUBM = SQRT(WUBOT**2 + WVBOT**2) 

     CALL E2N2D(WUBOT,WUBOT_N)
     CALL E2N2D(WVBOT,WVBOT_N)
  TAUBM_N = SQRT(WUBOT_N**2 + WVBOT_N**2) 
    IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TAUBM_N)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: BOTTOM_ROUGHNESS"

  RETURN
END SUBROUTINE BOTTOM_ROUGHNESS
!==============================================================================|
