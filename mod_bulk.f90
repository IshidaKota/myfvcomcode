










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

MODULE mod_bulk
  USE MOD_PREC
  USE ALL_VARS,ONLY : IPT
  USE ALL_VARS,ONLY : PAR, MYID, NPROCS
  USE MOD_PAR, ONLY : NC, EC, AEXCHANGE
  implicit none

  REAL, PARAMETER :: EPS = 1.0E-7_SP

CONTAINS


  !------------------------------------------------------------------
  ! Set the lower and upper limits of ocean roughness.
  ! ref: Davis et al. (2008)
  !
  ! input : ocean roughness (Z0, on cell)
  ! output: ocean roughness (Z0, on cell)
  !
  ! Siqi Li, 2021-01-27
  !------------------------------------------------------------------
  SUBROUTINE Z0_limits_Davis2008
    !
    USE ALL_VARS,  only : N, Z0
    !
    IMPLICIT NONE
    !
    REAL(SP) :: Z0_min, Z0_max
    INTEGER  :: I
    !
    Z0_min = 0.125E-6_SP
    Z0_max = 2.850E-3_SP
    !
    DO I = 1, N

      Z0(I) = MAX(Z0(I), Z0_min)
      Z0(I) = MIN(Z0(I), Z0_max)

    END DO 
    !
  END SUBROUTINE Z0_limits_Davis2008

  !------------------------------------------------------------------
  ! Add the ocean roughness for a smooth surface to Z0
  ! ref: Businger (1973)
  ! 
  ! input : air temperature (T_AIR, on node), with 1
  !         friction velocity (U_star on cell), with 1
  ! output: ocean roughness (Z0, on cell)
  !
  ! Siqi Li, 2021-01-27
  !------------------------------------------------------------------
  SUBROUTINE Z0_Businger
    !
    USE ALL_VARS,  only : N, NV, T_AIR, U_star, Z0
    !
    IMPLICIT NONE
    !
    REAL(SP) :: Ta, Vair
    INTEGER  :: I
    !
    DO I = 1, N

      Ta = ( T_AIR(NV(I,1))+T_AIR(NV(I,2))+T_AIR(NV(I,3)) ) / 3.0_SP
      Vair = 1.326E-5_SP * ( 1.0_SP + 6.542E-3_SP*Ta       +   &
                                      8.301E-6_SP*Ta*Ta    -   &
                                      4.840E-9_SP*Ta*Ta*Ta )

      Z0(I) = Z0(I) + 0.11_SP * Vair / max(U_star(I),EPS) 

    END Do
    !
  END SUBROUTINE Z0_Businger

  !------------------------------------------------------------------
  ! Calculate drag coefficient based on ocean roughness.
  ! ref: Charnock (1955)
  !
  ! input : wind height (ZUU)
  !         ocean roughness (Z0, on cell)
  ! output: drag coefficient (Cd, on cell)
  !
  ! Siqi Li, 2021-01-27
  !------------------------------------------------------------------
  SUBROUTINE Cd_on_Z0
    !
    USE ALL_VARS,  only : N, Z0, Cd, ZUU
    !
    IMPLICIT NONE
    !
    REAL(SP), PARAMETER :: KAPPA=0.4_SP
    INTEGER  :: I
    !
    DO I = 1, N

      Cd(I) = ( KAPPA / max(log(ZUU/Z0(I)),EPS) )**2.0_SP

    END DO
    !
  END SUBROUTINE Cd_on_Z0

  !------------------------------------------------------------------
  ! Calculate drag coefficient using LP1981 method.
  ! ref: Large and Pond (1981)
  !
  ! input : u-wind (UUWIND, on cell)
  !         v-wind (VVWIND, on cell)
  ! output: drag coefficient (Cd, on cell)
  !
  ! Siqi Li, 2021-01-27
  !------------------------------------------------------------------
  SUBROUTINE Cd_LP1981
    !
    USE ALL_VARS,  only : N, UUWIND, VVWIND, Cd
    !
    IMPLICIT NONE
    !
    REAL(SP) :: Cd_min, Cd_max, WSD
    INTEGER  :: I
    !
    Cd_min = 1.205E-3_SP
    Cd_max = 2.115E-3_SP
    !
    DO I = 1, N
      
      WSD = SQRT(UUWIND(I)**2_SP + VVWIND(I)**2_SP)
      IF (WSD<11) THEN
        Cd(I) = Cd_min
      ELSEIF (WSD>25) THEN
        Cd(I) = Cd_max
      ELSE
        Cd(I) = (0.49_SP  + 0.065_SP*WSD) * 1E-3_SP
      END IF

    END DO

    IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, Cd)

    !
  END SUBROUTINE Cd_LP1981




  !------------------------------------------------------------------
  ! Calculate wind stress.
  !
  ! input : u-wind (UUWIND, on cell)
  !         v-wind (VVWIND, on cell)
  !         drag coefficient (Cd, on cell)
  ! output: u-wind stress (WUSURF2, on cell)
  !         v-wind stress (WUSURF2, on cell)
  !
  ! Siqi Li, 2021-01-27
  !------------------------------------------------------------------
  SUBROUTINE TAU_on_Cd_WIND(MODE)
    !
    USE ALL_VARS,  only : N, UUWIND, VVWIND, Cd
    USE ALL_VARS,  only : WUSURF, WVSURF, WUSURF2, WVSURF2
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=3), INTENT(IN) :: MODE
    REAL(SP), PARAMETER :: RHO_air = 1.2_SP
    REAL(SP)  :: SPD 
    INTEGER   :: I
    
    SELECT CASE (MODE)
    CASE ('INT')

      DO I = 1, N

        SPD = SQRT(UUWIND(I)**2.0_SP +VVWIND(I)**2.0_SP)

        WUSURF(I) = RHO_air * Cd(I) * SPD * UUWIND(I)
        WVSURF(I) = RHO_air * Cd(I) * SPD * VVWIND(I)

      END DO

      IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, WUSURF)
      IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, WVSURF)

    CASE ('EXT')
        
      DO I = 1, N
         
        SPD = SQRT(UUWIND(I)**2.0_SP +VVWIND(I)**2.0_SP)

        WUSURF2(I) = RHO_air * Cd(I) * SPD * UUWIND(I)
        WVSURF2(I) = RHO_air * Cd(I) * SPD * VVWIND(I)

      END DO

      IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, WUSURF2)
      IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, WVSURF2)

    END SELECT

  END SUBROUTINE TAU_on_Cd_WIND

  !------------------------------------------------------------------
  ! Calculate TAU on x and y direction
  !
  ! Siqi Li, 2021-01-27
  !------------------------------------------------------------------
  SUBROUTINE TAU_to_XY(MODE)
    !
    USE ALL_VARS,  only : N, TAU_WIND, UUWIND, VVWIND
    USE ALL_VARS,  only : WUSURF, WVSURF, WUSURF2, WVSURF2
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=3), INTENT(IN) :: MODE
    REAL(SP)  :: stress1, stress2
    INTEGER   :: I

!write(IPT,*) TAU_WIND(154), WUSURF(154), WVSURF(154), WUSURF2(154), WVSURF2(154)
    SELECT CASE (MODE)
    CASE ('INT')

      DO I = 1, N

        stress1 = SQRT(UUWIND(I)**2.0_SP + VVWIND(I)**2.0_SP)
        stress2 = TAU_WIND(I)

        WUSURF(I) = UUWIND(I) *stress2/max(stress1,EPS)
        WVSURF(I) = VVWIND(I) *stress2/max(stress1,EPS)

      END DO

        IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, WUSURF)
        IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, WVSURF)

    CASE ('EXT')

      DO I = 1, N

        stress1 = SQRT(UUWIND(I)**2.0_SP + VVWIND(I)**2.0_SP)
        stress2 = TAU_WIND(I)

        WUSURF2(I) = UUWIND(I) *stress2/max(stress1,EPS)
        WVSURF2(I) = VVWIND(I) *stress2/max(stress1,EPS)

      END DO

       IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, WUSURF2)
       IF(PAR) CALL AEXCHANGE(EC, MYID, NPROCS, WVSURF2)

    END SELECT
!write(IPT,*) TAU_WIND(154), WUSURF(154), WVSURF(154), WUSURF2(154), WVSURF2(154)
    !
  END SUBROUTINE TAU_to_XY

  !------------------------------------------------------------------
  ! Update wind stress based on the following five methods:
  !   -> LP1981
  !   -> COARE
  !   -> TY2001
  !   -> OOST
  !   -> DGHQ
  !
  ! Siqi Li, 2021-01-27
  !------------------------------------------------------------------
  SUBROUTINE UPDATE_WINDSTRESS(MODE)
    !
    USE ALL_VARS,  only : WIND_STRESS_METHOD
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=3), INTENT(IN) :: MODE

    SELECT CASE (WIND_STRESS_METHOD)

    CASE ('LP1981')
      CALL Cd_LP1981
      CALL TAU_on_Cd_WIND(MODE)

    CASE ('COARE')
      CALL TAU_to_XY(MODE)


    CASE DEFAULT
      ! It is SAFE to do nothing here, because we have already
      ! checked the WIND_STRESS_METHOD in mod_force.F

    END SELECT
    !
  END SUBROUTINE UPDATE_WINDSTRESS

  SUBROUTINE ASIMPLE_DRAG(spdx,spdy,strx,stry)
    IMPLICIT NONE
    REAL(SP),ALLOCATABLE, TARGET, INTENT(IN)  :: SPDX(:),SPDY(:)
    REAL(SP),ALLOCATABLE, TARGET, INTENT(INOUT) :: STRX(:),STRY(:)

    REAL(SP), POINTER :: SPDXP(:),SPDYP(:)
    REAL(SP), POINTER :: STRXP(:),STRYP(:)
    
    SPDXP => SPDX
    SPDYP => SPDY

    STRXP => STRX
    STRYP => STRY
    
    CALL PSIMPLE_DRAG(SPDXP,SPDYP,STRXP,STRYP)


  END SUBROUTINE ASIMPLE_DRAG



  SUBROUTINE PSIMPLE_DRAG(spdx,spdy,strx,stry)
    IMPLICIT NONE
    REAL(SP), POINTER,INTENT(IN)  :: SPDX(:),SPDY(:)
    REAL(SP), POINTER, INTENT(INOUT) :: STRX(:),STRY(:)
    INTEGER :: I, N
    REAL(SP) :: CD, WDS, TX, TY


    IF(.not.Associated(SPDX)) WRITE(6,*) "SIMPLE DRAG: SPDX is not associated"
    IF(.not.Associated(SPDY)) WRITE(6,*) "SIMPLE DRAG: SPDY is not associated"
    IF(.not.Associated(STRX)) WRITE(6,*) "SIMPLE DRAG: STRX is not associated"
    IF(.not.Associated(STRY)) WRITE(6,*) "SIMPLE DRAG: STRY is not associated"

    N = UBOUND(SPDX,1)


    IF(N /= UBOUND(SPDY,1)) WRITE(6,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"
    IF(N /= UBOUND(STRY,1)) WRITE(6,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"
    IF(N /= UBOUND(STRX,1)) WRITE(6,*) "SIMPLE DRAG: MIS-MATCHED DIMENSIONS"

    DO I=1,N
       TX = SPDX(I)
       TY = SPDY(I)
       WDS=SQRT(TX*TX+TY*TY)
       CD=1.2E-3
       IF (WDS >= 11.0_SP) CD=(0.49_SP+0.065_SP*WDS)*1.E-3_SP
       IF (WDS >= 25.0_SP) CD=(0.49_SP+0.065_SP*25.0_SP)*1.E-3_SP

       STRX(I) = 1.2_SP*CD*TX*WDS
       STRY(I) = 1.2_SP*CD*TY*WDS

    END DO

  END SUBROUTINE PSIMPLE_DRAG
  



END MODULE mod_bulk
