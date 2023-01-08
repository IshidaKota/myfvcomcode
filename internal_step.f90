










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

SUBROUTINE INTERNAL_STEP
  USE MOD_NESTING
  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_OBCS
  USE MOD_TIME
  USE EQS_OF_STATE
  USE MOD_WD
  USE MOD_ASSIM
  USE MOD_PAR
  USE MOD_ICE
  USE MOD_ICE2D
  USE MOD_NORTHPOLE
! lwang added for SST mld assimilation
  USE MOD_MLD_RHO
! lwang
  USE MOD_DYE




  IMPLICIT NONE
  INTEGER :: I,J,K
  REAL(SP) :: UTMP,VTMP
  integer :: i1,i2

  REAL(SP), ALLOCATABLE :: TB1(:,:),SB1(:,:),Q2B(:,:),Q2LB(:,:),UB(:,:),VB(:,:)
  REAL(SP), ALLOCATABLE :: DYEB(:,:)
  REAL(SP) :: RCOFT
  INTEGER :: ITER
!------------------------------------------------------------------------------|

  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "Start: internal_step"

!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
  RAMP = 0.0_SP
  IF(IRAMP /= 0) THEN
     RAMP=TANH(real(IINT,sp)/real(IRAMP,sp))
  ELSE
     RAMP = 1.0_SP
  END IF

!----OFFLINE SEDIMENT UPDATE
!----loop end for offline sediment



!----SET UP WATER QUALITY MODEL COEFFICIENTS-----------------------------------!


  

  !----ADJUST CONSISTENCY BETWEEN 3-D VELOCITY AND VERT-AVERAGED VELOCITIES------!
    CALL ADJUST2D3D(1)

  !----SPECIFY THE SOLID BOUNDARY CONDITIONS OF U&V INTERNAL MODES---------------!
    CALL BCOND_GCN(5,0)

    IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
! end defined semi-implicit && nh

  !----SPECIFY THE SURFACE FORCING OF INTERNAL MODES-----------------------------!
  CALL BCOND_GCN(8,0)

  ! New Open Boundary Condition ----3

  !end !defined (TWO_D_MODEL)
  
  !----SPECIFY THE BOTTOM ROUGHNESS AND CALCULATE THE BOTTOM STRESSES------------!
  CALL BOTTOM_ROUGHNESS


!==============================================================================!
!  CALCULATE DISPERSION (GX/GY) AND BAROCLINIC PRESSURE GRADIENT TERMS         !
!==============================================================================!

    CALL ADVECTION_EDGE_GCN(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !

  IF(RECALCULATE_RHO_MEAN) THEN
     IF(RECALC_RHO_MEAN .LT. IntTime)THEN
        RECALC_RHO_MEAN = RECALC_RHO_MEAN + DELT_RHO_MEAN
        CALL RHO_PMEAN    
     END IF
  END IF


  IF(.NOT. BAROTROPIC)THEN                    !Barotropic Flow ?            !
     SELECT CASE(BAROCLINIC_PRESSURE_GRADIENT)
     CASE ("sigma levels")
        CALL BAROPG      !Sigma Level Pressure Gradient!
     CASE('z coordinates')
        CALL PHY_BAROPG  !Z Level Pressure Gradient    !
     CASE DEFAULT
        CALL FATAL_ERROR("UNKNOW BAROCLINIC PRESURE GRADIENT TYPE",&
             & TRIM(BAROCLINIC_PRESSURE_GRADIENT))
     END SELECT

  END IF                                      !                             !


  ADX2D = 0.0_SP ; ADY2D = 0.0_SP             !Initialize GX/GY Terms       !
  DRX2D = 0.0_SP ; DRY2D = 0.0_SP             !Initialize BCPG for Ext Mode !

  DO K=1,KBM1
     DO I=1, N
        ADX2D(I)=ADX2D(I)+ADVX(I,K)   !*DZ1(I,K)
        ADY2D(I)=ADY2D(I)+ADVY(I,K)   !*DZ1(I,K)
        DRX2D(I)=DRX2D(I)+DRHOX(I,K)  !*DZ1(I,K)
        DRY2D(I)=DRY2D(I)+DRHOY(I,K)  !*DZ1(I,K)
     END DO
  END DO

  CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
  ADX2D = ADX2D - ADVUA                       !Subtract to Form GX
  ADY2D = ADY2D - ADVVA                       !Subtract to Form GY

  !----INITIALIZE ARRAYS USED TO CALCULATE AVERAGE UA/E  OVER EXTERNAL STEPS-----!
  UARD = 0.0_SP
  VARD = 0.0_SP
  EGF  = 0.0_SP


!!# if defined (WET_DRY)
!!  UARDS = 0.0_SP
!!  VARDS = 0.0_SP
!!# endif

  IF(IOBCN > 0) THEN
     UARD_OBCN(1:IOBCN)=0.0_SP
  END IF
! end defined semi-implicit
! end defined (TWO_D_MODEL)


  ! New Open Boundary Condition ----4


!JQIGOTO# if defined (WAVE_ONLY)
!JQIGOTO  goto 5501
!JQIGOTO# endif

!JQIGOTO
!JQIGOTO
  !==============================================================================!
  !  LOOP OVER EXTERNAL TIME STEPS                                               !
  !==============================================================================!
  DO IEXT=1,ISPLIT

     IF (DBG_SET(DBG_SBRIO)) WRITE(IPT,*) "/// EXT SETP: ",IEXT

     ExtTime = ExtTime + IMDTE

     CALL EXTERNAL_STEP
     
     
  END DO


!  new update  !  jqi ggao 0730/2007  for E-P

!  new update  !  jqi ggao 0730/2007

! end defined semi-implicit

!==============================================================================!
!==============================================================================!
!                     BEGIN THREE D ADJUSTMENTS
!==============================================================================!
!==============================================================================!

  !==============================================================================!
  !    ADJUST INTERNAL VELOCITY FIELD TO CORRESPOND TO EXTERNAL                  !
  !==============================================================================!

    CALL ADJUST2D3D(2)

  ! New Open Boundary Condition ----9
  
  !==============================================================================!
  !     CALCULATE INTERNAL VELOCITY FLUXES                                       |
  !==============================================================================!
  !                                                    !
!!

!end !defined (TWO_D_MODEL)

!==============================================================================!
!     CALCULATE INTERNAL VELOCITY FLUXES                                       |
!==============================================================================!


    CALL VERTVL_EDGE     ! Calculate/Update Sigma Vertical Velocity (Omega)   !


    CALL VISCOF_H        ! Calculate horizontal diffusion coefficient scalars !
  
  IF(.NOT. RK_3D_ON)THEN
    CALL ADV_UV_EDGE_GCN ! Horizontal Advect/Diff + Vertical Advection        !

     CALL VDIF_UV      ! Implicit Integration of Vertical Diffusion of U/V  !

    IF(ADCOR_ON) THEN
      CALL ADCOR
     CALL VDIF_UV      ! Implicit Integration of Vertical Diffusion of U/V  !
!      CALL VDIF_UV   ! Implicit Integration of Vertical Diffusion of U/V  !

    ENDIF


    CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input           !




  ELSE
   ALLOCATE(UB(0:NT,KB))      ; UB = 0.0_SP
   ALLOCATE(VB(0:NT,KB))      ; VB = 0.0_SP
   UB = U
   VB = V
   RCOFT = 0.5_SP
   DO ITER = 1,2
    CALL ADV_UV_EDGE_GCN_RK(UB,VB) ! Horizontal Advect/Diff + Vertical Advection        !

     CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !

    IF(ADCOR_ON) THEN
      CALL ADCOR
      CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !
    ENDIF


    CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input

    U = UF
    V = VF
    IF(ITER /= 2)THEN
      U = RCOFT*U+(1.0_SP-RCOFT)*UB
      V = RCOFT*V+(1.0_SP-RCOFT)*VB
    END IF
   END DO
   DEALLOCATE(UB,VB)
  END IF

    IF(NESTING_ON )THEN
      CALL SET_VAR(intTime,U=UF)
      CALL SET_VAR(intTime,V=VF)
    END IF  


!if !defined (SEMI_IMPLICIT)

           
  CALL WREAL           ! Calculate True Vertical Velocity (W)               !

!  CALL REPORT("before VISCOF_H")

  !==============================================================================!
  !    TURBULENCE MODEL SECTION                                                  |
  !==============================================================================!
!     IF(PAR)CALL EXCHANGE(EC,NT,1,MYID,NPROCS,WUSURF,WVSURF)
!     IF(PAR)CALL EXCHANGE(EC,NT,1,MYID,NPROCS,WUBOT,WVBOT)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUSURF,WVSURF)
     IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WUBOT,WVBOT)

  
  SELECT CASE(VERTICAL_MIXING_TYPE)
  CASE('closure')
     !=================General Ocean Turbulence Model==========================!
     !===================Original FVCOM MY-2.5/Galperin 1988 Model=============!
   IF(.NOT. RK_3D_ON)THEN  

     CALL ADV_Q(Q2,Q2F)       !!Advection of Q2 

     CALL ADV_Q(Q2L,Q2LF) 

     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2F)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2LF)


     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2             !Conservation Correction   !
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2L            !Conservation Correction   !
!end !defined (ONE_D_MODEL)



     CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
!     IF(PAR)CALL EXCHANGE(NC,MT,KB,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange   !
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange   !
      Q2  = Q2F
      Q2L = Q2LF

  ELSE
    ALLOCATE(Q2B(0:MT,KB))      ;  Q2B = 0.0_SP
    ALLOCATE(Q2LB(0:MT,KB))     ;  Q2LB = 0.0_SP
    Q2B  = Q2
    Q2LB = Q2L
    RCOFT = 0.5_SP
    DO ITER = 1,2
     CALL ADV_Q_RK(Q2,Q2B,Q2F)       !!Advection of Q2 
     CALL ADV_Q_RK(Q2L,Q2LB,Q2LF)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2F)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,Q2LF)

     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2             !Conservation Correction   !
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_Q2L            !Conservation Correction   !

     CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange !
      Q2  = Q2F
      Q2L = Q2LF

      IF(ITER /= 2)THEN
        Q2 = RCOFT*Q2+(1.0_SP-RCOFT)*Q2B
        Q2L = RCOFT*Q2L+(1.0_SP-RCOFT)*Q2LB
      END IF
    END DO
    DEALLOCATE(Q2B,Q2LB) 
  END IF
! end if defined(GOTM)
  CASE('constant')
     KM = UMOL
     KH = UMOL*VPRNU
  END SELECT

  
!     IF(PAR)CALL EXCHANGE(NC,MT,KB,MYID,NPROCS,KM,KQ,KH)
  IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,KM,KQ,KH)
  CALL N2E3D(KM,KM1)




!==============================================================================!
!    SEDIMENT MODEL SECTION  
!      Advance Sed Model (Erode/Deposit/Advect/Diffuse)      
!      Change bathymetry (if MORPHO_MODEL=T)  
!        Note:  morph array returned from sed:  Deposition:  morph > 0
!        0 < morpho_factor < 1  
!      morpho_model and morpho_factor are stored and set in mod_sed.F              
!==============================================================================!
  

  !==============================================================================!
  !    UPDATE TEMPERATURE IN NON-BAROTROPIC CASE                                 !
  !==============================================================================!
  IF(TEMPERATURE_ACTIVE)THEN
  IF(.NOT. RK_3D_ON)THEN
     
     CALL ADV_T                                     !Advection                 !

     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,TF1)


     !#                                                   if !defined (1)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_T            !Conservation Correction   !
     !#                                                   endif
        !end !defined (ONE_D_MODEL)                                    !                          !

     IF(CASENAME(1:3) == 'gom')THEN
        CALL VDIF_TS_GOM(1,TF1)
     ELSE  
        CALL VDIF_TS(1,TF1)                            !Vertical Diffusion        !
     END IF

    
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TF1) !Interprocessor Exchange   !
     CALL BCOND_TS(1)                               !Boundary Conditions       !

  ELSE
    ALLOCATE(TB1(0:MT,KB))        ; TB1 = 0.0_SP
    TB1 = T1
    RCOFT = 0.5_SP
    DO ITER = 1,2
     CALL ADV_T_RK(TB1)                                     !Advection                 !
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,TF1)

     !#                                                   if !defined !(1)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_T            !Conservation Correction   !
     !#                                                   endif

     IF(CASENAME(1:3) == 'gom')THEN
        CALL VDIF_TS_GOM(1,TF1)
     ELSE
        CALL VDIF_TS(1,TF1)                            !Vertical Diffusion
     END IF

     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,TF1) !Interprocessor Exchange   !
     CALL BCOND_TS(1)                               !Boundary Conditions       !
     
     T1 = TF1
     IF(ITER /= 2)THEN
       T1 = RCOFT*T1+(1.0_SP-RCOFT)*TB1
     END IF
    END DO
    DEALLOCATE(TB1)
  ENDIF

     IF(NESTING_ON )THEN
       CALL SET_VAR(intTime,T1=TF1)
     END IF

!!$!QXU{
!!$#  if defined (SST_GRID_ASSIM)
!!$      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM_GRD)CALL TEMP_NUDGING
!!$      IF(ASSIM_FLAG==1 .AND. SST_ASSIM_GRD)CALL TEMP_NUDGING
!!$#  else
!!$      IF(ASSIM_FLAG==0 .AND. .NOT. SST_ASSIM)CALL TEMP_NUDGING
!!$      IF(ASSIM_FLAG==1 .AND. SST_ASSIM)CALL TEMP_NUDGING
!!$#  endif
!!$!QXU}
     




!     IF(NESTING_ON )THEN
!       CALL SET_VAR(intTime,T1=TF1)
!     END IF

!J. Ge for tracer advection
     IF(BACKWARD_ADVECTION .EQV. .TRUE.)THEN
       IF(BACKWARD_STEP==1)THEN
        T0 = T1
       ELSEIF(BACKWARD_STEP==2)THEN
         T2 = T0
         T0 = T1
       ENDIF
     ENDIF
!J. Ge for tracer advection
     T1 = TF1                                       !Update to new time level  !

     CALL N2E3D(T1,T)                               !Shift to Elements         !


  END IF                                         !                          !

  !==============================================================================!
  !    UPDATE SALINITY IN NON-BAROTROPIC CASE                                    !
  !==============================================================================!

  IF(SALINITY_ACTIVE)THEN
  IF(.NOT. RK_3D_ON)THEN


     CALL ADV_S                                     !Advection                 !

     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SF1)
     

     !#                                                   if !defined (1)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_S       !Conservation Correction   !
     !#                                                   endif
     !end !defined (ONE_D_MODEL)                                    !                          !

     IF(CASENAME(1:3) == 'gom')THEN
        CALL VDIF_TS_GOM(2,SF1)
     ELSE  
        CALL VDIF_TS(2,SF1)                            !Vertical Diffusion        !
     END IF

     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SF1) !Interprocessor Exchange   !

     CALL BCOND_TS(2)                               !Boundary Conditions       !

  ELSE
    ALLOCATE(SB1(0:MT,KB))         ; SB1 = 0.0_SP
    SB1 = S1
    RCOFT = 0.5_SP
    DO ITER = 1,2
     CALL ADV_S_RK(SB1)                                     !Advection                 !
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SF1)

     !#                                                   if !defined     !(1)
     IF(SCALAR_POSITIVITY_CONTROL) CALL FCT_S       !Conservation Correction   !
     !#                                                   endif

     IF(CASENAME(1:3) == 'gom')THEN
        CALL VDIF_TS_GOM(2,SF1)
     ELSE
        CALL VDIF_TS(2,SF1)                            !Vertical Diffusion
     END IF

     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SF1) !Interprocessor Exchange   !

     CALL BCOND_TS(2)                               !Boundary Conditions       !

     S1 = SF1
     IF(ITER /= 2)THEN
       S1 = RCOFT*S1+(1.0_SP-RCOFT)*SB1
     END IF
    END DO
    DEALLOCATE(SB1)
  END IF

     IF(NESTING_ON )THEN
       CALL SET_VAR(intTime,S1=SF1)
     END IF




!     IF(NESTING_ON )THEN
!       CALL SET_VAR(intTime,S1=SF1)
!     END IF

!J. Ge for tracer advection
     IF(BACKWARD_ADVECTION .EQV. .TRUE.)THEN
       IF(BACKWARD_STEP==1)THEN
        S0 = S1
       ELSEIF(BACKWARD_STEP==2)THEN
         S2 = S0
         S0 = S1
       ENDIF
     ENDIF
!J. Ge for tracer advection
     S1 = SF1                                     !Update to new time level  !

     CALL N2E3D(S1,S)                               !Shift to Elements         !


  END IF                      

  !==============================================================================!
  
  !==============================================================================!
  !    UPDATE DYE IN NON-BAROTROPIC CASE                                         !
  !==============================================================================!
  !     IF(DYE_ON.AND.IINT.GE.IINT_SPE_DYE_B) THEN     !                          !                          !
  IF(DYE_ON) THEN     !
   IF(.NOT. RK_3D_ON)THEN                          !                          !
     CALL ADV_DYE                                   !Advection                 !
     
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,DYEF)
     
     
     
     CALL VDIF_DYE(DYEF)                            !Vertical Diffusion        !
     
     !check
     !!     DYE = DYEF                                     !Update to new time level  !
     !!     IF(IINT.GE.IINT_SPE_DYE_B) CALL ARCHIVE
     !check    
     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,DYEF) !Interprocessor Exchange   !
     CALL BCOND_DYE                                 !Boundary Conditions       !
     DYE = DYEF                                     !Update to new time level  !
     

      !!     IF(MSR) WRITE(IPT,*) 'CALL Dye_on--iint=',iint,IINT_SPE_DYE_B
   ELSE
     ALLOCATE(DYEB(0:MT,KB))     ; DYEB=0.0_SP
     DYEB = DYE
     RCOFT = 0.5_SP
     DO ITER = 1,2
     CALL ADV_DYE_RK(DYEB)                                   !Advection                 !
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,DYEF)

     CALL VDIF_DYE(DYEF)                            !Vertical Diffusion        !

     IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,DYEF) !Interprocessor Exchange   !
     CALL BCOND_DYE                                 !Boundary Conditions       !
     DYE = DYEF                                     !Update to new time level  !

     IF(ITER /= 2)THEN
       DYE = RCOFT*DYE+(1.0_SP-RCOFT)*DYEB
     ENDIF
    END DO
    DEALLOCATE(DYEB)
   END IF
  END IF                                         !                          !
  !==================================================================================!
  !endif defined (DYE)

!==================================================================================!
!    ADJUST TEMPERATURE AND SALINITY AT RIVER MOUTHS
!==================================================================================!
  




  
  
  !==============================================================================!
  !     UPDATE THE DENSITY IN NON-BAROTROPIC CASE                                |
  !==============================================================================!
  IF(.NOT.BAROTROPIC)THEN
   SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
   CASE(SW_DENS1)
      CALL DENS1
   CASE(SW_DENS2)
      CALL DENS2
   CASE(SW_DENS3)
      CALL DENS3
   CASE DEFAULT
      CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
           & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
   END SELECT
  END IF
  !==============================================================================!
  !     MIMIC CONVECTIVE OVERTURNING TO STABILIZE VERTICAL DENSITY PROFILE       |
  !==============================================================================!
  
  IF(CONVECTIVE_OVERTURNING)THEN
     CALL CONV_OVER
     IF(.NOT.BAROTROPIC)THEN
        SELECT CASE(SEA_WATER_DENSITY_FUNCTION)
        CASE(SW_DENS1)
           CALL DENS1
        CASE(SW_DENS2)
           CALL DENS2
        CASE(SW_DENS3)
           CALL DENS3
        CASE DEFAULT
           CALL FATAL_ERROR("INVALID DENSITY FUNCTION SELECTED:",&
                & "   "//TRIM(SEA_WATER_DENSITY_FUNCTION) )
        END SELECT
        
     END IF
  END IF
!==============================================================================!
!==============================================================================!
  ! end if !defined (TWO_D_MODEL)
!==============================================================================!
!==============================================================================!
  
  

  !==============================================================================!
  !     UPDATE VELOCITY FIELD (NEEDED TO WAIT FOR SALINITY/TEMP/TURB/TRACER)     |
  !==============================================================================!
  U = UF
  V = VF


  !==============================================================================!
  !    PERFORM DATA EXCHANGE FOR ELEMENT BASED INFORMATION AT PROC BNDRIES       |
  !==============================================================================!
  
  IF(PAR)THEN
     CALL AEXCHANGE(EC,MYID,NPROCS,U,V)
     CALL AEXCHANGE(NC,MYID,NPROCS,Q2,Q2L)
     CALL AEXCHANGE(EC,MYID,NPROCS,RHO,T,S)
     CALL AEXCHANGE(NC,MYID,NPROCS,S1,T1,RHO1)
     
     CALL AEXCHANGE(NC,MYID,NPROCS,DYE)
     
     CALL AEXCHANGE(EC,MYID,NPROCS,VISCOFM)
     CALL AEXCHANGE(NC,MYID,NPROCS,VISCOFH)
  END IF

!---> Siqi Li, @20210809  
  IF (NC_MLD) THEN
    IF (abs(NC_DAT%FTIME%NEXT_IO -IntTime)<0.1_SP*IMDTI) THEN
      CALL MLD_RHO_CAL(M, KBM1, RHO1(1:M,:), MLD_OUT(1:M)) 
      IF (PAR) CALL AEXCHANGE(NC,MYID,NPROCS,MLD_OUT) 
    END IF
  END IF
!<--- 

  !==============================================================================!
  !     PERFORM DATA EXCHANGE FOR WATER QUALITY VARIABLES                        |
  !==============================================================================!
  ! end if !defined (TWO_D_MODEL)
  


  !
  !----SHIFT SEA SURFACE ELEVATION AND DEPTH TO CURRENT TIME LEVEL---------------!
  !
  ET  = EL  
  DT  = D 
  ET1 = EL1
  DT1 = D1
  
  IF(WETTING_DRYING_ON) CALL WD_UPDATE(3)
  
  ! New Open Boundary Condition ----10

!JQIGOTO# if defined (WAVE_ONLY)
!JQIGOTO   5501 continue
!JQIGOTO# endif
!JQIGOTO
!JQIGOTO

!!$  CALL DUMP_PROBE_DATA 


  if(dbg_set(dbg_sbr)) write(ipt,*)&
       & "End: internal_step"

END SUBROUTINE INTERNAL_STEP
