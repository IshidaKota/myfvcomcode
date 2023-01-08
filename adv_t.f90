










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
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|

SUBROUTINE ADV_T               

  !------------------------------------------------------------------------------|

  USE ALL_VARS
  USE MOD_UTILS
  USE MOD_PAR
  USE BCS
  USE MOD_OBCS
  USE MOD_WD 
  USE MOD_SPHERICAL
  USE MOD_NORTHPOLE


  IMPLICIT NONE
  REAL(SP), DIMENSION(0:MT,KB)      :: XFLUX,XFLUX_ADV,RF
  REAL(SP), DIMENSION(0:MT)         :: PUPX,PUPY,PVPX,PVPY  
  REAL(SP), DIMENSION(0:MT)         :: PTPX,PTPY,PTPXD,PTPYD,VISCOFF
  REAL(SP), DIMENSION(3*(NT),KBM1)  :: DTIJ 
  REAL(SP), DIMENSION(3*(NT),KBM1)  :: UVN
  REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1 !,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
  REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN,TTIME,ZDEP
  REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT,STPOINT1,STPOINT2
  REAL(SP) :: FACT,FM1
  INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
  REAL(SP) :: T1MIN, T1MAX, T2MIN, T2MAX


!!$#  if defined (SPHERICAL)
!!$  REAL(DP) TY,TXPI,TYPI
!!$  REAL(DP) :: XTMP1,XTMP
!!$  REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
!!$  REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
!!$  REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
!!$  REAL(DP) :: TXPI_TMP,TYPI_TMP
!!$#  endif

  REAL(SP) :: TMIN,TMAX,XXXX
  REAL(SP), DIMENSION(0:MT,KB)     :: T1_S      !! temporary temperature in modified upwind
  REAL(SP), DIMENSION(0:MT,KB)     :: T1_SF     !! temporary temperature in modified upwind
  REAL(SP), DIMENSION(0:MT,KB)     :: WWWS     
  REAL(SP), DIMENSION(0:MT,KB)     :: WWWSF   
  REAL(SP), DIMENSION(0:MT)        :: DTWWWS  
  REAL(SP), DIMENSION(0:MT,KB)     :: ZZZFLUX   !! temporary total flux in corrected part
  REAL(SP), DIMENSION(0:MT,KB)     :: BETA      !! temporary beta coefficient in corrected part
  REAL(SP), DIMENSION(0:MT,KB)     :: BETAIN    !! temporary beta coefficient in corrected part
  REAL(SP), DIMENSION(0:MT,KB)     :: BETAOUT   !! temporary beta coefficient in corrected part
  REAL(SP), DIMENSION(0:MT,KB)     :: T1_FRESH  !! for source term which also bring mass volume
  REAL(SP), DIMENSION(0:MT,KB)     :: OFFS      !! Offset to temperature field to avoid values less than zero.
  INTEGER ITERA, NTERA

   REAL(SP)  CONV_T(1:KB), DISS_T(1:KB)
   REAL(SP)  SL_H(0:KB), T_TEMP(0:KB)
   REAL(SP)  SL_U, SL_F


  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"Start: adv_t"
  !------------------------------------------------------------------------------!

  SELECT CASE(HORIZONTAL_MIXING_TYPE)
  CASE ('closure')
     FACT = 1.0_SP
     FM1  = 0.0_SP
  CASE('constant')
     FACT = 0.0_SP
     FM1  = 1.0_SP
  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOW HORIZONTAL MIXING TYPE:",&
          & TRIM(HORIZONTAL_MIXING_TYPE) )
  END SELECT

  OFFS = 5.0_SP

  !
  !--Initialize Fluxes-----------------------------------------------------------!
  !
  XFLUX     = 0.0_SP
  XFLUX_ADV = 0.0_SP

  !
  !--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
  !
!!#  if !defined (WET_DRY)
  DO I=1,NCV
     I1=NTRG(I)
     !     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
        DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
        ! USE U,V
        UVN(I,K)  = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I) 
     END DO
  END DO
!!#  else
!!  DO I=1,NCV
!!     I1=NTRG(I)
!!     !     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!        DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
!!        ! USE US,VS
!!        UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!#      if defined (SEMI_IMPLICIT)
!!       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
!!       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
!!#      endif
!!     END DO
!!  END DO
!!#  endif

  !
  !--Add the Shortwave Radiation Body Force--------------------------------------!
  !
  RF = 0.0_sp
  IF(HEATING_ON) THEN
   !<-----ishid debug
 !  write(ipt,*) "now we are on body case"
   !>-----ishid debug
     SELECT CASE(HEATING_TYPE)
     CASE('body')

        ! REMOVED START TIME FOR BODY HEATING
        !     TTIME=THOUR
        !     if(ttime < thour_hs) then
        !       RF=0.0_SP
        !     else

        DO  K=1,KBM1
           DO  I=1,M
              !         IF(TTIME < THOUR_HS)THEN
              !           RF(I,K)=0.0_SP
              !         ELSE
              ! REMOVED DEPTH CHECK FROM ECOM-SI, FVCOM DOES NOT NEED IT!
              !          ZEDP = 0.0_SP
              !           IF(DT(I) > 0.0_SP) THEN
              ZDEP=0.5_SP*(Z(I,K)+Z(I,K+1))*DT(I)
              !           END IF
!JQI              RF(I,K)=-SWRAD(I)*((RHEAT/ZETA1)*EXP(ZDEP/ZETA1) &
!JQI                   +((1-RHEAT)/ZETA2)*EXP(ZDEP/ZETA2))*DT(I)
              RF(I,K)=-SWRAD(I)*((RHEAT/ZETA1)*EXP(ZDEP/ZETA1) &
                   +((1-RHEAT)/ZETA2)*EXP(ZDEP/ZETA2))*DT(I)*DZ(I,K)
!CHEN              RF(I,K)=-SWRAD(I)*((RHEAT/(ZETA1*ZETA1))*EXP(ZDEP/ZETA1) &
!CHEN                   +((1-RHEAT)/(ZETA2*ZETA2))*EXP(ZDEP/ZETA2))*DT(I)*DZ(I,K)
              !         END IF
           END DO
        END DO
        !     endif

     CASE('flux')
        RF = 0.0_SP
     CASE('surface')
        RF = 0.0_SP
     CASE DEFAULT
        CALL FATAL_ERROR('The surface heating type is set incorrectly:',&
             & TRIM(HEATING_TYPE))
     END SELECT
  END IF

!--ADJUST VOLUME'S HEAT CONTENT FOR EVAPORATION AND PRECIPITATION ------------------!
  IF (PRECIPITATION_ON) THEN
     RF(:,1)=RF(:,1)+ROFVROS*(QEVAP+QPREC)*T1(:,1)
  END IF

    ! Adding offset to avoid less-than-zero problems
      T1 = T1+OFFS
      TMEAN1 = TMEAN1+OFFS

  !
  !--Calculate the Advection and Horizontal Diffusion Terms----------------------!
  !

  DO K=1,KBM1
     PTPX  = 0.0_SP
     PTPY  = 0.0_SP
     PTPXD = 0.0_SP
     PTPYD = 0.0_SP
     DO I=1,M
        DO J=1,NTSN(I)-1
           I1=NBSN(I,J)
           I2=NBSN(I,J+1)
	 
!J. Ge for tracer advection	 
       IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN
         FFD=0.5_SP*(T1(I1,K)+T1(I2,K)-TMEAN1(I1,K)-TMEAN1(I2,K))
         FF1=0.5_SP*(T1(I1,K)+T1(I2,K))
       ELSE
         IF(BACKWARD_STEP==1)THEN
           FFD=0.5_SP*((T0(I1,K)+T1(I1,K))*0.5+(T0(I2,K)+T1(I2,K))*0.5-TMEAN1(I1,K)-TMEAN1(I2,K))
           FF1=0.5_SP*((T0(I1,K)+T1(I1,K))*0.5+(T0(I2,K)+T1(I2,K))*0.5)
         ELSEIF(BACKWARD_STEP==2)THEN
           FFD=0.5_SP*((T2(I1,K)+T0(I1,K)+T1(I1,K))/3.0_SP+(T2(I2,K)+T0(I2,K)+T1(I2,K))/3.0_SP-TMEAN1(I1,K)-TMEAN1(I2,K))
           FF1=0.5_SP*((T2(I1,K)+T0(I1,K)+T1(I1,K))/3.0_SP+(T2(I2,K)+T0(I2,K)+T1(I2,K))/3.0_SP)
         ENDIF
       ENDIF
!J. Ge for tracer advection
	 
!!$#        if defined (SPHERICAL)
!!$           XTMP  = VX(I2)*TPI-VX(I1)*TPI
!!$           XTMP1 = VX(I2)-VX(I1)
!!$           IF(XTMP1 >  180.0_SP)THEN
!!$              XTMP = -360.0_SP*TPI+XTMP
!!$           ELSE IF(XTMP1 < -180.0_SP)THEN
!!$              XTMP =  360.0_SP*TPI+XTMP
!!$           END IF
!!$           TXPI=XTMP*COS(DEG2RAD*VY(I))
!!$           TYPI=(VY(I1)-VY(I2))*tpi
!!$           ! ERROR HERE
!!$           !#    if defined (NORTHPOLE)
!!$           IF(NODE_NORTHAREA(I) == 1)THEN
!!$              VX1_TMP = REARTH * COS(VY(I1)*DEG2RAD) * COS(VX(I1)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+SIN(VY(I1)*DEG2RAD))
!!$              VY1_TMP = REARTH * COS(VY(I1)*DEG2RAD) * SIN(VX(I1)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+SIN(VY(I1)*DEG2RAD))
!!$
!!$              VX2_TMP = REARTH * COS(VY(I2)*DEG2RAD) * COS(VX(I2)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+SIN(VY(I2)*DEG2RAD))
!!$              VY2_TMP = REARTH * COS(VY(I2)*DEG2RAD) * SIN(VX(I2)*DEG2RAD) &
!!$                   * 2._SP /(1._SP+SIN(VY(I2)*DEG2RAD))
!!$
!!$              TXPI = (VX2_TMP-VX1_TMP)/(2._SP /(1._SP+SIN(VY(I)*DEG2RAD)))
!!$              TYPI = (VY1_TMP-VY2_TMP)/(2._SP /(1._SP+SIN(VY(I)*DEG2RAD)))
!!$              IF(I /= NODE_NORTHPOLE)THEN
!!$                 TXPI_TMP = TYPI*COS(VX(I)*DEG2RAD)-TXPI*SIN(VX(I)*DEG2RAD)
!!$                 TYPI_TMP = TXPI*COS(VX(I)*DEG2RAD)+TYPI*SIN(VX(I)*DEG2RAD)
!!$                 TYPI_TMP = -TYPI_TMP
!!$
!!$                 TXPI = TXPI_TMP
!!$                 TYPI = TYPI_TMP
!!$              END IF
!!$           END IF
!!$           ! END ERROR
!!$           
!!$           PTPX(I)=PTPX(I)+FF1*TYPI
!!$           PTPY(I)=PTPY(I)+FF1*TXPI
!!$           PTPXD(I)=PTPXD(I)+FFD*TYPI
!!$           PTPYD(I)=PTPYD(I)+FFD*TXPI
!!$#        else
!!$           PTPX(I)=PTPX(I)+FF1*(VY(I1)-VY(I2))
!!$           PTPY(I)=PTPY(I)+FF1*(VX(I2)-VX(I1))
!!$           PTPXD(I)=PTPXD(I)+FFD*(VY(I1)-VY(I2))
!!$           PTPYD(I)=PTPYD(I)+FFD*(VX(I2)-VX(I1))
!!$#        endif

           
           PTPX(I)=PTPX(I)+FF1*DLTYTRIE(i,j)
           PTPY(I)=PTPY(I)+FF1*DLTXTRIE(i,j)
           PTPXD(I)=PTPXD(I)+FFD*DLTYTRIE(i,j)
           PTPYD(I)=PTPYD(I)+FFD*DLTXTRIE(i,j)

        END DO
! gather all neighboring control volumes connecting at dam node 
        PTPX(I)=PTPX(I)/ART2(I)
        PTPY(I)=PTPY(I)/ART2(I)
        PTPXD(I)=PTPXD(I)/ART2(I)
        PTPYD(I)=PTPYD(I)/ART2(I)

     END DO

     IF(K == KBM1)THEN
        DO I=1,M
           PFPXB(I) = PTPX(I)
           PFPYB(I) = PTPY(I)
        END DO
     END IF

     DO I=1,M

        VISCOFF(I) = VISCOFH(I,K)

     END DO
     IF(K == KBM1) THEN
        AH_BOTTOM(1:M) = (FACT*VISCOFF(1:M) + FM1) * NN_HVC(1:M)
     END IF


     DO I=1,NCV_I
        IA=NIEC(I,1)
        IB=NIEC(I,2)


!!$        XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
!!$        YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
!!$#      if defined (SPHERICAL)
!!$        X1_DP=XIJE(I,1)
!!$        Y1_DP=YIJE(I,1)
!!$        X2_DP=XIJE(I,2)
!!$        Y2_DP=YIJE(I,2)
!!$        CALL ARCC(X2_DP,Y2_DP,X1_DP,Y1_DP,XII,YII)
!!$        XI=XII		
!!$        XTMP  = XI*TPI-VX(IA)*TPI
!!$        XTMP1 = XI-VX(IA)
!!$        IF(XTMP1 >  180.0_SP)THEN
!!$           XTMP = -360.0_SP*TPI+XTMP
!!$        ELSE IF(XTMP1 < -180.0_SP)THEN
!!$           XTMP =  360.0_SP*TPI+XTMP
!!$        END IF
!!$        DXA=XTMP*COS(DEG2RAD*VY(IA))  
!!$        DYA=(YI-VY(IA))*TPI
!!$
!!$        XTMP  = XI*TPI-VX(IB)*TPI
!!$        XTMP1 = XI-VX(IB)
!!$        IF(XTMP1 >  180.0_SP)THEN
!!$           XTMP = -360.0_SP*TPI+XTMP
!!$        ELSE IF(XTMP1 < -180.0_SP)THEN
!!$           XTMP =  360.0_SP*TPI+XTMP
!!$        END IF
!!$
!!$        DXB=XTMP*COS(DEG2RAD*VY(IB))
!!$        DYB=(YI-VY(IB))*TPI
!!$#      else
!!$        DXA=XI-VX(IA)
!!$        DYA=YI-VY(IA)
!!$        DXB=XI-VX(IB)
!!$        DYB=YI-VY(IB)
!!$#      endif
!!$        FIJ1=T1(IA,K)+DXA*PTPX(IA)+DYA*PTPY(IA)
!!$        FIJ2=T1(IB,K)+DXB*PTPX(IB)+DYB*PTPY(IB)

!J. Ge for tracer advection
        IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN
          FIJ1=T1(IA,K)+DLTXNCVE(I,1)*PTPX(IA)+DLTYNCVE(I,1)*PTPY(IA)
          FIJ2=T1(IB,K)+DLTXNCVE(I,2)*PTPX(IB)+DLTYNCVE(I,2)*PTPY(IB)
        ELSE
          IF(BACKWARD_STEP==1)THEN
            FIJ1=(T0(IA,K)+T1(IA,K))*0.5+DLTXNCVE(I,1)*PTPX(IA)+DLTYNCVE(I,1)*PTPY(IA)
            FIJ2=(T0(IB,K)+T1(IB,K))*0.5+DLTXNCVE(I,2)*PTPX(IB)+DLTYNCVE(I,2)*PTPY(IB)
          ELSEIF(BACKWARD_STEP==2)THEN
            FIJ1=(T2(IA,K)+T0(IA,K)+T1(IA,K))/3.0_SP+DLTXNCVE(I,1)*PTPX(IA)+DLTYNCVE(I,1)*PTPY(IA)
            FIJ2=(T2(IA,K)+T0(IB,K)+T1(IB,K))/3.0_SP+DLTXNCVE(I,2)*PTPX(IB)+DLTYNCVE(I,2)*PTPY(IB)
          ENDIF
        ENDIF
!J. Ge for tracer advection

        T1MIN=MINVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
        T1MIN=MIN(T1MIN, T1(IA,K))
        T1MAX=MAXVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
        T1MAX=MAX(T1MAX, T1(IA,K))
        T2MIN=MINVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
        T2MIN=MIN(T2MIN, T1(IB,K))
        T2MAX=MAXVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
        T2MAX=MAX(T2MAX, T1(IB,K))
        IF(FIJ1 < T1MIN) FIJ1=T1MIN
        IF(FIJ1 > T1MAX) FIJ1=T1MAX
        IF(FIJ2 < T2MIN) FIJ2=T2MIN
        IF(FIJ2 > T2MAX) FIJ2=T2MAX

        UN=UVN(I,K)

!        VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)
        ! David moved HPRNU and added HVC
        VISCOF=(FACT*0.5_SP*(VISCOFF(IA)*NN_HVC(IA)+VISCOFF(IB)*NN_HVC(IB)) + FM1*0.5_SP*(NN_HVC(IA)+NN_HVC(IB)))

!JQI NOV2021
        TXX=0.5_SP*(PTPXD(IA)+PTPXD(IB))*VISCOF
        TYY=0.5_SP*(PTPYD(IA)+PTPYD(IB))*VISCOF
!JQI NOV2021

        FXX=-DTIJ(I,K)*TXX*DLTYE(I)
        FYY= DTIJ(I,K)*TYY*DLTXE(I)

        EXFLUX=-UN*DTIJ(I,K)* &
             ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

        XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
        XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

        XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
        XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)


     END DO


  END DO !! K LOOP

  IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,XFLUX,XFLUX_ADV)

  DO K=1,KBM1
     IF(IOBCN > 0) THEN
        DO I=1,IOBCN
           I1=I_OBC_N(I)
           XFLUX_OBC(I,K)=XFLUX_ADV(I1,K)
        END DO
     END IF
  END DO



  !--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
  !
  !--------------------------------------------------------------------------------
  !   S. HU
  !   Using smolarkiewicz, P. K; A fully multidimensional positive definite advection
  !   transport algorithm with small implicit diffusion, Journal of Computational
  !   Physics, 54, 325-362, 1984
  !--------------------------------------------------------------------------------

  !-----combine all the horizontal flux first-----------------------------------

  !-------fresh water part-----------------------------------------------------
  T1_FRESH=T1

  IF(RIVER_TS_SETTING == 'calculated') THEN
     IF(RIVER_INFLOW_LOCATION == 'node') THEN
        IF(NUMQBC >= 1) THEN
           DO J=1,NUMQBC
              JJ=INODEQ(J)
              STPOINT=TDIS(J)+OFFS(1,1)
              DO K=1,KBM1
              !   STPOINT    = T1(JJ,K)
                 T1_FRESH(JJ,K)=TDIS(J)+OFFS(1,1)            ! for gauge guess in the following
                 XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT     !/DZ(JJ,K)
              END DO
           END DO
        END IF
     ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
        IF(NUMQBC > 0) THEN
           DO J=1,NUMQBC
              J1=N_ICELLQ(J,1)
              J2=N_ICELLQ(J,2)
              STPOINT=TDIS(J)+OFFS(1,1)
              DO K=1,KBM1
              !   STPOINT1 = T1(J1,K)
              !   STPOINT2 = T1(J2,K)
                 T1_FRESH(J1,K)=TDIS(J)+OFFS(1,1)
                 T1_FRESH(J2,K)=TDIS(J)+OFFS(1,1)
                 XFLUX(J1,K)=XFLUX(J1,K)-  &
                      QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT    !1     !/DZ1(J1,K)
                 XFLUX(J2,K)=XFLUX(J2,K)-  &      
                      QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT    !2     !/DZ1(J2,K)
              END DO
           END DO
        END IF
     END IF
  END IF

  ! We remove the offset at the end of the smolarkiewicz loop

  !------------------------------------------------------------------------------

  ! The horizontal term of advection is neglected here
  DO K=1,KBM1
     DO I=1,M
        IF(ISONB(I) == 2) THEN
           XFLUX(I,K)=0.
        ENDIF
     END DO
  END DO

  ! Initialize variables of 1
  T1_S=0._SP
  T1_SF=0._SP
  WWWS=0._SP
  WWWSF=0._SP
  DTWWWS=0._SP
  ZZZFLUX=0._SP
  BETA=0._SP
  BETAIN=0._SP
  BETAOUT=0._SP

  !!   first loop for vertical upwind
  !!   flux including horizontal and vertical upwind
  DO K=1,KBM1
     DO I=1,M
           IF(K == 1) THEN
              TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*T1(I,K)   &
                   -(WTS(I,K+1)+ABS(WTS(I,K+1)))*T1(I,K+1) &
                   +(WTS(I,K)+ABS(WTS(I,K)))*T1(I,K)    
           ELSE IF(K == KBM1) THEN
              TEMP = +(WTS(I,K)-ABS(WTS(I,K)))*T1(I,K-1)     &
                   +(WTS(I,K)+ABS(WTS(I,K)))*T1(I,K)
           ELSE
              TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*T1(I,K)   &
                   -(WTS(I,K+1)+ABS(WTS(I,K+1)))*T1(I,K+1) &
                   +(WTS(I,K)-ABS(WTS(I,K)))*T1(I,K-1)     &
                   +(WTS(I,K)+ABS(WTS(I,K)))*T1(I,K)
           END IF
           TEMP = 0.5_SP*TEMP 

           IF(K == 1)THEN
              TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
              TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
              TMAX = MAX(TMAX,T1(I,K+1),T1(I,K),T1_FRESH(I,K))
              TMIN = MIN(TMIN,T1(I,K+1),T1(I,K),T1_FRESH(I,K))
           ELSEIF(K == KBM1) THEN
              TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
              TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
              TMAX = MAX(TMAX,T1(I,K-1),T1(I,K),T1_FRESH(I,K))
              TMIN = MIN(TMIN,T1(I,K-1),T1(I,K),T1_FRESH(I,K))
           ELSE
              TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
              TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
              TMAX = MAX(TMAX,T1(I,K+1),T1(I,K-1),T1(I,K),T1_FRESH(I,K))
              TMIN = MIN(TMIN,T1(I,K+1),T1(I,K-1),T1(I,K),T1_FRESH(I,K))
           END IF

         ! Total (horizontal + vertical) flux to the cell
           ZZZFLUX(I,K) = TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K)/ART1(I)*(DTI/DT(I))/DZ(I,K)

         ! Updated salinity without limiter minus current salinity
           XXXX = ZZZFLUX(I,K)*DT(I)/DTFA(I)+T1(I,K)-T1(I,K)*DT(I)/DTFA(I)
           
           ! Added by Akvaplan 2018 to avoid single precision-crash
         IF((ABS(XXXX).LT.ABS(TMAX-T1(I,K))).AND.(XXXX.LT.0.0_SP)) THEN
            T1_SF(I,K) = T1(I,K)-XXXX
         ELSE IF((ABS(XXXX).LT.ABS(TMIN-T1(I,K))).AND.(XXXX.GT.0.0_SP)) THEN 
            T1_SF(I,K) = T1(I,K)-XXXX
         ELSE IF(XXXX.EQ.0) THEN
            T1_SF(I,K) = T1(I,K)
         ELSE
            BETA(I,K)=0.5*(1.-SIGN(1._SP,XXXX)) * (TMAX-T1(I,K))/(ABS(XXXX)+1.E-10) &                               
                 +0.5*(1.-SIGN(1._SP,-XXXX)) * (T1(I,K)-TMIN)/(ABS(XXXX)+1.E-10)                                     
            T1_SF(I,K)=T1(I,K)-MIN(1._SP,BETA(I,K))*XXXX
         END IF

!           ZZZFLUX(I,K) = TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K)/ART1(I)*(DTI/DT(I))/DZ(I,K) 
!           XXXX = ZZZFLUX(I,K)*DT(I)/DTFA(I)+T1(I,K)-T1(I,K)*DT(I)/DTFA(I) 

!           BETA(I,K)=0.5*(1.-SIGN(1._SP,XXXX)) * (TMAX-T1(I,K))/(ABS(XXXX)+1.E-10) &
!                +0.5*(1.-SIGN(1._SP,-XXXX)) * (T1(I,K)-TMIN)/(ABS(XXXX)+1.E-10)

!           T1_SF(I,K)=T1(I,K)-MIN(1.,BETA(I,K))*XXXX

     END DO
  END DO  !! SIGMA LOOP

  !----------------------------------------------------------------------------------------
  NTERA = 4
  DO ITERA=1,NTERA   !! Smolaricizw Loop 
     IF(ITERA == 1)THEN
        WWWSF  = WTS
        T1_S   = T1_SF
        DTWWWS = DT
     ELSE
        WWWSF  = WWWS
        T1_S   = T1_SF
        DTWWWS = DTFA
     END IF
     DO K=2,KBM1
        DO I=1,M
           TEMP=ABS(WWWSF(I,K))-DTI*(WWWSF(I,K))*(WWWSF(I,K))/DZ(I,K)/DTWWWS(I)
           WWWS(I,K)=TEMP*(T1_S(I,K-1)-T1_S(I,K))/(ABS(T1_S(I,K-1))+ABS(T1_S(I,K))+1.E-14)

           IF(TEMP < 0.0_SP .OR. T1_S(I,K) == 0.0_SP)THEN
              WWWS(I,K)=0. 
           END IF
        END DO
     END DO
     DO I=1,M
        WWWS(I,1)=0.
     END DO

     DO I=1,M
        TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),1))
        TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),1))
        TMAX = MAX(TMAX,T1(I,2),T1(I,1),T1_FRESH(I,1))
        TMIN = MIN(TMIN,T1(I,2),T1(I,1),T1_FRESH(I,1))

        TEMP=0.5*((WWWS(I,2)+ABS(WWWS(I,2)))*T1_S(I,2))*(DTI/DTFA(I))/DZ(I,1)
        BETAIN(I,1)=(TMAX-T1_S(I,1))/(TEMP+1.E-10)

        TEMP=0.5*((WWWS(I,1)+ABS(WWWS(I,1)))*T1_S(I,1)-        &
             (WWWS(I,2)-ABS(WWWS(I,2)))*T1_S(I,1))*(DTI/DTFA(I))/DZ(I,1)
        BETAOUT(I,1)=(T1_S(I,1)-TMIN)/(TEMP+1.E-10)

        WWWSF(I,1)=0.5*MIN(1.,BETAOUT(I,1))*(WWWS(I,1)+ABS(WWWS(I,1))) + &
             0.5*MIN(1.,BETAIN(I,1))*(WWWS(I,1)-ABS(WWWS(I,1)))
     END DO

     DO K=2,KBM1-1
        DO I=1,M
           TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
           TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
           TMAX = MAX(TMAX,T1(I,K+1),T1(I,K-1),T1(I,K),T1_FRESH(I,K))
           TMIN = MIN(TMIN,T1(I,K+1),T1(I,K-1),T1(I,K),T1_FRESH(I,K))

           TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*T1_S(I,K+1)-  &
                (WWWS(I,K)-ABS(WWWS(I,K)))*T1_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
           BETAIN(I,K)=(TMAX-T1_S(I,K))/(TEMP+1.E-10)

           TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)-        &
                (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*T1_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
           BETAOUT(I,K)=(T1_S(I,K)-TMIN)/(TEMP+1.E-10)

           WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
                0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
        END DO
     END DO


     K=KBM1
     DO I=1,M
        TMAX = MAXVAL(T1(NBSN(I,1:NTSN(I)),K))
        TMIN = MINVAL(T1(NBSN(I,1:NTSN(I)),K))
        TMAX = MAX(TMAX,T1(I,K-1),T1(I,K),T1_FRESH(I,K))
        TMIN = MIN(TMIN,T1(I,K-1),T1(I,K),T1_FRESH(I,K))

        TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*T1_S(I,K+1)-  &
             (WWWS(I,K)-ABS(WWWS(I,K)))*T1_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
        BETAIN(I,K)=(TMAX-T1_S(I,K))/(TEMP+1.E-10)

        TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)-        &
             (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*T1_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
        BETAOUT(I,K)=(T1_S(I,K)-TMIN)/(TEMP+1.E-10)

        WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + &
             0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
     END DO



     WWWS=WWWSF 

     DO K=1,KBM1
        DO I=1,M
              IF(K == 1) THEN
                 TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*T1_S(I,K)   &
                      -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*T1_S(I,K+1) &
                      +(WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)
              ELSE IF(K == KBM1) THEN
                 TEMP = +(WWWS(I,K)-ABS(WWWS(I,K)))*T1_S(I,K-1)     &
                      +(WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)
              ELSE
                 TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*T1_S(I,K)   &
                      -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*T1_S(I,K+1) &
                      +(WWWS(I,K)-ABS(WWWS(I,K)))*T1_S(I,K-1)     &
                      +(WWWS(I,K)+ABS(WWWS(I,K)))*T1_S(I,K)
              END IF
              TEMP = 0.5_SP*TEMP
              T1_SF(I,K)=(T1_S(I,K)-TEMP*(DTI/DTFA(I))/DZ(I,K))
        END DO
     END DO  !! SIGMA LOOP
  END DO  !! Smolarvizw Loop
  !--------------------------------------------------------------------------
  ! End of smolarkiewicz upwind loop
  !--------------------------------------------------------------------------
      T1 = T1-OFFS
      T1_SF = T1_SF-OFFS
      TMEAN1 = TMEAN1-OFFS
     




  ! APPLY GROUND WATER TEMPERATURE FORCING
  IF(GROUNDWATER_ON .and. GROUNDWATER_TEMP_ON)THEN
     DO I=1,M
        XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*BFWTMP(I)
     END DO
  ELSEIF(GROUNDWATER_ON) THEN
     DO I=1,M
 !J. Ge for tracer advection
        IF(BACKWARD_ADVECTION.eqv..FALSE.)THEN
          XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*T1(I,KBM1)
        ELSE
          IF(BACKWARD_STEP==1)THEN
            XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*(T0(I,KBM1)+T1(I,KBM1))*0.5
          ELSEIF(BACKWARD_STEP==2)THEN
            XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS(I)*(T2(I,KBM1)+T0(I,KBM1)+T1(I,KBM1))/3.0_SP
          ENDIF
        ENDIF
!J. Ge for tracer advection
     END DO
  END IF


  !
  !--Update Temperature----------------------------------------------------------!
  !


  DO I=1,M
     DO K=1,KBM1
        XFLUX(I,K) = XFLUX(I,K) - RF(I,K)*ART1(I)    !/DZ(I,K)
        TF1(I,K)=T1_SF(I,K) 
        !	TF1(I,K)=T1_SF(I,K) - RF(I,K)*DTI/DTFA(I)/DZ(I,K)  ! considering about heat flux source term, not tested yet
     END DO
  END DO



  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*)"End: adv_t"

END SUBROUTINE ADV_T
!==============================================================================|
