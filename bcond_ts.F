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
!   Set Boundary Conditions on Temperature and Salinity                        |
!    NCON2 = 1:  SET CONDITIONS SPECIFIC TO TEMPERATURE                        |
!    NCON2 = 2:  SET CONDITIONS SPECIFIC TO SALINITY                           |
!==============================================================================|

SUBROUTINE BCOND_TS(NCON2)     

  !------------------------------------------------------------------------------|
  USE ALL_VARS
  USE BCS
  USE MOD_UTILS
  USE MOD_OBCS
  USE MOD_FORCE

  IMPLICIT NONE
  REAL(SP) :: S2D,S2D_NEXT,S2D_OBC,T2D,T2D_NEXT,T2D_OBC,XFLUX2D,TMP,RAMP_TS


  INTEGER  :: I,J,K,J1,J11,J22,NCON2
  REAL(SP), ALLOCATABLE :: TTMP(:,:),STMP(:,:)
  !<----------ishid debug
   REAL(SP) :: dx,dy,dj1x,dj1y,dp1,dp2,dp12,flag_i,flag_o,INTMP,OUTTMP

   INTEGER :: EJ,J_NEXT
  !>----------ishid debug
  REAL(SP) ::TMAX,TMIN,SMAX,SMIN

  !------------------------------------------------------------------------------|

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "Start: bcond_ts"

  !
  !--SET CONDITIONS FOR FRESH WATER INFLOW---------------------------------------|
  !
  IF(RIVER_TS_SETTING == 'specified') THEN
     IF(NUMQBC > 0) THEN
        IF(RIVER_INFLOW_LOCATION == 'node') THEN
           DO I=1,NUMQBC
              J11=INODEQ(I)
              DO K=1,KBM1
                 TF1(J11,K)=TDIS(I)
                 SF1(J11,K)=SDIS(I)
              END DO
           END DO
        ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
           DO I=1,NUMQBC
              J11=N_ICELLQ(I,1)
              J22=N_ICELLQ(I,2)
              DO K=1,KBM1
                 TF1(J11,K)=TDIS(I)
                 SF1(J11,K)=SDIS(I)
                 TF1(J22,K)=TDIS(I)
                 SF1(J22,K)=SDIS(I)
              END DO
           END DO
        END IF
     END IF
  END IF


  IF(.NOT. OBC_ON) RETURN


  !
  !  SET TEMPERATURE CONDITIONS ON OUTER BOUNDARY
  !
  IF(NCON2 == 1) THEN

     IF(OBC_TEMP_NUDGING) CALL UPDATE_OBC_TEMP(IntTime,TEMP_OBC)

     ALLOCATE(TTMP(IOBCN,KBM1));  TTMP = 0.0_SP
     DO I=1,IOBCN
        J=I_OBC_N(I)
        J1=NEXT_OBC(I)
        T2D=0.0_SP
        T2D_NEXT=0.0_SP
        XFLUX2D=0.0_SP
        DO K=1,KBM1
           T2D=T2D+T1(J,K)*DZ(J,K)
           T2D_NEXT=T2D_NEXT+TF1(J1,K)*DZ(J1,K)
           XFLUX2D=XFLUX2D+XFLUX_OBC(I,K)           !*DZ(J,K)
        END DO

     
!<--------ishid debug 20230107
        !REAL(SP) :: dx,dy
        !INTEGER :: EJ
        !REAL(SP),ALLOCATABLE :: res(:),

        TMP=XFLUX2D+T2D*UARD_OBCN(I)
        T2D_OBC=(T2D*DT(J)-TMP*DTI/ART1(J))/D(J)

        CALL BCOND_T_PERTURBATION(T2D_NEXT,T2D,TTMP,I,J,J1)


        IF (I .le. IOBCN-1)  THEN
         J_NEXT = I_OBC_N(I+1)
         dx = vx(J_NEXT) - vx(J) !delta x in meters for obc edge
         dy = vy(J_NEXT) - vy(J) !delta y in meters for obc edge  
         dj1x = vx(J1) -vx(J)
         dj1y = vy(J1) -vy(J)
         IF (dx >0) THEN
            EJ = LISBCE_2(I)
         ELSE
            EJ = LISBCE_2(IOBCN-I)
         ENDIF
        !ESLE IF (I == IOBCN) THEN
         !???????????????????????????????????????????????????
        ENDIF
        !write(IPT,*) "I=", I ,"J=",J,"EJ=",EJ



        


        DO K = 1,KBM1
         dp1 = - dy * U(EJ,K)   +   dx * V(EJ,K) !????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????90???????????????????????????????????????????????????cell????????????????????????????????????????????????
         dp2 = - dy * dj1x     +   dx * dj1y !???????????????????????????????????????????????? NEXT_OBC() ????????????????????????????????????????????????

         dp12 = dp1 * dp2
         flag_i = MAX(SIGN(1.0_SP,  dp12), 0.0_SP) ! flag_i = 1??????????????????=0???????????????
         flag_o = MAX(SIGN(1.0_SP, -dp12), 0.0_SP) ! flag_o = 1??????????????????=0???????????????
         !INTMP !?????????????????????    !OUTTMP !?????????????????????
           ! ??????
           IF(OBC_TEMP_NUDGING) THEN
              INTMP = T1(J,K) - OBC_TEMP_NUDGING_TIMESCALE*RAMP*(T1(J,K)&
              &-TEMP_OBC(I,K))
           ELSE
              INTMP = T1(J,K)
           ENDIF

            !??????
           OUTTMP=T2D_OBC+TTMP(I,K)

           TMAX = MAXVAL(T1(NBSN(J,1:NTSN(J)),K))
           TMIN = MINVAL(T1(NBSN(J,1:NTSN(J)),K))

            IF(K == 1)THEN
               TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/  &
               (DZ(J,K)+DZ(J,K+1)))
               TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/  &
               (DZ(J,K)+DZ(J,K+1)))
            ELSE IF(K == KBM1)THEN
               TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
                  (DZ(J,K)+DZ(J,K-1)))
               TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  & 
                  (DZ(J,K)+DZ(J,K-1)))
            ELSE
               TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
                  (DZ(J,K)+DZ(J,K-1)),                             &
                  (T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/           &
                  (DZ(J,K)+DZ(J,K+1)))
               TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
                  (DZ(J,K)+DZ(J,K-1)),                             &
                  (T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/           &
                  (DZ(J,K)+DZ(J,K+1)))
            END IF

            IF(TMIN-OUTTMP > 0.0_SP)OUTTMP = TMIN
            IF(OUTTMP-TMAX > 0.0_SP)OUTTMP = TMAX

            TF1(J,K) = flag_i * INTMP + flag_o * (1.0_SP-flag_i) * OUTTMP
            !IF (K==4 .and. I==7) write(IPT,*) "TEMP:INTMP=",INTMP," OUTTMP=, " ,OUTTMP ,"flag_i=",flag_i,"flag_o=",flag_o ,"TF1=",TF1(J,K)
         ENDDO
         
        !...
        !
!            ! IF THE FLOW IS OUT OF THE DOMAIN
!        IF(UARD_OBCN(I) > 0.0_SP) THEN

!           TMP=XFLUX2D+T2D*UARD_OBCN(I)
!           T2D_OBC=(T2D*DT(J)-TMP*DTI/ART1(J))/D(J)

!           CALL BCOND_T_PERTURBATION(T2D_NEXT,T2D,TTMP,I,J,J1)

!           DO K=1,KBM1
!              TF1(J,K)=T2D_OBC+TTMP(I,K)
              !           TF1(J,K)=T2D_OBC+(TF1(J1,K)-T2D_NEXT)
!           END DO

!           DO K=1,KBM1
!              TMAX = MAXVAL(T1(NBSN(J,1:NTSN(J)),K))
!              TMIN = MINVAL(T1(NBSN(J,1:NTSN(J)),K))

!              IF(K == 1)THEN
!                 TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K+1)))
!                 TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K+1)))
!              ELSE IF(K == KBM1)THEN
!                 TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K-1)))
!                 TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  & 
!                      (DZ(J,K)+DZ(J,K-1)))
!              ELSE
!                 TMAX = MAX(TMAX,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K-1)),                             &
!                      (T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/           &
!                      (DZ(J,K)+DZ(J,K+1)))
!                 TMIN = MIN(TMIN,(T1(J,K)*DZ(J,K-1)+T1(J,K-1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K-1)),                             &
!                      (T1(J,K)*DZ(J,K+1)+T1(J,K+1)*DZ(J,K))/           &
!                      (DZ(J,K)+DZ(J,K+1)))
!              END IF

!              IF(TMIN-TF1(J,K) > 0.0_SP)TF1(J,K) = TMIN
!              IF(TF1(J,K)-TMAX > 0.0_SP)TF1(J,K) = TMAX

!           END DO

!        ELSE! IF THE FLOW IS INTO THE DOMAIN

!           IF(OBC_TEMP_NUDGING) THEN
!              DO K=1,KBM1
!                 TF1(J,K) = T1(J,K) - OBC_TEMP_NUDGING_TIMESCALE*RAMP*(T1(J,K)&
!                      &-TEMP_OBC(I,K))
!              END DO
!           ELSE
!              DO K=1,KBM1
!              TF1(J,K) = T1(J,K)
!            END DO
!         END IF

!       END IF

     END DO
     DEALLOCATE(TTMP)

!>----------------ishid debug 20230107

     !
     !  SET SALINITY CONDITIONS ON OUTER BOUNDARY
     !
  ELSE IF(NCON2 == 2) THEN

     IF (OBC_SALT_NUDGING) CALL UPDATE_OBC_SALT(IntTime,SALT_OBC)

     ALLOCATE(STMP(IOBCN,KBM1));  STMP = 0.0_SP

     DO I=1,IOBCN
        J=I_OBC_N(I)
        J1=NEXT_OBC(I)
        S2D=0.0_SP
        S2D_NEXT=0.0_SP
        XFLUX2D=0.0_SP
        DO K=1,KBM1
           S2D=S2D+S1(J,K)*DZ(J,K)
           S2D_NEXT=S2D_NEXT+SF1(J1,K)*DZ(J1,K)
           XFLUX2D=XFLUX2D+XFLUX_OBC(I,K)             !*DZ(J,K)
        END DO
!<-------------------------ishid debug
        TMP=XFLUX2D+S2D*UARD_OBCN(I)
        S2D_OBC=(S2D*DT(J)-TMP*DTI/ART1(J))/D(J)

        CALL BCOND_S_PERTURBATION(S2D_NEXT,S2D,STMP,I,J,J1)
        IF (I .le. IOBCN-1)  THEN
         J_NEXT = I_OBC_N(I+1)
         dx = vx(J_NEXT) - vx(J) !delta x in meters for obc edge
         dy = vy(J_NEXT) - vy(J) !delta y in meters for obc edge  
         dj1x = vx(J1) -vx(J)
         dj1y = vy(J1) -vy(J)
         IF (dx >0) THEN
            EJ = LISBCE_2(I)
         ELSE
            EJ = LISBCE_2(IOBCN-I)
         ENDIF
        !ESLE IF (I == IOBCN) THEN
         !???????????????????????????????????????????????????
        ENDIF
        !write(IPT,*) "I=", I ,"J=",J,"EJ=",EJ

      DO K = 1,KBM1
         dp1 = - dy * U(EJ,K)   +   dx * V(EJ,K) !????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????90???????????????????????????????????????????????????cell????????????????????????????????????????????????
         dp2 = - dy * dj1x     +   dx * dj1y !???????????????????????????????????????????????? NEXT_OBC() ????????????????????????????????????????????????

         dp12 = dp1 * dp2
         flag_i = MAX(SIGN(1.0_SP,  dp12), 0.0_SP) ! flag_i = 1??????????????????=0???????????????
         flag_o = MAX(SIGN(1.0_SP, -dp12), 0.0_SP) ! flag_o = 1??????????????????=0???????????????

           ! ??????
           IF(OBC_SALT_NUDGING) THEN
              INTMP = S1(J,K) - OBC_SALT_NUDGING_TIMESCALE*RAMP*(S1(J,K)&
              &-SALT_OBC(I,K))
           ELSE
              INTMP = S1(J,K)
           ENDIF
   
            !??????
           OUTTMP=S2D_OBC+STMP(I,K)

           SMAX = MAXVAL(S1(NBSN(J,1:NTSN(J)),K))
           SMIN = MINVAL(S1(NBSN(J,1:NTSN(J)),K))

            IF(K == 1)THEN
               SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/  &
               (DZ(J,K)+DZ(J,K+1)))
               SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/  &
               (DZ(J,K)+DZ(J,K+1)))
            ELSE IF(K == KBM1)THEN
               SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
                  (DZ(J,K)+DZ(J,K-1)))
               SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  & 
                  (DZ(J,K)+DZ(J,K-1)))
            ELSE
               SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
                  (DZ(J,K)+DZ(J,K-1)),                             &
                  (S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/           &
                  (DZ(J,K)+DZ(J,K+1)))
               SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
                  (DZ(J,K)+DZ(J,K-1)),                             &
                  (S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/           &
                  (DZ(J,K)+DZ(J,K+1)))
            END IF

            IF(SMIN-OUTTMP > 0.0_SP)OUTTMP = SMIN
            IF(OUTTMP-SMAX > 0.0_SP)OUTTMP = SMAX
         SF1(J,K) = flag_i * INTMP + flag_o * (1.0_SP-flag_i) * OUTTMP
         !IF (K==4 .and. I==7) write(IPT,*) "sal INTMP=",INTMP," OUTTMP=, " ,OUTTMP ,"flag_i=",flag_i,"flag_o=",flag_o,"SF1=",SF1(J,K)
      ENDDO
!        ! IF THE FLOW IS OUT OF THE DOMAIN
!        IF(UARD_OBCN(I) > 0.0_SP) THEN
!           TMP=XFLUX2D+S2D*UARD_OBCN(I)
!           S2D_OBC=(S2D*DT(J)-TMP*DTI/ART1(J))/D(J)
!
!           CALL BCOND_S_PERTURBATION(S2D_NEXT,S2D,STMP,I,J,J1)
!
!           DO K=1,KBM1
!              SF1(J,K)=S2D_OBC+STMP(I,K)  
!              !           SF1(J,K)=S2D_OBC+(SF1(J1,K)-S2D_NEXT)  
!           END DO

!           DO K=1,KBM1
!              SMAX = MAXVAL(S1(NBSN(J,1:NTSN(J)),K))
!              SMIN = MINVAL(S1(NBSN(J,1:NTSN(J)),K))

!              IF(K == 1)THEN
!                 SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K+1)))
!                 SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K+1)))
!              ELSE IF(K == KBM1)THEN
!                 SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K-1)))
!                 SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K-1)))
!              ELSE
!                 SMAX = MAX(SMAX,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K-1)),                             &
!                      (S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/           &
!                      (DZ(J,K)+DZ(J,K+1)))
!                 SMIN = MIN(SMIN,(S1(J,K)*DZ(J,K-1)+S1(J,K-1)*DZ(J,K))/  &
!                      (DZ(J,K)+DZ(J,K-1)),                             &
!                      (S1(J,K)*DZ(J,K+1)+S1(J,K+1)*DZ(J,K))/           &
!                      (DZ(J,K)+DZ(J,K+1)))
!              END IF
!
!              IF(SMIN-SF1(J,K) > 0.0_SP) SF1(J,K) = SMIN
!              IF(SF1(J,K)-SMAX > 0.0_SP) SF1(J,K) = SMAX
!
!           END DO
!        ELSE ! IF THE FLOW IS INTO THE DOMAIN
!
!           IF (OBC_SALT_NUDGING) THEN
!              DO K=1,KBM1
!                 SF1(J,K) = S1(J,K) - OBC_SALT_NUDGING_TIMESCALE*RAMP*(S1(J,K)&
!                      &-SALT_OBC(I,K))
!              END DO
!           ELSE
!              DO K=1,KBM1
!                 SF1(J,K) = S1(J,K)
!              END DO
!           END IF
!
  !      END IF
!>--------------ishid debug 20230107
     END DO
     DEALLOCATE(STMP)
!>----------------ishid debug 20230107
  ELSE
     PRINT*, 'NCON2 NOT IN THE LIST'
     PRINT*, 'MUST BE 1 OR 2'
     CALL PSTOP
  END IF



  !
  !--SET BOUNDARY CONDITIONS-----------------------------------------------------|
  !
  DO K=1,KBM1
     T(0,K)=0.0_SP
     S(0,K)=0.0_SP
  END DO

  IF(DBG_SET(DBG_SBR)) WRITE(IPT,*) "End: bcond_ts"

  RETURN
END SUBROUTINE BCOND_TS
!==============================================================================|
