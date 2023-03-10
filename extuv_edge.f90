










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
!     ACCUMLATE FLUXES FOR EXTERNAL MODE                                       |
!==============================================================================|

   SUBROUTINE EXTUV_EDGE(K)       

!==============================================================================|
   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_WD

   USE MOD_NORTHPOLE




   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:NT) :: RESX,RESY,TMP
   REAL(SP) :: UAFT,VAFT
   INTEGER  :: I

!==============================================================================|

   if(dbg_set(dbg_sbr)) write(ipt,*) "Start: extuv_edge.F"

!
!--ACCUMULATE RESIDUALS FOR EXTERNAL MODE EQUATIONS----------------------------|
!
   UAFT = UAF(0)
   VAFT = VAF(0)

   ! THIS APPEARS TO BE TO PREVENT DIVISION BY ZERO, BUT IT IS A
   ! STRANGE WAY TO DO IT!
   H1(0)= H1(1)

!!#  if defined (WET_DRY)
!!   IF(K == 3)THEN

!!#    if !defined (NH)
!!     RESX = ADX2D + ADVUA + DRX2D + PSTX - COR*VA*D1*ART  &
!!            -(WUSURF2 + WUBOT)*ART
!!     RESY = ADY2D + ADVVA + DRY2D + PSTY + COR*UA*D1*ART  &
!!            -(WVSURF2 + WVBOT)*ART
!!#    else
!!     RESX = ADX2D + ADVUA + DRX2D + PSTX - COR*VA*D1*ART  &
!!            -(WUSURF2 + WUBOT)*ART + NHQ2DX
!!     RESY = ADY2D + ADVVA + DRY2D + PSTY + COR*UA*D1*ART  &
!!            -(WVSURF2 + WVBOT)*ART + NHQ2DY
!!#    endif

!!#  if defined (SPHERICAL)
!!     RESX = RESX -UA*VA/REARTH*TAN(DEG2RAD*YC)*D1*ART
!!     RESY = RESY +UA*UA/REARTH*TAN(DEG2RAD*YC)*D1*ART
!!#  endif

!!!
!!!--UPDATE----------------------------------------------------------------------|
!!!

!!     UAF = (UARK*(H1+ELRK1)-ALPHA_RK(K)*DTE*RESX/ART)/(H1+ELF1)
!!     VAF = (VARK*(H1+ELRK1)-ALPHA_RK(K)*DTE*RESY/ART)/(H1+ELF1)
!!     UAS = UAF
!!     VAS = VAF
!!   END IF
!!#  endif

   DO I=1,NT

       RESX(I) = ADX2D(I)+ADVUA(I)+DRX2D(I)+PSTX(I)-COR(I)*VA(I)*D1(I)*ART(I)  &
                 -(WUSURF2(I)+WUBOT(I))*ART(I)
       RESY(I) = ADY2D(I)+ADVVA(I)+DRY2D(I)+PSTY(I)+COR(I)*UA(I)*D1(I)*ART(I)  &
                 -(WVSURF2(I)+WVBOT(I))*ART(I)



!
!--UPDATE----------------------------------------------------------------------|
!

       UAF(I)  =  (UARK(I)*(H1(I)+ELRK1(I))-ALPHA_RK(K)*DTE*RESX(I)/ART(I))/(H1(I)+ELF1(I))
       VAF(I)  =  (VARK(I)*(H1(I)+ELRK1(I))-ALPHA_RK(K)*DTE*RESY(I)/ART(I))/(H1(I)+ELF1(I))
   END DO

   
   VAF(0) = VAFT
   UAF(0) = UAFT

!
!--ADJUST EXTERNAL VELOCITY IN SPONGE REGION-----------------------------------|
!
!old:   UAF = UAF-CC_SPONGE*UAF
!old:   VAF = VAF-CC_SPONGE*VAF
! ---- new: Karsten Lettmann: 2012.06.25 -------
   UAF = UAF/(1.0_SP+CC_SPONGE*UAF**2.0_SP)
   VAF = VAF/(1.0_SP+CC_SPONGE*VAF**2.0_SP)
! ------- end new -------------------------------


!
!--STORE VARIABLES FOR MOMENTUM BALANCE CHECK----------------------------------|
!

   if(dbg_set(dbg_sbr)) write(ipt,*) "End: extuv_edge.F"

   END SUBROUTINE EXTUV_EDGE
!==============================================================================|
