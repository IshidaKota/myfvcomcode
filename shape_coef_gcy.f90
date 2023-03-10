










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

   SUBROUTINE SHAPE_COEF_GCY

!----------------------------------------------------------------------!
!  This subroutine is used to calculate the coefficient for a linear   !
!  function on the x-y plane, i.e.:                                    !
!                                                                      !
!                     r(x,y;phai)=phai_c+cofa1*x+cofa2*y               !
!                                                                      !
!  This subroutine is used for ghost cell boundary condition cases     !
!----------------------------------------------------------------------!

   USE ALL_VARS
   USE MOD_UTILS
   IMPLICIT NONE
   REAL(DP) X1,X2,X3,Y1,Y2,Y3,DELT,AI1,AI2,AI3,BI1,BI2,BI3,CI1,CI2,CI3
   REAL(DP) DELTX,DELTY,TEMP1,ANG1,ANG2,B1,B2,ANGLE
   INTEGER  I,II,J,JJ,J1,J2
   REAL(SP) AA1,AA2,BB1,BB2,CC1,CC2,XTMP1,YTMP1,XTMP2,YTMP2
!
!---------------interior and boundary cells------------------------------------!
!

   IF(DBG_SET(DBG_LOG)) THEN
      WRITE(IPT,*) "!"
      WRITE(IPT,*) "!           SETTING UP LINEAR INTEROPOLATION COEFFICIENTS"
   END IF

   DO I = 1, N
     IF(ISBCE(I) == 0) THEN
       Y1 = YC(NBE(I,1))-YC(I)
       Y2 = YC(NBE(I,2))-YC(I)
       Y3 = YC(NBE(I,3))-YC(I)
       X1 = XC(NBE(I,1))-XC(I)
       X2 = XC(NBE(I,2))-XC(I)
       X3 = XC(NBE(I,3))-XC(I)
     ELSE IF(ISBCE(I) == 1) THEN
       DO J = 1, 3
         IF(NBE(I,J) == 0) JJ = J
       END DO
       J1 = JJ+1-INT((JJ+1)/4)*3
       J2 = JJ+2-INT((JJ+2)/4)*3



       IF(JJ == 1)THEN
         X1 = XTMP1-XC(I)
         Y1 = YTMP1-YC(I)
         X2 = XC(NBE(I,J1))-XC(I)
         Y2 = YC(NBE(I,J1))-YC(I)
         X3 = XC(NBE(I,J2))-XC(I)
         Y3 = YC(NBE(I,J2))-YC(I)
       ELSE IF(JJ == 2)THEN
         X1 = XC(NBE(I,J2))-XC(I)
         Y1 = YC(NBE(I,J2))-YC(I)
         X2 = XTMP1-XC(I)
         Y2 = YTMP1-YC(I)
         X3 = XC(NBE(I,J1))-XC(I)
         Y3 = YC(NBE(I,J1))-YC(I)
       ELSE
         X1 = XC(NBE(I,J1))-XC(I)
         y1 = YC(NBE(I,J1))-YC(I)
         X2 = XC(NBE(I,J2))-XC(I)
         Y2 = YC(NBE(I,J2))-YC(I)
         X3 = XTMP1-XC(I)
         Y3 = YTMP1-YC(I)
       END IF
!     else if(isbce(i) == 2) then
!       do j=1,4
!         a1u(i,j)=0.0_SP
!         a2u(i,j)=0.0_SP
!       end do
     ELSE IF(ISBCE(I) == 2) THEN
       DO J = 1, 3
         IF(NBE(I,J) == 0) JJ = J
       END DO
       J1 = JJ+1-INT((JJ+1)/4)*3
       J2 = JJ+2-INT((JJ+2)/4)*3

       DELTX = VX(NV(I,J1))-VX(NV(I,J2))
       DELTY = VY(NV(I,J1))-VY(NV(I,J2))

       ALPHA(I) = ATAN2(DELTY,DELTX)
       ALPHA(I) = ALPHA(I)-3.1415926_SP/2.0_SP

       AA1 = -DELTY
       BB1 = DELTX
       CC1 = -AA1*VX(NV(I,J1))-BB1*VY(NV(I,J1))

       AA2 = BB1
       BB2 = -AA1
       CC2 = -AA2*XC(I)-BB2*YC(I)

       XTMP1 = -(CC1*BB2-CC2*BB1)/(AA1*BB2-AA2*BB1)
       YTMP1 = -(CC1*AA2-CC2*AA1)/(BB1*AA2-BB2*AA1)

       IF(JJ == 1)THEN
         X1 = (XTMP1-XC(I))*2.0_SP
         Y1 = (YTMP1-YC(I))*2.0_SP
         X2 = XC(NBE(I,J1))-XC(I)
         Y2 = YC(NBE(I,J1))-YC(I)
         X3 = XC(NBE(I,J2))-XC(I)
         Y3 = YC(NBE(I,J2))-YC(I)
       ELSE IF(JJ == 2)THEN
         X1 = XC(NBE(I,J2))-XC(I)
         Y1 = YC(NBE(I,J2))-YC(I)
         X2 = (XTMP1-XC(I))*2.0_SP
         Y2 = (YTMP1-YC(I))*2.0_SP
         X3 = XC(NBE(I,J1))-XC(I)
         Y3 = YC(NBE(I,J1))-YC(I)
       ELSE
         X1 = XC(NBE(I,J1))-XC(I)
         y1 = YC(NBE(I,J1))-YC(I)
         X2 = XC(NBE(I,J2))-XC(I)
         Y2 = YC(NBE(I,J2))-YC(I)
         X3 = (XTMP1-XC(I))*2.0_SP
         Y3 = (YTMP1-YC(I))*2.0_SP
       END IF
     ELSE IF(ISBCE(I) == 3)THEN
       DO J = 1, 3
         IF(NBE(I,J) /= 0) JJ = J
       END DO
       J1 = JJ+1-INT((JJ+1)/4)*3
       J2 = JJ+2-INT((JJ+2)/4)*3



       IF(JJ == 1)THEN
         X1 = XC(NBE(I,JJ))-XC(I)
         Y1 = YC(NBE(I,JJ))-YC(I)
         X2 = XTMP1-XC(I)
         Y2 = YTMP1-YC(I)
         X3 = XTMP2-XC(I)
         Y3 = YTMP2-YC(I)
       ELSE IF(JJ == 2)THEN
         X1 = XTMP2-XC(I)
         Y1 = YTMP2-YC(I)
         X2 = XC(NBE(I,JJ))-XC(I)
         Y2 = YC(NBE(I,JJ))-YC(I)
         X3 = XTMP1-XC(I)
         Y3 = YTMP1-YC(I)
       ELSE
         X1 = XTMP1-XC(I)
         Y1 = YTMP1-YC(I)
         X2 = XTMP2-XC(I)
         Y2 = YTMP2-YC(I)
         X3 = XC(NBE(I,JJ))-XC(I)
         Y3 = YC(NBE(I,JJ))-YC(I)
       END IF
     END IF

     X1 = X1/1000.0_SP
     X2 = X2/1000.0_SP
     X3 = X3/1000.0_SP
     Y1 = Y1/1000.0_SP
     Y2 = Y2/1000.0_SP
     Y3 = Y3/1000.0_SP

     DELT = (X1*Y2-X2*Y1)**2+(X1*Y3-X3*Y1)**2+(X2*Y3-X3*Y2)**2
     DELT = DELT*1000._SP

     A1U(I,1) = (Y1+Y2+Y3)*(X1*Y1+X2*Y2+X3*Y3)- &
                (X1+X2+X3)*(Y1**2+Y2**2+Y3**2)
     A1U(I,1) = A1U(I,1)/DELT
     A1U(I,2) = (Y1**2+Y2**2+Y3**2)*X1-(X1*Y1+X2*Y2+X3*Y3)*Y1
     A1U(I,2) = A1U(I,2)/DELT
     A1U(I,3) = (Y1**2+Y2**2+Y3**2)*X2-(X1*Y1+X2*Y2+X3*Y3)*Y2
     A1U(I,3) = A1U(I,3)/DELT
     A1U(I,4) = (Y1**2+Y2**2+Y3**2)*X3-(X1*Y1+X2*Y2+X3*Y3)*Y3
     A1U(I,4) = A1U(I,4)/DELT

     A2U(I,1) = (X1+X2+X3)*(X1*X1+X2*X2+X3*X3)- &
                (Y1+Y2+Y3)*(X1**2+X2**2+X3**2)
     A2U(I,1) = A2U(I,1)/DELT
     A2U(I,2) = (X1**2+X2**2+X3**2)*Y1-(X1*Y1+X2*Y2+X3*Y3)*X1
     A2U(I,2) = A2U(I,2)/DELT
     A2U(I,3) = (X1**2+X2**2+X3**2)*Y2-(X1*Y1+X2*Y2+X3*Y3)*X2
     A2U(I,3) = A2U(I,3)/DELT
     A2U(I,4) = (X1**2+X2**2+X3**2)*Y3-(X1*Y1+X2*Y2+X3*Y3)*X3
     A2U(I,4) = A2U(I,4)/DELT
!     end if

     X1 = VX(NV(I,1))-XC(I)
     X2 = VX(NV(I,2))-XC(I)
     X3 = VX(NV(I,3))-XC(I)
     Y1 = VY(NV(I,1))-YC(I)
     Y2 = VY(NV(I,2))-YC(I)
     Y3 = VY(NV(I,3))-YC(I)


     AI1 = Y2-Y3
     AI2 = Y3-Y1
     AI3 = Y1-Y2
     BI1 = X3-X2
     BI2 = X1-X3
     BI3 = X2-X1
     CI1 = X2*Y3-X3*Y2
     CI2 = X3*Y1-X1*Y3
     CI3 = X1*Y2-X2*Y1

     AW0(I,1) = -CI1/2./ART(I)
     AW0(I,2) = -CI2/2./ART(I)
     AW0(I,3) = -CI3/2./ART(I)
     AWX(I,1) = -AI1/2./ART(I)
     AWX(I,2) = -AI2/2./ART(I)
     AWX(I,3) = -AI3/2./ART(I)
     AWY(I,1) = -BI1/2./ART(I)
     AWY(I,2) = -BI2/2./ART(I)
     AWY(I,3) = -BI3/2./ART(I)
   END DO

!<This part may be not useful again but keep it for verification
   ang1=359.9_SP/180.0_SP*3.1415926_SP
   ang2=-0.1_SP/180.0_SP*3.1415926_SP

   do i=1,m
     if((isonb(i).eq.1).and.(ntve(i).gt.2)) then
       angle=alpha(nbve(i,ntve(i)))-alpha(nbve(i,1))
       if(angle.gt.ang1) then
         angle=100000.0_SP
       else if(angle.gt.3.1415926_SP) then
         angle=angle-2.0_SP*3.1415926_SP
       else if(angle.lt.-3.1415926_SP) then
         angle=angle+2.0_SP*3.1415926_SP
       else if(angle.lt.ang2) then
         angle=100000.0_SP
       end if
       do j=2,ntve(i)-1
         ii=nbve(i,j)
         if(isbce(ii).ne.1) then
           alpha(ii)=alpha(nbve(i,1))+ &
                     angle/float(ntve(i)-1)*float(j-1)
         end if
       end do
     end if
   end do
!end>

   IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*) "!"
       WRITE(IPT,*) "!  INTERP COEFFICIENTS   :    COMPLETE"
    END IF

   RETURN
   END SUBROUTINE SHAPE_COEF_GCY

