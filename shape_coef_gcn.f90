










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

   SUBROUTINE SHAPE_COEF_GCN

!----------------------------------------------------------------------!
!  This subrountine is used to calculate the coefficient for a linear  !
!  function on the x-y plane, i.e.:                                    !
!                     r(x,y;phai)=phai_c+cofa1*x+cofa2*y               !
!     innc(i)=0    cells on the boundary                               !
!     innc(i)=1    cells in the interior                               !
!----------------------------------------------------------------------!
     
   USE ALL_VARS
   USE MOD_UTILS
   IMPLICIT NONE
   REAL(DP) X1,X2,X3,Y1,Y2,Y3,DELT,AI1,AI2,AI3,BI1,BI2,BI3,CI1,CI2,CI3
   REAL(DP) DELTX,DELTY,TEMP1,ANG1,ANG2,B1,B2,ANGLE
   INTEGER  I,II,J,JJ,J1,J2
!
!---------------interior cells-----------------------------------------!
!
   IF(DBG_SET(DBG_LOG)) THEN
      WRITE(IPT,*) "!"
      WRITE(IPT,*) "!           SETTING UP LINEAR INTEROPOLATION COEFFICIENTS"
   END IF
 
   DO I=1,N
     IF(ISBCE(I) == 0)THEN
       Y1 = YC(NBE(I,1))-YC(I)
       Y2 = YC(NBE(I,2))-YC(I)
       Y3 = YC(NBE(I,3))-YC(I)
       X1=XC(NBE(I,1))-XC(I)
       X2=XC(NBE(I,2))-XC(I)
       X3=XC(NBE(I,3))-XC(I)

       X1=X1/1000.0_SP
       X2=X2/1000.0_SP
       X3=X3/1000.0_SP
       Y1=Y1/1000.0_SP
       Y2=Y2/1000.0_SP
       Y3=Y3/1000.0_SP

       delt=(x1*y2-x2*y1)**2+(x1*y3-x3*y1)**2+(x2*y3-x3*y2)**2
       delt=delt*1000.

       a1u(i,1)=(y1+y2+y3)*(x1*y1+x2*y2+x3*y3)- &
                (x1+x2+x3)*(y1**2+y2**2+y3**2)
       a1u(i,1)=a1u(i,1)/delt
       a1u(i,2)=(y1**2+y2**2+y3**2)*x1-(x1*y1+x2*y2+x3*y3)*y1
       a1u(i,2)=a1u(i,2)/delt
       a1u(i,3)=(y1**2+y2**2+y3**2)*x2-(x1*y1+x2*y2+x3*y3)*y2
       a1u(i,3)=a1u(i,3)/delt
       a1u(i,4)=(y1**2+y2**2+y3**2)*x3-(x1*y1+x2*y2+x3*y3)*y3
       a1u(i,4)=a1u(i,4)/delt

       a2u(i,1)=(x1+x2+x3)*(x1*y1+x2*y2+x3*y3)- &
                (y1+y2+y3)*(x1**2+x2**2+x3**2)
       a2u(i,1)=a2u(i,1)/delt
       a2u(i,2)=(x1**2+x2**2+x3**2)*y1-(x1*y1+x2*y2+x3*y3)*x1
       a2u(i,2)=a2u(i,2)/delt
       a2u(i,3)=(x1**2+x2**2+x3**2)*y2-(x1*y1+x2*y2+x3*y3)*x2
       a2u(i,3)=a2u(i,3)/delt
       a2u(i,4)=(x1**2+x2**2+x3**2)*y3-(x1*y1+x2*y2+x3*y3)*x3
       a2u(i,4)=a2u(i,4)/delt
     end if

     x1=vx(nv(i,1))-xc(i)
     x2=vx(nv(i,2))-xc(i)
     x3=vx(nv(i,3))-xc(i)
     y1=vy(nv(i,1))-yc(i)
     y2=vy(nv(i,2))-yc(i)
     y3=vy(nv(i,3))-yc(i)


     ai1=y2-y3
     ai2=y3-y1
     ai3=y1-y2
     bi1=x3-x2
     bi2=x1-x3
     bi3=x2-x1
     ci1=x2*y3-x3*y2
     ci2=x3*y1-x1*y3
     ci3=x1*y2-x2*y1

     aw0(i,1)=-ci1/2./art(i)
     aw0(i,2)=-ci2/2./art(i)
     aw0(i,3)=-ci3/2./art(i)
     awx(i,1)=-ai1/2./art(i)
     awx(i,2)=-ai2/2./art(i)
     awx(i,3)=-ai3/2./art(i)
     awy(i,1)=-bi1/2./art(i)
     awy(i,2)=-bi2/2./art(i)
     awy(i,3)=-bi3/2./art(i)
   end do

!
!--------boundary cells------------------------------------------------!
!
   do i=1,n
     if(isbce(i) > 1) then
       do j=1,4
         a1u(i,j)=0.0_SP
         a2u(i,j)=0.0_SP
       end do
     else if(isbce(i) == 1) then
       do j=1,3
         if(nbe(i,j) == 0) jj=j
       end do
       j1=jj+1-int((jj+1)/4)*3
       j2=jj+2-int((jj+2)/4)*3
       x1=vx(nv(i,j1))-xc(i)
       x2=vx(nv(i,j2))-xc(i)
       y1=vy(nv(i,j1))-yc(i)
       y2=vy(nv(i,j2))-yc(i)

       delt=x1*y2-x2*y1
       b1=(y2-y1)/delt
       b2=(x1-x2)/delt
       deltx=vx(nv(i,j1))-vx(nv(i,j2))
       delty=vy(nv(i,j1))-vy(nv(i,j2))

       alpha(i)=atan2(delty,deltx)
       alpha(i)=alpha(i)-3.1415926_SP/2.0_SP
       x1=xc(nbe(i,j1))-xc(i)
       x2=xc(nbe(i,j2))-xc(i)
       y1=yc(nbe(i,j1))-yc(i)
       y2=yc(nbe(i,j2))-yc(i)
       temp1=x1*y2-x2*y1

       if(abs(temp1).lt.1.e-6_SP)  then
         print*, 'shape_f of solid b. c. temp1=0'
         print*, 'i,jj,j1,j2,x1,x2,y1,y2'
         print*, i,jj,j1,j2,x1,x2,y1,y2
         print*, 'x1*y2==',x1*y2
         print*, 'x2*y1==',x2*y1
         call pstop
       end if

       a1u(i,1)=0.0_SP
       a1u(i,jj+1)=0.0_SP
       a1u(i,j1+1)=0.0_SP
       a1u(i,j2+1)=0.0_SP

       a2u(i,1)=0.0_SP
       a2u(i,jj+1)=0.0_SP
       a2u(i,j1+1)=0.0_SP
       a2u(i,j2+1)=0.0_SP
     end if
   end do


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


   IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*) "!"
       WRITE(IPT,*) "!  INTERP COEFFICIENTS   :    COMPLETE"
    END IF

   return
   end subroutine shape_coef_gcn
