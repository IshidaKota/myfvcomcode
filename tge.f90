










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

   SUBROUTINE TRIANGLE_GRID_EDGE

! NEEDS UPDATE TO FIX ERROR MESSAGES

!==============================================================================!
!  This program is used to define the non-overlapped, unstructured             !
!  triangular meshes used for flux computations. The mesh could be             !
!  created using the commerical software called "sms8.0" or other              !
!  mesh generation programs. The mesh file generated by sms8.0 can             !
!  be directly used for this subroutine, while the mesh file                   !
!  generated using other programs must be coverted the data format             !
!  to meet the required format used here.                                      !
!==============================================================================!
!     variable list:							       !
!  vx(m)    :: vx(i) = x-coordinate of node i (input from mesh)	               !
!  vy(m)    :: vy(i) = y-coordinate of node i (input from mesh)	               !
!  nv(n,3)  :: nv(i:1-3) = 3 node numbers of element i                         !
!  xc(n)    :: xc(i) = x-coordinate of element i (calculated from vx)          !
!  yc(n)    :: yc(i) = y-coordinate of element i (calculated from vy)          !
!  cor(n)   :: cor(i) = f plane coriolis at element i                          !
!                                                                              !
!  nbe(n,3) :: nbe(i,1->3) = element index of 1->3 neighbors of element i      !
!  isbce(n) :: flag if element is on the boundary, see below for values        !
!  isonb(m) :: flag is node is on the boundary, see below for values           !
!                                                                              !
!  ntve(m)  :: the number of neighboring elements of node m                    !
!  nbve(m,ntve(m)) :: nbve(i,1->ntve(i)) = ntve elements containing node i     !
!  nbvt(m,ntve(m)) :: nbvt(i,j) = the node number of node i in element         ! 
!                     nbve(i,j) (has a value of 1,2,or 3)                      ! 
!                                                                              !
!   ne       ::   number of unique element edges                               !
!  iec(ne,2) :: nbe(i,1->2) cell number of cell(s) connected to edge i         !
!   isbc(ne) :: flag marking edge property                                     !
!              isbc(i) = 0:  element edge i is in the interior                 !
!              isbc(i) = 1:  element edge i is on the boundary                 !
!ienode(ne,2):: ienode(i,1->2) node numbers at each end of element edge i      ! 
!  xijc(ne)  :: xijc(i) = x-coordinate of mid point of element edge i          ! 
!  yijc(ne)  :: yijc(i) = y-coordinate of mid point of element edge i          ! 
!  dltxyc(ne):: dltxyc(i) = length of element edge i                           !
!  dltxc(ne) :: dltxc(i) = deltax (x-projection) of element edge i             !
!  dltyc(ne) :: dltyc(i) = deltay (y-projection) of element edge i             !
!  sitac(ne) :: sitac(i) =  arctg(dltyc,dltxc) (angle of inclination of edge)  !
!                                                                              !
!==============================================================================!
!     classification of the triangles nodes, and edges                         !
!                                                                              !
!     isonb(i)=0:  node in the interior computational domain                   !
!     isonb(i)=1:  node on the solid boundary                                  !
!     isonb(i)=2:  node on the open boundary                                   !
!                                                                              !
!     isbce(i)=0:  element in the interior computational domain                !
!     isbce(i)=1:  element on the solid boundary                               !
!     isbce(i)=2:  element on the open boundary                                !
!     isbce(i)=3:  element with 2 solid boundary edges                         !
!                                                                              !
!      isbc(i)=0:  element edge in interior                                    !
!      isbc(i)=1:  element edge on boundary                                    !
!==============================================================================!


!==============================================================================|
!   FIND NEIGHBORING ELEMENTS, MARK BOUNDARY NODES AND ELEMENTS                |
!									       |
!   NBE(N,3) :  NBE(I,J) = IDENTITY OF NBOR ELMNT TO TRI I ON EDGE J           |
!   IBCE(N)  :  DESCRIBED IN SUBROUTINE HEADING			               |	
!   ISONB(M):  DESCRIBED IN SUBROUTINE HEADING			               |	
!==============================================================================|
   USE ALL_VARS
   USE MOD_SPHERICAL


   USE MOD_PAR

   USE MOD_OBCS
   IMPLICIT NONE
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TEMP,TEMP2,NB_TMP,ISET
   INTEGER I,J,DIF1,DIF2,DIF3,II,JJ,NTMP,NCNT,INEY,NFLAG
   INTEGER ITMP1,ITMP2,ITMP3,JN,JJB,IBCETMP,NCTMP,NCETMP,NPT
   INTEGER, ALLOCATABLE :: CELLS(:,:),CELLCNT(:),NBET(:,:)
   REAL(SP)  DTMP
   INTEGER N1,N2,N3,J1,J2,J3,SBUF, IERR
   
   ! STORAGE FOR ISONB, ISBCE TO RUN EXCHANGE AND SET HALO'S
   REAL(SP), ALLOCATABLE :: FTEMP(:)

   Character(len=8)  :: tmpstr
   Character(len=200):: errstr


   REAL(SP) DELTX,DELTY,ALPHA1,ALPHA2

!----------------------------REPORT--------------------------------------------!
   IF (DBG_SET(DBG_LOG))THEN
      WRITE(IPT,*  )'!'
      WRITE(IPT,*)'!           SETTING UP TRIS/ELEMENTS/EDGES/CVS          '
      WRITE(IPT,*  )'!'
   END IF
!----------------------------INITIALIZE----------------------------------------!
   
   ISBCE = 0
   ISONB = 0
   NBE   = 0

!
!----DETERMINE NBE(i=1:n,j=1:3): INDEX OF 1 to 3 NEIGHBORING ELEMENTS----------!
!
   ALLOCATE(NBET(NT,3)) ; NBET = 0
   ALLOCATE(CELLS(MT,50)) ; CELLS = 0
   ALLOCATE(CELLCNT(MT))  ; CELLCNT = 0
   DO I=1,NT
     N1 = NV(I,1) ; CELLCNT(N1) = CELLCNT(N1)+1
     N2 = NV(I,2) ; CELLCNT(N2) = CELLCNT(N2)+1
     N3 = NV(I,3) ; CELLCNT(N3) = CELLCNT(N3)+1
     CELLS(NV(I,1),CELLCNT(N1)) = I
     CELLS(NV(I,2),CELLCNT(N2)) = I
     CELLS(NV(I,3),CELLCNT(N3)) = I
   END DO
   if(maxval(cellcnt) > 50)write(ipt,*)'bad',maxval(cellcnt)
   DO I=1,NT
     N1 = NV(I,1)
     N2 = NV(I,2)
     N3 = NV(I,3)
     DO J1 = 1,CELLCNT(N1) 
     DO J2 = 1,CELLCNT(N2) 
       IF((CELLS(N1,J1) == CELLS(N2,J2)).AND. CELLS(N1,J1) /= I)NBE(I,3) = CELLS(N1,J1)
     END DO
     END DO
     DO J2 = 1,CELLCNT(N2) 
     DO J3 = 1,CELLCNT(N3) 
       IF((CELLS(N2,J2) == CELLS(N3,J3)).AND. CELLS(N2,J2) /= I)NBE(I,1) = CELLS(N2,J2)
     END DO
     END DO
     DO J1 = 1,CELLCNT(N1) 
     DO J3 = 1,CELLCNT(N3) 
       IF((CELLS(N1,J1) == CELLS(N3,J3)).AND. CELLS(N1,J1) /= I)NBE(I,2) = CELLS(N3,J3)
     END DO
     END DO
   END DO
   DEALLOCATE(CELLS,CELLCNT)

   ! OLD - USED TO MAKE GLOBAL INDEX ARRAYS TO SAVE PARALLEL OUTPUT
   ! REPLACED IN MOD_NCDIO WITH MORE GENERIC TYPE GRID DATA!
!   IF (PAR) THEN
!      ALLOCATE(NBEGL(0:NT,3)); NBEGL      = 0
!      DO I = 1,NT
!         NBEGL(I,:) = EGID_X(NBE(I,:))
!      END DO
!   ELSE
!      NBEGL => NBE
!   END IF

   IF (DBG_SET(DBG_LOG))WRITE(IPT,*)  '!  NEIGHBOR FINDING      :    COMPLETE'
!
!--ENSURE ALL ELEMENTS HAVE AT LEAST ONE NEIGHBOR------------------------------!
!
   NFLAG = 0
   DO I=1,NT
     IF(SUM(NBE(I,1:3))==0)THEN 
       NFLAG = 1
       WRITE(IPT,*)'ELEMENT ',I,' AT ',XC(I),YC(I),' HAS NO NEIGHBORS'
       CALL PSTOP
     END IF
   END DO
   IF(NFLAG == 1) CALL PSTOP
     
!
!----IF ELEMENT ON BOUNDARY SET ISBCE(I)=1 AND ISONB(J)=1 FOR BOUNDARY NODES J-!
!

   DO I=1,NT 
     IF(MIN(NBE(I,1),NBE(I,2),NBE(I,3))==0)THEN    !!ELEMENT ON BOUNDARY
       ISBCE(I) = 1
       IF(NBE(I,1) == 0)THEN 
         ISONB(NV(I,2)) = 1 ; ISONB(NV(I,3)) = 1
       END IF
       IF(NBE(I,2) ==0) THEN
         ISONB(NV(I,1)) = 1 ; ISONB(NV(I,3)) = 1
       END IF
       IF(NBE(I,3) ==0) THEN
         ISONB(NV(I,1)) = 1 ; ISONB(NV(I,2)) = 1
       END IF
     END IF
   END DO
   IF (DBG_SET(DBG_LOG)) WRITE(IPT,*)  '!  ISONB SETTING         :    COMPLETE'
        
!==============================================================================|
!             DEFINE NTVE, NBVE, NBVT                                          !
!                                                                              !
! ntve(1:m):           total number of the surrounding triangles               !
!                      connected to the given node                             !
! nbve(1:m, 1:ntve+1): the identification number of surrounding                !
!                      triangles with a common node (counted clockwise)        !
! nbvt(1:m,ntve(1:m)): the idenfication number of a given node over            !
!                      each individual surrounding triangle(counted            !
!                      clockwise)                                              !
! ntsn(1:m):           total number of surrounding nodes                       !
! nbsn(1:m, ntsn):     the identification number of surrounding nodes          !
!                      (counted clockwise)                                     !
! nbse(1:m,2*ntsn):    the identification number of control volume s           !
!                      edges between two neighbor nodes                        !
!==============================================================================|

!
!----DETERMINE MAX NUMBER OF SURROUNDING ELEMENTS------------------------------!
!
   MX_NBR_ELEM = 0
   DO I=1,M
     NCNT = 0
     DO J=1,NT
!       IF( (NV(J,1)-I)*(NV(J,2)-I)*(NV(J,3)-I)==0) NCNT = NCNT + 1
!       IF( (NV(J,1)-I) ==0 .OR. (NV(J,2)-I) ==0 .OR. (NV(J,3)-I)==0) NCNT = NCNT + 1
       IF( FLOAT(NV(J,1)-I)*FLOAT(NV(J,2)-I)*FLOAT(NV(J,3)-I) == 0.0_SP) &
         NCNT = NCNT + 1
     END DO
     MX_NBR_ELEM = MAX(MX_NBR_ELEM,NCNT)
   END DO
!   WRITE(IPT,*) 'MAXIMUM NUMBER OF NEIGHBOR ELEMENTS',MX_NBR_ELEM

   IF(PAR) THEN
      SBUF = MX_NBR_ELEM
      CALL MPI_ALLREDUCE(SBUF,MX_NBR_ELEM,1,MPI_INTEGER,MPI_MAX,MPI_FVCOM_GROUP,IERR)
   END IF

!
!----ALLOCATE ARRAYS BASED ON MX_NBR_ELEM--------------------------------------!
! 
   ALLOCATE(NBVE(M,MX_NBR_ELEM+1)); NBVE = 0
   ALLOCATE(NBVT(M,MX_NBR_ELEM+1)); NBVT = 0
!   ALLOCATE(NBSN(M,MX_NBR_ELEM+2))
   ALLOCATE(NBSN(M,MX_NBR_ELEM+3)); NBSN = 0
!
!--DETERMINE NUMBER OF SURROUNDING ELEMENTS FOR NODE I = NTVE(I)---------------!
!--DETERMINE NBVE - INDICES OF NEIGHBORING ELEMENTS OF NODE I------------------!
!--DETERMINE NBVT - INDEX (1,2, or 3) OF NODE I IN NEIGHBORING ELEMENT---------!
!
       
   DO I=1,M
     NCNT=0
     DO J=1,NT
!       IF( (NV(J,1)-I) == 0 .OR.  (NV(J,2)-I) == 0 .OR. (NV(J,3)-I) == 0)THEN 
        IF (FLOAT(NV(J,1)-I)*FLOAT(NV(J,2)-I)*FLOAT(NV(J,3)-I) == 0.0_SP)THEN
         NCNT = NCNT+1
         NBVE(I,NCNT)=J
         IF((NV(J,1)-I) == 0) NBVT(I,NCNT)=1
         IF((NV(J,2)-I) == 0) NBVT(I,NCNT)=2
         IF((NV(J,3)-I) == 0) NBVT(I,NCNT)=3
       END IF
     ENDDO
     NTVE(I)=NCNT
   ENDDO
!
!--Reorder Order Elements Surrounding a Node to Go in a Cyclical Procession----!
!--Determine NTSN  = Number of Nodes Surrounding a Node (+1)-------------------!
!--Determine NBSN  = Node Numbers of Nodes Surrounding a Node------------------!
!
!Lettmann&JQI   ALLOCATE(NB_TMP(M,MX_NBR_ELEM+1))
   ALLOCATE(NB_TMP(MX_NBR_ELEM+1,2))
   DO I=1,M
     IF(ISONB(I) == 0) THEN
       NB_TMP(1,1)=NBVE(I,1)
       NB_TMP(1,2)=NBVT(I,1)
       DO J=2,NTVE(I)+1
         II=NB_TMP(J-1,1)
         JJ=NB_TMP(J-1,2)
         NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
         JJ=NB_TMP(J,1)
         IF((NV(JJ,1)-I) == 0) NB_TMP(J,2)=1
         IF((NV(JJ,2)-I) == 0) NB_TMP(J,2)=2
         IF((NV(JJ,3)-I) == 0) NB_TMP(J,2)=3
       ENDDO

       DO J=2,NTVE(I)+1
         NBVE(I,J)=NB_TMP(J,1)
       ENDDO

       DO J=2,NTVE(I)+1
         NBVT(I,J)=NB_TMP(J,2)
       ENDDO

       NTMP=NTVE(I)+1
       IF(NBVE(I,1) /= NBVE(I,NTMP)) THEN
          PRINT*, MYID,ngid(I),nTMP,'NBVE(I,nTMP) NOT CORRECT!!'
          PRINT*, "NBVE(I,:)=",EGID(NBVE(I,:))
          CALL PSTOP
       ENDIF

       IF(NBVT(I,1) /= NBVT(I,NTMP)) THEN
          PRINT*, ngid(I),'NBVT(I) NOT CORRECT!!'
          PRINT*, "NBVT(I,:)=",NBVT(I,:)
          CALL PSTOP
       END IF

       NTSN(I)=NTVE(I)

       DO J=1,NTSN(I)
         II=NBVE(I,J)
         JJ=NBVT(I,J)
         NBSN(I,J)=NV(II,JJ+1-INT((JJ+1)/4)*3)
       ENDDO

       NTSN(I)=NTSN(I)+1
       NBSN(I,NTSN(I))=NBSN(I,1)

     ELSE 
           JJB=0

       DO J=1,NTVE(I)
         JJ=NBVT(I,J)
         IF(NBE(NBVE(I,J),JJ+2-INT((JJ+2)/4)*3) == 0) THEN
           JJB=JJB+1
           NB_TMP(JJB,1)=NBVE(I,J)
           NB_TMP(JJB,2)=NBVT(I,J)
         END IF
       ENDDO

       IF(JJB /= 1) THEN
         WRITE(IPT,*) 'ERROR IN ISONB !,I,J', I,J
         CALL FATAL_ERROR("ERROR IN TGE DETERMINING ISONB")
       END IF

       DO J=2,NTVE(I)
         II=NB_TMP(J-1,1)
         JJ=NB_TMP(J-1,2)
         NB_TMP(J,1)=NBE(II,JJ+1-INT((JJ+1)/4)*3)
         JJ=NB_TMP(J,1)
         IF((NV(JJ,1)-I) == 0) NB_TMP(J,2)=1
         IF((NV(JJ,2)-I) == 0) NB_TMP(J,2)=2
         IF((NV(JJ,3)-I) == 0) NB_TMP(J,2)=3
       ENDDO

       DO J=1,NTVE(I)
         NBVE(I,J)=NB_TMP(J,1)
         NBVT(I,J)=NB_TMP(J,2)
       ENDDO

       NBVE(I,NTVE(I)+1)=0
       NTSN(I)=NTVE(I)+1
       NBSN(I,1)=I

       DO J=1,NTSN(I)-1
         II=NBVE(I,J)
         JJ=NBVT(I,J)
         NBSN(I,J+1)=NV(II,JJ+1-INT((JJ+1)/4)*3)
       ENDDO

       J=NTSN(I)
       II=NBVE(I,J-1)
       JJ=NBVT(I,J-1)
       NBSN(I,J+1)=NV(II,JJ+2-INT((JJ+2)/4)*3)
       NTSN(I)=NTSN(I)+2
       NBSN(I,NTSN(I))=I
     END IF
   END DO
   DEALLOCATE(NB_TMP)
   IF(MX_NBR_ELEM+3 -MAXVAL(NTSN) < 0)THEN
      WRITE(IPT,*)'CHECK NTSN/NBSN',MAXVAL(NTSN),MX_NBR_ELEM+3
      CALL PSTOP
   END IF
  

   ! OLD - USED TO MAKE GLOBAL INDEX ARRAYS TO SAVE PARALLEL OUTPUT
   ! REPLACED IN MOD_NCDIO WITH MORE GENERIC TYPE GRID DATA!
!   IF(PAR) THEN
!      ALLOCATE(NBSNGL(M,MX_NBR_ELEM+3))
!      ALLOCATE(NBVEGL(M,MX_NBR_ELEM+1))
!      DO I = 1,M
!         NBVEGL(I,:)=EGID_X(NBVE(I,:))
!
!         NBSNGL(I,:) = NGID_X(NBSN(I,:))
!      END DO
!   ELSE
!      NBVEGL => NBVE
!
!      NBSNGL =>NBSN
!   END IF

   IF (DBG_SET(DBG_LOG))WRITE(IPT,*)  '!  NBVE/NBVT             :    COMPLETE'
     

!==============================================================================!
!  Define the parameters of each triangular edge                               !
!                                                                              !
!  ne           :    number of unique element edges                            !
!  iec(1:ne,1:2):    counting number identifying two connected cells           !
!  isbc(1:ne):       0: triangle s edge in the interior                        !
!                    1: triangle s edge on the boundary                        !
!  ienode(1:ne,1:2): the identification number of two end points of a          !
!                    edge                                                      !
!  xijc(1:ne):       the x-coordinate location of the middle points            !
!                    of a edge                                                 !
!  yijc(1:ne):       the y-coordinate location of the middle points            !
!                    of a edge                                                 !
!  dltxyc(1:ne):     length of the edge                                        !
!  dltxc(1:ne):      vx(ienode(i,2))-vx(idnode(i,1))                           !
!  dltyc(1:ne):      vy(ienode(i,2))-vy(idnode(i,1))                           !
!  sitac(1:ne):      arctg(dltyc,dltxc)                                        !
!==============================================================================!

   ALLOCATE(ISET(NT,3),TEMP((NT)*3,2),TEMP2((NT)*3,2))
   ISET = 0
   NE = 0
   TEMP = 0 
   TEMP2 = 0
   DO I=1,NT
     DO J=1,3
       IF(ISET(I,J) == 0)THEN
         NE   = NE + 1
         INEY = NBE(I,J)
         ISET(I,J) = 1
         DO JN=1,3
           IF(I == NBE(INEY,JN)) ISET(INEY,JN) = 1
         END DO
         TEMP(NE,1) = I ; TEMP(NE,2) = INEY
         TEMP2(NE,1) = NV(I,J+1-INT((J+1)/4)*3)
         TEMP2(NE,2) = NV(I,J+2-INT((J+2)/4)*3)
       END IF
     END DO
   END DO
   DEALLOCATE(ISET)
!
!--ALLOCATE ARRAYS REQUIRING NUMBER OF EDGES-----------------------------------!
!
   ALLOCATE(IEC(NE,2))
   ALLOCATE(IENODE(NE,2))
   ALLOCATE(XIJC(NE))
   ALLOCATE(YIJC(NE))
   ALLOCATE(DLTXYC(NE))
   ALLOCATE(DLTXC(NE))
   ALLOCATE(DLTYC(NE))
   ALLOCATE(SITAC(NE))
   ALLOCATE(ISBC(NE))


   IEC(:,1) = TEMP(1:NE,1)
   IEC(:,2) = TEMP(1:NE,2)
   IENODE(:,1) = TEMP2(1:NE,1)
   IENODE(:,2) = TEMP2(1:NE,2)


   DEALLOCATE(TEMP,TEMP2)
         
         
!
!------MARK ELEMENT EDGES THAT ARE ON THE BOUNDARY-----------------------------!
!
   ISBC = 0
   DO I=1,NE
     IF((IEC(I,1) == 0) .OR. (IEC(I,2) == 0)) ISBC(I) = 1 
   END DO

!
!------CALCULATE ELEMENT EDGE METRICS------------------------------------------!
!
   DO I=1,NE
     DLTXC(I) =  VX(IENODE(I,2))-VX(IENODE(I,1))
     DLTYC(I) =  VY(IENODE(I,2))-VY(IENODE(I,1))
     XIJC(I)  = (VX(IENODE(I,1))+VX(IENODE(I,2)))/2.0_SP
     YIJC(I)  = (VY(IENODE(I,1))+VY(IENODE(I,2)))/2.0_SP
     DLTXYC(I)= SQRT(DLTXC(I)**2+DLTYC(I)**2)
     SITAC(I) = ATAN2(DLTYC(I),DLTXC(I))
   END DO
    IF (DBG_SET(DBG_LOG))WRITE(IPT,*)  '!  EDGE SETUP            :    COMPLETE'


!==============================================================================!
!  read triangular mesh parameters on open boundary :                          !
!  iobce:   number of open boundary cells.                                     !
!  isbcn:   number of open boundary nodes.                                     !
!  i_obc_e: counter number of open boundary cells                              !
!  i_obc_n: counter number of open boundary nodes                              !
!==============================================================================!

!
!----TRAVERSE  BOUNDARY NODE NUMBERS AND SET ISONB(NODE)=2---------------------!
!
   DO I=1,IOBCN
     ISONB(I_OBC_N(I))=2
   ENDDO

!
!---- SET HALO VALUES IF PAR
! 

   IF(PAR) THEN
      CALL AEXCHANGE(NC,MYID,NPROCS,ISONB)
   END IF

!
!----DETERMINE IF ELEMENT IS ON OPEN BOUNDARY (CONTAINS EDGE ON OPEN BOUNDARY)-!
!
   IBCETMP=0
   DO I=1,N
     ITMP1=ISONB(NV(I,1))
     ITMP2=ISONB(NV(I,2))
     ITMP3=ISONB(NV(I,3))

     IF(SUM(ISONB(NV(I,1:3))) == 4) THEN
       ISBCE(I)=2
       IBCETMP =IBCETMP+1
     ELSE IF(SUM(ISONB(NV(I,1:3))) > 4) THEN
       PRINT*,'SORRY, THE BOUNDARY CELL',EGID(I),'IS NOT GOOD FOR MODEL.'
       PRINT*,'IT HAS EITHER TWO SIDES OF OPEN BOUNDARY OR ONE OPEN BOUNDARY'
       PRINT*,'AND ONE SOLID BOUNDARY. PLEASE CHECK AND MODIFIED IT.'
       PRINT*,'THIS MESSAGE IS IN SUBROUTINE TRIANGLE_GRID_EDGE (TGE.F)'
       PRINT*,'STOP RUNNING...'
       CALL PSTOP

     END IF
   END DO

   DO I=1,NT
     IF((NBE(I,1)+NBE(I,2)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,1)+NBE(I,2) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,2)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
     IF((NBE(I,1)+NBE(I,3) == 0).AND.(ISBCE(I) /= 2)) ISBCE(I)=3
   ENDDO

! SET HALO VALUES IF PAR
   IF(PAR) THEN
      CALL AEXCHANGE(EC,MYID,NPROCS,ISBCE)
   END IF


!==============================================================================!
!  xije(1:nc,1:2):  the x coordinate locations of starting and ending          !
!                   points of the control volumes edges                        !
!  yije(1:nc,1:2):  the y coordinate locations of starting and ending          !
!                   points of the control volumes edges                        !
!  niec(1:nc,1:2):  the counting number of left and right nodes                !
!                   conected to this control volumes edge from                 !
!                   starting point to ending point                             !
!  dltxe(1:nc):     the x distance of individual edges                         !
!  dltye(1:nc)      the y distance of individual edges                         !
!  dltxye(1:nc):    the length of individual edges                             !
!  ntrg(1:nc)  :    element associated with this control volume edge           !
!==============================================================================!
   NCTMP  = 0
   NCETMP = 0 

   DO I=1,NE
     IF(ISBC(I) == 0) THEN
       IF(IEC(I,1) <= N)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*N)
       END IF
       XIJE(NPT,1) = XC(IEC(I,1))
       YIJE(NPT,1) = YC(IEC(I,1))
       XIJE(NPT,2) = XIJC(I)
       YIJE(NPT,2) = YIJC(I)
       NIEC(NPT,1) = IENODE(I,1)
       NIEC(NPT,2) = IENODE(I,2)
       NTRG(NPT)   = IEC(I,1)
       DLTXE(NPT)  = XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)  = YIJE(NPT,2)-YIJE(NPT,1)
       DTMP        = DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT) = SQRT(DTMP)
       SITAE(NPT)  = ATAN2(DLTYE(NPT),DLTXE(NPT))

       IF(IEC(I,2) <= N)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*N)
       END IF
       XIJE(NPT,1)=XC(IEC(I,2))
       YIJE(NPT,1)=YC(IEC(I,2))
       XIJE(NPT,2)=XIJC(I)
       YIJE(NPT,2)=YIJC(I)
       NIEC(NPT,1)=IENODE(I,2)
       NIEC(NPT,2)=IENODE(I,1)
       NTRG(NPT)=IEC(I,2)
       DLTXE(NPT)=XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)=YIJE(NPT,2)-YIJE(NPT,1)
       DTMP=DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT)=SQRT(DTMP)
       SITAE(NPT)=ATAN2(DLTYE(NPT),DLTXE(NPT))
     ELSE IF(ISBC(I) == 1) THEN
       IF(IEC(I,1) <= N)THEN
         NCTMP=NCTMP+1
         NPT  =NCTMP
       ELSE
         NCETMP = NCETMP + 1
         NPT    = NCETMP+(3*N)
       END IF
       IF(IEC(I,1) == 0) THEN
         PRINT*, I,'IEC(I,1)===0'
         CALL PSTOP
       END IF
       XIJE(NPT,1)=XC(IEC(I,1))
       YIJE(NPT,1)=YC(IEC(I,1))
       XIJE(NPT,2)=XIJC(I)
       YIJE(NPT,2)=YIJC(I)
       NIEC(NPT,1)=IENODE(I,1)
       NIEC(NPT,2)=IENODE(I,2)
       NTRG(NPT)=IEC(I,1)
       DLTXE(NPT)=XIJE(NPT,2)-XIJE(NPT,1)
       DLTYE(NPT)=YIJE(NPT,2)-YIJE(NPT,1)
       DTMP=DLTXE(NPT)*DLTXE(NPT)+DLTYE(NPT)*DLTYE(NPT)
       DLTXYE(NPT)=SQRT(DTMP)
       SITAE(NPT)=ATAN2(DLTYE(NPT),DLTXE(NPT))
     ELSE
       WRITE(IPT,*) 'ISBC(I) NOT CORRECT, I==',I
       CALL PSTOP
     END IF
   ENDDO

   NCV_I  = NCTMP
   NCV    = NCETMP+NCTMP
      
   IF(NCV /= 3*(NT)) THEN
     PRINT*,'NCV IS NOT CORRECT, PLEASE CHECK THE SETUP'
     CALL PSTOP
   END IF
   IF(NCV_I /= 3*N) THEN
     PRINT*,'NCV_I IS NOT CORRECT, PLEASE CHECK THE SETUP'
      CALL PSTOP
   END IF

   DO I=1,NCV_I
     IF(NIEC(I,1) > M .OR. NIEC(I,2) > M)THEN
       write(ipt,*)'problemas',niec(i,1),niec(i,2),m
       CALL PSTOP
     END IF
   END DO

!==============================================================================!
!  nisbce_1/nisbce_2/nisbce_3:  number of elements with isbce of 1,2,3         !
!  lisbce_1/lisbce_2/lisbce_3:  list of elements with isbce of 1,2,3           !
!  epor                      :  element porosity (=0 if isbce = 2)             !
!==============================================================================!

!
!--COUNT NUMBER OF ELEMENTS OF EACH TYPE (ISBCE=1,2,3)-------------------------!
!
   NISBCE_1 = 0 ; NISBCE_2 = 0 ; NISBCE_3 = 0
   DO I=1,N 
     IF(ISBCE(I) == 1) NISBCE_1 = NISBCE_1 + 1
     IF(ISBCE(I) == 2) NISBCE_2 = NISBCE_2 + 1
     IF(ISBCE(I) == 3) NISBCE_3 = NISBCE_3 + 1
   END DO

!
!--ALLOCATE ELEMENT TYPE ARRAYS LISBCE_1,LISBCE_2,LISBCE_3---------------------!
!
   IF(NISBCE_1 > 0)THEN
     ALLOCATE( LISBCE_1(NISBCE_1) )
   ELSE
!     WRITE(IPT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=1'
   END IF
     
   IF(NISBCE_2 > 0)THEN
     ALLOCATE( LISBCE_2(NISBCE_2) )
   ELSE
!     WRITE(IPT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=2'
   END IF
     
   IF(NISBCE_3 > 0)THEN
     ALLOCATE( LISBCE_3(NISBCE_3) )
   ELSE
!     WRITE(IPT,*)  '!  WARNING               :    NO ELEMENTS WITH ISBCE=3'
   END IF

!
!--LOAD ELEMENT TYPE ARRAYS LISBCE_1,LISBCE_2,LISBCE_3--------------------------!
!
   NISBCE_1 = 0 ; NISBCE_2 = 0 ; NISBCE_3 = 0
   DO I=1,N 
     IF(ISBCE(I) == 1) THEN
       NISBCE_1 = NISBCE_1 + 1
       LISBCE_1(NISBCE_1) = I
     END IF
     IF(ISBCE(I) == 2) THEN
       NISBCE_2 = NISBCE_2 + 1
       LISBCE_2(NISBCE_2) = I
     END IF
     IF(ISBCE(I) == 3) THEN
       NISBCE_3 = NISBCE_3 + 1
       LISBCE_3(NISBCE_3) = I
     END IF
   END DO

!
!--SET ELEMENT POROSITY---------------------------------------------------------!
!
   ALLOCATE(EPOR(0:NT)) ; EPOR = 1.0_SP
   DO I=1,N
     IF(ISBCE(I) == 2) EPOR(I) = 0.0_SP
   END DO
     
   IF (DBG_SET(DBG_LOG))WRITE(IPT,*)  '!  NISBCE/LISBCE/EPOR    :    COMPLETE'
   IF (DBG_SET(DBG_LOG))WRITE(IPT,*)  '!  TRIS/EDGES/CVOLS      :    COMPLETE'

   RETURN
   END SUBROUTINE TRIANGLE_GRID_EDGE
!==============================================================================!
