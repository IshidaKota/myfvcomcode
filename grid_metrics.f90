










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

SUBROUTINE GRID_METRICS
  USE MOD_PAR
  USE MOD_OBCS
  USE ALL_VARS
  USE MOD_NORTHPOLE
  USE MOD_SETUP
  USE MOD_STATION_TIMESERIES
  USE MOD_SPARSE_TIMESERIES
  USE MOD_NESTING, ONLY: NESTING_ON


  IMPLICIT NONE
  
  INTEGER:: STATUS
  INTEGER:: I

  ! THESE SUBROUTINES SHOULD BE FURTHER BROKEN DOWN INTO THEIR
  ! COMPONENT OPERATIONS AND THE WHOLE THING SHOULD BE MODULAR
  
  
!============================================
!Set up fluxes and control Volumes
!============================================
  CALL TRIANGLE_GRID_EDGE
  IF(OUT_STATION_TIMESERIES_ON) CALL TRIANGLE_GRID_EDGE_GL    
!============================================
!Calculate Element and Control Volume Areas
!============================================
  CALL CELL_AREA            


!====================================================
! Calculate Shape Coefficients for Flux Construction
!====================================================
  CALL SHAPE_COEF_GCN       


!============================================
!Calculate Node Control Volume Edge Lengths
!============================================
  CALL EDGE_LEN


!====================================================
!  EXCHANGE SHAPE FACTOR INFORMATION
!====================================================
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,A1U,A2U) 
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,AWX,AWY,AW0) 
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,ALPHA) 
  IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,ART)

  ! THESE MUST BE SET CORRECTLY IN THE HALO FOR OUTPUT!
  IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,ART1,ART2)
  IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,NTSN,NTVE)

!====================================================
! Calculate the horizontal gradient of the bathymetry
!====================================================
 CALL DEPTH_GRADIENT


!====================================================
! SETUP BOUNDARY NEIGHBORS AND SUCH
!====================================================
  CALL SETUP_OBC

!-------------------Read dam info--------------------

!====================================================
! SETUP THE LONGSHORE FLOW INDEX AND ANGLES
!====================================================
  IF(OBC_LONGSHORE_FLOW_ON) CALL longshore_flow

END SUBROUTINE GRID_METRICS

