










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

!==============================================================================!
!  VERSION 4.4.2
!==============================================================================!

PROGRAM FVCOM

  !================================================================================!
  !                                                                                !
  !                           USG-FVCOM                                            !
  !    The Unstructured Grid Finite Volume Coastal Ocean Model                     !
  !                                                                                !
  !    The USG-FVCOM (publicly called FVCOM) was developed by Drs. Changsheng Chen !
  !  and Hedong Liu at the Marine Ecosystem Dynamics Modeling Laboratory at the    !
  !  School of Marine Science and Technology (SMAST), University of Massachusetts- !
  !  Dartmouth (UMASSD) and Dr. Robert C. Beardsley at the Department of Physical  !
  !  Oceanography, Woods Hole Oceanographic Institution (WHOI). This code was      !
  !  rewritten in FORTRAN 90/2K, modularized, rendered somewhat understandable,    !
  !  and parallelized by Geoff Cowles at SMAST/UMASSD. FVCOM is being upgraded by  !
  !  the UMASSD-WHOI joint team led by Dr. Chen and Dr. Beardsley.                 ! 
  !                                                                                !
  !    The Development was initially supported by the Georgia Sea Grant College    !
  !  Program for the study of the complex dynamic system in Georgia estuaries. The !
  !  code improvement has been supported by Dr. Chen's research grants received    !
  !  from NSF and NOAA Coastal Ocean research grants and SMAST fishery research    !
  !  grants.                                                                       !
  !                                                                                !
  !    FVCOM is a prognostic, unstructured grid, Finite-Volume free-surface three- !
  !  dimensional hydrostatic/non-hydrostatic primitive equations, Coastal Ocean    !
  !  circulation Model (Chen et al., 2003; 2006, 2007). The model computes the     !
  !  momentum, continuity, temperature salinity, and density equations and is      !
  !  closed physically and mathematically using user-specified turbulent parameters!
  !  or turbulent closure submodel. The irregular bottom slope is represented using!
  !  a generalized terrain-following coordinate transformation with spatially      !
  !  variable vertical distribution (Pietrzak et al., 2002). The horizontal grids  !
  !  comprise unstructured triangular cells. The spatial fluxes of momentum are    !
  !  discretized using a second-order accurate finite-volume method (Kobayashi et  !
  !  al. 1999) in conjunction with a vertical velocity adjustment enforce exact    !
  !  volume conservation in individual control volumes. The flux of scalars (e. g. !
  !  temperature, salinity) is calculated by the second-order accurate upwind      !
  !  scheme. The finite-volume method (FVM) used in this model combines the        !
  !  advantages of the finite-element method (FEM) for geometric flexibility and   !
  !  the finite-difference method (FDM) for simple discrete computation. Current,  !
  !  temperature, and salinity in the model are computed in the integral form of   !
  !  the equations, which provides a better representation of the conservation laws!
  !  for mass, momentum, and heat in the coastal region with complex geometry.     !
  !                                                                                !
  !    FVCOM was originally developed for the coastal and estuarine applications.  !
  !  With implementation of spherical coordinates, this model is capable of using  !
  !  for the global ocean.                                                         !
  !                                                                                !
  !    All users should read this agreement carefully.  A user, who receives any   !
  !  version of the source code of FVCOM, must accept all the terms and conditions !
  !  of this agreement and also agree that this agreement is like any written      !
  !  negotiated agreement signed by you. You may be required to have another       !
  !  written agreement directly with Dr. Changsheng Chen at SMAST/UMASSD that      !
  !  supplements all or portions of this agreement. Dr. Changsheng Chen, leader of !
  !  the FVCOM development team, owns all intellectual property rights to the      !
  !  software. The University of Massachusetts-Dartmouth owns the copyright of the !
  !  software. All copyrights are reserved. Unauthorized reproduction and re-      !
  !  distribution of this program are expressly prohibited. This program is only   !
  !  permitted for use in non-commercial academic research and education.          !
  !  Commercial use must be approved by Dr. Chen. Registration is required for all !
  !  new users. Users should realize that this model software is a research product!
  !  without any warranty. Users must take full responsibility for any mistakes    !
  !  that might be caused by any inappropriate modification of the source code.    !
  !  Modification is not encouraged for users who do not have a deep understanding !
  !  of the finite-volume numerical methods used in FVCOM. Contributions made to   !
  !  correcting and modifying the programs will be credited, but will not affect   !
  !  copyrights. No duplicate configurations of FVCOM are allowed in the same      !
  !  geographical region, particularly in the regions where FVCOM has been already !
  !  been applied. Users who want to use FVCOM in a region that the SMAST/UMASS    !
  !  Marine Ecosystem Dynamics Modeling (MEDM) group (led by Dr. Chen) is working  !
  !  on must request permission from Dr. Chen. No competition is allowed in the    !
  !  same region using FVCOM, especially with Dr. Chen's group. FVCOM has been     !
  !  validated for many standard model test cases.  Users are welcome to do any    !
  !  further model validation experiments. These experiments shall be carried out  !
  !  in collaboration with the SMAST/UMASSD model development team. To avoid or    !
  !  reduce deriving any incorrect conclusions due to an inappropriate use of      !
  !  FVCOM, users are required to contact the scientific leaders of the FVCOM      !
  !  development team (Dr. Chen at SMAST/UMASSD and Dr. Beardsley at WHOI) before  !
  !  any formal publications are prepared for model validation.                    !
  !                                                                                !
  !    For public use, all users should name this model as "FVCOM". In any         !
  !  publications with the use of FVCOM, acknowledgement must be included. The     !
  !  rationale behind this FVCOM distribution policy is straightforward.  New      !
  !  researchers and educators who want to use FVCOM and agree to the above        !
  !  requirements get free access to the latest version of FVCOM and the collective!
  !  wisdom and experience of the FVCOM development team and existing users.       !
  !  Problems arising in new FVCOM applications, both related to conceptual as well!
  !  as numerical and coding issues, can be shared with the development team and   !
  !  other users who can work together on physics and code improvements that over  !
  !  time will lead to a better FVCOM.                                             !
  !                                                                                !
  !    FVCOM has been developed to date with state and federal funding with the    !
  !  idea that FVCOM will become a community model that new users start to use the !
  !  model and its scientific usefulness and numerical accuracy and efficiency     !
  !  continue to improve. The FVCOM distribution policy is designed to encourage   !
  !  this transition while maintaining a central core group responsible for overall!
  !  FVCOM development direction, implementing official code improvements, and     !
  !  maintaining well tested and documented updated code versions.                 !
  !                                                                                !
  !                                                                                !
  !  External forces used to drive this model:                                     !
  !                                                                                !
  !  1) Tidal amplitudes and phases at open boundaries (initial designs include 8  !
  !         tidal constituents, more can be added as needed) plus astronomical     !
  !         forcing implemented via gradients of tidal-generating potential for    !
  !         each tidal constituent and various corrections due to the Earth tidal  !
  !         and ocean loading.                                                     !
  !  2) Wind stress [3 ways: a) uniform wind speed and direction, b) spatially     !
  !         distributed wind velocity field, and c) the weather model output wind  !
  !         fields (NCEP, MM5, WRF, ECMWF).                                        !
  !  3) Air pressure gradients [2 ways: weather model output and Hurricane or      !
  !         typhoon model predicted.                                               !
  !  4) Surface heat flux [3 ways: a) uniform heat flux, b) spatially distributed  !
  !         heat flux, and c) the weather model-output heat flux fields. All the   !
  !         surface net heat flux and short-wave radiation are needed in the input !
  !         file.                                                                  !
  !  5) Precipitation via evaporation.                                             !
  !  6) River discharges: specify the location and discharge volume, temperature,  !
  !         and salinity.                                                          !
  !  7) Groundwater input: currently diffused bottom flux only.                    !
  !                                                                                !
  !  Note: FVCOM can be run directly by coupling with the weather model (WRF or    !
  !        MM5). This is used in the US Northeast Coastal Ocean Forecast System    !
  !        (NECOFS) (http://fvcom.smast.umassd.edu). Modes use for air-sea coupling!
  !        is not included in this release.                                        !
  !                                                                                !
  !  Hydrodynamics initial conditions:                                             !
  !                                                                                !
  !  The model can be prognostically run for both barotropic and baroclinic cases. !
  !                                                                                !
  !  Tidal forcing can be added into the system with zero velocity field at initial!
  !  or specified the 3-D tidal initial velocity field using the model-predicted   !
  !  harmonic tidal currents.                                                      !
  !                                                                                !
  !  Initial fields of temperature and salinity needed to be specified by using    !
  !  either limatological field, real-time observed field or idealized functions.  !
  !  The model has included Gregorian time for the time simulation for tidal       !
  !  currents.                                                                     !
  !                                                                                !
  !  For the purpose of interdisciplinary studies, biological, chemical, and       !
  !  sediment suspension models are available for FVCOM.  These submodels are      !
  !  directly driven by the FVCOM physical model. These submodels are called       !
  !  offline FVCOM modes.  A description of these submodels follows.               !
  !                                                                                !
  !  Generalized Biological Model (GBM)-a software platform that allows users to   !
  !        build his own biological model in FVCOM (developed by a team including  !
  !        Drs. Chen, Tian and Qi)                                                 !
  !                                                                                !
  !  NPZ model-built using GBM for 3 component nutrient-phytoplankton-zooplankton. !
  !                                                                                ! 
  !  NPZD model-an 4 and 8 component nutrient-phytoplankton-zooplankton-detritus   !
  !        model constructed by GBM                                                !
  !                                                                                !
  !  NPZDB-model-a 9 phosphorus-controlled component nutrient-phytoplankton-       !
  !        zooplankton-detritus-bacteria model constructed by GBM                  !
  !                                                                                !
  !  FVCOM-WQM: The water quality model modified from US-EPA WASP with inclusion of!
  !        the benthic process. This code was converted from the structured grid   !
  !        WQM (developed by Zheng and Chen) to the unstructured grid version.     !
  !                                                                                !
  !  UG-RCA: The unstructured grid Row/column version of the water quality model   !
  !        ESOP (RCA). UG-RCA developed by the FVCOM development team (Drs. Chen,  !
  !        J. Qi and R. Tian) on converting HydroQual's structured grid RCA to     !
  !        unstructured grid finite-volume version.                                !
  !                                                                                !
  !  UG-CE-QUAL-ICM: A unstructured grid version of the army corp water quality    !
  !        model. The code was written by Drs. Qi and Chen at UMASSD and modified  !
  !        and tested by Drs. Kim, Labiosa, Khanhaonkar and Yang at PNL.           !
  !                                                                                !
  !  FVCOM-SED: The 3-D sediment module (developed by Cowles based on the U.S.G.S. !
  !        national community sediment transport model and modified by the FVCOM   !
  !        development team staff-Ge, Wu, Chen and FVCOM users to include cohesive !
  !        process). This code can run online coupled with FVCOM or offline driven !
  !        by the FVCOM output                                                     ! 
  !                                                                                ! 
  !  FVCOM-SWAVE: The unstructured-grid (UG) version of the Simulating Wave        !
  !        Nearshore (SWAN) surface wave model. The code was developed by Drs. Qi  !
  !        and Chen at SMAST/UMASSD (Qi et al. 2008)                               !
  !                                                                                !
  !  Lagrangian particle tracking:                                                 !
  !                                                                                !
  !  FVCOM-LAG: A bilinear interpolation scheme is used to determine the particle  !
  !        velocity for the Lagrangian particle tracking. A random walk process    !
  !        also could be included with a specified function related to horizontal  !
  !        and vertical diffusion coefficients                                     !
  !                                                                                !
  !  Key references:                                                               !
  !                                                                                !
  !  Chen, C., H. Liu, and R. C. Beardsley, 2003. An unstructured grid, finite-    !
  !      volume, three-dimensional, primitive equations ocean model: application to!
  !      coastal ocean and estuaries, Journal of Atmospheric and Oceanic           !
  !      Technology,  20, 159-186.                                                 !
  !                                                                                !
  !  Chen, C, R. C. Beardsley and G. Cowles, 2006. An unstructured grid, finite-   !
  !      volume coastal ocean model (FVCOM) system. Special Issue entitled Advance !
  !      in Computational Oceanography, Oceanography, 19(1), 78-89.                !
  !                                                                                !
  !  Chen, C. H. Huang, R. C. Beardsley, H. Liu, Q. Xu and G. Cowles, 2007. A      !
  !      finite-volume numerical approach for coastal ocean circulation studies:   !
  !      comparisons with finite difference models. J. Geophys. Res. 112, C03018,  !
  !      doi:10.1029/2006JC003485.                                                 !
  !                                                                                !
  !  Chen, C., P. Malanotte-Rizzoli, J. Wei, R. C. Beardsley, Z. Lai, P. Xue, S.   !
  !      Lyu, Q. Xu, J. Qi, and G. W. Cowles, 2009. Application and comparison of  !
  !      Kalman filters for coastal ocean problems: An experiment with FVCOM, J.   !
  !      Geophys. Res., 114, C05011, doi:10.1029/2007JC004548.                     !
  !                                                                                !
  !  Cowles, G., 2008. Parallelization of the FVCOM Coastal Ocean Model,           !
  !      International Journal of High Performance Computing Applications, 22(2),  !
  !      177-193.                                                                  !
  !                                                                                !
  !  Huang, H. C. Chen, G. Cowles, C. D. Winant, R. C. Beardsley, K. Hedstrom, D.  !
  !      B. Haidvogel, 2008. FVCOM validation experiments: comparisons with ROMS   !
  !      for three idealized test problems. J. Geophys. Res, 113, C07042, doi:     !
  !      10.1029/2007JC004557.                                                     !
  !                                                                                !
  !  Lai, Z, C. Chen, G. Cowles and R. C. Beardsley, 2009. A non-hydrostatic       !
  !      version of FVCOM, part I: validation experiments. J. Geophys. Res.,       !
  !      in revision.                                                              !
  !                                                                                !
  !  Lai, Z., C. Chen, G. Cowles and R. C. Beardsley, 2009. A non-hydrostatic      !
  !      version of FVCOM, part II: mechanistic study of tidally generated         !
  !      nonlinear internal waves in Massachusetts Bay. J. Geophys. Res.,          !
  !      submitted.                                                                !
  !                                                                                !
  !  Qi, J. C. Chen, R. C. Beardsley, W. Perrie G. Cowles and Z. Lai, 2009. An     !
  !      unstructured-grid finite-volume surface wave model (FVCOM-SWAVE):         !
  !      implementation, validations and applications. Ocean Modelling, 28,        !
  !      153-166. doi:10.1016/j.ocemod.2009.01.007.                                !
  !                                                                                !
  !  Please direct criticisms and suggestions to                                   !
  !                                                                                !
  !               Changsheng Chen                                                  !
  !               School for Marine Science and Technology                         ! 
  !               University of Massachusetts-Dartmouth                            !
  !               New Bedford, MA 02742                                            !
  !               Phone: 508-910-6388, Fax: 508-910-6371                           !
  !               E-mail: c1chen@umassd.edu, cchen@whoi.edu                        !
  !               Web: http://fvcom.smast.umassd.edu                               !
  !                                                                                !
  !  What are new for version 3.1?                                                 !
  !                                                                                ! 
  !  1) Fully sea-ice coupling. The sea ice model was developed by converting the  !
  !     structured grid CICE into the unstructured grid version called UG-CICE).   !
  !     This is one component of G. Gao's Ph.D. thesis work supervised by C. Chen, !
  !     R. C. Beardsley and A. Proshutinsky. The code was written by G. Gao with   !
  !     assistance of J. Qi and C. Chen and tested by G. Gao and C. Chen.          !
  !                                                                                !
  !  2) Non-hydrostatic version. The non-hydrostatic dynamics is added into FVCOM  !
  !     as one component of Z. Lai's Ph.D. thesis work supervised by C. Chen, G.   !
  !     Cowles and R. C. Beardsley.  The code was written by Z. Lai with assistance!
  !     of C. Chen and G. Cowles, and tested by Z. Lai and C. Chen.                !
  !                                                                                ! 
  !  3) Fully current-wave interaction. An unstructured grid surface wave model    !
  !     (called FVCOM-SWAVE) was developed by Qi et al. (2008). This wave model    !
  !     is modified from SWAN by converting it to unstructured grid version using  !
  !     the second-order accurate FVCOM algorithms. The implementing FVCOM-SWAVE   !
  !     into FVCOM was initialized by A. Wu at SMAST/UMASSD as one component of his!
  !     Ph.D. thesis work. The work is supervised by C. Chen at SMAST/UMASSD. The  !
  !     code is re-organized, modified and re-debugged to fit more flexible and    !
  !     general setup of FVCOM by J. Qi and C. Chen.                               !
  !                                                                                !
  !  4) Fully current-wave-sediment interaction. Coupling of wave-sediment was     !
  !     initialized by A. Wu and J. Ge at SMAST/UMASSD as a visiting student       !
  !     project. The code was re-organized, modified and re-debugged by J. Qi and  !
  !     C. Chen when it is implemented into version 3.1.                           !
  !                                                                                !
  !  5) NetCDF Input and output. Version 3.1 uses the NetCDF input and output      !
  !     implemented into FVCOM by David Stuebe at SMAST/UMASSD. The code is also   !
  !     modified to improve the efficiency of inter-node data exchange and model   !
  !     output writing.                                                            !
  !                                                                                !
  !  6) Semi-implicit solver. In addition to the mode-split version, a semi-       !
  !     implicit solver was implemented into FVCOM by Z. Lai and C. Chen. Using    !
  !     this option requires the installation of PETCs.  With this implementation, !
  !     FVCOM can be run by choosing either mode-split scheme or semi-implicit     !
  !     scheme. The semi-implicit version of FVCOM can significantly increase the  !
  !     computational efficiency for regional and global scale application.        !
  !                                                                                !
  !  7) Nesting software module.  This module is built in FVCOM to allow multi-sub-!
  !     regional FVCOM domains to be run simultaneously through nesting approach.  !
  !     The module was originally written by P. Xue at SMAST/UMASSD when he and    !
  !     Chen applied FVCOM Kalman Filter to conduct The Observing System Simulation!
  !     Experiments (OSSEs) in Nantucket Sound. This module was modified by Dave   !
  !     Stuebe at SMAST/UMASSD and implemented into the version 3.0 of FVCOM, and  !
  !     graded by J. Qi and C. Chen at SMAST/UMASSD for current-wave interaction   !
  !     modules in in version 3.1.                                                 !
  !                                                                                !
  !  8) Dike-Groyne module. This module is implemented into FVCOM to deal with the !
  !     vertical wall under or above the sea surface. The algorithm was derived by ! 
  !     C. Chen at SMAST/UMASSD, and the parallelized code was written and tested  !
  !     by J. Qi, J. Ge and C. Chen at SMAST/UMASSD.                               !
  !                                                                                !
  !  9) SST/SSH assimilation module.  A SSH assimilation module is added to include!
  !     the satellite-altimeter data into FVCOM. The algorithm is the same as that !
  !     used for SST assimilation. The SST/SSH is merged together by Q. Xu, Z. Lai !
  !     and C. Chen at SMAST/UMASSD and implemented into version 3.1. C. Chen      !
  !     improved the SST algorithm to avoid the unreal thin layer caused by SST    !
  !     nudging.                                                                   !
  !                                                                                ! 
  !  10) Kalman Filter Assimilation modules. The Kalman Filter package is upgraded !
  !     by adding a Signular Evolutive Interpolated Kalman filter (SEIK). The      !
  !     source code of SEIK code is provided by Dr. Lars Nerger at AWI, Germany.   !
  !     In addition to Reduced Rank Kalman Filter (RRKF), Ensemble Kalman Filter   !
  !     (EnKF), we also reconstruct Ensemble Transform Kalman Filter (EnTKF) for   !
  !     the OSSEs application. The Kalman Filter assimilation modules were         !
  !     originally developed by a team effort led by C. Chen at SMAST/UMASSD and P.!
  !     Rizzoli/MIT. The team members include P. Xue, Z. Lai, Q. Xu and J. Qi at   !
  !     SMAST/UMASSD, J. Wei and S. Lyu at MIT and R. C. Beardsley at WHOI. The    !
  !     package is upgraded and implemented into version 3.1 by P. Xue, J.Qi and   !
  !     C.Chen. Dr. Lars Nerger owns the copyright of SEIK method used in FVCOM.   !                                   !
  !                                                                                !
  !  11) Optimal Interpolation Assimilation Module. This module is implemented into!
  !     FVCOM to integrate vertical profiles of hydrographic/current data into the !
  !     assimilation. The code was written by Q. Xu, J. Qi and C. Chen at          !
  !     SMAST/UMASSD.                                                              !
  !                                                                                ! 
  !  12) The pre-process programs to convert the input from the old version of     !
  !     FVCOM to version 3.0 or up.                                                !
  !                                                                                !
  !  13) Bugs in the wet/dry treatments on heat flux at the wet/dry boundary are   !
  !     corrected. The bug was detected by T. Kim and fixed by C. Chen and J. Qi at!
  !     SMAST/UMASSD. The second-order limiter for advection terms in the wet/dry  !
  !     treatment was originally coded in Fortran 77 version of FVCOM. This limiter!
  !     is added back to updated version 2.6 or up of FVCOM by C. Chen and J. Qi at!
  !     SMAST/UMASSD.                                                              !
  !                                                                                !
  !     Bugs in selecting General Turbulence Module is corrected by J. Qi and C.   !
  !  Chen at SMAST/UMASSD.                                                         !
  !                                                                                !
  !     Non-general definitions and commands in arrays, variables and NetCDF       !
  !  modules are changed by G. Cowles and Z. Lai to allow to use the version 3.1 of!
  !  FVCOM on IBM super computer. These problems only occur in version 3.0 or up.  !
  !  The bugs in PETCs were detected by Z. Lai, C. Chen and J. Qi at SMAST/UMASSD  !
  !  and corrected by the PETCs development team.                                  ! 
  !                                                                                ! 
  !    Bugs in the Lagrangian-tracking in the generalized terrain-following        !
  !  coordinates are corrected. Many minor bugs are also corrected.                !
  !                                                                                !
  !  Enjoy!      C. Chen at SMAST/UMASSD at Dec. 8 2009.                           !
  !                                                                                !
  !  Note for the grid input file in the FVCOM:                                    !
  !                                                                                !
  !  Some users used their own Matlab to create the grid input "*_grd.dat". If     !
  !  doing it, please check the way how nodes on each triangular cell are          !
  !  identified in the code. In FVCOM (see chapter 3 of users' manual), on each    !
  !  triangular cell, the three nodes are identified using integral numbers that   !
  !  are counted clockwise from 1 to 3. When the FVCOM was coded, it was tested    !
  !  using the grid that was created using the SMS software. Three nodes on each   !
  !  triangular cell in the SMS-generated grid are identified counter-clockwise,   !
  !  which are opposite to what are coded in FVCOM. For this reason, we coded the  !
  !  program in FVCOM to covert the SMS-generated grid. We have checked most of the!
  !  unstructured grid generation software, they all have the same definition as   !
  !  SMS. If one uses his/her own Matlab program to generate grid with the same    !
  !  definition of node number described in the users' manual, this conversion     !
  !  should be unnecessary.                                                        !
  !  We have added a program in FVCOM v4.3  or later to check the identification of!
  !  node numbers on each triangular cell when the grid file is input. This program!
  !  will stop the model run if the grid file does not meet the required format.   !
  !                                                                                !
  !                          C. Chen at SMAST/UMASSD at Apr. 23 2018.              !    
  !================================================================================!

  !==============================================================================!
  !  INCLUDE MODULES                                                             !
  !==============================================================================!

  USE MOD_UTILS
  USE CONTROL

  USE MOD_PAR  

  USE MOD_STARTUP

  USE MOD_TIME
  USE MOD_CLOCK

  USE MOD_INPUT
  USE MOD_NCDIO
  USE MOD_NCLL

  USE MOD_SETUP
  USE MOD_SET_TIME

  USE MOD_FORCE
  USE MOD_OBCS

  USE MOD_NESTING

  USE MOD_REPORT
  USE PROBES
  USE MOD_BOUNDSCHK !bounds checking
  USE MOD_DYE


  



! Added by researchers at Akvaplan-niva 2018, idealized tests give promising results.



! for periodic lateral boundary conditions






  USE MOD_STATION_TIMESERIES 
  USE MOD_SPARSE_TIMESERIES


  !------------------------------------------------------------------------------|
  IMPLICIT NONE

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string

  INTEGER :: IERR

  type(watch) Timer
  
  type(TIME) ::GET_BEGIN
  integer status

  !==============================================================================!
  ! INITIALIZE ALL CONTROL VARIABLES
  !==============================================================================!
  CALL INITIALIZE_CONTROL("FVCOM")

  ! INTIALIZE MPI CONTROL VARIABLES
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_COMM_FVCOM  = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
  MPI_FVCOM_GROUP = MPI_COMM_FVCOM ! FOR NOW MAKE THEM EQUAL
  MPI_IO_GROUP    = MPI_COMM_FVCOM ! FOR NOW MAKE THEM EQUAL

  !==============================================================================!
  !   INITIALIZE A STOP WATCH TIMER FOR TESTING SUBROUTINE EFFICENCY             !
  !==============================================================================!
  CALL WATCH_INIT(TIMER)

  !==============================================================================!
  !   IMPORT CASENAME AND COMMAND LINE ARGUMENTS AND START LOG FILE              !
  !==============================================================================!
  CALL COMMANDLINEIO(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)       
  if(DBG_SET(dbg_log)) Call WRITE_BANNER(PAR,NPROCS,MYID)

  !==============================================================================!
  ! SET DEFAULT VALUES AND READ NAME LISTS                                            
  !==============================================================================!

  CALL NAMELIST

  !==============================================================================!
  !   SET MODEL CONTROL PARAMTERS BASED ON NAME LIST HERE                        !
  !==============================================================================!
  CALL CNTRL_PRMTRS

  !==============================================================================!
  !   SET THE STARTUP TYPE TO BE USED!                                           !
  !==============================================================================!
  CALL SET_STARTUP_TYPE ! see: startup_type.F

  !==============================================================================!
  !   OPEN ALL FILES NEEDED BASED ON THE RUN PARAMETERS                          !
  !==============================================================================!
  CALL OPEN_ALL

  !==============================================================================!
  !   SET MODEL TIME BASED ON THE NAMELIST TIME STRINGS OR RESTART FILE          !
  !==============================================================================!
  CALL SETUP_TIME

  !==============================================================================!
  !   LOAD GRID CONNECTIVITY AND OBC LIST FOR METIS DECOMPOSITION                !
  !==============================================================================!
  CALL LOAD_GRID

  !==============================================================================!
  !   SETUP THE DOMAIN FOR PARALLEL OR SERIAL RUNNING                            !
  !==============================================================================!
  CALL SETUP_DOMAIN

  !==============================================================================!
  !   ALLOCATE ALL DOMAIN SIZE VARIABLES HERE                                    !
  !==============================================================================!
  CALL ALLOCATE_ALL
  CALL ALLOC_VARS_DYE
  !==============================================================================!
  !   LOAD/SETUP PHYSICAL QUANTITIES (CORIOLIS, GRAVITY, SPONGE LAYER, XY/LATLON)!
  !==============================================================================!
  CALL COORDS_N_CONST

  !==============================================================================!
  ! CALCULATE GRID METRICS - NEIGHBORS, GRADIENTS, CELL AREA, INTERP COEFF'S     !
  !==============================================================================!
  CALL GRID_METRICS

  !==============================================================================!
  ! SETUP THE SEDIMENT MODEL (MUST COME BEFORE SETUP_FORCING)                    ! 
  !==============================================================================!

   !SETUP TVD ADVECTION

!JQI  !==============================================================================!
!JQI  !  SETUP THE MODEL FORCING                                                     !
!JQI  !==============================================================================!
!JQI  CALL SETUP_FORCING

  !==============================================================================!
  !  GET THE PARAMETERS OF BIOLOGICAL MODEL                                      !
  !==============================================================================!

  !==============================================================================!
  !  SETUP THE MODEL FORCING                                                     !
  !==============================================================================!
  CALL SETUP_FORCING

  !==============================================================================!
  !  SETUP OTHER TOOLS, MODELS AND DATA ASSIMILATION                             !
  !==============================================================================!


  !  SETUP PETSc FOR SEMI_IMPLICIT AND NON-HYDROSTATIC MODULE





  ! SETUP DATA ASSIMILATION MODE
! New Open Boundary Condition ----2

  !==============================================================================!
  !  SET THE INITIAL CONDITIONS FOR THE MODEL                                    !
  !==============================================================================!
  CALL STARTUP

!qxu{ test function READ_DATETIME
!  GET_BEGIN = READ_DATETIME('2009-01-01T00:00:00.0','ymd','UTC',status)
!  CALL PRINT_TIME(GET_BEGIN,IPT,'2009-01-01T00:00:00.0')
!qxu}

  !==============================================================================!
  !  CALL ARCHIVE TO SETUP THE OUTPUT AND DUMP CONSTANT VALUES                   !
  !==============================================================================!
  CALL BCOND_GCN(8,0)
 
  
  CALL ARCHIVE
  ! ORDER MATTERS - ARCHIVE_NEST MUST GO AFTER ARCHIVE DURING SETUP!
  CALL ARCHIVE_NEST

  CALL SET_PROBES(PROBES_ON,PROBES_NUMBER,PROBES_FILE)
  IF(OUT_STATION_TIMESERIES_ON)CALL READ_STATION_FILE

  ! Setup Bounds checking (shutdown if variables exceed threshold)
  CALL SETUP_BOUNDSCHK !bounds checking


  
  IF(OUT_STATION_TIMESERIES_ON)THEN
    CALL GET_OUTPUT_FILE_INTERVAL(TRIM(OUT_INTERVAL),INTERVAL_TIME_SERIES)
    CALL OUT_STATION_TIMESERIES
    TIME_SERIES = STARTTIME + INTERVAL_TIME_SERIES
  END IF
  !==============================================================================!
  !  SELECT THE RUN MODE AND EXECUTE THE MAIN LOOP
  !==============================================================================!
  SELECT CASE(FVCOM_RUN_MODE)
     ! RUN MODE SET IN mod_assim.F(set_assim_param)

     ! =============================================================================!
     ! == PURE SIMULATION MODE - Instantanious data assimilation only ==============!
  CASE(FVCOM_PURE_SIM)
     ! =============================================================================!


     !==============================================================================!
     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
     !==============================================================================!
     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
        write(IPT,*) "!============== STARTING MAIN LOOP AT:==============="
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        write(IPT,*) "!===================================================="
     end if

     CALL REPORT('INITIAL CONDITIONS')

     if(DBG_SET(dbg_log)) THEN
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
        write(IPT,*) "!===================================================="
     end if

     !////////////////////////// MAIN LOOP //////////////////////////////////////////
     DO IINT=ISTART,IEND

        IntTime=IntTime + IMDTI

        CALL INTERNAL_STEP

        !==============================================================================!
        !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
        !==============================================================================!
        if(DBG_SET(dbg_log)) &
             & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
        
        IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

        !==============================================================================!
        !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
        !==============================================================================!

!for Eulerian velocity output

        CALL ARCHIVE


        CALL DUMP_PROBE_DATA
        IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
        !==============================================================================!
        !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
        !==============================================================================!
        CALL SHUTDOWN_CHECK(D1)

        !==============================================================================!
        !  CALL BOUNDS CHECK TO SEE IF VARIABLES EXCEED USER-DEFINED THRESHOLDS 
        !==============================================================================!
        CALL BOUNDSCHK  !bounds checking

        !==============================================================================!
        !    LAGRANGIAN PARTICLE TRACKING                                              |
        !==============================================================================!

        !==============================================================================!
        !    NESTING OUTPUT                                                            |
        !==============================================================================!
        IF(NCNEST_ON)      CALL ARCHIVE_NEST

     END DO
     !////////////////////////// END MAIN LOOP //////////////////////////////////////

  ! ================================================================================!
  ! == THIS MAIN LOOP USES NUDGING OR OI METHODS TO ASSIMILATE =====================!
  CASE(FVCOM_NUDGE_OI_ASSIM)
  ! ================================================================================!

!  ! ================================================================================!
!  ! == THIS MAIN LOOP USES OPTIMAL INTERPOLATION  METHODS TO ASSIMILATE ============!
!  CASE(FVCOM_OI_ASSIM)
!  ! ================================================================================!
!
!     CALL SET_OIASSIM_INTERVALS
!     CALL ALLOC_BUFFER_OI
!
!     ! SET THE ASSIMILATION/SIMULATION RESET TIME
!     ASSIM_RESTART_TIME = StartTime
!
!     ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
!     INT_START = ISTART
!     INT_COUNT = ISTART
!     IF(SST_OIASSIM.OR.SSTGRD_OIASSIM)THEN
!       INT_END = ISTART + 2*(IEND-ISTART)
!     ELSEIF(TS_OIASSIM.OR.CUR_OIASSIM)THEN
!       INT_END = IEND
!     END IF
!     !==============================================================================!
!     !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
!     !==============================================================================!
!     if(DBG_SET(dbg_log)) THEN
!        write(IPT,*) "===================================================="
!        write(IPT,*) "===================================================="
!        write(IPT,*) "======STARTING MAIN LOOP OIASSIMILATION MODE:======="
!        Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)
!        write(IPT,*) "===================================================="
!     end if
!
!     CALL REPORT('INITIAL CONDITIONS')
!     
!     if(DBG_SET(dbg_log)) THEN
!        write(IPT,*) "===================================================="
!        write(IPT,*) "===================================================="
!        write(IPT,*) "===================================================="
!     end if

!     DO WHILE(IntTime < EndTime)
!        IF(SST_OIASSIM.OR.SSTGRD_OIASSIM)THEN
!          ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY
!        ELSEIF(TS_OIASSIM.OR.CUR_OIASSIM)THEN
!          ASSIM_RESTART_TIME = EndTime
!        END IF
!        SST_SAVE_INDEX = 0
!        SST_SAVE_TIME  = IntTime + sst_save_interval
!
!        CALL OI_SAVE_STATE
!
!
!        if(DBG_SET(dbg_log)) THEN
!           Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)
!           write(IPT,*) "======= Start 1 Day Simulation  ===================="
!        end if
!        !==============================================================================!
!        !    RUN PURE SIMULATION MODE:
!        !==============================================================================!
!
!        OIASSIM_FLAG = 0
!        DO WHILE(IntTime < ASSIM_RESTART_TIME)
!
!           IINT = IINT + 1
!           IntTime=IntTime + IMDTI
!
!           CALL INTERNAL_STEP
!
!           !==============================================================================!
!           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
!           !==============================================================================!
!           if(DBG_SET(dbg_log)) &
!                & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)
!
!           IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')
!
!           !==============================================================================!
!           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
!           !==============================================================================!
!           IF(TS_OIASSIM.OR.CUR_OIASSIM)THEN
!
!             CALL ARCHIVE
!
!             CALL DUMP_PROBE_DATA
!
!           END IF
!
!           CALL SHUTDOWN_CHECK(D1)
!
!        END DO
!
!        IF(SST_OIASSIM.OR.SSTGRD_OIASSIM)THEN
!        !==============================================================================!
!        !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
!        !==============================================================================!
!          CALL OI_RESTORE_STATE
!
!          if(DBG_SET(dbg_log)) THEN
!             Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
!             write(IPT,*) "======= Start 1 Day Assimilation  ===================="
!          end if
!
!        !==============================================================================!
!        !    RUN DATA ASSIMILATION MODE:
!        !==============================================================================!
!          OIASSIM_FLAG = 1
!          DO WHILE(IntTime < ASSIM_RESTART_TIME)
!
!             IINT = IINT + 1
!             IntTime=IntTime + IMDTI
!
!             CALL INTERNAL_STEP
!
!           !==============================================================================!
!           !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
!           !==============================================================================!
!             if(DBG_SET(dbg_log)) &
!                  & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)
!
!             IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')
!
!           !==============================================================================!
!           !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
!           !==============================================================================!
!             CALL ARCHIVE
!
!             CALL DUMP_PROBE_DATA
!           !==============================================================================!
!           !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
!           !==============================================================================!
!             CALL SHUTDOWN_CHECK(D1)
!
!
!           !==============================================================================!
!           !    LAGRANGIAN PARTICLE TRACKING                                              |
!           !==============================================================================!
!# if defined (LAG_PARTICLE)
!             CALL LAG_UPDATE
!# endif
!
!          END DO
!        END IF
!        ! RUN THE NEXT DAY
!
!     END DO
    
  CASE(FVCOM_RRKF_WITHOUT_SSA)
  CASE(FVCOM_RRKF_WITH_SSA)
    
!============================================END RRKF_WITH_SSA================================| 





  CASE(FVCOM_ENKF_WITHOUT_SSA)
    
!============================================END ENKF_WITHOUT_SSA================================|     

  CASE(FVCOM_ENKF_WITH_SSA)
    
!============================================END ENKF_WITH_SSA================================|     
  CASE DEFAULT
     CALL FATAL_ERROR("UNKNOWN FVCOM_RUN_MODE :'"//TRIM(FVCOM_RUN_MODE),&
          & "Options are the following: '"//TRIM(FVCOM_PURE_SIM)//"&
          &' OR '"//TRIM(FVCOM_NUDGE_OI_ASSIM) )
  END SELECT





  
  if(DBG_SET(dbg_log)) write(IPT,*)"TADA!"
  CALL PSHUTDOWN

END PROGRAM FVCOM
