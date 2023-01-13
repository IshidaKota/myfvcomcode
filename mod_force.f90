










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

MODULE MOD_FORCE
  USE ALL_VARS
  USE MOD_INTERP
  USE BCS
  USE MOD_TIME
  USE MOD_NCTOOLS
  USE MOD_NCLL
  USE MOD_UTILS
  USE MOD_SPHERICAL
  USE MOD_PAR
  USE MOD_INPUT
  USE MOD_HEATFLUX

  IMPLICIT NONE

  SAVE
  PRIVATE

  ! COMMON FILE TYPE SHARED BY SEVERAL TYPES OF FORCING
!  CHARACTER(LEN=80),PUBLIC, PARAMETER :: WRF2FVCOM_SOURCE = &
!       & "wrf2fvcom version 0.14 (2007-07-19) (Bulk method: COARE 2.6Z)" 

  CHARACTER(LEN=80),PUBLIC, PARAMETER :: WRF2FVCOM_SOURCE = &
       & "wrf2fvcom version"

  CHARACTER(LEN=80),PUBLIC, PARAMETER :: fvcom_grid_SOURCE = &
       & "fvcom grid (unstructured) surface forcing" 

  CHARACTER(LEN=80),PUBLIC, PARAMETER :: fvcom_cap_grid_SOURCE = &
       & "FVCOM grid (unstructured) surface forcing" 

  CHARACTER(LEN=80),PUBLIC, PARAMETER :: wrf_grid_SOURCE = &
       & "wrf grid (structured) surface forcing" 

  CHARACTER(LEN=80),PUBLIC, PARAMETER :: surf_forcing_pt_SOURCE = &
       & "single-point time-dependent surface forcing"


  ! TIDAL FORCING VARIABLES FOR UPDATE AND SETUP
  INTEGER, PUBLIC :: TIDE_FORCING_TYPE
  INTEGER, PARAMETER, PUBLIC :: TIDE_FORCING_SPECTRAL = 1
  INTEGER, PARAMETER, PUBLIC :: TIDE_FORCING_TIMESERIES = 2
  TYPE(NCFILE), POINTER :: TIDE_FILE
  TYPE(NCVAR), POINTER ::  TIDE_ELV_N, TIDE_ELV_P
  CHARACTER(LEN=Char_max_attlen), PUBLIC,ALLOCATABLE :: TIDE_FORCING_COMMENTS(:)


  ! RIVER FORCING VARIABLES FOR UPDATE AND SETUP
  ! NOTE RIVERS ARE A PAIN - EACH PROCESSOR HAS TO FETCH ITS OWN DATA!
  CHARACTER(LEN=Char_max_attlen), PUBLIC,ALLOCATABLE :: RIVER_FORCING_COMMENTS(:)
  TYPE A_RIVER_FILE
     TYPE(NCFILE), POINTER :: NCF
     INTEGER RIVERS_IN_FILE
     
     TYPE(TIME) :: RIVER_PERIOD
     INTEGER, ALLOCATABLE :: RIV_FILE2LOC(:)
     ! USAGE :   DO I = 1, RIVERS_IN_FILE
     !              J = RIV_FILE2LOC(I)
     !              IF (J/=0) QDIS(J) = FILE_DIS(I)   
     TYPE(NCVAR), POINTER :: FLUX_N, FLUX_P
     TYPE(NCVAR), POINTER :: TEMP_N, TEMP_P
     TYPE(NCVAR), POINTER :: SALT_N, SALT_P
     ! ADD MORE HERE!
  END TYPE A_RIVER_FILE
  TYPE(A_RIVER_FILE), ALLOCATABLE :: RIVER_FORCING(:)



  ! =================================================================
  ! GROUND WATER FORCING VARIABLES FOR UPDATE AND SETUP
  TYPE(NCFILE), POINTER :: GWATER_FILE

  INTEGER :: GWATER_FORCING_TYPE
  INTEGER,PARAMETER :: GWATER_IS_XXX     = 0
  INTEGER,PARAMETER :: GWATER_IS_FVCOMGRID   = 1

  INTEGER :: GWATER_UNITS
  INTEGER, PARAMETER :: GWATER_M3S_1=1
  INTEGER, PARAMETER :: GWATER_MS_1=2

  TYPE(TIME) :: GWATER_PERIOD

  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE :: GWATER_FORCING_COMMENTS(:)
  TYPE(INTERP_WEIGHTS),POINTER :: GWATER_INTP_N
  TYPE(INTERP_WEIGHTS),POINTER :: GWATER_INTP_C
  TYPE(NCVAR), POINTER :: GWATER_FLUX_N, GWATER_FLUX_P    ! Discharge
  TYPE(NCVAR), POINTER :: GWATER_TEMP_N, GWATER_TEMP_P    ! Temperature
  TYPE(NCVAR), POINTER :: GWATER_SALT_N, GWATER_SALT_P    ! Salinity

  ! =================================================================
  ! OPEN BOUNDARY CONDITION SALINITY VARIABLES FOR UPDATE AND SETUP
  TYPE(NCFILE), POINTER :: OBC_S_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC :: OBC_S_COMMENTS
  INTEGER :: OBC_S_TYPE
  INTEGER,PARAMETER :: OBC_S_SIGMA = 1
  TYPE(NCVAR), POINTER :: OBC_S_N, OBC_S_P  

  ! =================================================================
  ! OPEN BOUNDARY CONDITION TEMPERATURE VARIABLES FOR UPDATE AND SETUP
  TYPE(NCFILE), POINTER :: OBC_T_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC :: OBC_T_COMMENTS
  INTEGER :: OBC_T_TYPE
  INTEGER,PARAMETER :: OBC_T_SIGMA = 1
  TYPE(NCVAR), POINTER :: OBC_T_N, OBC_T_P   

  ! =================================================================

  ! =================================================================
  ! =================================================================

  ! =================================================================
  ! SURFACE HEAT FORCING FILE DATA 
  INTEGER :: HEAT_FORCING_TYPE

  INTEGER,PARAMETER :: HEAT_IS_WRFGRID     = 0
  INTEGER,PARAMETER :: HEAT_IS_FVCOMGRID   = 1

  TYPE(TIME) :: HEAT_PERIOD

  TYPE(NCFILE), POINTER :: HEAT_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE :: HEAT_FORCING_COMMENTS(:)
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE :: HEAT_CALCULATE_COMMENTS(:)
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE :: HEAT_SOLAR_COMMENTS(:)
  TYPE(INTERP_WEIGHTS),POINTER :: HEAT_INTP_N
  TYPE(INTERP_WEIGHTS),POINTER :: HEAT_INTP_C
  TYPE(NCVAR), POINTER :: HEAT_SWV_N, HEAT_SWV_P   ! SHORT WAVE
  !  TYPE(NCVAR), POINTER :: HEAT_LWV_N, HEAT_LWV_P   ! LONG WAVE
  !  TYPE(NCVAR), POINTER :: HEAT_LTNT_N, HEAT_LTNT_P ! LATENT
  !  TYPE(NCVAR), POINTER :: HEAT_SNS_N, HEAT_SNS_P   ! SENSIBLE
  TYPE(NCVAR), POINTER :: HEAT_NET_N, HEAT_NET_P   ! NET HEAT FLUX
  TYPE(NCVAR), POINTER :: T_AIR_N, T_AIR_P   
  TYPE(NCVAR), POINTER :: RH_AIR_N, RH_AIR_P   
  TYPE(NCVAR), POINTER :: PA_AIR_N, PA_AIR_P   
  TYPE(NCVAR), POINTER :: DLW_AIR_N, DLW_AIR_P   
  TYPE(NCVAR), POINTER :: DSW_AIR_N, DSW_AIR_P   
  ! =================================================================
  ! SURFACE WIND STRESS FILE DATA 
  INTEGER :: WINDS_FORCING_TYPE

  INTEGER, PARAMETER :: WINDS_ARE_WRFGRID   = 0
  INTEGER, PARAMETER :: WINDS_ARE_FVCOMGRID = 1
  INTEGER, PARAMETER :: WINDS_ARE_PT_SOURCE = 2
  

  TYPE(TIME) :: WINDS_PERIOD

  TYPE(NCFILE), POINTER :: WINDS_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE    :: WINDS_FORCING_COMMENTS(:)
  TYPE(INTERP_WEIGHTS),POINTER :: WINDS_INTP_N
  TYPE(INTERP_WEIGHTS),POINTER :: WINDS_INTP_C
  TYPE(NCVAR), POINTER :: WINDS_STRX_N, WINDS_STRX_P   ! STRESS IN THE X DIRECTION
  TYPE(NCVAR), POINTER :: WINDS_STRY_N, WINDS_STRY_P   ! STRESS IN THE Y DIRECTION








!Jadon
  ! =================================================================
  ! SURFACE WAVE STRESS FILE DATA 
  INTEGER :: WAVES_FORCING_TYPE

  INTEGER, PARAMETER :: WAVES_ARE_WRFGRID   = 0
  INTEGER, PARAMETER :: WAVES_ARE_FVCOMGRID = 1

  TYPE(TIME) :: WAVES_PERIOD

  TYPE(NCFILE), POINTER :: WAVES_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE    :: WAVES_FORCING_COMMENTS(:)
  TYPE(INTERP_WEIGHTS),POINTER :: WAVES_INTP_N
  TYPE(INTERP_WEIGHTS),POINTER :: WAVES_INTP_C
  TYPE(NCVAR), POINTER :: WAVES_HEIGHT_N,    WAVES_HEIGHT_P    ! WAVE HEIGHT
  TYPE(NCVAR), POINTER :: WAVES_LENGTH_N,    WAVES_LENGTH_P    ! WAVE LENGTH
  TYPE(NCVAR), POINTER :: WAVES_DIRECTION_N, WAVES_DIRECTION_P ! WAVE DIRECTION
  TYPE(NCVAR), POINTER :: WAVES_PERIOD_N,    WAVES_PERIOD_P    ! WAVE PERIOD
  TYPE(NCVAR), POINTER :: WAVES_PER_BOT_N,   WAVES_PER_BOT_P   ! BOTTOM PERIOD
  TYPE(NCVAR), POINTER :: WAVES_UB_BOT_N,    WAVES_UB_BOT_P    ! BOTTOM VELOCITY
  ! =================================================================
  ! SURFACE PRECIPTATION DATA
  INTEGER :: PRECIP_FORCING_TYPE

  INTEGER,PARAMETER :: PRECIP_IS_WRFGRID     = 0
  INTEGER,PARAMETER :: PRECIP_IS_FVCOMGRID      = 1

  TYPE(TIME) :: PRECIP_PERIOD

  TYPE(NCFILE), POINTER :: PRECIP_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE:: PRECIP_FORCING_COMMENTS(:)
  TYPE(INTERP_WEIGHTS),POINTER :: PRECIP_INTP_N
  TYPE(INTERP_WEIGHTS),POINTER :: PRECIP_INTP_C
  TYPE(NCVAR), POINTER :: PRECIP_PRE_N, PRECIP_PRE_P   ! PRECIPITATION
  TYPE(NCVAR), POINTER :: PRECIP_EVP_N, PRECIP_EVP_P   ! EVAPORATION
  !<------ishid debug 20230104
  ! IF (PRECIPITATION_BULK_ON) THEN
   !IF
   !  TYPE(NCVAR), POINTER :: PRECIP_UWND_N, PRECIP_UWND_P   ! PRECIPITATION
   !  TYPE(NCVAR), POINTER :: PRECIP_VWND_N, PRECIP_VWND_P   ! EVAPORATION
     TYPE(NCVAR), POINTER :: PRECIP_PAIR_N, PRECIP_PAIR_P   ! PRECIPITATION
     TYPE(NCVAR), POINTER :: PRECIP_RH_N, PRECIP_RH_P   ! relative humidity
  ! ENDIF
  !>------ishid debug 20230104
  ! =================================================================
  ! AIR PRESSURE FILE DATA 
  INTEGER :: AIRPRESSURE_FORCING_TYPE

  INTEGER, PARAMETER :: AIRPRESSURE_IS_WRFGRID   = 0
  INTEGER, PARAMETER :: AIRPRESSURE_IS_FVCOMGRID = 1

  TYPE(TIME) :: AIRPRESSURE_PERIOD

  TYPE(NCFILE), POINTER :: AIRPRESSURE_P_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE    :: AIRPRESSURE_FORCING_COMMENTS(:)
  TYPE(INTERP_WEIGHTS),POINTER :: AIRPRESSURE_INTP_N
  TYPE(INTERP_WEIGHTS),POINTER :: AIRPRESSURE_INTP_C
  TYPE(NCVAR), POINTER :: AIR_PRESSURE_N, AIR_PRESSURE_P   


  ! =================================================================
  ! ICE MODEL DATA
  INTEGER :: ICE_FORCING_TYPE

  INTEGER,PARAMETER :: ICE_IS_WRFGRID     = 0
  INTEGER,PARAMETER :: ICE_IS_FVCOMGRID      = 1

  TYPE(TIME) :: ICE_PERIOD

  TYPE(NCFILE), POINTER :: ICE_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE:: ICE_FORCING_COMMENTS(:)
  TYPE(INTERP_WEIGHTS),POINTER :: ICE_INTP_N
  TYPE(INTERP_WEIGHTS),POINTER :: ICE_INTP_C
  TYPE(NCVAR), POINTER :: ICE_SWV_N, ICE_SWV_P   ! SHORT WAVE
  TYPE(NCVAR), POINTER :: ICE_SAT_N, ICE_SAT_P   ! SEA LEVEL AIR TEMPERATURE
  TYPE(NCVAR), POINTER :: ICE_SPQ_N, ICE_SPQ_P   ! SPECFIC HUMIDITY
  TYPE(NCVAR), POINTER :: ICE_CLD_N, ICE_CLD_P   ! CLOUD COVER

  ! =================================================================
  ! ICING MODEL DATA
  INTEGER :: ICING_FORCING_TYPE

  INTEGER,PARAMETER :: ICING_IS_WRFGRID     = 0
  INTEGER,PARAMETER :: ICING_IS_FVCOMGRID      = 1

  TYPE(TIME) :: ICING_PERIOD

  TYPE(NCFILE), POINTER :: ICING_FILE
  CHARACTER(LEN=Char_max_attlen), PUBLIC, ALLOCATABLE:: ICING_FORCING_COMMENTS(:)
  TYPE(INTERP_WEIGHTS),POINTER :: ICING_INTP_N
  TYPE(INTERP_WEIGHTS),POINTER :: ICING_INTP_C
  TYPE(NCVAR), POINTER :: ICING_SAT_N, ICING_SAT_P   ! SEA LEVEL AIR PRESSURE
  TYPE(NCVAR), POINTER :: ICING_WSPX_N, ICING_WSPX_P   ! SEA LEVEL AIR TEMPERATURE
  TYPE(NCVAR), POINTER :: ICING_WSPY_N, ICING_WSPY_P   ! SPECFIC HUMIDITY

  PUBLIC :: SETUP_FORCING
  PUBLIC :: UPDATE_GROUNDWATER
  PUBLIC :: UPDATE_WIND
  PUBLIC :: UPDATE_WAVE
  
  PUBLIC :: UPDATE_HEAT_CALCULATED



  PUBLIC :: UPDATE_PRECIPITATION
  PUBLIC :: UPDATE_AIRPRESSURE
  PUBLIC :: UPDATE_TIDE
  PUBLIC :: UPDATE_RIVERS
  PUBLIC :: UPDATE_OBC_TEMP
  PUBLIC :: UPDATE_OBC_SALT
  PUBLIC :: UPDATE_ICE
  PUBLIC :: UPDATE_ICING

CONTAINS

  SUBROUTINE SETUP_FORCING
    IMPLICIT NONE
    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START SETUP_FORCING"

    IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*  )'!'
       WRITE(IPT,*  )'!           SETTING UP PRESCRIBED BOUNDARY CONDITIONS   '
       WRITE(IPT,*  )'!'
    END IF


    ! NULLIFY EVERYTHING
    NULLIFY(TIDE_FILE, TIDE_ELV_N, TIDE_ELV_P)

    NULLIFY(GWATER_FILE)

    NULLIFY(HEAT_FILE,HEAT_INTP_N, HEAT_INTP_C, T_AIR_P, T_AIR_N,  &
         &  RH_AIR_P, RH_AIR_N, PA_AIR_P, PA_AIR_N, DLW_AIR_P,     &
	 &  DLW_AIR_N, DSW_AIR_P, DSW_AIR_N)
    !<----------ishid debug 20230105
    !NULLIFY(PRECIP_FILE,PRECIP_PRE_P,PRECIP_PRE_N,PRECIP_UWND_N,PRECIP_UWND_P,PRECIP_VWND_N,PRECIP_VWND_P, &
    !      & PRECIP_PAIR_N,PRECIP_PAIR_N,PRECIP_RH_N,PRECIP_RH_P)
     NULLIFY(PRECIP_FILE,PRECIP_PAIR_N,PRECIP_PAIR_N,PRECIP_RH_N,PRECIP_RH_P)
    !>---------ishid debug 20230105
    NULLIFY(WINDS_FILE,WINDS_INTP_N,WINDS_INTP_C, WINDS_STRX_N,&
         & WINDS_STRX_P, WINDS_STRY_N, WINDS_STRY_P)

    NULLIFY(AIRPRESSURE_P_FILE,AIRPRESSURE_INTP_N,AIRPRESSURE_INTP_C, AIR_PRESSURE_N,&
         & AIR_PRESSURE_P)

    NULLIFY(WAVES_FILE,WAVES_INTP_N,WAVES_INTP_C, &
         & WAVES_HEIGHT_N,    WAVES_HEIGHT_P,     &
         & WAVES_LENGTH_N,    WAVES_LENGTH_P,     &
         & WAVES_DIRECTION_N, WAVES_DIRECTION_P,  & 
         & WAVES_PERIOD_N,    WAVES_PERIOD_P,     &
         & WAVES_PER_BOT_N,   WAVES_PER_BOT_P,    &
         & WAVES_UB_BOT_N,    WAVES_UB_BOT_P        )





    CALL TIDAL_ELEVATION
    CALL OBC_TEMPERATURE
    CALL OBC_SALINITY
    CALL RIVER_DISCHARGE
    CALL SURFACE_WINDSTRESS
    CALL SURFACE_PRECIPITATION
    CALL SURFACE_AIRPRESSURE
    CALL SURFACE_HEATING_CALCULATED




    CALL GROUND_WATER

    ! ORDER IS IMPORTANT! ICE AND ICING MAY USE THE SAME POINTERS TO
    ! REFERENCE FILE! THEY MUST BE SET UP LAST.
    CALL ICING_FORCING

    CALL ICE_MODEL_FORCING

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END SETUP_FORCING"
  END SUBROUTINE SETUP_FORCING
  !================================================================
  !================================================================
  SUBROUTINE GROUND_WATER
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND

    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    CHARACTER(len=60) :: tempstrng, flowstrng, saltstrng
    TYPE(TIME) :: TIMETEST

    INTEGER :: LATS, LONS, I, Ntimes

    INTEGER :: STATUS

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START GROUND_WATER"

    IF (.NOT. GROUNDWATER_ON ) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "!  GROUND WATER FORCING IS OFF!"
       ALLOCATE(GWATER_FORCING_COMMENTS(1))
       GWATER_FORCING_COMMENTS="GROUND WATER FORCING IS OFF!"
       RETURN
    END IF



    ! DETERMINE HOW TO LOAD THE DATA
    SELECT CASE(GROUNDWATER_KIND)
    CASE (CNSTNT)
       
       write(flowstrng,'(f8.4)') groundwater_flow
       write(tempstrng,'(f8.4)') groundwater_temp
       write(saltstrng,'(f8.4)') groundwater_salt
       
       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*)"! SETTING UP CONSTANT GROUNDWATER FORCING: "
          WRITE(IPT,*)"      Flow Rate: "//trim(flowstrng)
          WRITE(IPT,*)"           Temp: "//trim(tempstrng)
          WRITE(IPT,*)"           Salt: "//trim(saltstrng)
       END IF
       
       ALLOCATE(GWATER_FORCING_COMMENTS(4))
       GWATER_FORCING_COMMENTS(1) = "Using constant groundwater forcing from run file:"
       GWATER_FORCING_COMMENTS(2) = "Flow Rate:"//trim(flowstrng)
       IF(GROUNDWATER_TEMP_ON) THEN
          GWATER_FORCING_COMMENTS(3) = "Temperature (specified):"//trim(tempstrng)
       ELSE
          GWATER_FORCING_COMMENTS(3) = "Temperature is calculated"
       END IF

       IF(GROUNDWATER_SALT_ON) THEN
          GWATER_FORCING_COMMENTS(4) = "Salinity (specified):"//trim(saltstrng)
       ELSE
          GWATER_FORCING_COMMENTS(4) = "Salinity is calculated"
       END IF
       RETURN
       
    CASE(STTC)
       
       CALL FATAL_ERROR("STATIC GROUNDWATER Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDENT GROUNDWATER Not Set Up Yet")

    CASE(PRDC)
       
       GWATER_FILE => FIND_FILE(FILEHEAD,trim(GROUNDWATER_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE
       ATT => FIND_ATT(GWATER_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(GWATER_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")

       IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          GWATER_FORCING_TYPE = GWATER_IS_FVCOMGRID

       ELSEIF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          GWATER_FORCING_TYPE = GWATER_IS_FVCOMGRID

       ELSE
          CALL PRINT_FILE(GWATER_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE GROUNDWATER FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(GWATER_FORCING_COMMENTS(5))
       GWATER_FORCING_COMMENTS(1) = "FVCOM periodic GroundWater forcing:"
       GWATER_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(GroundWater_FILE)
       GWATER_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))
       IF(GROUNDWATER_TEMP_ON) THEN
          GWATER_FORCING_COMMENTS(4) = "Temperature is specified"
       ELSE
          GWATER_FORCING_COMMENTS(4) = "Temperature is calculated"
       END IF

       IF(GROUNDWATER_SALT_ON) THEN
          GWATER_FORCING_COMMENTS(5) = "Salinity is specified"
       ELSE
          GWATER_FORCING_COMMENTS(5) = "Salinity is calculated"
       END IF

       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(GWATER_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       GWATER_PERIOD = get_file_time(GWATER_FILE,ntimes)


       IF (ZEROTIME /= get_file_time(GWATER_FILE,1)) THEN

          CALL PRINT_REAL_TIME(get_file_time(GWATER_FILE,1),IPT,"FIRST FILE TIME",TIMEZONE)
          CALL PRINT_REAL_TIME(get_file_time(GWATER_FILE,ntimes),IPT,"LAST FILE TIME",TIMEZONE)

          CALL FATAL_ERROR&
               &("TO USE PERIODIC FORCING THE FILE TIME MUST COUNT FROM ZERO",&
               & "THE DIFFERENCE BETWEEN THE CURRENT MODEL TIME AND THE START TIME,",&
               & "MODULO THE FORCING PERIOD, DETERMINES THE CURRENT FORCING")
       END IF

       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! USING PERIODIC GroundWater FORCING:"
          CALL PRINT_TIME(GWATER_PERIOD,IPT,"PERIOD")
       END IF

    CASE(VRBL)
       
       GWATER_FILE => FIND_FILE(FILEHEAD,trim(GROUNDWATER_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE
       ATT => FIND_ATT(GWATER_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(GWATER_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")

       IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          GWATER_FORCING_TYPE = GWATER_IS_FVCOMGRID

       ELSEIF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          GWATER_FORCING_TYPE = GWATER_IS_FVCOMGRID

       ELSE
          CALL PRINT_FILE(GWATER_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE GROUNDWATER FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(GWATER_FORCING_COMMENTS(5))
       GWATER_FORCING_COMMENTS(1) = "FVCOM variable GroundWater forcing:"
       GWATER_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(GroundWater_FILE)
       GWATER_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       IF(GROUNDWATER_TEMP_ON) THEN
          GWATER_FORCING_COMMENTS(4) = "Temperature is specified"
       ELSE
          GWATER_FORCING_COMMENTS(4) = "Temperature is calculated"
       END IF

       IF(GROUNDWATER_SALT_ON) THEN
          GWATER_FORCING_COMMENTS(5) = "Salinity is specified"
       ELSE
          GWATER_FORCING_COMMENTS(5) = "Salinity is calculated"
       END IF

       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(GWATER_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")

       NTIMES = DIM%DIM
       
       ! CHECK THE FILE TIME AND COMPARE WITH MODEL RUN TIME
       TIMETEST = get_file_time(GWATER_FILE,1)
       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN THE GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE FORCING TIME SERIES")

       TIMETEST = get_file_time(GWATER_FILE,ntimes)
       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN THE GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            &"THE MODEL RUN ENDS AFTER THE FORCING TIME SERIES")

    CASE DEFAULT
       CALL FATAL_ERROR("GROUND_WATER: UNKNOWN GROUND WATER KIND?")

    END SELECT


    !==================================================================
    SELECT CASE(GWATER_FORCING_TYPE)
    !==================================================================
    CASE(GWATER_IS_FVCOMGRID)
    !==================================================================

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP GROUND WATER FORCING FROM A 'fvcom grid' FILE"

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_DIM(GWATER_FILE,'node',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN THE GROUND WATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            & "COULD NOT FIND DIMENSION 'node'")

       if (mgl /= dim%dim) CALL FATAL_ERROR&
            &("GROUNDWATER: the number of nodes in the file does not match the fvcom grid?")


       DIM => FIND_DIM(GWATER_FILE,'nele',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN THE GROUND WATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            & "COULD NOT FIND DIMENSION 'nele'")

       if (ngl /= dim%dim) CALL FATAL_ERROR&
            &("GROUNDWATER: the number of elements in the file does not match the fvcom grid?")

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       ! GROUND WATER FLUX DATA
       VAR => FIND_VAR(GWATER_FILE,"groundwater_flux",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN THE GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            & "COULD NOT FIND VARIABLE 'groundwater_flux'")

       ATT => FIND_ATT(VAR,"units",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN THE GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            & "COULD NOT FIND THE UNITS FOR THE VARIABLE 'groundwater_flux'")
       
       IF (ATT%CHR(1)(1:len_trim("m3 s-1")) == "m3 s-1") THEN
          GWATER_UNITS = GWATER_M3S_1
       ELSEIF (ATT%CHR(1)(1:len_trim("m s-1")) == "m s-1") THEN
          GWATER_UNITS = GWATER_MS_1
       ELSE
          CALL FATAL_ERROR &
            & ("IN THE GROUNDWATER FILE OBJECT",&
            & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
            & "UNKNOWN UNITS FOR THE VARIABLE 'groundwater_flux'")
       END IF

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       GWATER_FLUX_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN GROUNDWATER")
       CALL NC_CONNECT_PVAR(GWATER_FLUX_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       GWATER_FLUX_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN GROUNDWATER")
       CALL NC_CONNECT_PVAR(GWATER_FLUX_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! GROUNDWATER INFLOW TEMPERATURE
       IF(GROUNDWATER_TEMP_ON)THEN
          VAR => FIND_VAR(GWATER_FILE,"groundwater_temp",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN THE GROUNDWATER FILE OBJECT",&
               & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
               & "COULD NOT FIND VARIABLE 'groundwater_temp'")
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          GWATER_TEMP_N => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN GROUNDWATER")
          CALL NC_CONNECT_PVAR(GWATER_TEMP_N,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
          
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          GWATER_TEMP_P => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN GROUNDWATER")
          CALL NC_CONNECT_PVAR(GWATER_TEMP_P,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
       END IF

       ! GROUNDWATER INFLOW SALINITY
       IF(GROUNDWATER_SALT_ON)THEN
          VAR => FIND_VAR(GWATER_FILE,"groundwater_salt",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN THE GROUNDWATER FILE OBJECT",&
               & "FILE NAME: "//TRIM(GROUNDWATER_FILE),&
               & "COULD NOT FIND VARIABLE 'groundwater_salt'")
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          GWATER_SALT_N => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN GROUNDWATER")
          CALL NC_CONNECT_PVAR(GWATER_SALT_N,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
          
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          GWATER_SALT_P => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN GROUNDWATER")
          CALL NC_CONNECT_PVAR(GWATER_SALT_P,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
       END IF


    !==================================================================
    CASE DEFAULT
    !==================================================================
       CALL FATAL_ERROR("CAN NOT RECOGNIZE GROUND WATER FILE TYPE!")
    !==================================================================
    END SELECT
    !==================================================================

    ! ---------- new: 2016 , april, Karsten Lettmann after Hint by Qi and ayumi.fujisaki@noaa.gov------
    ! Initialize some variables 
    GWATER_FLUX_P%curr_stkcnt = 0 ; GWATER_FLUX_N%curr_stkcnt = 0
    GWATER_TEMP_P%curr_stkcnt = 0 ; GWATER_TEMP_N%curr_stkcnt = 0
    GWATER_SALT_P%curr_stkcnt = 0 ; GWATER_SALT_N%curr_stkcnt = 0
    ! --------- end new ----------------------------------------------------------------




    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END GROUND_WATER"
  END SUBROUTINE GROUND_WATER
  !================================================================
  !================================================================
  SUBROUTINE TIDAL_ELEVATION
    IMPLICIT NONE

    ! VARIABLES TO CHECK OBC NODE LIST
    INTEGER MYNOBC
    INTEGER, ALLOCATABLE :: MYOBCLIST(:)

    REAL(SP), POINTER :: STORAGE_VEC(:)
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR

    ! SOME HANDY VARIABLES TO PLAY WITH
    LOGICAL FOUND, VALID
    INTEGER NTIMES
    TYPE(TIME) :: TIMETEST
    real(SP) rbuf,float_time
    integer status, I, J

    ! SOME TEST STUFF FOR BRACKET
    Character(len=80):: dstring
    Character(len=80) :: dformat, tzone

    REAL(SP), ALLOCATABLE :: MYPERIOD(:)

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START TIDAL_ELEVATION"


    ! ONLY ESCAPE EARLY IF EQUI_TIDE IS OFF. OTHERWISE WE STILL NEED
    ! THE TIDAL FORCING FILE
    IF (.NOT. OBC_ELEVATION_FORCING_ON ) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) &
        "! TIDAL ELEVATION FORCING IS EITHER IN NESTING OR OFF!"
       ALLOCATE(TIDE_FORCING_COMMENTS(1))
       TIDE_FORCING_COMMENTS="TIDAL ELEVATION FORCING IS EITHER IN NESTING OR OFF!"
       RETURN
    END IF


    ! BOTH ASCII AND NETCDF NON JULIAN DATA FILES HAVE A NCFILE
    ! POINTER. THE DUMMY POINTER WAS CREATED FOR THE ASCII FILE AS A
    ! WAY TO 'TRICK' THE CODE AND CONTAIN THE NUMBER OF CONTROL VARIABLES

    ! FIND THE TIDAL FORCING FILE OBJECT
    TIDE_FILE => FIND_FILE(FILEHEAD,trim(OBC_ELEVATION_FILE),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND OPEN BOUNDARY CONDITION ELEVATION FILE OBJECT",&
         & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE))

    ATT => FIND_ATT(TIDE_FILE,"type",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("IN OPEN BOUNDARY CONDITION ELEVATION FILE OBJECT",&
         & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
         &"COULD NOT FIND GLOBAL ATTRIBURE: 'type'")


    SELECT CASE(TRIM(ATT%CHR(1)))
       !=================================
       ! NON JULIAN ELEVATION FORCING DATA
    CASE("FVCOM NON JULIAN ELEVATION FORCING FILE",&
         & "FVCOM SPECTRAL ELEVATION FORCING FILE")
       !================================== 
       ATT => FIND_ATT(TIDE_FILE,"components",FOUND)
       IF(FOUND) THEN
          ALLOCATE(TIDE_FORCING_COMMENTS(size(ATT%CHR)+1))
          TIDE_FORCING_COMMENTS(1)= "Spectral Tidal Forcing Components:"
          TIDE_FORCING_COMMENTS(2:)= ATT%CHR(:)
       ELSE
          CALL WARNING("ATTRIBUTE 'components' IS MISSING IN THE TIDAL FORCING FILE!")
          ALLOCATE(TIDE_FORCING_COMMENTS(1))
          TIDE_FORCING_COMMENTS = "Spectral Tidal Forcing Components&
               &: UNKNOWN"
       END IF


       DIM => FIND_DIM(TIDE_FILE,'tidal_components',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND DIMENSION 'tidal_components'")

       NTIDECOMPS = DIM%DIM

       ! LOAD TIDAL PERIOD DATA
       ALLOCATE(PERIOD(NTIDECOMPS),stat=status)
       IF (0 /= status) CALL FATAL_ERROR("TIDAL_ELEVATION COULD NOT &
            &ALLOCATE 'NTIDECOMPS'")

       VAR => FIND_VAR(TIDE_FILE,'tide_period',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND THE VARIABLE 'tide_period'")

       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND PERIOD VARIRIABLE'S ATTRIBUTE 'units'")

       if(trim(ATT%CHR(1)) .NE. 'seconds')  CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"PERIOD VARIRIABLE ATTRIBUTE 'units' SHOULD BE 'seconds'")

       CALL NC_CONNECT_AVAR(VAR,PERIOD)
       CALL NC_READ_VAR(VAR)

       ! LOAD THE Time Origin data if present
       VAR => FIND_VAR(TIDE_FILE,'time_origin', FOUND)
       IF (FOUND) THEN

          IF( IS_VALID_DATETIME( VAR,tzone)) THEN

             CALL NC_CONNECT_AVAR(VAR,dstring)
             CALL NC_READ_VAR(VAR)

             SPECTIME = READ_DATETIME(dstring,DATE_FORMAT,tzone,status)
             IF(Status == 0) CALL FATAL_ERROR&
                  & ("Could not read date in 'time_origin' attribute of spectral forcing file")

          ELSE IF(IS_VALID_FLOAT_DAYS( VAR,tzone)) THEN

             CALL NC_CONNECT_AVAR(VAR,float_time)
             CALL NC_READ_VAR(VAR)

             SPECTIME = days2time(float_time)

          ELSE IF(IS_VALID_FLOAT_SECONDS( VAR,tzone)) THEN

             CALL NC_CONNECT_AVAR(VAR,float_time)
             CALL NC_READ_VAR(VAR)

             SPECTIME = seconds2time(float_time)

          ELSE
             CALL FATAL_ERROR("SPECTRAL TIDAL FORCING TIME ORIGIN VA&
                  &RIABLE MUST BE A CHARACTER STRING Date or a float&
                  &ing point time???")
          END IF

       ELSE
          CALL WARNING("Setting Spectral Tidal Phase Time Orgin to 0.0 MJD")
          SPECTIME%MJD =0
          SPECTIME%MUSOD =0
       END IF



  if(IOBCN_GL>0) then
       ! LOOK FOR THE DIMENSIONS 'nobc' and 'tidal_compontents'
       DIM => FIND_DIM(TIDE_FILE,'nobc',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND DIMENSION 'nobc'")

       IF(IOBCN_GL /= DIM%DIM) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"THE DIMENSION 'nobc' MUST MATCH THE NUMBER OF OBC NODES")


       ! LOAD GLOBAL OPEN BOUNDARY NODE NUMBER DATA AND COMPARE WITH
       ! OBC.DAT/RESTART FILE INPUT
       ALLOCATE(MYOBCLIST(IOBCN))
       VAR => FIND_VAR(TIDE_FILE,'obc_nodes',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND VARIABLE 'obc_nodes'")
       CALL NC_CONNECT_AVAR(VAR,MYOBCLIST)
       CALL NC_READ_VAR(VAR)

       DO I = 1, IOBCN

          IF(I_OBC_N(I) /= NLID(MYOBCLIST(I))) THEN
             write(IPT,*) "NLID(MYOBCLIST)= ", NLID(MYOBCLIST(I)), "; I=",I
             write(IPT,*) "I_OBC_N= ", I_OBC_N(I), "; I=",I
             CALL FATAL_ERROR &
                  & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
                  & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
                  &"THE LIST OF BOUNDARY NODES DOES NOT MATCH")
          END IF
       END DO

       ! LOAD THE ELEVATION REFERENCE LEVEL DATA
       ALLOCATE(EMEAN(IOBCN))
       VAR => FIND_VAR(TIDE_FILE,'tide_Eref',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND VARIABLE 'tide_Eref'")

       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND ELEVATION REFERENCE VARIRIABLE'S ATTRIBUTE 'units'")

       if(trim(ATT%CHR(1)) .NE. 'meters')  CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"ELEVATION REFERENCE VARIRIABLE ATTRIBUTE 'units' SHOULD BE 'meters'")

       CALL NC_CONNECT_AVAR(VAR,EMEAN)
       CALL NC_READ_VAR(VAR)

       ! LOAD THE AMPLITUDE DATA
       ALLOCATE(APT(IOBCN,NTIDECOMPS))
       VAR => FIND_VAR(TIDE_FILE,'tide_Eamp',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND VARIABLE 'tide_Eamp'")

       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND AMPLITUDE VARIRIABLE'S ATTRIBUTE 'units'")

       if(trim(ATT%CHR(1)) .NE. 'meters')  CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"AMPLITUDE VARIRIABLE ATTRIBUTE 'units' SHOULD BE 'meters'")


       CALL NC_CONNECT_AVAR(VAR,APT)
       CALL NC_READ_VAR(VAR)

       ! LOAD THE PHASE DATA
       ALLOCATE(PHAI(IOBCN,NTIDECOMPS))
       VAR => FIND_VAR(TIDE_FILE,'tide_Ephase',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SPECTRAL ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND VARIABLE 'tide_Ephase'")

       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND PHASE VARIRIABLE'S ATTRIBUTE 'units'")

       if(ATT%CHR(1)(1:7) .NE. 'degrees')  CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"PHASE VARIRIABLE ATTRIBUTE 'units' SHOULD BE 'degrees'")

       CALL NC_CONNECT_AVAR(VAR,PHAI)
       CALL NC_READ_VAR(VAR)


       PHAI = MOD(PHAI,360.0_SP)

       !--REPORT RESULTS--------------------------------------------------------------!

       RBUF = MAXVAL(APT)
       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*)'!'
          WRITE(IPT,*  )'!  SPECTRAL TIDE         :    SET'
          WRITE(IPT,101)'!  MAX TIDE AMPLITUDE    : ',RBUF
          CALL PRINT_REAL_TIME(SpecTime,IPT,"Tide Time Origin")
       END IF

    endif

       ! SET THE TIDE FORCING TYPE FOR USE IN UPDATE
       TIDE_FORCING_TYPE = TIDE_FORCING_SPECTRAL


       !=========================
       ! NON JULIAN FORCING DATA IN AN ASCII FILE 
    CASE("ASCII FILE DUMMY ATTRIBUTE")
       !=========================

       
       CALL LOAD_JULIAN_OBC(NTIDECOMPS,TIDE_FORCING_COMMENTS&
            &,PERIOD,APT_EQI, BETA_EQI, TIDE_TYPE,APT,PHAI,EMEAN,SPECTIME)

        PHAI = MOD(PHAI,360.0_SP)

       !--REPORT RESULTS--------------------------------------------------------------!

       RBUF = MAXVAL(APT)
       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*)'!'
          WRITE(IPT,*  )'!  SPECTRAL TIDE         :    SET'
          WRITE(IPT,101)'!  MAX TIDE AMPLITUDE    : ',RBUF
          CALL PRINT_REAL_TIME(SpecTime,IPT,"Tide Time Origin")
       END IF

       ! SET THE TIDE FORCING TYPE FOR USE IN UPDATE
       TIDE_FORCING_TYPE = TIDE_FORCING_SPECTRAL


       !=========================
       ! TIME SERIES FORCING DATA
    CASE("FVCOM JULIAN TIME SERIES ELEVATION FORCING FILE", &
         & "FVCOM TIME SERIES ELEVATION FORCING FILE")
       !=========================

       ALLOCATE(TIDE_FORCING_COMMENTS(1))
       ATT => FIND_ATT(TIDE_FILE,"title",FOUND)
       IF(FOUND) THEN
          TIDE_FORCING_COMMENTS = "Tidal Forcing Time Series Title: "&
               &//TRIM(ATT%CHR(1))
       ELSE
          CALL WARNING("ATTRIBUTE 'title' IS MISSING IN THE TIDAL FORCING FILE!")
          TIDE_FORCING_COMMENTS = "Tidal Forcing Time Series Title: UNKNOWN"
       END IF



       ! LOOK FOR THE DIMENSIONS 'nobc' and 'time'
       DIM => FIND_DIM(TIDE_FILE,'time',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND DIMENSION 'time'")

       NTIMES = DIM%DIM

       DIM => FIND_DIM(TIDE_FILE,'nobc',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND DIMENSION 'nobc'")

       if(IOBCN_GL>0) then
       IF(IOBCN_GL /= DIM%DIM) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"THE DIMENSION 'nobc' MUST MATCH THE NUMBER OF OBC NODES")
       endif

       ! lOAD GLOBAL OPEN BOUNDARY NOD NUMBER DATA AND COMPARE WITH
       ! OBC.DAT/RESTART FILE INPUT
       ALLOCATE(MYOBCLIST(IOBCN))
       VAR => FIND_VAR(TIDE_FILE,'obc_nodes',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND VARIABLE 'obc_nodes'")
       CALL NC_CONNECT_AVAR(VAR,MYOBCLIST)
       CALL NC_READ_VAR(VAR)

       DO I = 1, IOBCN

          IF(SERIAL) THEN
             IF(I_OBC_N(I) /= MYOBCLIST(I)) THEN
                write(IPT,*) "NLID(MYOBCLIST)= ", MYOBCLIST(I), "; I=",I
                write(IPT,*) "I_OBC_N= ", I_OBC_N(I), "; I=",I
                CALL FATAL_ERROR &
                     & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
                     & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
                     &"THE LIST OF BOUNDARY NODES DOES NOT MATCH")
             END IF
          ELSE
          END IF
       END DO

       ! LOAD TIME AND CHECK TO MAKE SURE THE TIME RANGE IS VALID

       TIMETEST = GET_FILE_TIME(TIDE_FILE,1)

       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE ELVATION TIME SERIES")

       TIMETEST = GET_FILE_TIME(TIDE_FILE,NTIMES)


       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"THE MODEL RUN ENDS AFTER THE ELVATION TIME SERIES")

       VAR => FIND_VAR(TIDE_FILE,'elevation',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND VARIABLE 'elevation'")

       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"COULD NOT FIND ELEVATION VARIRIABLE'S ATTRIBUTE 'units'")

       if(TRIM(ATT%CHR(1)) .NE. 'meters')  CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TIME SERIES ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"ELEVATION VARIRIABLE ATTRIBUTE 'units' SHOULD BE 'meters'")


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_VEC(IOBCN), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN TIDAL_ELEVATION")
       TIDE_ELV_N => reference_var(var)
       CALL NC_CONNECT_PVAR(TIDE_ELV_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_VEC(IOBCN), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN TIDAL_ELEVATION")
       TIDE_ELV_P => reference_var(var)
       CALL NC_CONNECT_PVAR(TIDE_ELV_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       tide_elv_p%curr_stkcnt = 0; tide_elv_n%curr_stkcnt = 0  ! Siqi Li, 2021-01-27


       ! SINCE NO DATA HAS BEEN LOADED THERE IS NOT MUCH TO REPORT
       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*)'!'
          WRITE(IPT,*  )'!  TIME SERIES TIDE      :    SET'
       END IF

       ! SET THE TIDE FORCING TYPE FOR USE IN UPDATE
       TIDE_FORCING_TYPE = TIDE_FORCING_TIMESERIES


       ! NOT MUCH ELSE TO REPORT SINCE WE DON'T LOAD ANY DATA NOW

       !=====================================
       ! DEFAULT CASE IF GLOBAL ATTRIBUTES OF FILE ARE INCORRECT
    CASE DEFAULT
       !=====================================
       CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION ELEVATION FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_ELEVATION_FILE),&
            &"THE GLOBAL ATTRIBURE 'type' RETURNED UNKNOWN TYPE:",&
            & TRIM(ATT%CHR(1)))
    END SELECT

    ! Removed by Siqi Li, 2021-01-27
    ! afm 20151112 & EJA 20160921
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3. 
!    tide_elv_p%curr_stkcnt = 0; tide_elv_n%curr_stkcnt = 0

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END TIDAL_ELEVATION"

101 FORMAT(1X,A26,F10.4)  
  END SUBROUTINE TIDAL_ELEVATION
  !================================================================
  !================================================================

  ! Rules about river names: River names must be unique
  ! examples: r1, r2, r3
  !         : Mississippi, Connecticut, St.Lawrence
  ! BAD EXAMPLES: Miss, Mississippi, Missouri
  !             : R1, R2, R

  SUBROUTINE RIVER_DISCHARGE
    IMPLICIT NONE
    INTEGER :: I, J,K, fcnt,rcnt,status, nfiles,nrs,IOS,NS
    TYPE(A_RIVER_FILE) DUMMY

    TYPE(NCFILE),POINTER :: NCF
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    TYPE(NCVAR), POINTER :: DUM_P

    REAL(SP), POINTER :: STORAGE_VEC(:)
    LOGICAL :: FOUND, MINE
    CHARACTER(LEN=7) :: chr
    character(len=20), allocatable :: dist_strings(:)

    REAL(SP) :: MYDIST(KBM1)

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START RIVER_DISCHARGE"

    NULLIFY(STORAGE_VEC)

    ! TRANSLATE TO THE OLD NAME FOR TOTAL NUMBER OF RIVERS
    NUMQBC_GL = river_number

    IF (river_number == 0 ) THEN
       NUMQBC = 0
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! THERE ARE NO RIVERS IN THIS MODEL"
       ALLOCATE(RIVER_FORCING_COMMENTS(1))
       RIVER_FORCING_COMMENTS(1) = "THERE ARE NO RIVERS IN T&
            &HIS MODEL"
       RETURN
    END IF


    IF (DBG_SET(DBG_SBRIO))THEN
       WRITE(IPT,*) "Total Number Of Rivers = ",river_number 
       WRITE(IPT,*) "RIVER_TS_SETTING = "//TRIM(RIVER_TS_SETTING)
       WRITE(IPT,*) "RIVER_INFLOW_LOCATION = "//trim(RIVER_INFLOW_LOCATION)
       WRITE(IPT,*) "RIVER_KIND = "//trim(RIVER_KIND)

       WRITE(IPT,*)"============================="
       DO I =1,river_number 
          write(IPT,*) "River number:",I
          WRITE(IPT,*) "River File ="//TRIM(RIVERS(i)%FILE)
          WRITE(IPT,*) "River Name ="//trim(RIVERS(i)%NAME)
          WRITE(IPT,*) "River Location =",RIVERS(i)%LOCATION

          WRITE(IPT,*) "River Distribution ="//TRIM(RIVERS(i)%distribution)

          WRITE(IPT,*)"============================="
       END DO
    END IF

    IF(TRIM(RIVER_TS_SETTING) /= 'calculated' .AND. TRIM(RIVER_TS_SETTING) /= 'specified') THEN
       CALL FATAL_ERROR("RIVER_TS_SETTING NOT CORRECT IN NAMELIST","SHOULD BE 'calculated' or 'specified'")
    END IF

    ALLOCATE(RIVER_FORCING_COMMENTS(river_number + 3))
    WRITE(CHR,'(I7)')river_number
    RIVER_FORCING_COMMENTS(1) = "THERE ARE "//TRIM(adjustl(CHR))//" RIVERS IN THIS MODEL."

    RIVER_FORCING_COMMENTS(2) = "RIVER INFLOW IS ON THE "//TRIM(RIVER_INFLOW_LOCATION)&
         &//"s WHERE TEMPERATURE AND SALINITY ARE "//TRIM(RIVER_TS_SETTING)//" IN THE MODEL."

    RIVER_FORCING_COMMENTS(3) = "THE FOLLOWING RIVER NAMES ARE USED:"
    DO I =1,river_number
       RIVER_FORCING_COMMENTS(3+I) = TRIM(RIVERS(i)%NAME)
    END DO

    !    SELECT CASE(trim(RIVER_INFLOW_LOCATION))

    ! CHECK TO MAKE SURE THE LOCATION IS VALID AND ADD NAMES TO
    ! COMMENT STRING
    nfiles = 0
    NUMQBC = 0
    DO  I =1,river_number 


       ! CHECK THE INFLOW LOCATION OF EACH RIVER
       ! AND
       ! COUNT THE NUMBER THAT BELONG TO EACH PROCESSOR
       SELECT CASE(trim(RIVER_INFLOW_LOCATION))
       CASE('node')

          IF(RIVERS(i)%LOCATION > MGL .or. RIVERS(i)%LOCATION < 1)THEN
             write(chr,'(I7)') RIVERS(i)%LOCATION
             CALL FATAL_ERROR ("RIVER_DISCHARGE: FOR THE RIVER NAMED: "&
                  &//trim(RIVERS(i)%NAME),"THE RIVER GRID LOCATION IN&
                  & THE NAME LIST IS NOT IN THE GLOBAL DOMAIN",&
                  & "YOU SPECIFIED NODE NUMBER: "//chr)
          END IF

          ! COUNT THE NUMBER OF RIVERS OWNED BY THIS PROC
          IF (NLID(RIVERS(I)%LOCATION) .GT. 0) NUMQBC = NUMQBC + 1


       CASE('edge') 
          IF(RIVERS(i)%LOCATION > NGL .or. RIVERS(i)%LOCATION < 1)THEN
             write(chr,'(I7)') RIVERS(i)%LOCATION
             CALL FATAL_ERROR ("RIVER_DISCHARGE: FOR THE RIVER NAMED: "&
                  &//trim(RIVERS(i)%NAME),"THE RIVER GRID LOCATION IN&
                  & THE NAME LIST IS NOT IN THE GLOBAL DOMAIN",&
                  & "YOU SPECIFIED CELL NUMBER: "//chr)
          END IF

          ! COUNT THE NUMBER OF RIVERS OWNED BY THIS PROC
          IF (ELID(RIVERS(I)%LOCATION) .GT. 0) NUMQBC = NUMQBC + 1




       CASE DEFAULT
          CALL FATAL_ERROR("RIVER_INFLOW_LOCATION: NOT CORRECT IN NAMELIST",&
               & "SHOULD BE 'node' or 'edge'")
       END SELECT


       ! COUNT THE NUMBER OF FILES
       NCF => FIND_FILE(FILEHEAD,TRIM(Rivers(I)%FILE),FOUND)
       IF (.NOT.FOUND) CALL FATAL_ERROR &
            & ("RIVER_DISCHARGE: COULD NOT FIND RIVER FILE OBJECT NAMED: &
            &"//TRIM(RIVERS(I)%FILE))
       
       IF(.NOT. associated(NCF%FTIME)) CALL FATAL_ERROR&
            &("RIVER FILE DID NOT LOAD PROPERLY",& 
            & "File name:"//trim(NCF%FNAME),&
            &"Please check the time format!")


       IF (NCF%FTIME%PREV_STKCNT /= 999) THEN
          nfiles = nfiles +1
          NCF%FTIME%PREV_STKCNT = 999
       END IF

    END DO


    ! ALLOCATE THE SPACE FOR THE RIVER FILES AND NULLIFY POINTERS
    ALLOCATE(RIVER_FORCING(NFILES))
    DO I =1,NFILES
       NULLIFY(RIVER_FORCING(I)%NCF)
       NULLIFY(RIVER_FORCING(I)%FLUX_N)
       NULLIFY(RIVER_FORCING(I)%FLUX_P)
       NULLIFY(RIVER_FORCING(I)%TEMP_N)
       NULLIFY(RIVER_FORCING(I)%TEMP_P)
       NULLIFY(RIVER_FORCING(I)%SALT_N)
       NULLIFY(RIVER_FORCING(I)%SALT_P)
    END DO


    ! ALLOCATE THE SPACE FOR THE RIVER DATA
    ALLOCATE(RIV_GL2LOC(NUMQBC))
    ALLOCATE(INODEQ(NUMQBC))
    ALLOCATE(ICELLQ(NUMQBC))  ! THE CELL INDEX
    ALLOCATE(N_ICELLQ(NUMQBC,2)) ! THE NODES BOUNDING THE EDGE
    ALLOCATE(VQDIST(NUMQBC,KBM1));  VQDIST = 0.0_SP
    ALLOCATE(QDIS(NUMQBC));         QDIS   = 0.0_SP
    ALLOCATE(QDIS2(NUMQBC));        QDIS2  = 0.0_SP
    ALLOCATE(TDIS(NUMQBC));         TDIS   = 0.0_SP
    ALLOCATE(SDIS(NUMQBC));         SDIS   = 0.0_SP
    ALLOCATE(QAREA(NUMQBC));        QAREA  = 0.0_SP
    ALLOCATE(ANGLEQ(NUMQBC));       ANGLEQ = 0.0_SP
    ALLOCATE(VLCTYQ(NUMQBC));       VLCTYQ = 0.0_SP
    ALLOCATE(RDISQ(NUMQBC,2));      RDISQ  = 0.0_SP

    FCNT = 0
    RCNT = 0
    DO  I =1,river_number 
       ! SET MINE TO FALSE
       MINE = .false. 
       ! Mine is set true if this river name belongs
       ! to this processor



       ! MAKE THE LOCAL AND GLOBAL INDEX
       SELECT CASE(trim(RIVER_INFLOW_LOCATION))
       CASE('node')

          IF (NLID(RIVERS(I)%LOCATION) .GT. 0) THEN
             rcnt = rcnt + 1
             INODEQ(rcnt)       = NLID(RIVERS(I)%LOCATION)
             RIV_GL2LOC(rcnt)   = RIVERS(I)%LOCATION
             MINE = .TRUE.
          END IF




       CASE('edge')

          IF (ELID(RIVERS(I)%LOCATION) .GT. 0) THEN
             rcnt = rcnt + 1
             ICELLQ(rcnt)       = ELID(RIVERS(I)%LOCATION)
             RIV_GL2LOC(rcnt)   = RIVERS(I)%LOCATION
             MINE = .true.
          END IF

       END SELECT


       ! MAKE VQDIST
       IF(MINE) THEN

!--------------------------------------------------------------------
!--------------------------------------------------------------------

          CALL SET_DISTRIBUTION(RIVERS(I)%DISTRIBUTION,RIVER_INFLOW_LOCATION,RIVERS(I)%LOCATION,MYDIST)
          VQDIST(RCNT,1:kbm1)=MYDIST
          
!--------------------------------------------------------------------
!--------------------------------------------------------------------

       END IF

       ! NOW PUT THE RIVERS IN THE RIVER_FORCING TYPE
       NCF => FIND_FILE(FILEHEAD,TRIM(Rivers(I)%FILE),FOUND)

       IF (NCF%FTIME%PREV_STKCNT == 999) THEN
          FCNT = FCNT + 1
          RIVER_FORCING(FCNT)%NCF => NCF
          NCF%FTIME%PREV_STKCNT = 0

          DIM => FIND_DIM(NCF,'rivers',FOUND)
          IF(.not.FOUND) CALL FATAL_ERROR &
               & ("COULD NOT FIND DIMENSION 'rivers'",&
               & "In the file: "//trim(NCF%FNAME) )

          RIVER_FORCING(FCNT)%RIVERS_IN_FILE=DIM%DIM

          ALLOCATE(RIVER_FORCING(FCNT)%RIV_FILE2LOC(DIM%DIM))
          RIVER_FORCING(FCNT)%RIV_FILE2LOC = 0

       END IF

       ! FIGURE OUT WHICH RIVER IT IS WHICH FILE AND SET THE INDEX 
       DO J = 1, NFILES
          ! IN THIS FILE FIND THE NAME
          IF (associated(NCF,RIVER_FORCING(J)%NCF)) THEN
             K=SEARCH_NAME(NCF,Rivers(I)%NAME)

             ! THIS IS MY RIVER SO SET THE INDEX
             IF (MINE) RIVER_FORCING(J)%RIV_FILE2LOC(K)=RCNT

          END IF
       END DO


    END DO

    IF (FCNT /= NFILES) CALL FATAL_ERROR&
         ("RIVER_DISCHARGE: WE LOST A RIVER FILE IN THE MIDDLE OF NOWHERE!")

    IF (RCNT /= NUMQBC) CALL FATAL_ERROR&
         ("RIVER_DISCHARGE: WE LOST A RIVER IN THE MIDDLE OF NOWHERE!")


    DO I = 1, NFILES

       NCF => RIVER_FORCING(I)%NCF

       SELECT CASE (RIVER_KIND)
       CASE(PRDC)
          CALL CHECK_RIVER_FILE(NCF, RIVER_FORCING(I)%river_period)
       CASE(VRBL)          
          CALL CHECK_RIVER_FILE(NCF)
       CASE DEFAULT
          CALL FATAL_ERROR("Invalid RIVER_KIND in namelist runfile:",&
               & " Options are: "//TRIM(PRDC)//" or "//TRIM(VRBL))
       END SELECT

       DIM => FIND_DIM(NCF,'rivers',FOUND)

       nrs = DIM%dim

       ! GET THE FLUX VARIABLE
       VAR => FIND_VAR(NCF,'river_flux',FOUND)

       ALLOCATE(STORAGE_VEC(nrs), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN RIVER_DISCHARGE")
       RIVER_FORCING(I)%FLUX_N => reference_var(var)
       CALL NC_CONNECT_PVAR(RIVER_FORCING(I)%FLUX_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ALLOCATE(STORAGE_VEC(nrs), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN RIVER_DISCHARGE")
       RIVER_FORCING(I)%FLUX_P => reference_var(var)
       CALL NC_CONNECT_PVAR(RIVER_FORCING(I)%FLUX_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! GET THE TEMPERATURE VARIABLE
       VAR => FIND_VAR(NCF,'river_temp',FOUND)

       ALLOCATE(STORAGE_VEC(nrs), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN RIVER_DISCHARGE")
       RIVER_FORCING(I)%TEMP_N => reference_var(var)
       CALL NC_CONNECT_PVAR(RIVER_FORCING(I)%TEMP_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ALLOCATE(STORAGE_VEC(nrs), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN RIVER_DISCHARGE")
       RIVER_FORCING(I)%TEMP_P => reference_var(var)
       CALL NC_CONNECT_PVAR(RIVER_FORCING(I)%TEMP_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! GET THE SALINITY VARIABLE
       VAR => FIND_VAR(NCF,'river_salt',FOUND)

       ALLOCATE(STORAGE_VEC(nrs), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN RIVER_DISCHARGE")
       RIVER_FORCING(I)%SALT_N => reference_var(var)
       CALL NC_CONNECT_PVAR(RIVER_FORCING(I)%SALT_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ALLOCATE(STORAGE_VEC(nrs), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN RIVER_DISCHARGE")
       RIVER_FORCING(I)%SALT_P => reference_var(var)
       CALL NC_CONNECT_PVAR(RIVER_FORCING(I)%SALT_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! afm 20160516 & EJA 20160921
       RIVER_FORCING(I)%FLUX_N%curr_stkcnt = 0; RIVER_FORCING(I)%FLUX_P%curr_stkcnt = 0
       RIVER_FORCING(I)%TEMP_N%curr_stkcnt = 0; RIVER_FORCING(I)%TEMP_P%curr_stkcnt = 0
       RIVER_FORCING(I)%SALT_N%curr_stkcnt = 0; RIVER_FORCING(I)%SALT_P%curr_stkcnt = 0 



    END DO

    ! SET THE RIVER BNDRY METRICS - USED TO BE SET_BNDRY
    CALL SET_RIVER_BNDRY_METRICS


    IF(DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*)"!"
       WRITE(IPT,*)"!  RIVER FORCING ON"
       WRITE(IPT,*)'!  GLOBAL NUMBER OF RIVERS :',river_number
       WRITE(IPT,*)'!  NUMBER OF RIVER FILES   :', nfiles
    END IF

    IF(DBG_SET(DBG_SCL)) THEN
       WRITE(IPT,*)"/////////=============================///////////"
       WRITE(IPT,*)"   PRINTING RIVER FORCING DETAILS"

       WRITE(IPT,*)"  LOCAL NUMBER OF RIVERS   : ", numqbc

       WRITE(IPT,*)"============================="       
       DO I = 1,nfiles
          WRITE(IPT,*)" FILE NAME: "//TRIM(RIVER_FORCING(I)%NCF%FNAME)
          !             WRITE(IPT,*)" RIVER NAME: "//TRIM(RIVER_FORCING(I)%NAME)
          WRITE(IPT,*)" NUMBER IN FILE=",RIVER_FORCING(I)%RIVERS_IN_FILE
          WRITE(IPT,*)" RIV_FILE2LOC = ",RIVER_FORCING(I)%RIV_FILE2LOC
          WRITE(IPT,*)"============================="
       END DO

       WRITE(IPT,*)"/////////=============================///////////"
    END IF

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END RIVER_DISCHARGE"



  END SUBROUTINE RIVER_DISCHARGE
  !================================================================
  !================================================================
  FUNCTION SEARCH_NAME(NCF,NAME) RESULT(RES)
    ! OUTPUT is the TYPE we are trying to set
    ! Input is the River type we are searching from
    IMPLICIT NONE
    INTEGER :: RES
    TYPE(NCFILE), POINTER :: NCF
    CHARACTER(LEN=*) :: NAME

    INTEGER :: I, rvrs_in_file,strlen,status
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL FOUND

    IF (DBG_SET(DBG_IO)) THEN
       WRITE(IPT,*)"SEARCH_NAME (RIVERS)"
       WRITE(IPT,*)"============================="
       write(IPT,*)"LOOKING FOR: '"//TRIM(NAME)//"'"
       WRITE(IPT,*)"=========="
       CALL PRINT_FILE(NCF)
    END IF


    ! FIND THE RIVER NAME IN THE FILE
    DIM => FIND_DIM(NCF,'rivers',FOUND)
    rvrs_in_file = DIM%DIM

    VAR => FIND_VAR(NCF,'river_names',found)

    ! ATTEMPT TO ONLY READ NAMES LIST ONCE
    IF(.NOT. ASSOCIATED(VAR%VEC_CHR)) THEN

       ALLOCATE(VAR%VEC_CHR(DIM%DIM),STAT=STATUS)
       IF(STATUS/=0) CALL FATAL_ERROR("SEARCH_NAME: CAN NOT ALLOCATE TEMP!")

       CALL NC_READ_VAR(VAR)

    END IF

    DO I = 1,rvrs_in_file
       IF(VAR%VEC_CHR(I) .EQ. NAME)THEN
          RES = I
          RETURN
       END IF
    END DO

    

    ! SHOULD NOT BE HERE, SOMETHING IS WRONG

    IF (DBG_SET(DBG_LOG)) THEN
       WRITE(IPT,*)"============================="
       write(IPT,*)"LOOKING FOR: '"//TRIM(NAME)//"'; In File:"
       WRITE(IPT,*)"============================="
       CALL PRINT_FILE(NCF)
       DO I = 1,rvrs_in_file
          WRITE(IPT,*) "RIVER NAMES: "//TRIM(VAR%VEC_CHR(I))
       END DO
       WRITE(IPT,*)"============================="
       WRITE(IPT,*)"============================="
       WRITE(IPT,*)"============================="
    END IF
    

    CALL FATAL_ERROR("COULD NOT FIND CORRECT NAME IN RIVER FILE?")

  END FUNCTION SEARCH_NAME
  !================================================================
  !================================================================  
  SUBROUTINE CHECK_RIVER_FILE(NCF,PERIOD)
    ! CALL FATAL_ERROR IF THERE IS ANYTHING WRONG WITH THE FILE
    IMPLICIT NONE
    TYPE(NCFILE), POINTER ::NCF
    TYPE(TIME),OPTIONAL :: PERIOD
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR

    TYPE(TIME) :: fSTART, fEND
    INTEGER :: NTIMES,NS
    LOGICAL FOUND

    IF(.NOT. ASSOCIATED(NCF)) CALL FATAL_ERROR &
         & ("THE RIVER FILE OBJECT PASSED TO 'check_river_file' IS NOT ASSOCIATED")

    ! CHECK DIMENSIONS
    DIM => FIND_DIM(NCF,'time',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"IS MISSING THE 'time' DIMENSION")

    nTIMES = DIM%dim

    DIM => FIND_DIM(NCF,'namelen',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"IS MISSING THE 'namelen' DIMENSION")

    DIM => FIND_DIM(NCF,'rivers',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"IS MISSING THE 'rivers' DIMENSION")

    ! CHECK VARIABLES AND THEIR ATTS
    IF(.NOT. ASSOCIATED(NCF%FTIME)) CALL FATAL_ERROR &
         & ('THE RIVER FILE '//TRIM(NCF%FNAME),&
         'DOES NOT HAVE A RECONGIZED TIME VARIABLE')

    IF(PRESENT(PERIOD)) THEN

       ! CHECK START AND END TIME FOR THE FILE:
       FSTART = GET_FILE_TIME(NCF,1)
       IF(ZEROTIME /= FSTART) THEN
          CALL PRINT_TIME(FSTART,IPT,"River Data Start")
          CALL FATAL_ERROR &
               & ("Date of the first river data point must be 0.0 for periodic forcoing mode:",&
               & "The River File: "//TRIM(NCF%FNAME)//'; has a bad start date.')
       END IF
       
       PERIOD = GET_FILE_TIME(NCF,NTIMES)
       IF(PERIOD .LE. zerotime) THEN
          CALL PRINT_REAL_TIME(PERIOD,IPT,"River Data End")
          
          CALL FATAL_ERROR &
               & ("Date of the last river data point must be greater than or equal to zero for periodic forcing mode:",&
               & "The River File: "//TRIM(NCF%FNAME)//'; has a bad end date.')
       END IF


    ELSE
       ! CHECK START AND END TIME FOR THE FILE:
       FSTART = GET_FILE_TIME(NCF,1)
       IF(FSTART > StartTIME) THEN
          CALL PRINT_REAL_TIME(StartTime,IPT,"Model Start")
          CALL PRINT_REAL_TIME(FSTART,IPT,"River Data Start")
          CALL FATAL_ERROR &
               & ("Date of the first river data point must be less than or equal to the model start date:",&
               & "The River File: "//TRIM(NCF%FNAME)//'; has a bad start date.')
       END IF
       
       FEND = GET_FILE_TIME(NCF,NTIMES)
       IF(FEND < EndTime) THEN
          CALL PRINT_REAL_TIME(EndTime,IPT,"Model End")
          CALL PRINT_REAL_TIME(FEND,IPT,"River Data End")
          
          CALL FATAL_ERROR &
               & ("Date of the last river data point must be greater than or equal to the model end date:",&
               & "The River File: "//TRIM(NCF%FNAME)//'; has a bad end date.')
       END IF
    END IF

    VAR => FIND_VAR(NCF,'river_names',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"IS MISSING THE 'river_names' VARIABLE")
    ! HAS NO ATTRIBUTES

    VAR => FIND_VAR(NCF,'river_flux',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"IS MISSING THE 'river_flux' VARIABLE")

    ATT => FIND_ATT(VAR,'units',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"THE VARIABLE: "&
         &//TRIM(VAR%VARNAME), "IS MISSING THE ATTRIBUTE 'units'")

    IF(TRIM(ATT%CHR(1)) /= "m^3s^-1")CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"THE VARIABLE: "&
         &//TRIM(VAR%VARNAME), "THE ATTRIBUTE 'units' IS INCORRECT: EXPE&
         &CTING 'm^3s^-1'")

    VAR => FIND_VAR(NCF,'river_temp',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"IS MISSING THE 'river_temp' VARIABLE")

    ATT => FIND_ATT(VAR,'units',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"THE VARIABLE: "&
         &//TRIM(VAR%VARNAME), "IS MISSING THE ATTRIBUTE 'units'")

    IF(TRIM(ATT%CHR(1)) /= "Celsius")CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"THE VARIABLE: "&
         &//TRIM(VAR%VARNAME), "THE ATTRIBUTE 'units' IS INCORRECT: EXPE&
         &CTING 'Celsius'")

    VAR => FIND_VAR(NCF,'river_salt',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"IS MISSING THE 'river_salt' VARIABLE")

    ATT => FIND_ATT(VAR,'units',FOUND)
    IF(.NOT. FOUND)  CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"THE VARIABLE: "&
         &//TRIM(VAR%VARNAME), "IS MISSING THE ATTRIBUTE 'units'")

    IF(TRIM(ATT%CHR(1)) /= "PSU")CALL FATAL_ERROR &
         & ("THE RIVER FILE:"//TRIM(NCF%FNAME),"THE VARIABLE: "&
         &//TRIM(VAR%VARNAME), "THE ATTRIBUTE 'units' IS INCORRECT: EXPE&
         &CTING 'PSU'")

  END SUBROUTINE CHECK_RIVER_FILE
  !==============================================================================|
  !  SET METRICS FOR THE BOUNDARY CONDITIONS       			       |
  !==============================================================================|  
  SUBROUTINE SET_RIVER_BNDRY_METRICS           
    IMPLICIT NONE
    REAL(DP)  DX12,DY12,DX32,DY32,ATMP1,ATMP2,DXYTMP,HTMP,AREATMP
    REAL(DP)  XNORM,YNORM,XP,YP,XN,YN,XI,YI,FAC,XNEXT,YNEXT,MODNR
    INTEGER I,J,I1,I2,I3,J1,J2,II,ITMP,JTMP,INODE,JNODE,KNODE,NNORM
    CHARACTER(len=7) :: strng

    !------------------------------------------------------------------------------!

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START SET_RIVER_BNDRY_METRICS"

    SELECT CASE(trim(RIVER_INFLOW_LOCATION))
       !=====================================
       ! CASE NODE: if the river inflow is on the nodes
    CASE('node')
       !=====================================

       ! WHY DON'T WE REQUIRE THAT RIVERS ARE ON THE BOUNDARY FOR NODE?

       DO I=1,NUMQBC
          J=INODEQ(I)
          if (ISONB(J) /= 1) THEN
             write(strng,'(I6)') ngid(j)
             CALL FATAL_ERROR &
                  & ("You seem to be trying to put a river in the middle of the domain",&
                  & "The global node number you selected is: "//trim(adjustl(strng)),&
                  &"but this is not a solid boundary node?")
          END if
          I1=NBSN(J,2)
          I2=NBSN(J,1)
          I3=NBSN(J,NTSN(J)-1)
          DY12=VY(I1)-VY(I2)
          DY32=VY(I3)-VY(I2)
          DX12=VX(I1)-VX(I2)
          DX32=VX(I3)-VX(I2)

          ATMP1=ATAN2(DY12,DX12)
          ATMP2=ATAN2(DY32,DX32)
          !         IF(ATMP1 < 0.0_SP) ATMP1=ATMP1+2.0_SP*3.1415927_SP
          !         IF(ATMP2 < 0.0_SP) ATMP2=ATMP2+2.0_SP*3.1415927_SP
          IF(ATMP1 < ATMP2) ATMP1=ATMP1+2.0_SP*3.1415927_SP
          DXYTMP=SQRT(DX12**2+DY12**2)+SQRT(DX32**2+DY32**2)
          QAREA(I)=0.5_SP*DXYTMP*H(INODEQ(I))
          ANGLEQ(I)=(ATMP1-ATMP2)/2.+ATMP2
       END DO

       !=====================================
       ! CASE EDGE: if the river inflow is on the cell edge
    CASE('edge')
       !=====================================
       DO I=1,NUMQBC
          II=ICELLQ(I)
          IF(ISBCE(II) /= 1) THEN

             write(strng,'(I6)') egid(II)
             CALL FATAL_ERROR &
                  & ("You seem to be trying to put a river in the middle of the domain",&
                  & "The global cell number you selected is: "//trim(adjustl(strng)),&
                  &"but this is not a solid boundary node?")
          END IF
          ITMP=0
          DO J=1,3
             IF(NBE(II,J) == 0) THEN
                JTMP=J
                ITMP=ITMP+1
             END IF
          END DO
          IF(ITMP /= 1) THEN

             write(strng,'(I6)') egid(II)
             CALL FATAL_ERROR &
                  & ("You have selected an invalide cell for edge based river inflow.",&
                  & "The global cell number you selected is: "//trim(adjustl(strng)),&
                  & "This cell has the wrong number of solid boundaries!")
          END IF
          J1=JTMP+1-INT((JTMP+1)/4)*3
          J2=JTMP+2-INT((JTMP+2)/4)*3
          I1=NV(II,J1)
          I2=NV(II,J2)
          N_ICELLQ(I,1)=I1
          N_ICELLQ(I,2)=I2
          HTMP=0.5_SP*(H(I1)+H(I2))
          DY12=VY(I1)-VY(I2)
          DX12=VX(I1)-VX(I2)
          ATMP1=ATAN2(DY12,DX12)
          QAREA(I)=SQRT(DX12**2+DY12**2)*HTMP
          ANGLEQ(I)=ATMP1+3.1415927/2.0
          AREATMP=ART1(I1)+ART1(I2)
          RDISQ(I,1)=ART1(I1)/AREATMP
          RDISQ(I,2)=ART1(I2)/AREATMP
       END DO
       !=====================================
       ! DEFAULT CASE: if the name list has a bad value
    CASE DEFAULT
       !=====================================
       CALL FATAL_ERROR("RIVER_INFLOW_LOCATION: NOT CORRECT IN NAMELIST",&
            & "SHOULD BE 'node' or 'edge' - It passed River_Discharge: how?")
    END SELECT


    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END SET_RIVER_BNDRY_METRICS"

    RETURN
  END SUBROUTINE SET_RIVER_BNDRY_METRICS
  !================================================================
  !================================================================
  SUBROUTINE SET_DISTRIBUTION(NAME,TYPE,LOC,MYDIST)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: NAME,TYPE
    INTEGER, INTENT(IN) :: LOC
    REAL(SP), INTENT(OUT) :: MYDIST(KBM1)

    REAL(SP) :: MYZ(KBM1)
    REAL(SP) :: MYH,MYEL

    REAL(SP) :: TOTAL

    CHARACTER(LEN=12) :: IDX
    
    ! FOR GET VAL================
    INTEGER :: NLINE, NCHAR, INTVAL(150), NVAL
    REAL(DP) :: REALVAL(150)
    CHARACTER(LEN=40) :: VARNAME 
    CHARACTER(LEN=80) :: STRINGVAL(150)
    CHARACTER(LEN=7) :: VARTYPE
    LOGICAL :: LOGVAL
    !============================

    MYDIST = 0.0_SP
    
    IF (trim(TYPE)=='node') THEN
       MYDIST = DZ(NLID(LOC),1:kbm1)
       MYZ = ZZ(NLID(LOC),1:KBM1)
       MYH = H(NLID(LOC))
       MYEL = EL(NLID(LOC))

    ELSEIF (trim(TYPE)=='edge') THEN
       MYDIST = DZ1(ELID(LOC),1:kbm1)
       MYZ = ZZ1(ELID(LOC),1:KBM1)
       MYH = H1(ELID(LOC))
       MYEL = EL1(ELID(LOC))
    ELSE
       CALL FATAL_ERROR("BAD RIVER LOCATION (edge OR node) ?")
    END IF

    SELECT CASE(NAME(1:6))
       ! UNIFORM FUNCTION
    CASE ('unifor')
       
       ! ALREADY SET VALUES
       WRITE(IPT,*) "UNIFORM RIVER DISTRIBUTION",MYDIST
       ! HEAVISIDE FUNCTION
    CASE('heavis')
       
       NLINE=-1
       NCHAR = LEN_TRIM(NAME)
       CALL GET_VALUE(NLINE,NCHAR,NAME,VARNAME,VARTYPE,LOGVAL&
            &,STRINGVAL,REALVAL,INTVAL,NVAL)
       
       IF(VARTYPE /= "float") THEN
          WRITE(IDX,*) LOC
          CALL FATAL_ERROR&
            &("HEAVISIDE RIVER DISTRIBUTION MUST SET A FLOATING POINT VALUE",&
            &"River on "//TRIM(TYPE)//" number:"//TRIM(IDX))
       END IF

       IF(NVAL>1) CALL FATAL_ERROR&
            &("COULD NOT READ RIVER DISTRIBUTION STRING?",&
            & "BAD STRING:"//TRIM(NAME))

       IF(index(VARNAME,'depth')/=0)THEN
          

          MYZ = (MYH+MYEL)*MYZ+MYEL

          IF(MYZ(KBM1) > REALVAL(1) .OR. REALVAL(1) > MYZ(1)) THEN
             WRITE(IDX,*) LOC
             WRITE(IPT,*) "================================"
             WRITE(IPT,*) "HEAVISIDE CASE- depth",REALVAL(1)
             WRITE(IPT,*) "RIVER DEPTH = ",MYH
             WRITE(IPT,*) "RIVER SURFACE = ",MYEL
             
             CALL FATAL_ERROR("RIVER DISTRIBUTION: Depth value out of bounds!",&
                  & "River on "//TRIM(TYPE)//" number:"//TRIM(IDX))
          END IF

          WRITE(IPT,*) "DEPTH:",MYZ

          WHERE (MYZ<REALVAL(1))
             MYDIST = 0.0_SP
          END WHERE

          TOTAL = SUM(MYDIST)
          MYDIST = MYDIST/TOTAL

       ELSEIF(index(VARNAME,'sigma')/=0)THEN
          

          IF(-1.0_SP > REALVAL(1) .OR. REALVAL(1) >0.0_SP) THEN
             WRITE(IPT,*) "================================"
             WRITE(IPT,*) "HEAVISIDE CASE- sigma",REALVAL(1)
             WRITE(IDX,*) LOC
             CALL FATAL_ERROR&
                  &("RIVER DISTRIBUTION: Sigma value out of bounds!",&
                  & "River on "//TRIM(TYPE)//" number:"//TRIM(IDX))
          END IF

          WHERE (MYZ<REALVAL(1))
             MYDIST = 0.0_SP
          END WHERE

          TOTAL = SUM(MYDIST)
          MYDIST = MYDIST/TOTAL

       ELSE

          CALL FATAL_ERROR("RIVER DISTRIBUTION: UNKNOWN HEAVISIDE SETTING?",&
               & "BAD STRING:"//TRIM(NAME))

       END IF
       WRITE(IPT,*) "HEAVISIDE RIVER DISTRIBUTION",MYDIST

     CASE('linear')
       

       NLINE=-1
       NCHAR = LEN_TRIM(NAME)
       CALL GET_VALUE(NLINE,NCHAR,NAME,VARNAME,VARTYPE,LOGVAL&
            &,STRINGVAL,REALVAL,INTVAL,NVAL)
       
       IF(VARTYPE /= "float") THEN
          WRITE(IDX,*) LOC
          CALL FATAL_ERROR&
               &("LINEAR RIVER DISTRIBUTION MUST SET A FLOATING POINT VALUE",&
               &"River on "//TRIM(TYPE)//" number:"//TRIM(IDX))
       END IF

       IF(NVAL>1) CALL FATAL_ERROR&
            &("COULD NOT READ RIVER DISTRIBUTION STRING?",&
            & "BAD STRING:"//TRIM(NAME))
       
       IF(index(VARNAME,'slope')/=0)THEN
          
          IF(REALVAL(1) <0.0_SP) THEN
             WRITE(IPT,*) "================================"
             WRITE(IPT,*) "LINEAR CASE- slope",REALVAL(1)
             WRITE(IDX,*) LOC
             CALL FATAL_ERROR&
                  &("RIVER DISTRIBUTION: linear slope less than zero!",&
                  & "River on "//TRIM(TYPE)//" number:"//TRIM(IDX))
          END IF

          MYZ = (MYH+MYEL)*MYZ+MYEL

          MYZ = MYZ *MYDIST * REALVAL(1)

          DO WHILE(SUM(MYZ,MYZ>0.0_SP)<1.0_SP)
             MYZ=MYZ + MYDIST*0.01
          END DO

          WHERE (MYZ > 0.0_SP)
             MYDIST = MYZ
          ELSEWHERE
             MYDIST = 0.0_SP
          END WHERE

          TOTAL = SUM(MYDIST)
          MYDIST = MYDIST/TOTAL

       ELSE

          CALL FATAL_ERROR("RIVER DISTRIBUTION: UNKOWN LINEAR SETTING?",&
               & "BAD STRING:"//TRIM(NAME))

       END IF


       WRITE(IPT,*) "LINEAR RIVER DISTRIBUTION",MYDIST
       
    CASE DEFAULT
       
       CALL FATAL_ERROR("UNKNOWN RIVER DISTRIBUTION FUNCTION:"//TRIM(NAME),&
            &"SEE FVCOM MANUAL OR mod_force.F FOR OPTIONS!")
       
    END SELECT
    

  END SUBROUTINE SET_DISTRIBUTION

  !================================================================
  !================================================================
  SUBROUTINE OBC_TEMPERATURE
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND

    INTEGER MYNOBC
    INTEGER, ALLOCATABLE :: MYOBCLIST(:)
    INTEGER :: MYSIGLAY

    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    INTEGER :: NTIMES, I
    TYPE(TIME) :: TIMETEST

    INTEGER :: STATUS


    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START OBC_TEMPERATURE"

    IF (.NOT. OBC_TEMP_NUDGING) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! OPEN BOUNDARY TEMPERATURE NUDGING IS OFF!"
       OBC_T_COMMENTS="OPEN BOUNDARY TEMPERATURE NUDGING IS OFF!"
       RETURN
    END IF

    IF(DBG_SET(DBG_LOG)) write(IPT,*) "! OPEN BOUNDARY TEMPERATURE NUDGING IS ON!"

    OBC_T_COMMENTS="OPEN BOUNDARY TEMPERATURE NUDGING IS ON!"
    
    ! FIND THE TIDAL FORCING FILE OBJECT
    OBC_T_FILE => FIND_FILE(FILEHEAD,trim(OBC_TEMP_FILE),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND OPEN BOUNDARY CONDITION TEMPERATURE FILE OBJECT",&
         & "FILE NAME: "//TRIM(OBC_TEMP_FILE))

    ATT => FIND_ATT(OBC_T_FILE,"type",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("IN OPEN BOUNDARY CONDITION TEMPERATURE FILE OBJECT",&
         & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
         &"COULD NOT FIND GLOBAL ATTRIBURE: 'type'")


    SELECT CASE(TRIM(ATT%CHR(1)))
       !=================================
       ! TIME SERIES OBC TEMPERATURE NUDGING DATA
    CASE("FVCOM TIME SERIES OBC TS FILE")
       !==================================         

       OBC_T_TYPE = OBC_T_SIGMA

       ATT => FIND_ATT(OBC_T_FILE,"title",FOUND)
       IF(FOUND) THEN
          OBC_T_COMMENTS = "Open Boundary Temperature Data: "&
               &//TRIM(ATT%CHR(1))
       ELSE
          CALL WARNING("ATTRIBUTE 'title' IS MISSING IN THE TEMPERATURE NUDGING FILE!")
          OBC_T_COMMENTS = "Open Boundary Temperature Data: UNKOWN SOURCE"
       END IF

       ! LOOK FOR THE DIMENSIONS 'nobc' and 'time'
       DIM => FIND_DIM(OBC_T_FILE,'time',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"COULD NOT FIND DIMENSION 'time'")

       NTIMES = DIM%DIM

       DIM => FIND_DIM(OBC_T_FILE,'siglay',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"COULD NOT FIND DIMENSION 'siglay'")

       MYSIGLAY = DIM%DIM

       if(KBM1 /= MYSIGLAY) CALL FATAL_ERROR&
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"THE 'siglay' DIMENSION DOES NOT MATCH THE MODEL RUN!")

       DIM => FIND_DIM(OBC_T_FILE,'nobc',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"COULD NOT FIND DIMENSION 'nobc'")

       if(IOBCN_GL>0) then
       IF(IOBCN_GL /= DIM%DIM) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"THE DIMENSION 'nobc' MUST MATCH THE NUMBER OF OBC NODES")
       endif

       ! lOAD GLOBAL OPEN BOUNDARY NOD NUMBER DATA AND COMPARE WITH
       ! OBC.DAT/RESTART FILE INPUT
       ALLOCATE(MYOBCLIST(IOBCN))
       VAR => FIND_VAR(OBC_T_FILE,'obc_nodes',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"COULD NOT FIND VARIABLE 'obc_nodes'")
       CALL NC_CONNECT_AVAR(VAR,MYOBCLIST)
       CALL NC_READ_VAR(VAR)

       DO I = 1, IOBCN

          IF(SERIAL) THEN
             IF(I_OBC_N(I) /= MYOBCLIST(I)) THEN
                write(IPT,*) "NLID(MYOBCLIST)= ", MYOBCLIST(I), "; I=",I
                write(IPT,*) "I_OBC_N= ", I_OBC_N(I), "; I=",I
                CALL FATAL_ERROR &
                     & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
                     & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
                     &"THE LIST OF BOUNDARY NODES DOES NOT MATCH")
             END IF
          ELSE
          END IF
       END DO

       ! LOAD TIME AND CHECK TO MAKE SURE THE TIME RANGE IS VALID

       TIMETEST = GET_FILE_TIME(OBC_T_FILE,1)

       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE TEMPERATURE TIME SERIES")

       TIMETEST = GET_FILE_TIME(OBC_T_FILE,NTIMES)

       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"THE MODEL RUN ENDS AFTER THE TEMPERATURE TIME SERIES")

       VAR => FIND_VAR(OBC_T_FILE,'obc_temp',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"COULD NOT FIND VARIABLE 'obc_temp'")

       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"COULD NOT FIND TEMP VARIRIABLE'S ATTRIBUTE 'units'")

       if(TRIM(ATT%CHR(1)) .NE. 'Celsius' .and. TRIM(ATT%CHR(1)) .NE. 'Celcius')  CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"TEMP VARIRIABLE ATTRIBUTE 'units' SHOULD BE 'Celsius'")



       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(IOBCN,KBM1), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN OBC_TEMPERATURE")
       OBC_T_N => reference_var(var)
       CALL NC_CONNECT_PVAR(OBC_T_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(IOBCN,KBM1), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN OBC_TEMPERATURE")
       OBC_T_P => reference_var(var)
       CALL NC_CONNECT_PVAR(OBC_T_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)


       !=====================================
       ! DEFAULT CASE IF GLOBAL ATTRIBUTES OF FILE ARE INCORRECT
    CASE DEFAULT
       !=====================================
       CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION TEMPERATURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_TEMP_FILE),&
            &"THE GLOBAL ATTRIBURE 'type' RETURNED UNKNOWN TYPE:",&
            & TRIM(ATT%CHR(1)))
    END SELECT

    ! afm 20151112 & EJA 20160921
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3.     
    OBC_T_N%curr_stkcnt = 0; OBC_T_P%curr_stkcnt = 0

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END OBC_TEMPERATURE"

  END SUBROUTINE OBC_TEMPERATURE
  !================================================================
  !================================================================
  SUBROUTINE OBC_SALINITY
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND

    INTEGER MYNOBC
    INTEGER, ALLOCATABLE :: MYOBCLIST(:)
    INTEGER :: MYSIGLAY

    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    INTEGER :: NTIMES, I
    TYPE(TIME) :: TIMETEST

    INTEGER :: STATUS


    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START OBC_SALINITY"

    IF (.NOT. OBC_SALT_NUDGING) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! OPEN BOUNDARY SALINITY NUDGING IS OFF!"
       OBC_S_COMMENTS="OPEN BOUNDARY SALINITY NUDGING IS OFF!"
       RETURN
    END IF

    IF(DBG_SET(DBG_LOG)) write(IPT,*) "! OPEN BOUNDARY SALINITY NUDGING IS ON!"
    OBC_S_COMMENTS="OPEN BOUNDARY SALINITY NUDGING IS ON!"

    ! FIND THE TIDAL FORCING FILE OBJECT
    OBC_S_FILE => FIND_FILE(FILEHEAD,trim(OBC_SALT_FILE),FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("COULD NOT FIND OPEN BOUNDARY CONDITION SALINITY FILE OBJECT",&
         & "FILE NAME: "//TRIM(OBC_SALT_FILE))

    ATT => FIND_ATT(OBC_S_FILE,"type",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("IN OPEN BOUNDARY CONDITION SALINITY FILE OBJECT",&
         & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
         &"COULD NOT FIND GLOBAL ATTRIBURE: 'type'")


    SELECT CASE(TRIM(ATT%CHR(1)))
       !=================================
       ! TIME SERIES OBC SALINITY NUDGING DATA
    CASE("FVCOM TIME SERIES OBC TS FILE")
       !==================================         

       OBC_S_TYPE = OBC_S_SIGMA

       ATT => FIND_ATT(OBC_S_FILE,"title",FOUND)
       IF(FOUND) THEN
          OBC_S_COMMENTS = "Open Boundary Salinity Data: "&
               &//TRIM(ATT%CHR(1))
       ELSE
          CALL WARNING("ATTRIBUTE 'title' IS MISSING IN THE SALINITY NUDGING FILE!")
          OBC_S_COMMENTS = "Open Boundary Salinity Data: UNKOWN SOURCE"
       END IF

       ! LOOK FOR THE DIMENSIONS 'nobc' and 'time'
       DIM => FIND_DIM(OBC_S_FILE,'time',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"COULD NOT FIND DIMENSION 'time'")

       NTIMES = DIM%DIM

       DIM => FIND_DIM(OBC_S_FILE,'siglay',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"COULD NOT FIND DIMENSION 'siglay'")

       MYSIGLAY = DIM%DIM

       if(KBM1 /= MYSIGLAY) CALL FATAL_ERROR&
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"THE 'siglay' DIMENSION DOES NOT MATCH THE MODEL RUN!")

       DIM => FIND_DIM(OBC_S_FILE,'nobc',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"COULD NOT FIND DIMENSION 'nobc'")

       if(IOBCN_GL>0) then
       IF(IOBCN_GL /= DIM%DIM) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"THE DIMENSION 'nobc' MUST MATCH THE NUMBER OF OBC NODES")
       endif

       ! lOAD GLOBAL OPEN BOUNDARY NOD NUMBER DATA AND COMPARE WITH
       ! OBC.DAT/RESTART FILE INPUT
       ALLOCATE(MYOBCLIST(IOBCN))
       VAR => FIND_VAR(OBC_S_FILE,'obc_nodes',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"COULD NOT FIND VARIABLE 'obc_nodes'")
       CALL NC_CONNECT_AVAR(VAR,MYOBCLIST)
       CALL NC_READ_VAR(VAR)

       DO I = 1, IOBCN

          IF(SERIAL) THEN
             IF(I_OBC_N(I) /= MYOBCLIST(I)) THEN
                write(IPT,*) "NLID(MYOBCLIST)= ", MYOBCLIST(I), "; I=",I
                write(IPT,*) "I_OBC_N= ", I_OBC_N(I), "; I=",I
                CALL FATAL_ERROR &
                     & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
                     & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
                     &"THE LIST OF BOUNDARY NODES DOES NOT MATCH")
             END IF
          ELSE
          END IF
       END DO

       ! LOAD TIME AND CHECK TO MAKE SURE THE TIME RANGE IS VALID

       TIMETEST = GET_FILE_TIME(OBC_S_FILE,1)

       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE SALINITY TIME SERIES")

       TIMETEST = GET_FILE_TIME(OBC_S_FILE,NTIMES)


       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"THE MODEL RUN ENDS AFTER THE SALINITY TIME SERIES")

       VAR => FIND_VAR(OBC_S_FILE,'obc_salinity',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"COULD NOT FIND VARIABLE 'obc_salinity'")

       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"COULD NOT FIND TEMP VARIRIABLE'S ATTRIBUTE 'units'")

       if(TRIM(ATT%CHR(1)) .NE. 'PSU')  CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY NUDGING FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"TEMP VARIRIABLE ATTRIBUTE 'units' SHOULD BE 'PSU'")



       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(IOBCN,KBM1), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN OBC_SALINITY")
       OBC_S_N => reference_var(var)
       CALL NC_CONNECT_PVAR(OBC_S_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(IOBCN,KBM1), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN OBC_SALINITY")
       OBC_S_P => reference_var(var)
       CALL NC_CONNECT_PVAR(OBC_S_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)


       !=====================================
       ! DEFAULT CASE IF GLOBAL ATTRIBUTES OF FILE ARE INCORRECT
    CASE DEFAULT
       !=====================================
       CALL FATAL_ERROR &
            & ("IN OPEN BOUNDARY CONDITION SALINITY FILE OBJECT",&
            & "FILE NAME: "//TRIM(OBC_SALT_FILE),&
            &"THE GLOBAL ATTRIBURE 'type' RETURNED UNKNOWN TYPE:",&
            & TRIM(ATT%CHR(1)))
    END SELECT

    ! afm 20150914 & EJA 20160921
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3. 
    OBC_S_N%curr_stkcnt = 0; OBC_S_P%curr_stkcnt = 0

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END OBC_SALINITY"

  END SUBROUTINE OBC_SALINITY
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  !================================================================

  ! CURRENTELY NOT IN USE! REPLACEED BY INTERP_BILINEAR WHICH IS A
  ! GENERAL INTERPOLATION SCHEME FOR CURVILINEAR COORDINATES

!!$  SUBROUTINE SET_FILE_INTERP_QUAD(NCF,INTP_N,INTP_C)
!!$    IMPLICIT NONE
!!$    TYPE(NCFILE), POINTER :: NCF
!!$    TYPE(INTERP_WEIGHTS),POINTER :: INTP_N
!!$    TYPE(INTERP_WEIGHTS),POINTER :: INTP_C
!!$
!!$    TYPE(NCATT), POINTER :: ATT
!!$    TYPE(NCDIM), POINTER :: DIM
!!$    TYPE(NCVAR), POINTER :: VAR
!!$
!!$    INTEGER :: LATS, LONS, I, Ntimes
!!$    REAL(SP), POINTER :: XLON(:,:),XLAT(:,:)
!!$    REAL(SP), POINTER :: HEATX(:,:),HEATY(:,:)
!!$    REAL(SP), POINTER :: TMP1(:),TMP2(:)
!!$
!!$    real(sp) :: rzero
!!$
!!$    LOGICAL :: FOUND
!!$
!!$    IF(.NOT. ASSOCIATED(NCF))CALL FATAL_ERROR&
!!$         & ("SET_FILE_INTERP: FILE OBJECT ARGUMENT IS NOT ASSOCIATED!") 
!!$
!!$    ! EITHER BOTH WEIGHTS MUST ALREADY BE SET OR NONE 
!!$    IF(ASSOCIATED(NCF%INTERP_N)) THEN
!!$       IF(ASSOCIATED(NCF%INTERP_C)) THEN
!!$          INTP_N => NCF%INTERP_N
!!$          INTP_C => NCF%INTERP_C
!!$          RETURN
!!$       ELSE
!!$          CALL PRINT_FILE(NCF)
!!$          CALL FATAL_ERROR("ONLY ONE INTERP POINTER IS ASSOCAITED IN THIS FILE",&
!!$               & "SET_FILE_INTERP: IS NOT PREPARED TO HANDLE THIS.")
!!$       END IF
!!$    ELSE
!!$       IF(ASSOCIATED(NCF%INTERP_C))THEN
!!$          CALL PRINT_FILE(NCF)
!!$          CALL FATAL_ERROR("ONLY ONE INTERP POINTER IS ASSOCAITED IN THIS FILE",&
!!$               & "SET_FILE_INTERP: IS NOT PREPARED TO HANDLE THIS.")
!!$       END IF
!!$
!!$    END IF
!!$
!!$    
!!$    ATT => FIND_ATT(NCF,'DX',FOUND)
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ( "SET_FILE_INTERP:",&
!!$         & "FILE NAME: "//TRIM(NCF%FNAME),&
!!$         & "COULD NOT FIND ATTRIBUTE 'DX'")
!!$    
!!$    rzero = att%flt(1)
!!$    
!!$    DIM => FIND_DIM(NCF,'south_north',FOUND)  
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ("SET_FILE_INTERP:",&
!!$         & "FILE NAME: "//TRIM(NCF%FNAME),&
!!$         & "COULD NOT FIND DIMENSION 'south_north'")
!!$    
!!$    LATS = DIM%DIM
!!$    
!!$    DIM => FIND_DIM(NCF,'west_east',FOUND)
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ("SET_FILE_INTERP:",&
!!$         & "FILE NAME: "//TRIM(NCF%FNAME),&
!!$         & "COULD NOT FIND DIMENSION 'west_east'")
!!$    LONS = DIM%DIM
!!$    
!!$    
!!$    ! GET THE INTERPOLATION COEFFICIENTS
!!$    ALLOCATE(XLON(lons,lats))
!!$    ALLOCATE(XLAT(lons,lats))
!!$    
!!$    VAR => FIND_VAR(NCF,"XLAT",FOUND)
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ("SET_FILE_INTERP:",&
!!$         & "FILE NAME: "//TRIM(NCF%FNAME),&
!!$         & "COULD NOT FIND VARIABLE 'XLAT'")
!!$    
!!$    CALL NC_CONNECT_PVAR(VAR,XLAT)
!!$    CALL NC_READ_VAR(VAR)
!!$    
!!$    
!!$    VAR => FIND_VAR(NCF,"XLONG",FOUND)
!!$    IF(.not. FOUND) CALL FATAL_ERROR &
!!$         & ("SET_FILE_INTERP:",&
!!$         & "FILE NAME: "//TRIM(NCF%FNAME),&
!!$         & "COULD NOT FIND VARIABLE 'XLONG'")
!!$    
!!$    CALL NC_CONNECT_PVAR(VAR,XLON)
!!$    CALL NC_READ_VAR(VAR)
!!$    
!!$# if !defined(SPHERICAL)
!!$    ALLOCATE(HEATX(lons,lats))
!!$    ALLOCATE(HEATY(lons,lats))
!!$    
!!$    IF (.NOT. USE_PROJ) CALL FATAL_ERROR('PROJ IS NEEDED TO USE T&
!!$         &HIS TYPE OF FORCING FILE IN CARTESIAN MODE:',&
!!$         & ' RECOMPILE WITH projection 4')
!!$    
!!$    CALL DEGREES2METERS(XLON,XLAT,PROJECTION_REFERENCE,HEATX,HEATY,lons,lats)
!!$    
!!$    DEALLOCATE(XLAT,XLON)
!!$# else
!!$    HEATX => XLON
!!$    HEATY => XLAT
!!$    
!!$    NULLIFY(XLON)
!!$    NULLIFY(XLAT)
!!$# endif
!!$    TMP1 => XM
!!$    TMP2 => YM
!!$    ALLOCATE(INTP_N)
!!$    CALL SETUP_INTERP_QUAD_P(HEATX,HEATY,TMP1,TMP2,INTP_N,rzero)
!!$    
!!$    TMP1 => XMC
!!$    TMP2 => YMC
!!$    
!!$    ALLOCATE(INTP_C)
!!$    CALL SETUP_INTERP_QUAD_P(HEATX,HEATY,TMP1,TMP2,INTP_C,rzero)
!!$    
!!$    ! THIS SHOULD DEALLOCATE HEATX,HEATY FOR NON SPHERICAL
!!$    ! THIS SHOULD DEALLOCATE HEATX,HEATY WHICH ARE POINTED AT
!!$    ! XLONS AND XLATS IN THE SHPERICAL CASE
!!$    DEALLOCATE(HEATX, HEATY)
!!$    
!!$    
!!$    NCF%INTERP_N => INTP_N
!!$    NCF%INTERP_C => INTP_C
!!$    
!!$  END SUBROUTINE SET_FILE_INTERP_QUAD
  !================================================================
  SUBROUTINE SET_FILE_INTERP_BILINEAR(NCF,INTP_N,INTP_C,MASK_VAR_NAME)
    IMPLICIT NONE
    TYPE(NCFILE), POINTER :: NCF
    TYPE(INTERP_WEIGHTS),POINTER :: INTP_N
    TYPE(INTERP_WEIGHTS),POINTER :: INTP_C

    TYPE(NCATT), POINTER :: ATT
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR

    INTEGER :: LATS, LONS, I, Ntimes,j, ierr
    REAL(SP), POINTER :: XLON(:,:),XLAT(:,:)
    REAL(SP), POINTER :: HEATX(:,:),HEATY(:,:)
    REAL(SP), POINTER :: TMP1(:),TMP2(:)


    CHARACTER(LEN=80), OPTIONAL :: MASK_VAR_NAME
    REAL(SP), POINTER :: FMASK(:,:)
    INTEGER, POINTER :: MASK(:,:)


    LOGICAL :: FOUND

    IF(.NOT. ASSOCIATED(NCF))CALL FATAL_ERROR&
         & ("SET_FILE_INTERP: FILE OBJECT ARGUMENT IS NOT ASSOCIATED!") 

    ! EITHER BOTH WEIGHTS MUST ALREADY BE SET OR NONE 
    IF(ASSOCIATED(NCF%INTERP_N)) THEN
       IF(ASSOCIATED(NCF%INTERP_C)) THEN
          INTP_N => NCF%INTERP_N
          INTP_C => NCF%INTERP_C
          RETURN
       ELSE
          CALL PRINT_FILE(NCF)
          CALL FATAL_ERROR("ONLY ONE INTERP POINTER IS ASSOCAITED IN THIS FILE",&
               & "SET_FILE_INTERP: IS NOT PREPARED TO HANDLE THIS.")
       END IF
    ELSE
       IF(ASSOCIATED(NCF%INTERP_C))THEN
          CALL PRINT_FILE(NCF)
          CALL FATAL_ERROR("ONLY ONE INTERP POINTER IS ASSOCAITED IN THIS FILE",&
               & "SET_FILE_INTERP: IS NOT PREPARED TO HANDLE THIS.")
       END IF

    END IF

        
    DIM => FIND_DIM(NCF,'south_north',FOUND)  
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("SET_FILE_INTERP:",&
         & "FILE NAME: "//TRIM(NCF%FNAME),&
         & "COULD NOT FIND DIMENSION 'south_north'")
    
    LATS = DIM%DIM
    
    DIM => FIND_DIM(NCF,'west_east',FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("SET_FILE_INTERP:",&
         & "FILE NAME: "//TRIM(NCF%FNAME),&
         & "COULD NOT FIND DIMENSION 'west_east'")
    LONS = DIM%DIM
    
    
    ! GET THE INTERPOLATION COEFFICIENTS
    ALLOCATE(XLON(lons,lats))
    ALLOCATE(XLAT(lons,lats))
   
    VAR => FIND_VAR(NCF,"XLAT",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("SET_FILE_INTERP:",&
         & "FILE NAME: "//TRIM(NCF%FNAME),&
         & "COULD NOT FIND VARIABLE 'XLAT'")
    
    CALL NC_CONNECT_PVAR(VAR,XLAT)
    CALL NC_READ_VAR(VAR)
    
    
    VAR => FIND_VAR(NCF,"XLONG",FOUND)
    IF(.not. FOUND) CALL FATAL_ERROR &
         & ("SET_FILE_INTERP:",&
         & "FILE NAME: "//TRIM(NCF%FNAME),&
         & "COULD NOT FIND VARIABLE 'XLONG'")
    
    CALL NC_CONNECT_PVAR(VAR,XLON)
    CALL NC_READ_VAR(VAR)
    
    ALLOCATE(HEATX(lons,lats))
    ALLOCATE(HEATY(lons,lats))
    
    IF (.NOT. USE_PROJ) CALL FATAL_ERROR('PROJ IS NEEDED TO USE T&
         &HIS TYPE OF FORCING FILE IN CARTESIAN MODE:',&
         & ' RECOMPILE WITH projection 4')
    IF(MSR) CALL DEGREES2METERS(XLON,XLAT,PROJECTION_REFERENCE,HEATX,HEATY,lons,lats)
    
    DEALLOCATE(XLAT,XLON)

!!$    THIS IS VERY SLOW - LOAD DATA FROM A FILE IF NEEDED
!!$    ! MAKE A LAND MASK
!!$    ALLOCATE(MASK(lons,lats))
!!$     DO I = 1,lons
!!$       DO J = 1,lats
!!$          MASK(i,j) = FIND_ELEMENT_CONTAINING(HEATX(i,j)-vxmin,HEATY(i,j)-vymin)
!!$       END DO
!!$    END DO
!!$    WHERE (MASK >0)
!!$       MASK = 0
!!$    elsewhere
!!$       MASK = 1
!!$    END WHERE


    IF (PRESENT(MASK_VAR_NAME))THEN
       VAR => FIND_VAR(NCF,MASK_VAR_NAME,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("SET_FILE_INTERP:",&
            & "FILE NAME: "//TRIM(NCF%FNAME),&
            & "COULD NOT FIND VARIABLE 'XLONG'")

       
       select case(VAR%XTYPE)
       case(NF90_INT)
          ALLOCATE(MASK(lons,lats))
          CALL NC_CONNECT_PVAR(VAR,MASK)
          CALL NC_READ_VAR(VAR)

       case(NF90_FLOAT)
          ALLOCATE(MASK(lons,lats))
          ALLOCATE(FMASK(lons,lats))
          CALL NC_CONNECT_PVAR(VAR,FMASK)
          CALL NC_READ_VAR(VAR)
          MASK = anint(FMASK)
          deallocate(Fmask)

       case default
          call fatal_error("SET_FILE_INTERP_BILINEAR: Unknown mask variable xtype?")
       END select

       TMP1 => VX
       TMP2 => VY

       TMP1 = TMP1 + VXmin
       TMP2 = TMP2 + VYmin
              
       ALLOCATE(INTP_N)
       CALL SETUP_INTERP_BILINEAR_P(HEATX,HEATY,TMP1,TMP2,INTP_N,MASK)

       TMP1 = TMP1 - VXmin
       TMP2 = TMP2 - VYmin
       
       TMP1 => XC
       TMP2 => YC

       TMP1 = TMP1 + VXmin
       TMP2 = TMP2 + VYmin
       
       ALLOCATE(INTP_C)
       CALL SETUP_INTERP_BILINEAR_P(HEATX,HEATY,TMP1,TMP2,INTP_C,MASK)

       TMP1 = TMP1 - VXmin
       TMP2 = TMP2 - VYmin
       
    ELSE


       TMP1 => VX
       TMP2 => VY

       TMP1 = TMP1 + VXmin
       TMP2 = TMP2 + VYmin

       !    NO MASK
       ALLOCATE(INTP_N)
       CALL SETUP_INTERP_BILINEAR_P(HEATX,HEATY,TMP1,TMP2,INTP_N)
       
       TMP1 = TMP1 - VXmin
       TMP2 = TMP2 - VYmin

       TMP1 => XC
       TMP2 => YC

       TMP1 = TMP1 + VXmin
       TMP2 = TMP2 + VYmin
       
       ALLOCATE(INTP_C)
       CALL SETUP_INTERP_BILINEAR_P(HEATX,HEATY,TMP1,TMP2,INTP_C)
       
       TMP1 = TMP1 - VXmin
       TMP2 = TMP2 - VYmin
       
    END IF
    
    ! THIS SHOULD DEALLOCATE HEATX,HEATY FOR NON SPHERICAL
    ! THIS SHOULD DEALLOCATE HEATX,HEATY WHICH ARE POINTED AT
    ! XLONS AND XLATS IN THE SHPERICAL CASE
    DEALLOCATE(HEATX, HEATY)
    
    
    NCF%INTERP_N => INTP_N
    NCF%INTERP_C => INTP_C
    
  END SUBROUTINE SET_FILE_INTERP_BILINEAR
  !================================================================
  !================================================================
  !================================================================
  !================================================================
  SUBROUTINE SURFACE_HEATING_CALCULATED
    USE MOD_HEATFLUX
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND

    REAL(SP), POINTER :: STORAGE_ARR(:,:), STORAGE_VEC(:)
    CHARACTER(LEN=60) :: t_airstrng, rh_airstrng, pa_airstrng, dlw_airstrng, dsw_airstrng
    TYPE(TIME) :: TIMETEST

    INTEGER :: LATS, LONS, I, Ntimes

    INTEGER :: STATUS
    REAL(SP) :: TEMP

    CHARACTER(LEN=80)  :: ISTR
    ISTR = "./"//TRIM(INPUT_DIR)//"/"//trim(casename)
  
    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START SURFACE_HEATING_CALCULATED"

    NULLIFY(ATT,DIM,VAR,STORAGE_ARR,STORAGE_VEC)

    IF (.NOT. HEATING_CALCULATE_ON ) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! SURFACE HEAT FORCING IS OFF!"
       ALLOCATE(HEAT_CALCULATE_COMMENTS(1))
       HEAT_CALCULATE_COMMENTS(1) = "SURFACE HEAT FORCING IS OFF"
       RETURN
    END IF

!------------------------------------------------------------------------------|
!--------------READ IN LATITUDE------------------------------------------------!

   ALLOCATE(CORRG(0:MGL))  ; CORRG = 0.0_SP
   !<-------ishid dbg 20230104
   !CALL FOPEN(CORIOLISUNIT, TRIM(ISTR)//'_cor.dat',"cfr")
   !CALL FOPEN(CORIOLISUNIT, TRIM(casename)//'_cor.dat',"cfr")
   Call FOPEN(CORIOLISUNIT,trim(INPUT_DIR)//trim(CORIOLIS_FILE),'cfr')
   !>-------ishid dbg 20230104
   REWIND(CORIOLISUNIT)
   READ(CORIOLISUNIT,*)
   DO I=1,MGL
     READ(CORIOLISUNIT,*) TEMP,TEMP,CORRG(I)
   END DO
   CLOSE(CORIOLISUNIT)

!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
   ALLOCATE(CORR(0:MT)) ; CORR = 0.0_SP
   IF(SERIAL) CORR = CORRG

   DEALLOCATE(CORRG)
   
   IF(TRIM(COARE_VERSION) == 'COARE40VN')THEN
     ALLOCATE(RAIN(0:MT)) ; RAIN = 0.0_SP
     ALLOCATE(CP40(0:MT)) ; CP40 = 0.0_SP
     ALLOCATE(SIGH(0:MT)) ; SIGH = 0.0_SP
     ALLOCATE(ZI40(0:MT)) ; ZI40 = 0.0_SP
   END IF   
   
   IF(HEATING_FRESHWATER)THEN
     ALLOCATE(USRCOARE(0:MT)) ; USRCOARE = 0.0_SP
   END IF  
!----------------------------REPORT--------------------------------------------!
   IF(MSR)WRITE(IPT,*)'!'
   IF(MSR)WRITE(IPT,*)'!            SETTING UP PRESCRIBED BOUNDARY CONDITIONS  '
   IF(MSR)WRITE(IPT,*)'!'

!==============================================================================|
!   Input Meteorological Boundary Conditions for Calculating Heat Flux         |
!==============================================================================|
!    bulk air temperature at height 2m:   degree(C)    "t_air"                 |
!    relative humidity at height 2m:      (%)          "rh_air"                |
!---> Siqi li, 2021-01-27
!    surface pressure:                    pa           "pa_air"                |
!    surface pressure:                    mb           "pa_air"                |  (the unit was wrong for coare)  
!<--- Siqi Li, 2021-01-27
!    downward longwave radiation:         w/m^2        "dlw_air"               |
!    downward shortwave radiation:        w/m^2        "dsw_air"               |
!==============================================================================|

    ! DETERMINE HOW TO LOAD THE DATA
    SELECT CASE(HEATING_CALCULATE_KIND)
    CASE (CNSTNT)
       write(t_airstrng,'(f8.4)')   AIR_TEMPERATURE
       write(rh_airstrng,'(f8.4)')  RELATIVE_HUMIDITY
       write(pa_airstrng,'(f8.4)')  SURFACE_PRESSURE
       write(dlw_airstrng,'(f8.4)') LONGWAVE_RADIATION
       write(dsw_airstrng,'(f8.4)') SHORTWAVE_RADIATION

       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*)"! SETTING UP CONSTANT HEAT FORCING: "
          WRITE(IPT,*)"         Bulk Air Temperature: "//trim(t_airstrng)
          WRITE(IPT,*)"            Relative Hunidity: "//trim(rh_airstrng)
          WRITE(IPT,*)"             Surface Pressure: "//trim(pa_airstrng)
          WRITE(IPT,*)"  Downward Longwave Radiation: "//trim(dlw_airstrng)
          WRITE(IPT,*)" Downward shortwave Radiation: "//trim(dsw_airstrng)
       END IF

       ALLOCATE(HEAT_CALCULATE_COMMENTS(6))
       HEAT_CALCULATE_COMMENTS(1) = "Using constant heating from run file:"
       HEAT_CALCULATE_COMMENTS(2) = "Bulk Air Temperature:"//trim(t_airstrng)
       HEAT_CALCULATE_COMMENTS(3) = "Relative Humidity:"//trim(rh_airstrng)
       HEAT_CALCULATE_COMMENTS(4) = "Surface Pressure:"//trim(pa_airstrng)
       HEAT_CALCULATE_COMMENTS(5) = "Downward Longwave Radiation:"//trim(dlw_airstrng)
       HEAT_CALCULATE_COMMENTS(6) = "Downward Shortwave Radiation:"//trim(dsw_airstrng)

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC HEATING Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT HEATING Not Set Up Yet")

    CASE(PRDC)

       HEAT_FILE => FIND_FILE(FILEHEAD,trim(HEATING_CALCULATE_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE HEATING BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE
       ATT => FIND_ATT(HEAT_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(HEAT_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          HEAT_FORCING_TYPE = HEAT_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          HEAT_FORCING_TYPE = HEAT_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          HEAT_FORCING_TYPE = HEAT_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          HEAT_FORCING_TYPE = HEAT_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(HEAT_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE HEATING FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(HEAT_FORCING_COMMENTS(4))
       HEAT_FORCING_COMMENTS(1) = "FVCOM periodic surface heat forcing:"
       HEAT_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(HEATING_CALCULATE_FILE)

       HEAT_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(HEAT_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          HEAT_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          HEAT_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF

       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(HEAT_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       HEAT_PERIOD = get_file_time(HEAT_FILE,ntimes)

       HEAT_PERIOD = HEAT_PERIOD - get_file_time(HEAT_FILE,1)


       IF (HEAT_PERIOD /= get_file_time(HEAT_FILE,ntimes)) THEN

          CALL PRINT_REAL_TIME(get_file_time(HEAT_FILE,1),IPT,"FIRST FILE TIME",TIMEZONE)
          CALL PRINT_REAL_TIME(get_file_time(HEAT_FILE,ntimes),IPT,"LAST FILE TIME",TIMEZONE)

          CALL FATAL_ERROR&
               &("TO USE PERIODIC FORCING THE FILE TIME MUST COUNT FROM ZERO",&
               & "THE DIFFERENCE BETWEEN THE CURRENT MODEL TIME AND THE START TIME,",&
               & "MODULO THE FORCING PERIOD, DETERMINES THE CURRENT FORCING")
       END IF


       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! USING PERIODIC HEAT FORCING:"
          CALL PRINT_TIME(HEAT_PERIOD,IPT,"PERIOD")
       END IF

    CASE(VRBL)

       HEAT_FILE => FIND_FILE(FILEHEAD,trim(HEATING_CALCULATE_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE HEATING BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE
       ATT => FIND_ATT(HEAT_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(HEAT_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          HEAT_FORCING_TYPE = HEAT_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          HEAT_FORCING_TYPE = HEAT_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          HEAT_FORCING_TYPE = HEAT_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          HEAT_FORCING_TYPE = HEAT_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(HEAT_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE HEATING FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(HEAT_FORCING_COMMENTS(4))
       HEAT_FORCING_COMMENTS(1) = "FVCOM variable surface heat forcing file:"
 
       HEAT_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(HEATING_CALCULATE_FILE)
       HEAT_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(HEAT_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          HEAT_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          HEAT_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF

       ! CHECK DIMENSIONS
       DIM => FIND_UNLIMITED(HEAT_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            &"COULD NOT FIND UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       ! CHECK THE FILE TIME AND COMPARE WITH MODEL RUN TIME
       TIMETEST = get_file_time(HEAT_FILE,1)
       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE FORCING TIME SERIES")

       TIMETEST = get_file_time(HEAT_FILE,ntimes)
       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            &"THE MODEL RUN ENDS AFTER THE FORCING TIME SERIES")

    CASE DEFAULT
       CALL FATAL_ERROR("SURFACE_HEATING: UNKNOWN HEATING KIND?")
    END SELECT



    !==================================================================
    SELECT CASE(HEAT_FORCING_TYPE)
       !==================================================================
    CASE(HEAT_IS_WRFGRID)
       !==================================================================

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP HEAT FORCING FROM A 'wrf grid' FILE"

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_DIM(HEAT_FILE,'south_north',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND DIMENSION 'south_north'")

       LATS = DIM%DIM

       DIM => FIND_DIM(HEAT_FILE,'west_east',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND DIMENSION 'west_east'")
       LONS = DIM%DIM


       CALL SET_FILE_INTERP_BILINEAR(HEAT_FILE,HEAT_INTP_N,HEAT_INTP_C)

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

!
       ! BULK AIR TEMPERATURE DATA
       VAR => FIND_VAR(HEAT_FILE,"air_temperature",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(HEAT_FILE,"air_temperature",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'air_temperature'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       T_AIR_N => reference_var(var)
       CALL NC_CONNECT_PVAR(T_AIR_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(T_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       T_AIR_P => reference_var(var)
       CALL NC_CONNECT_PVAR(T_AIR_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(T_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!
       ! RELATIVE HUMIDITY
       VAR => FIND_VAR(HEAT_FILE,"relative_humidity",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(HEAT_FILE,"relative_humidity",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'relative_humidity'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       RH_AIR_N => reference_var(var)
       CALL NC_CONNECT_PVAR(RH_AIR_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(RH_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       RH_AIR_P => reference_var(var)
       CALL NC_CONNECT_PVAR(RH_AIR_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(RH_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!       
       ! SURFACE PRESSURE
       VAR => FIND_VAR(HEAT_FILE,"air_pressure",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(HEAT_FILE,"SLP",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'air_pressure' OR 'SLP'")

       !---> Siqi Li, 2021-01-27
       ! The unit of air pressure is really important. It can influence both
       ! water elevation and wind stress. In calculation, the unit of air 
       ! pressure is 'Pa'.
       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(VAR,'unit',FOUND)
       IF (FOUND) THEN
       call UPCASE(ATT%CHR(1))
       SELECT CASE (TRIM(ATT%CHR(1)))
       CASE ('MB', 'HPA')
         Pair_unit_factor = 100.0_SP
       CASE ('PA','PASCAL')
         Pair_unit_factor = 1.0_SP
       CASE DEFAULT
         CALL FATAL_ERROR("UNKNOWN UNIT of AIR PRESSURE: " &
                 & //TRIM(ATT%CHR(1)))
       END SELECT
       ELSE
         ! There is no unit attribute for the air_pressure variable
         ! We assume the unit is 'Pa'.
         Pair_unit_factor = 1.0_SP
       END IF
       !<--- Siqi Li, 2021-01-27

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       PA_AIR_N => reference_var(var)
       CALL NC_CONNECT_PVAR(PA_AIR_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(PA_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       PA_AIR_P => reference_var(var)
       CALL NC_CONNECT_PVAR(PA_AIR_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(PA_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!       
       ! DOWNWARD LONGWAVE RADIATION
       VAR => FIND_VAR(HEAT_FILE,"long_wave",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(HEAT_FILE,"Longwave",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'long_wave'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       DLW_AIR_N => reference_var(var)
       CALL NC_CONNECT_PVAR(DLW_AIR_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(DLW_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       DLW_AIR_P => reference_var(var)
       CALL NC_CONNECT_PVAR(DLW_AIR_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(DLW_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!       
       ! DOWNWARD SHORTWAVE RADIATION
       VAR => FIND_VAR(HEAT_FILE,"short_wave",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(HEAT_FILE,"Shortwave",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'short_wave'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       DSW_AIR_N => reference_var(var)
       CALL NC_CONNECT_PVAR(DSW_AIR_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(DSW_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       DSW_AIR_P => reference_var(var)
       CALL NC_CONNECT_PVAR(DSW_AIR_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(DSW_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!       

       !==================================================================
    CASE(HEAT_IS_FVCOMGRID)
       !==================================================================

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP HEAT FORCING FROM A 'fvcom grid' FILE"

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_DIM(HEAT_FILE,'node',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND DIMENSION 'node'")

       if (mgl /= dim%dim) CALL FATAL_ERROR&
            &("Surface Heating: the number of nodes in the file does not match the fvcom grid?")


       DIM => FIND_DIM(HEAT_FILE,'nele',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND DIMENSION 'nele'")

       if (ngl /= dim%dim) CALL FATAL_ERROR&
            &("Surface Heating: the number of elements in the file does not match the fvcom grid?")

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

!
       ! BULK AIR TEMPERATURE DATA
       VAR => FIND_VAR(HEAT_FILE,"air_temperature",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'air_temperature'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       T_AIR_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(T_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       T_AIR_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(T_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!
       ! RELATIVE HUMIDITY DATA
       VAR => FIND_VAR(HEAT_FILE,"relative_humidity",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'relative_humidity'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       RH_AIR_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(RH_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       RH_AIR_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(RH_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!
       ! SURFACE PRESSURE
       VAR => FIND_VAR(HEAT_FILE,"air_pressure",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(HEAT_FILE,"SLP",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'air_pressure' OR 'SLP'") ! Siqi Li, 2021-01-27
!            & "COULD NOT FIND VARIABLE 'air_pressure' OR "SLP")

       !---> Siqi Li, 2021-01-27
       ! The unit of air pressure is really important. It can influence both
       ! water elevation and wind stress. In calculation, the unit of air 
       ! pressure is 'Pa'.
       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(VAR,'unit',FOUND)
       IF (FOUND) THEN
       call UPCASE(ATT%CHR(1))
       SELECT CASE (TRIM(ATT%CHR(1)))
       CASE ('MB', 'HPA')
         Pair_unit_factor = 100.0_SP
       CASE ('PA','PASCAL')
         Pair_unit_factor = 1.0_SP
       CASE DEFAULT
         CALL FATAL_ERROR("UNKNOWN UNIT of AIR PRESSURE: " &
                 & //TRIM(ATT%CHR(1)))
       END SELECT
       ELSE
         ! There is no unit attribute for the air_pressure variable
         ! We assume the unit is 'Pa'.
         Pair_unit_factor = 1.0_SP
       END IF
       !<--- Siqi Li, 2021-01-27

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       PA_AIR_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(PA_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       PA_AIR_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(PA_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!
       ! DOWNWARD LONGWAVE RADIATION
       VAR => FIND_VAR(HEAT_FILE,"long_wave",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'long_wave'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       DLW_AIR_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(DLW_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       DLW_AIR_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(DLW_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!
       ! DOWNWARD SHORTWAVE RADIATION
       VAR => FIND_VAR(HEAT_FILE,"short_wave",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(HEATING_CALCULATE_FILE),&
            & "COULD NOT FIND VARIABLE 'short_wave'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       DSW_AIR_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(DSW_AIR_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       DSW_AIR_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE HEATING")
       CALL NC_CONNECT_PVAR(DSW_AIR_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
!
       !==================================================================
    CASE DEFAULT
       !==================================================================
       CALL FATAL_ERROR("CAN NOT RECOGNIZE HEATING FILE TYPE!")
       !==================================================================
    END SELECT
    !==================================================================

    ! afm 20151112 & EJA 20160921
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3. 
    t_air_n%curr_stkcnt   = 0;t_air_p%curr_stkcnt   = 0
    rh_air_n%curr_stkcnt  = 0;rh_air_p%curr_stkcnt  = 0
    pa_air_n%curr_stkcnt  = 0;pa_air_p%curr_stkcnt  = 0
    dlw_air_n%curr_stkcnt = 0;dlw_air_p%curr_stkcnt = 0
    dsw_air_n%curr_stkcnt = 0;dsw_air_p%curr_stkcnt = 0

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END SURFACE_HEATING_CALCULATED"
  END SUBROUTINE SURFACE_HEATING_CALCULATED
!========================================================================
!========================================================================
  !================================================================
  !================================================================
  SUBROUTINE ICE_MODEL_FORCING
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND

    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    CHARACTER(len=60) :: SATstrng, SLPstrng,SPQstrng,CLDstrng, SWVstrng
    TYPE(TIME) :: TIMETEST

    INTEGER :: LATS, LONS, I, Ntimes

    INTEGER :: STATUS

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START ICE MODEL FORCING"

    NULLIFY(ATT,DIM,VAR,STORAGE_ARR,STORAGE_VEC)

    IF (.NOT. ICE_MODEL) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! ICE MODEL IS OFF!"

       ALLOCATE(ICE_FORCING_COMMENTS(1))
       ICE_FORCING_COMMENTS(1) = "ICE MODEL FORCING IS OFF"
       RETURN
    END IF


    ! DETERMINE HOW TO LOAD THE DATA
    SELECT CASE(ICE_FORCING_KIND)
    CASE (CNSTNT)

       
       ALLOCATE(ICE_FORCING_COMMENTS(6))
       write(SATstrng,'(f8.4)') ICE_AIR_TEMP
       write(SPQstrng,'(f8.4)') ICE_SPEC_HUMIDITY
       write(CLDstrng,'(f8.4)') ICE_CLOUD_COVER
       write(SWVstrng,'(f8.4)') ICE_SHORTWAVE


       ICE_FORCING_COMMENTS(1) = "Using constant ice forcing:"
       ICE_FORCING_COMMENTS(2) = "Sea Leval Air Temp="//trim(SATstrng)
       ICE_FORCING_COMMENTS(4) = "Specific Humidity="//trim(SPQstrng)
       ICE_FORCING_COMMENTS(5) = "Cloud Cover="//trim(CLDstrng)
       ICE_FORCING_COMMENTS(6) = "Shortwave Radiation="//trim(SWVstrng)

       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*)"! SETTING UP CONSTANT ICE FORCING:"
          WRITE(IPT,*)"!    Sea Leval Air Temp="//trim(SATstrng)
          WRITE(IPT,*)"!    Specific Humidity="//trim(SPQstrng)
          WRITE(IPT,*)"!    Cloud Cover="//trim(CLDstrng)
          WRITE(IPT,*)"!    Shortwave Radiation="//trim(SWVstrng)
       END IF

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC ICE FORCING Not Set Up Yet")
       !      HEAT_FORCING_COMMENTS = "Using Static heating from file"

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT ICE FORCING Not Set Up Yet")
       !       HEAT_FORCING_COMMENTS = "Using TIME DEPENDENT heating from file"

    CASE(PRDC)

       ICE_FILE => FIND_FILE(FILEHEAD,trim(ICE_FORCING_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND ICE MODEL BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE))


       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE
       ATT => FIND_ATT(ICE_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(ICE_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          ICE_FORCING_TYPE = ICE_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          ICE_FORCING_TYPE = ICE_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          ICE_FORCING_TYPE = ICE_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          ICE_FORCING_TYPE = ICE_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(ICE_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE ICE FORCING FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(ICE_FORCING_COMMENTS(4))
       ICE_FORCING_COMMENTS(1) = "FVCOM periodic surface ice model forcing:"
       ICE_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(ICE_FORCING_FILE)

       ICE_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(ICE_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          ICE_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          ICE_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF

       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(ICE_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       ICE_PERIOD = get_file_time(ICE_FILE,ntimes)

       ICE_PERIOD = ICE_PERIOD - get_file_time(ICE_FILE,1)


       IF (ICE_PERIOD /= get_file_time(ICE_FILE,ntimes)) THEN

          CALL PRINT_TIME(get_file_time(ICE_FILE,1),IPT,"FIRST FILE TIME")
          CALL PRINT_TIME(get_file_time(ICE_FILE,ntimes),IPT,"LAST FILE TIME")

          CALL FATAL_ERROR&
               &("TO USE PERIODIC FORCING THE FILE TIME MUST COUNT FROM ZERO",&
               & "THE DIFFERENCE BETWEEN THE CURRENT MODEL TIME AND THE START TIME,",&
               & "MODULO THE FORCING PERIOD, DETERMINES THE CURRENT FORCING")
       END IF


       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! USING PERIODIC ICE FORCING:"
          CALL PRINT_TIME(ICE_PERIOD,IPT,"PERIOD")
       END IF

    CASE(VRBL)

       ICE_FILE => FIND_FILE(FILEHEAD,trim(ICE_FORCING_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND ICE FORCING BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE
       ATT => FIND_ATT(ICE_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(ICE_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          ICE_FORCING_TYPE = ICE_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          ICE_FORCING_TYPE = ICE_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          ICE_FORCING_TYPE = ICE_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          ICE_FORCING_TYPE = ICE_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(ICE_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE ICE FORCING FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE


       ALLOCATE(ICE_FORCING_COMMENTS(4))
       ICE_FORCING_COMMENTS(1) = "FVCOM variable surface ice model forcing:"
       ICE_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(ICE_FORCING_FILE)

       ICE_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(ICE_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          ICE_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          ICE_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF


       ! CHECK DIMENSIONS
       DIM => FIND_UNLIMITED(ICE_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            &"COULD NOT FIND UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       ! CHECK THE FILE TIME AND COMPARE WITH MODEL RUN TIME
       TIMETEST = get_file_time(ICE_FILE,1)
       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE FORCING TIME SERIES")

       TIMETEST = get_file_time(ICE_FILE,ntimes)
       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            &"THE MODEL RUN ENDS AFTER THE FORCING TIME SERIES")

    CASE DEFAULT
       CALL FATAL_ERROR("ICE FORCING: UNKNOWN ICE_FORCING KIND?")
    END SELECT



    !==================================================================
    SELECT CASE(ICE_FORCING_TYPE)
       !==================================================================
    CASE(ICE_IS_WRFGRID)
       !==================================================================

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP ICE FORCING FROM A 'wrf grid' FILE"


       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_DIM(ICE_FILE,'south_north',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND DIMENSION 'south_north'")

       LATS = DIM%DIM

       DIM => FIND_DIM(ICE_FILE,'west_east',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND DIMENSION 'west_east'")
       LONS = DIM%DIM


       CALL SET_FILE_INTERP_BILINEAR(ICE_FILE,ICE_INTP_N,ICE_INTP_C)
       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       !IF(ASSOCIATED(ICE_FILE,HEAT_FILE)) THEN 
       VAR => FIND_VAR(HEAT_FILE,"short_wave",FOUND)
       IF( FOUND ) THEN
          ! USE THE SAME MEMORY USED FOR OCEAN MODEL HEAT FLUX
          
          ICE_SWV_N => HEAT_SWV_N

          ICE_SWV_P => HEAT_SWV_P

       ELSE
          ! LOAD YOUR OWN DATA FOR THE ICE MODEL

          ! SHORT WAVE RADIATION DATA
          VAR => FIND_VAR(ICE_FILE,"short_wave",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(ICE_FILE,"Shortwave",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
               & "COULD NOT FIND VARIABLE 'short_wave'")
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
          ICE_SWV_N => reference_var(var)
          CALL NC_CONNECT_PVAR(ICE_SWV_N,STORAGE_ARR)
          NULLIFY(STORAGE_ARR)
          
          ! MAKE SPACE FOR THE INTERPOLATED DATA
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
          CALL NC_CONNECT_PVAR(ICE_SWV_N,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
          
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
          ICE_SWV_P => reference_var(var)
          CALL NC_CONNECT_PVAR(ICE_SWV_P,STORAGE_ARR)
          NULLIFY(STORAGE_ARR)
          
          ! MAKE SPACE FOR THE INTERPOLATED DATA
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
          CALL NC_CONNECT_PVAR(ICE_SWV_P,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)

       END IF


       ! Surface Air Temperature DATA
       VAR => FIND_VAR(ICE_FILE,"SAT",FOUND)
!       VAR => FIND_VAR(ICE_FILE,"T2",FOUND)
! afm 20150930 & EJA 20160921
       IF(.not. FOUND) VAR => FIND_VAR(ICE_FILE,"air_temperature",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(ICE_FILE,"T2",FOUND)

       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'T2','air_temperature' or 'SAT'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       ICE_SAT_N => reference_var(var)
       CALL NC_CONNECT_PVAR(ICE_SAT_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_SAT_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       ICE_SAT_P => reference_var(var)
       CALL NC_CONNECT_PVAR(ICE_SAT_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_SAT_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! Specific Humidity DATA
       VAR => FIND_VAR(ICE_FILE,"SPQ",FOUND)
!       VAR => FIND_VAR(ICE_FILE,"Q2",FOUND)
! afm 20160513 & EJA 20160921
       IF(.not. FOUND) VAR => FIND_VAR(ICE_FILE,"Q2",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'Q2' or 'SPQ'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       ICE_SPQ_N => reference_var(var)
       CALL NC_CONNECT_PVAR(ICE_SPQ_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_SPQ_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       ICE_SPQ_P => reference_var(var)
       CALL NC_CONNECT_PVAR(ICE_SPQ_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_SPQ_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! Cloud Cover DATA
       VAR => FIND_VAR(ICE_FILE,"cloud_cover",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'cloud_cover'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       ICE_CLD_N => reference_var(var)
       CALL NC_CONNECT_PVAR(ICE_CLD_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_CLD_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       ICE_CLD_P => reference_var(var)
       CALL NC_CONNECT_PVAR(ICE_CLD_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_CLD_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       !==================================================================
    CASE(ICE_IS_FVCOMGRID)
       !==================================================================

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP HEAT FORCING FROM A 'fvcom grid' FILE"

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_DIM(ICE_FILE,'node',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND DIMENSION 'node'")

       if (mgl /= dim%dim) CALL FATAL_ERROR&
            &("Ice Forcing: the number of nodes in the file does not match the fvcom grid?")


       DIM => FIND_DIM(ICE_FILE,'nele',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND DIMENSION 'nele'")

       if (ngl /= dim%dim) CALL FATAL_ERROR&
            &("Ice Forcing: the number of elements in the file does not match the fvcom grid?")

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       ! SHORT WAVE RADIATION DATA
       !IF(ASSOCIATED(ICE_FILE,HEAT_FILE)) THEN
! afm 20151112 &  & EJA 20160921
! commented out for Solar ----------
! With SOLAR, use SOLAR-derived shortwave
! Without SOLAR, use shortwave from forcing data  
       VAR => FIND_VAR(HEAT_FILE,"short_wave",FOUND)
       IF( FOUND ) THEN
          ! USE THE SAME MEMORY USED FOR OCEAN MODEL HEAT FLUX
          
          ICE_SWV_N => HEAT_SWV_N

          ICE_SWV_P => HEAT_SWV_P

       ELSE
          ! LOAD YOUR OWN DATA FOR THE ICE MODEL

          ! SHORT WAVE RADIATION DATA
          VAR => FIND_VAR(ICE_FILE,"short_wave",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(ICE_FILE,"Shortwave",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
               & "COULD NOT FIND VARIABLE 'short_wave'")
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          ICE_SWV_N => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
          CALL NC_CONNECT_PVAR(ICE_SWV_N,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
          
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          ICE_SWV_P => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
          CALL NC_CONNECT_PVAR(ICE_SWV_P,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)

       END IF



       ! Surface Air Temperature DATA
       VAR => FIND_VAR(ICE_FILE,"SAT",FOUND)
! afm 20151112 & EJA 20160921
       IF(.not. FOUND) VAR => FIND_VAR(ICE_FILE,"air_temperature",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(ICE_FILE,"T2",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'T2','air_temperature' or 'SAT'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ICE_SAT_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_SAT_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ICE_SAT_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_SAT_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! Specific Humidity DATA
       VAR => FIND_VAR(ICE_FILE,"SPQ",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(ICE_FILE,"Q2",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'Q2' or 'SPQ'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ICE_SPQ_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_SPQ_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ICE_SPQ_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_SPQ_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! Specific Humidity DATA
       VAR => FIND_VAR(ICE_FILE,"cloud_cover",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN ICE FORCING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICE_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'cloud_cover'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ICE_CLD_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_CLD_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ICE_CLD_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICE FORCING")
       CALL NC_CONNECT_PVAR(ICE_CLD_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       !==================================================================
    CASE DEFAULT
       !==================================================================
       CALL FATAL_ERROR("CAN NOT RECOGNIZE ICE FORCING FILE TYPE!")
       !==================================================================
    END SELECT
    !==================================================================

    ! afm 20151112 & EJA 20160921
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3. 
! afm 20180717
    ice_swv_n%curr_stkcnt   = 0;ice_swv_p%curr_stkcnt   = 0
    ice_sat_n%curr_stkcnt   = 0;ice_sat_p%curr_stkcnt   = 0
    ice_spq_n%curr_stkcnt   = 0;ice_spq_p%curr_stkcnt   = 0
    ice_cld_n%curr_stkcnt   = 0;ice_cld_p%curr_stkcnt   = 0

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END ICE MODEL FORCING"
  END SUBROUTINE ICE_MODEL_FORCING
  !================================================================
  !================================================================
  SUBROUTINE ICING_FORCING
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND

    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    CHARACTER(len=30) :: SATstrng, WSPDstrng
    TYPE(TIME) :: TIMETEST

    INTEGER :: LATS, LONS, I, Ntimes

    INTEGER :: STATUS

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START ICING_FORCING"

    NULLIFY(ATT,DIM,VAR,STORAGE_ARR,STORAGE_VEC)

    IF (.NOT. ICING_MODEL ) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! ICING MODEL IS OFF!"
       ALLOCATE(ICING_FORCING_COMMENTS(1))
       ICING_FORCING_COMMENTS(1) = "ICING MODEL IS OFF"
       RETURN
    END IF


    ! DETERMINE HOW TO LOAD THE DATA
    SELECT CASE(ICING_FORCING_KIND)
    CASE (CNSTNT)

       write(SATstrng,'(f8.4)')  ICING_AIR_TEMP
       write(WSPDstrng,'(f8.4)') ICING_WSPD

       ALLOCATE(ICING_FORCING_COMMENTS(3))
       ICING_FORCING_COMMENTS(1) = "Using Constant heating:"
       
       ICING_FORCING_COMMENTS(2) = "Sea Level Air Temperature:"//trim(SATstrng)
       ICING_FORCING_COMMENTS(3) = "Wind Speed:"//trim(wspdstrng)
       
       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! SETTING UP CONSTANT ICING: "
          WRITE(IPT,*) "!    Sea Level Air Temperature:"//trim(SATstrng)
          WRITE(IPT,*) "!    Wind Speed:"//trim(wspdstrng)
       END IF
       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC HEATING Not Set Up Yet")
       !      HEAT_FORCING_COMMENTS = "Using Static heating from file"

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT HEATING Not Set Up Yet")
       !       HEAT_FORCING_COMMENTS = "Using TIME DEPENDENT heating from file"

    CASE(PRDC)

       ICING_FILE => FIND_FILE(FILEHEAD,trim(ICING_FORCING_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE ICING BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE
       ATT => FIND_ATT(ICING_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(ICING_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          ICING_FORCING_TYPE = ICING_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          ICING_FORCING_TYPE = ICING_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          ICING_FORCING_TYPE = ICING_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          ICING_FORCING_TYPE = ICING_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(ICING_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE ICING FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(ICING_FORCING_COMMENTS(4))
       ICING_FORCING_COMMENTS(1) = "FVCOM periodic surface icing forcing:"
       ICING_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(ICING_FORCING_FILE)

       ICING_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(ICING_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          ICING_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          ICING_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF

       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(ICING_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       ICING_PERIOD = get_file_time(ICING_FILE,ntimes)

       ICING_PERIOD = ICING_PERIOD - get_file_time(ICING_FILE,1)


       IF (ICING_PERIOD /= get_file_time(ICING_FILE,ntimes)) THEN

          CALL PRINT_REAL_TIME(get_file_time(ICING_FILE,1),IPT,"FIRST FILE TIME",TIMEZONE)
          CALL PRINT_REAL_TIME(get_file_time(ICING_FILE,ntimes),IPT,"LAST FILE TIME",TIMEZONE)

          CALL FATAL_ERROR&
               &("TO USE PERIODIC FORCING THE FILE TIME MUST COUNT FROM ZERO",&
               & "THE DIFFERENCE BETWEEN THE CURRENT MODEL TIME AND THE START TIME,",&
               & "MODULO THE FORCING PERIOD, DETERMINES THE CURRENT FORCING")
       END IF


       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! USING PERIODIC ICING FORCING:"
          CALL PRINT_TIME(ICING_PERIOD,IPT,"PERIOD")
       END IF

    CASE(VRBL)

       ICING_FILE => FIND_FILE(FILEHEAD,trim(ICING_FORCING_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE ICING BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE
       ATT => FIND_ATT(ICING_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(ICING_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")

       ICING_FORCING_COMMENTS ="VARIABLE ICING: "//TRIM(ATT%CHR(1))

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          ICING_FORCING_TYPE = ICING_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          ICING_FORCING_TYPE = ICING_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          ICING_FORCING_TYPE = ICING_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          ICING_FORCING_TYPE = ICING_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(ICING_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE ICING FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(ICING_FORCING_COMMENTS(4))
       ICING_FORCING_COMMENTS(1) = "FVCOM variable surface icing forcing:"
       ICING_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(ICING_FORCING_FILE)

       ICING_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(ICING_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          ICING_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          ICING_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF


       ! CHECK DIMENSIONS
       DIM => FIND_UNLIMITED(ICING_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            &"COULD NOT FIND UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       ! CHECK THE FILE TIME AND COMPARE WITH MODEL RUN TIME
       TIMETEST = get_file_time(ICING_FILE,1)
       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE FORCING TIME SERIES")

       TIMETEST = get_file_time(ICING_FILE,ntimes)
       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            &"THE MODEL RUN ENDS AFTER THE FORCING TIME SERIES")

    CASE DEFAULT
       CALL FATAL_ERROR("SURFACE_ICING: UNKNOWN ICING KIND?")
    END SELECT



    !==================================================================
    SELECT CASE(ICING_FORCING_TYPE)
       !==================================================================
    CASE(ICING_IS_WRFGRID)
       !==================================================================

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP ICING FORCING FROM A 'wrf grid' FILE"

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_DIM(ICING_FILE,'south_north',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            & "COULD NOT FIND DIMENSION 'south_north'")

       LATS = DIM%DIM

       DIM => FIND_DIM(ICING_FILE,'west_east',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            & "COULD NOT FIND DIMENSION 'west_east'")
       LONS = DIM%DIM


       CALL SET_FILE_INTERP_BILINEAR(ICING_FILE,ICING_INTP_N,ICING_INTP_C)

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       ! SEA LEVEL AIR TEMPERATURE
       IF(ASSOCIATED(ICE_FILE,ICING_FILE)) THEN
          ! USE THE ALREADY LOADED DATA FROM THE ICE MODEL
          ICING_SAT_N => ICE_SAT_N
          ICING_SAT_P => ICE_SAT_P

       ELSE

          VAR => FIND_VAR(ICING_FILE,"T2",FOUND)
          !       IF(.not. FOUND) VAR => FIND_VAR(ICING_FILE,"T2",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE HEATING BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
               & "COULD NOT FIND VARIABLE 'T2'")
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
          ICING_SAT_N => reference_var(var)
          CALL NC_CONNECT_PVAR(ICING_SAT_N,STORAGE_ARR)
          NULLIFY(STORAGE_ARR)
          
          ! MAKE SPACE FOR THE INTERPOLATED DATA
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
          CALL NC_CONNECT_PVAR(ICING_SAT_N,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
          
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
          ICING_SAT_P => reference_var(var)
          CALL NC_CONNECT_PVAR(ICING_SAT_P,STORAGE_ARR)
          NULLIFY(STORAGE_ARR)
          
          ! MAKE SPACE FOR THE INTERPOLATED DATA
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
          CALL NC_CONNECT_PVAR(ICING_SAT_P,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)

       END IF

       ! NET WIND SPEED X
       VAR => FIND_VAR(HEAT_FILE,"U10",FOUND)
!       IF(.not. FOUND) VAR => FIND_VAR(HEAT_FILE,"U10",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'U10'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       ICING_WSPX_N => reference_var(var)
       CALL NC_CONNECT_PVAR(ICING_WSPX_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       CALL NC_CONNECT_PVAR(ICING_WSPX_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       ICING_WSPX_P => reference_var(var)
       CALL NC_CONNECT_PVAR(ICING_WSPX_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       CALL NC_CONNECT_PVAR(ICING_WSPX_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! NET WIND SPEED Y
       VAR => FIND_VAR(HEAT_FILE,"V10",FOUND)
!       IF(.not. FOUND) VAR => FIND_VAR(HEAT_FILE,"U10",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'V10'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       ICING_WSPY_N => reference_var(var)
       CALL NC_CONNECT_PVAR(ICING_WSPY_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       CALL NC_CONNECT_PVAR(ICING_WSPY_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       ICING_WSPY_P => reference_var(var)
       CALL NC_CONNECT_PVAR(ICING_WSPY_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       CALL NC_CONNECT_PVAR(ICING_WSPY_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       !==================================================================
    CASE(ICING_IS_FVCOMGRID)
       !==================================================================

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP ICING FORCING FROM A 'fvcom grid' FILE"

       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_DIM(ICING_FILE,'node',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            & "COULD NOT FIND DIMENSION 'node'")

       if (mgl /= dim%dim) CALL FATAL_ERROR&
            &("Surface ICing: the number of nodes in the file does not match the fvcom grid?")


       DIM => FIND_DIM(HEAT_FILE,'nele',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            & "COULD NOT FIND DIMENSION 'nele'")

       if (ngl /= dim%dim) CALL FATAL_ERROR&
            &("Surface Icing: the number of elements in the file does not match the fvcom grid?")

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       ! Sea Surface Air Temperature
       IF(ASSOCIATED(ICE_FILE,ICING_FILE)) THEN
          ! USE THE SAME MEMORY USED FOR OCEAN MODEL HEAT FLUX
          
          ICING_SAT_N => ICE_SAT_N

          ICING_SAT_P => ICE_SAT_P

       ELSE       
          VAR => FIND_VAR(ICING_FILE,"T2",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
               & "COULD NOT FIND VARIABLE 'T2'")
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          ICING_SAT_N => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
          CALL NC_CONNECT_PVAR(ICING_SAT_N,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
          
          
          ! MAKE SPACE FOR THE DATA FROM THE FILE
          ICING_SAT_P => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
          CALL NC_CONNECT_PVAR(ICING_SAT_P,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)
          
       END IF
       
       ! Wind Speed X
       VAR => FIND_VAR(HEAT_FILE,"U10",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE ICING BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(ICING_FORCING_FILE),&
            & "COULD NOT FIND VARIABLE 'U10'")
       
       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ICING_WSPX_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       CALL NC_CONNECT_PVAR(ICING_WSPX_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ICING_WSPX_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN ICING_FORCING")
       CALL NC_CONNECT_PVAR(ICING_WSPX_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       !==================================================================
    CASE DEFAULT
       !==================================================================
       CALL FATAL_ERROR("CAN NOT RECOGNIZE ICING FILE TYPE!")
       !==================================================================
    END SELECT
    !==================================================================


    ! ---------- new: 2016 , april, Karsten Lettmann after Hint by Qi and ayumi.fujisaki@noaa.gov------
    ! Initialize some variables 
    ! afm 20150914
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3. 
    ICING_SAT_P%CURR_STKCNT = 0;  ICING_SAT_N%CURR_STKCNT = 0
    ICING_WSPX_P%CURR_STKCNT = 0; ICING_WSPX_N%CURR_STKCNT = 0
    ICING_WSPY_P%CURR_STKCNT = 0; ICING_WSPY_N%CURR_STKCNT = 0
    ! --------------------- end new ----------------------------------------------------



    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END ICING FORCING"
  END SUBROUTINE ICING_FORCING
!================================================================
!================================================================
  SUBROUTINE SURFACE_WINDSTRESS
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND
    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    TYPE(TIME) :: TIMETEST
    INTEGER :: LATS, LONS, I, Ntimes
    INTEGER :: STATUS
    CHARACTER(len=60) :: xstr, ystr

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START SURFACE_WINDSTRESS"

    NULLIFY(ATT,DIM,VAR,STORAGE_ARR,STORAGE_VEC)


    IF (.NOT. WIND_ON ) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! SURFACE WIND FORCING IS OFF!"
       ALLOCATE(WINDS_FORCING_COMMENTS(1))
       WINDS_FORCING_COMMENTS(1) = "SURFACE WIND FORCING IS OFF"
       RETURN
    END IF

    IF (WIND_TYPE /= SPEED .and.WIND_TYPE /= STRESS) CALL FATAL_ERROR&
         &("YOU MUST SELECT A WIND TYPE IN THE RUNFILE: '"&
         &//TRIM(SPEED)//", or '"//TRIM(STRESS)//"'")

!---> Siqi Li, 2021-01-27
    IF (WIND_ON) THEN

      SELECT CASE (TRIM(WIND_STRESS_METHOD))

        CASE ('LP1981')
          IF(DBG_SET(DBG_LOG)) write(IPT,*) "! WIND_STRESS_METHOD : LP1981!"

        CASE ('COARE')
          IF(DBG_SET(DBG_LOG)) write(IPT,*) "! WIND_STRESS_METHOD : COARE!"

        CASE ('TY2001')
          CALL FATAL_ERROR("TO USE TY2001 FOR WIND_STRESS_METHOD, &
                     &      RECOMPILE FVCOM WITH WAVE_ONLY or WAVE_CURRENT_INTERACTION")

        CASE ('OOST')
          CALL FATAL_ERROR("TO USE OOST FOR WIND_STRESS_METHOD, &
                     &      RECOMPILE FVCOM WITH 1 and &
                     &       WAVE_ONLY / WAVE_CURRENT_INTERACTION")

        CASE ('DGHQ')
          CALL FATAL_ERROR("TO USE DGHQ FOR WIND_STRESS_METHOD, &
                     &      RECOMPILE FVCOM WITH 1 and &
                     &       WAVE_ONLY / WAVE_CURRENT_INTERACTION")


        CASE DEFAULT
          CALL FATAL_ERROR               &
               &   ("WRONG WIND_STRESS_METHOD OPTIONS.")

      END SELECT


    END IF
!<--- Siqi Li, 2021-01-27

! DETERMINE HOW TO LOAD THE DATA
    SELECT CASE(WIND_KIND)
    CASE (CNSTNT)

       write(xstr,'(f8.4)') WIND_X
       write(ystr,'(f8.4)') WIND_Y

       IF (WIND_TYPE == SPEED)THEN
          
          IF(DBG_SET(DBG_LOG)) THEN
             WRITE(IPT,*)"! SETTING UP CONSTANT WIND SPEED FORCING: " 
             WRITE(IPT,*)"      Xspeed: "//trim(xstr)
             WRITE(IPT,*)"      Yspeed: "//trim(ystr)
          END IF
          
          ALLOCATE(WINDS_FORCING_COMMENTS(3))
          WINDS_FORCING_COMMENTS(1) = "Using constant wind speed from run file:"
          WINDS_FORCING_COMMENTS(2) = "Xspeed:"//trim(xstr)
          WINDS_FORCING_COMMENTS(3) = "Yspeed:"//trim(ystr)
          RETURN
       ELSEIF(WIND_TYPE == STRESS)THEN
          
          IF(DBG_SET(DBG_LOG)) THEN
             WRITE(IPT,*)"! SETTING UP CONSTANT WIND STRESS FORCING: " 
             WRITE(IPT,*)"      Xstress: "//trim(xstr)
             WRITE(IPT,*)"      Ystress: "//trim(ystr)
          END IF
          
          ALLOCATE(WINDS_FORCING_COMMENTS(3))
          WINDS_FORCING_COMMENTS(1) = "Using constant wind stress from run file:"
          WINDS_FORCING_COMMENTS(2) = "Xstress:"//trim(xstr)
          WINDS_FORCING_COMMENTS(3) = "Ystress:"//trim(ystr)
          RETURN

       END IF

   CASE(STTC)

       CALL FATAL_ERROR("STATIC WIND Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT WIND Not Set Up Yet")

    CASE(PRDC)
    
       WINDS_FILE => FIND_FILE(FILEHEAD,trim(WIND_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE WIND BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE   
       ATT => FIND_ATT(WINDS_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(WINDS_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
       

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_WRFGRID
       ELSE
          CALL PRINT_FILE(WINDS_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE WIND FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE
       
       ALLOCATE(WINDS_FORCING_COMMENTS(4))
       WINDS_FORCING_COMMENTS(1) = "FVCOM periodic surface Wind forcing:"
       WINDS_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(WIND_FILE)

       WINDS_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(WINDS_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          WINDS_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          WINDS_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF


       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING
       
       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(WINDS_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            &"COULD NOT FIND THE UNLMITIED DIMENSION")

       NTIMES = DIM%DIM

       WINDS_PERIOD = get_file_time(WINDS_FILE,ntimes)

       WINDS_PERIOD = WINDS_PERIOD - get_file_time(WINDS_FILE,1)

       IF (WINDS_PERIOD /= get_file_time(WINDS_FILE,ntimes)) THEN

          CALL PRINT_REAL_TIME(get_file_time(WINDS_FILE,1),IPT,"FIRST FILE TIME",TIMEZONE)
          CALL PRINT_REAL_TIME(get_file_time(WINDS_FILE,ntimes),IPT,"LAST FILE TIME",TIMEZONE)

          CALL FATAL_ERROR&
               &("TO USE PERIODIC FORCING THE FILE TIME MUST COUNT FROM ZERO",&
               & "THE DIFFERENCE BETWEEN THE CURRENT MODEL TIME AND THE START TIME,",&
               & "MODULO THE FORCING PERIOD, DETERMINES THE CURRENT FORCING")
       END IF


       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! USING PERIODIC WIND FORCING:"
          CALL PRINT_TIME(WINDS_PERIOD,IPT,"PERIOD")
       END IF


    CASE(VRBL)

       WINDS_FILE => FIND_FILE(FILEHEAD,trim(WIND_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE WIND BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE   
       ATT => FIND_ATT(WINDS_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(WINDS_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
       
       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_WRFGRID
       ELSE IF (ATT%CHR(1)(1:len_trim(surf_forcing_pt_source)) ==&
            & TRIM(surf_forcing_pt_source)) THEN
          WINDS_FORCING_TYPE = WINDS_ARE_PT_SOURCE
       ELSE
          CALL PRINT_FILE(WINDS_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE WIND FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(WINDS_FORCING_COMMENTS(4))
       WINDS_FORCING_COMMENTS(1) = "FVCOM variable surface Wind forcing:"
       WINDS_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(WIND_FILE)
       
       WINDS_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))
       
       ATT_DATE => FIND_ATT(WINDS_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          WINDS_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          WINDS_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF
              
       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(WINDS_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")
       
       NTIMES = DIM%DIM
       
       ! CHECK THE FILE TIME AND COMPARE WITH MODEL RUN TIME
       TIMETEST = get_file_time(WINDS_FILE,1)
       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE FORCING TIME SERIES")
       
       TIMETEST = get_file_time(WINDS_FILE,ntimes)
       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            &"THE MODEL RUN ENDS AFTER THE FORCING TIME SERIES")
       
    CASE DEFAULT
       CALL FATAL_ERROR("SURFACE_WINDSTRESS: UNKNOWN WIND KIND?")
    END SELECT

!==============================================================
    SELECT CASE(WINDS_FORCING_TYPE)
!==============================================================
    CASE(WINDS_ARE_WRFGRID)
!==============================================================


       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP WIND STRESS FORCING FROM A 'wrf grid' FILE"

       DIM => FIND_DIM(WINDS_FILE,'south_north',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            & "COULD NOT FIND DIMENSION 'south_north'")

       LATS = DIM%DIM

       DIM => FIND_DIM(WINDS_FILE,'west_east',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            & "COULD NOT FIND DIMENSION 'west_east'")
       LONS = DIM%DIM

       CALL SET_FILE_INTERP_bilinear(WINDS_FILE,WINDS_INTP_N,WINDS_INTP_C)

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       IF (WIND_TYPE == SPEED)THEN

          ! WIND SPEED IN THE X or EAST-WEST DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"uwind_speed",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(WINDS_FILE,"U10",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'uwind_speed' or 'U10'")

       ELSEIF(WIND_TYPE == STRESS)THEN
          ! WIND STRESS IN THE X or EAST-WEST DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"uwind_stress",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(WINDS_FILE,"Stress_U",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'uwind_stress' or 'Stress_U'")
       END IF

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       WINDS_STRX_N => reference_var(var)
       CALL NC_CONNECT_PVAR(WINDS_STRX_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:NT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRX_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       WINDS_STRX_P => reference_var(var)
       CALL NC_CONNECT_PVAR(WINDS_STRX_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:NT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRX_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
       
       IF (WIND_TYPE == SPEED)THEN

          ! WIND SPEED IN THE Y or NORTH SOUTH DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"vwind_speed",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(WINDS_FILE,"V10",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'vwind_speed' or 'V10'")

       ELSEIF(WIND_TYPE == STRESS)THEN
          ! WIND STRESS IN THE Y or NORTH SOUTH DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"vwind_stress",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(WINDS_FILE,"Stress_V",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'vwind_stress' or 'Stress_V'")
       END IF
          
       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       WINDS_STRY_N => reference_var(var)
       CALL NC_CONNECT_PVAR(WINDS_STRY_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:NT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRY_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       WINDS_STRY_P => reference_var(var)
       CALL NC_CONNECT_PVAR(WINDS_STRY_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:NT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRY_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


!==============================================================
    CASE(WINDS_ARE_FVCOMGRID)
!==============================================================
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP WIND STRESS FORCING FROM A 'FVCOM GRID' FILE"

       DIM => FIND_DIM(WINDS_FILE,'node',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            & "COULD NOT FIND DIMENSION 'node'")

       if (mgl /= dim%dim) CALL FATAL_ERROR&
            &("Surface Windstress: the number of nodes in the file does not match the fvcom grid?")

       DIM => FIND_DIM(WINDS_FILE,'nele',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            & "COULD NOT FIND DIMENSION 'nele'")

       if (ngl /= dim%dim) CALL FATAL_ERROR&
            &("Surface Windstress: the number of elements in the file does not match the fvcom grid?")


       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       IF (WIND_TYPE == SPEED)THEN
       
          ! WIND SPEED IN THE X or EAST-WEST DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"uwind_speed",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(WINDS_FILE,"U10",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'uwind_speed' or 'U10'")

       ELSEIF(WIND_TYPE == STRESS)THEN
          
          ! WIND STRESS IN THE X or EAST-WEST DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"uwind_stress",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'uwind_stress'")
       
       END IF

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WINDS_STRX_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:NT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRX_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WINDS_STRX_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:NT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRX_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       IF (WIND_TYPE == SPEED)THEN

          ! WIND SPEED IN THE Y or NORTH SOUTH DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"vwind_speed",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(WINDS_FILE,"V10",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'vwind_speed' or 'V10'")
       ELSEIF(WIND_TYPE == STRESS)THEN
          
          ! WIND STRESS IN THE Y or NORTH SOUTH DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"vwind_stress",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'vwind_stress'")
       END IF

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WINDS_STRY_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:NT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRY_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WINDS_STRY_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:NT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRY_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
       
!==============================================================
    CASE(WINDS_ARE_PT_SOURCE)
!==============================================================
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP WIND STRESS FORCING FROM A 'FVCOM GRID' FILE"


       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       IF (WIND_TYPE == SPEED)THEN
       
          ! WIND SPEED IN THE X or EAST-WEST DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"uwind_speed",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(WINDS_FILE,"U10",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'uwind_speed' or 'U10'")

       ELSEIF(WIND_TYPE == STRESS)THEN
          
          ! WIND STRESS IN THE X or EAST-WEST DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"uwind_stress",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'uwind_stress'")
       
       END IF

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WINDS_STRX_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(1), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRX_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WINDS_STRX_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(1), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRX_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       IF (WIND_TYPE == SPEED)THEN

          ! WIND SPEED IN THE Y or NORTH SOUTH DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"vwind_speed",FOUND)
          IF(.not. FOUND) VAR => FIND_VAR(WINDS_FILE,"V10",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'vwind_speed' or 'V10'")
       ELSEIF(WIND_TYPE == STRESS)THEN
          
          ! WIND STRESS IN THE Y or NORTH SOUTH DIRECTION
          VAR => FIND_VAR(WINDS_FILE,"vwind_stress",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(WIND_FILE),&
               & "COULD NOT FIND VARIABLE 'vwind_stress'")
       END IF

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WINDS_STRY_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(1), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRY_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WINDS_STRY_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(1), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE WINDSTRESS")
       CALL NC_CONNECT_PVAR(WINDS_STRY_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
       
!==============================================================
    CASE DEFAULT
!==============================================================
       CALL FATAL_ERROR("CAN NOT RECOGNIZE WIND FILE TYPE!")
!==============================================================
    END SELECT
!==============================================================

    ! afm 20151112 & EJA 20160921
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3. 
    winds_strx_n%curr_stkcnt = 0; winds_strx_p%curr_stkcnt = 0
    winds_stry_n%curr_stkcnt = 0; winds_stry_p%curr_stkcnt = 0

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END SURFACE_WINDSTRESS"
  END SUBROUTINE SURFACE_WINDSTRESS

!================================================================













!================================================================
  SUBROUTINE SURFACE_WAVE
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND
    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    TYPE(TIME) :: TIMETEST
    INTEGER :: LATS, LONS, I, Ntimes
    INTEGER :: STATUS
    CHARACTER(len=60) :: w_hs, w_len,w_dir,w_per,w_per_bot,w_ub_bot

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START SURFACE_WAVE"

    NULLIFY(ATT,DIM,VAR,STORAGE_ARR,STORAGE_VEC)


    IF (.NOT. WAVE_ON ) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! SURFACE WAVE FORCING IS OFF!"
       ALLOCATE(WAVES_FORCING_COMMENTS(1))
       WAVES_FORCING_COMMENTS(1) = "SURFACE WAVE FORCING IS OFF"
       RETURN
    END IF


! DETERMINE HOW TO LOAD THE DATA
    SELECT CASE(WAVE_KIND)
    CASE (CNSTNT)

       write(w_hs,     '(f8.4)') WAVE_HEIGHT
       write(w_len,    '(f8.4)') WAVE_LENGTH
       write(w_dir,    '(f8.4)') WAVE_DIRECTION
       write(w_per,    '(f8.4)') WAVE_PERIOD
       write(w_per_bot,'(f8.4)') WAVE_PER_BOT
       write(w_ub_bot, '(f8.4)') WAVE_UB_BOT

          
       IF(DBG_SET(DBG_LOG)) THEN
         WRITE(IPT,*)"! SETTING UP CONSTANT SURFACE WAVE FORCING: " 
         WRITE(IPT,*)"  wave height   : "//trim(w_hs)
         WRITE(IPT,*)"  wave length   : "//trim(w_len)
         WRITE(IPT,*)"  wave direction: "//trim(w_dir)
         WRITE(IPT,*)"  wave period   : "//trim(w_per)
         WRITE(IPT,*)"  wave per_bot  : "//trim(w_per_bot)
         WRITE(IPT,*)"  wave ub_bot   : "//trim(w_ub_bot)
       END IF
          
       ALLOCATE(WAVES_FORCING_COMMENTS(7))
       WAVES_FORCING_COMMENTS(1) = "Using constant surface wave from run file:"
       WAVES_FORCING_COMMENTS(2) = "  wave height   : "//trim(w_hs)
       WAVES_FORCING_COMMENTS(3) = "  wave length   : "//trim(w_len)
       WAVES_FORCING_COMMENTS(4) = "  wave direction: "//trim(w_dir)
       WAVES_FORCING_COMMENTS(5) = "  wave period   : "//trim(w_per)
       WAVES_FORCING_COMMENTS(6) = "  wave per_bot  : "//trim(w_per_bot)
       WAVES_FORCING_COMMENTS(7) = "  wave ub_bot   : "//trim(w_ub_bot)
       RETURN


   CASE(STTC)

       CALL FATAL_ERROR("STATIC WAVE Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT WAVE Not Set Up Yet")

    CASE(PRDC)
    
       WAVES_FILE => FIND_FILE(FILEHEAD,trim(WAVE_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE WAVE BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE   
       ATT => FIND_ATT(WAVES_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(WAVES_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
       

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          WAVES_FORCING_TYPE = WAVES_ARE_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          WAVES_FORCING_TYPE = WAVES_ARE_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          WAVES_FORCING_TYPE = WAVES_ARE_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          WAVES_FORCING_TYPE = WAVES_ARE_WRFGRID

       ELSE
          CALL PRINT_FILE(WAVES_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE WAVE FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE
       
       ALLOCATE(WAVES_FORCING_COMMENTS(4))
       WAVES_FORCING_COMMENTS(1) = "FVCOM periodic surface wave forcing:"
       WAVES_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(WAVE_FILE)

       WAVES_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(WAVES_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          WAVES_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          WAVES_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF


       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING
       
       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(WAVES_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            &"COULD NOT FIND THE UNLMITIED DIMENSION")

       NTIMES = DIM%DIM

       WAVES_PERIOD = get_file_time(WAVES_FILE,ntimes)

       WAVES_PERIOD = WAVES_PERIOD - get_file_time(WAVES_FILE,1)

       IF (WAVES_PERIOD /= get_file_time(WAVES_FILE,ntimes)) THEN

          CALL PRINT_REAL_TIME(get_file_time(WAVES_FILE,1),IPT,"FIRST FILE TIME",TIMEZONE)
          CALL PRINT_REAL_TIME(get_file_time(WAVES_FILE,ntimes),IPT,"LAST FILE TIME",TIMEZONE)

          CALL FATAL_ERROR&
               &("TO USE PERIODIC FORCING THE FILE TIME MUST COUNT FROM ZERO",&
               & "THE DIFFERENCE BETWEEN THE CURRENT MODEL TIME AND THE START TIME,",&
               & "MODULO THE FORCING PERIOD, DETERMINES THE CURRENT FORCING")
       END IF


       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! USING PERIODIC WAVE FORCING:"
          CALL PRINT_TIME(WAVES_PERIOD,IPT,"PERIOD")
       END IF


    CASE(VRBL)

       WAVES_FILE => FIND_FILE(FILEHEAD,trim(WAVE_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE WAVE BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE   
       ATT => FIND_ATT(WAVES_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(WAVES_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
       
       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          WAVES_FORCING_TYPE = WAVES_ARE_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source).or.(ATT%CHR(1)(1:5)=='fvcom')) THEN
          WAVES_FORCING_TYPE = WAVES_ARE_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source).or.(ATT%CHR(1)(1:5)=='FVCOM')) THEN
          WAVES_FORCING_TYPE = WAVES_ARE_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          WAVES_FORCING_TYPE = WAVES_ARE_WRFGRID

       ELSE
          CALL PRINT_FILE(WAVES_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE WAVE FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(WAVES_FORCING_COMMENTS(4))
       WAVES_FORCING_COMMENTS(1) = "FVCOM variable surface wave forcing:"
       WAVES_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(WAVE_FILE)
       
       WAVES_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))
       
       ATT_DATE => FIND_ATT(WAVES_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          WAVES_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          WAVES_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF
              
       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(WAVES_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")
       
       NTIMES = DIM%DIM
       
       ! CHECK THE FILE TIME AND COMPARE WITH MODEL RUN TIME
       TIMETEST = get_file_time(WAVES_FILE,1)
       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE FORCING TIME SERIES")
       
       TIMETEST = get_file_time(WAVES_FILE,ntimes)
       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            &"THE MODEL RUN ENDS AFTER THE FORCING TIME SERIES")
       
    CASE DEFAULT
       CALL FATAL_ERROR("SURFACE_WAVE: UNKNOWN WAVE KIND?")
    END SELECT

!==============================================================
    SELECT CASE(WAVES_FORCING_TYPE)
!==============================================================
    CASE(WAVES_ARE_WRFGRID)
!==============================================================

       CALL FATAL_ERROR("WAVE based on WRF grid Not Set Up Yet")

!==============================================================
    CASE(WAVES_ARE_FVCOMGRID)
!==============================================================
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP WAVE FORCING FROM A 'FVCOM GRID' FILE"

       DIM => FIND_DIM(WAVES_FILE,'node',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            & "COULD NOT FIND DIMENSION 'node'")

       if (mgl /= dim%dim) CALL FATAL_ERROR&
            &("Surface Wave: the number of nodes in the file does not match the fvcom grid?")


       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!


       
       ! WAVE HEIGHT
       VAR => FIND_VAR(WAVES_FILE,"hs",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            & "COULD NOT FIND VARIABLE 'hs' ")


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_HEIGHT_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN WAVE HEIGHT")
       CALL NC_CONNECT_PVAR(WAVES_HEIGHT_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_HEIGHT_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN WAVE HEIGHT")
       CALL NC_CONNECT_PVAR(WAVES_HEIGHT_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! WAVE LENGTH
       VAR => FIND_VAR(WAVES_FILE,"wlen",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            & "COULD NOT FIND VARIABLE 'wlen' ")


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_LENGTH_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN WAVE LENGTH")
       CALL NC_CONNECT_PVAR(WAVES_LENGTH_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_LENGTH_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN WAVE LENGTH")
       CALL NC_CONNECT_PVAR(WAVES_LENGTH_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

      ! WAVE DIRECTION
       VAR => FIND_VAR(WAVES_FILE,"dirm",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(WAVES_FILE,"wdir",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            & "COULD NOT FIND VARIABLE 'dirm' or 'wdir' ")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_DIRECTION_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN WAVE DIRECTION")
       CALL NC_CONNECT_PVAR(WAVES_DIRECTION_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_DIRECTION_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN WAVE DIRECTION")
       CALL NC_CONNECT_PVAR(WAVES_DIRECTION_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
 
      ! WAVE PERIOD
       VAR => FIND_VAR(WAVES_FILE,"tpeak",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            & "COULD NOT FIND VARIABLE 'tpeak' ")


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_PERIOD_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN WAVE PERIOD")
       CALL NC_CONNECT_PVAR(WAVES_PERIOD_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_PERIOD_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN WAVE PERIOD")
       CALL NC_CONNECT_PVAR(WAVES_PERIOD_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


      ! BOTTOM WAVE PERIOD
       VAR => FIND_VAR(WAVES_FILE,"pwave_bot",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(WAVES_FILE,"tmbot",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            & "COULD NOT FIND VARIABLE 'pwave_bot' or 'tmbot' ")


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_PER_BOT_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN BOTTOM WAVE PERIOD")
       CALL NC_CONNECT_PVAR(WAVES_PER_BOT_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_PER_BOT_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN BOTTOM WAVE PERIOD")
       CALL NC_CONNECT_PVAR(WAVES_PER_BOT_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


      ! BOTTOM ORBITAL VELOCITY
       VAR => FIND_VAR(WAVES_FILE,"ub_bot",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(WAVES_FILE,"ubot",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WAVE BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WAVE_FILE),&
            & "COULD NOT FIND VARIABLE 'ub_bot' or 'ubot' ")


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_UB_BOT_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN BOTTOM WAVE ORBITAL VELOCITY")
       CALL NC_CONNECT_PVAR(WAVES_UB_BOT_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       WAVES_UB_BOT_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN BOTTOM WAVE ORBITAL VELOCITY")
       CALL NC_CONNECT_PVAR(WAVES_UB_BOT_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       

!==============================================================
    CASE DEFAULT
!==============================================================
       CALL FATAL_ERROR("CAN NOT RECOGNIZE WAVE FILE TYPE!")
!==============================================================
    END SELECT
!==============================================================

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END SURFACE_WAVE"
  END SUBROUTINE SURFACE_WAVE

!================================================================



!================================================================
  SUBROUTINE SURFACE_AIRPRESSURE
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND
    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    TYPE(TIME) :: TIMETEST
    INTEGER :: LATS, LONS, I, Ntimes
    INTEGER :: STATUS
    CHARACTER(len=60) :: airpressurestr

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START SURFACE_AIRPRESSURE"

    NULLIFY(ATT,DIM,VAR,STORAGE_ARR,STORAGE_VEC)


    IF (.NOT. AIRPRESSURE_ON ) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! SURFACE AIR PRESSURE FORCING IS OFF!"
       ALLOCATE(AIRPRESSURE_FORCING_COMMENTS(1))
       AIRPRESSURE_FORCING_COMMENTS(1) = "SURFACE AIR PRESSURE FORCING IS OFF"
       RETURN
    END IF


! DETERMINE HOW TO LOAD THE DATA
    SELECT CASE(AIRPRESSURE_KIND)
    CASE (CNSTNT)

       write(airpressurestr,'(f8.4)') AIRPRESSURE_VALUE


       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*)"! SETTING UP CONSTANT AIR PRESSURE FORCING: " 
          WRITE(IPT,*)"      Air pressure: "//trim(airpressurestr)
       END IF
       
       ALLOCATE(AIRPRESSURE_FORCING_COMMENTS(3))
       AIRPRESSURE_FORCING_COMMENTS(1) = "Using constant air pressure from run file:"
       AIRPRESSURE_FORCING_COMMENTS(2) = "Air pressure:"//trim(airpressurestr)
       RETURN

   CASE(STTC)

       CALL FATAL_ERROR("STATIC AIR PRESSURE Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT AIR PRESSURE Not Set Up Yet")

    CASE(PRDC)
    
       AIRPRESSURE_P_FILE => FIND_FILE(FILEHEAD,trim(AIRPRESSURE_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE   
       ATT => FIND_ATT(AIRPRESSURE_P_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(AIRPRESSURE_P_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
       

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          AIRPRESSURE_FORCING_TYPE = AIRPRESSURE_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          AIRPRESSURE_FORCING_TYPE = AIRPRESSURE_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          AIRPRESSURE_FORCING_TYPE = AIRPRESSURE_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          AIRPRESSURE_FORCING_TYPE = AIRPRESSURE_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(AIRPRESSURE_P_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE AIR PRESSURE FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE
       
       ALLOCATE(AIRPRESSURE_FORCING_COMMENTS(4))
       AIRPRESSURE_FORCING_COMMENTS(1) = "FVCOM periodic surface Air Pressure forcing:"
       AIRPRESSURE_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(AIRPRESSURE_FILE)

       AIRPRESSURE_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(AIRPRESSURE_P_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          AIRPRESSURE_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          AIRPRESSURE_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF


       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING
       
       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(AIRPRESSURE_P_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            &"COULD NOT FIND THE UNLMITIED DIMENSION")

       NTIMES = DIM%DIM

       AIRPRESSURE_PERIOD = get_file_time(AIRPRESSURE_P_FILE,ntimes)

       AIRPRESSURE_PERIOD = AIRPRESSURE_PERIOD - get_file_time(AIRPRESSURE_P_FILE,1)

       IF (AIRPRESSURE_PERIOD /= get_file_time(AIRPRESSURE_P_FILE,ntimes)) THEN

          CALL PRINT_REAL_TIME(get_file_time(AIRPRESSURE_P_FILE,1),IPT,"FIRST FILE TIME",TIMEZONE)
          CALL PRINT_REAL_TIME(get_file_time(AIRPRESSURE_P_FILE,ntimes),IPT,"LAST FILE TIME",TIMEZONE)

          CALL FATAL_ERROR&
               &("TO USE PERIODIC FORCING THE FILE TIME MUST COUNT FROM ZERO",&
               & "THE DIFFERENCE BETWEEN THE CURRENT MODEL TIME AND THE START TIME,",&
               & "MODULO THE FORCING PERIOD, DETERMINES THE CURRENT FORCING")
       END IF


       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! USING PERIODIC AIR PRESSURE FORCING:"
          CALL PRINT_TIME(AIRPRESSURE_PERIOD,IPT,"PERIOD")
       END IF


    CASE(VRBL)

       AIRPRESSURE_P_FILE => FIND_FILE(FILEHEAD,trim(AIRPRESSURE_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE   
       ATT => FIND_ATT(AIRPRESSURE_P_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(AIRPRESSURE_P_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
       
       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          AIRPRESSURE_FORCING_TYPE = AIRPRESSURE_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          AIRPRESSURE_FORCING_TYPE = AIRPRESSURE_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          AIRPRESSURE_FORCING_TYPE = AIRPRESSURE_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          AIRPRESSURE_FORCING_TYPE = AIRPRESSURE_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(AIRPRESSURE_P_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE AIR PRESSURE FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(AIRPRESSURE_FORCING_COMMENTS(4))
       AIRPRESSURE_FORCING_COMMENTS(1) = "FVCOM variable surface Air Pressure:"
       AIRPRESSURE_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(AIRPRESSURE_FILE)
       
       AIRPRESSURE_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))
       
       ATT_DATE => FIND_ATT(AIRPRESSURE_P_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          AIRPRESSURE_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          AIRPRESSURE_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF
              
       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(AIRPRESSURE_P_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")
       
       NTIMES = DIM%DIM
       
       ! CHECK THE FILE TIME AND COMPARE WITH MODEL RUN TIME
       TIMETEST = get_file_time(AIRPRESSURE_P_FILE,1)
       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE FORCING TIME SERIES")
       
       TIMETEST = get_file_time(AIRPRESSURE_P_FILE,ntimes)
       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            &"THE MODEL RUN ENDS AFTER THE FORCING TIME SERIES")
       
    CASE DEFAULT
       CALL FATAL_ERROR("SURFACE_AIRPRESSURE: UNKNOWN ARE PRESSURE KIND?")
    END SELECT

!==============================================================
    SELECT CASE(AIRPRESSURE_FORCING_TYPE)
!==============================================================
    CASE(AIRPRESSURE_IS_WRFGRID)
!==============================================================


       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP AIR PRESSURE FORCING FROM A 'wrf grid' FILE"

       DIM => FIND_DIM(AIRPRESSURE_P_FILE,'south_north',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            & "COULD NOT FIND DIMENSION 'south_north'")

       LATS = DIM%DIM

       DIM => FIND_DIM(AIRPRESSURE_P_FILE,'west_east',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            & "COULD NOT FIND DIMENSION 'west_east'")
       LONS = DIM%DIM

       CALL SET_FILE_INTERP_bilinear(AIRPRESSURE_P_FILE,AIRPRESSURE_INTP_N,AIRPRESSURE_INTP_C)

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       ! AIR PRESSURE
       VAR => FIND_VAR(AIRPRESSURE_P_FILE,"air_pressure",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(AIRPRESSURE_P_FILE,"pressure_air",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(AIRPRESSURE_P_FILE,"SLP",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            & "COULD NOT FIND VARIABLE 'air_pressure' OR 'SLP'")

       !---> Siqi Li, 2021-01-27
       ! The unit of air pressure is really important. It can influence both
       ! water elevation and wind stress. In calculation, the unit of air 
       ! pressure is 'Pa'.
       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(VAR,'unit',FOUND)
       IF (FOUND) THEN
       call UPCASE(ATT%CHR(1))
       SELECT CASE (TRIM(ATT%CHR(1)))
       CASE ('MB', 'HPA')
         Pair_unit_factor = 100.0_SP
       CASE ('PA','PASCAL')
         Pair_unit_factor = 1.0_SP
       CASE DEFAULT
         CALL FATAL_ERROR("UNKNOWN UNIT of AIR PRESSURE: " &
                 & //TRIM(ATT%CHR(1)))
       END SELECT
       ELSE
         ! There is no unit attribute for the air_pressure variable
         ! We assume the unit is 'Pa'.
         Pair_unit_factor = 1.0_SP
       END IF
       !<--- Siqi Li, 2021-01-27

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE AIR PRESSURE")
       AIR_PRESSURE_N => reference_var(var)
       CALL NC_CONNECT_PVAR(AIR_PRESSURE_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE AIR PRESSURE")
       CALL NC_CONNECT_PVAR(AIR_PRESSURE_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE AIR PRESSURE")
       AIR_PRESSURE_P => reference_var(var)
       CALL NC_CONNECT_PVAR(AIR_PRESSURE_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE AIR PRESSURE")
       CALL NC_CONNECT_PVAR(AIR_PRESSURE_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

!==============================================================
    CASE(AIRPRESSURE_IS_FVCOMGRID)
!==============================================================
       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP AIR PRESSURE FORCING FROM A 'FVCOM GRID' FILE"

       DIM => FIND_DIM(AIRPRESSURE_P_FILE,'node',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            & "COULD NOT FIND DIMENSION 'node'")

       if (mgl /= dim%dim) CALL FATAL_ERROR&
            &("Surface Air Pressure: the number of nodes in the file does not match the fvcom grid?")

       DIM => FIND_DIM(AIRPRESSURE_P_FILE,'nele',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            & "COULD NOT FIND DIMENSION 'nele'")

       if (ngl /= dim%dim) CALL FATAL_ERROR&
            &("Surface Air Pressure: the number of elements in the file does not match the fvcom grid?")


       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       ! AIR PRESSURE
       VAR => FIND_VAR(AIRPRESSURE_P_FILE,"air_pressure",FOUND)
       IF(.not. FOUND) VAR => FIND_VAR(AIRPRESSURE_P_FILE,"SLP",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE AIR PRESSURE FILE OBJECT",&
            & "FILE NAME: "//TRIM(AIRPRESSURE_FILE),&
            & "COULD NOT FIND VARIABLE 'air_pressure' or 'SLP'")

       !---> Siqi Li, 2021-01-27
       ! The unit of air pressure is really important. It can influence both
       ! water elevation and wind stress. In calculation, the unit of air 
       ! pressure is 'Pa'.
       ATT => FIND_ATT(VAR,'units',FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(VAR,'unit',FOUND)
       IF (FOUND) THEN
       call UPCASE(ATT%CHR(1))
       SELECT CASE (TRIM(ATT%CHR(1)))
       CASE ('MB', 'HPA')
         Pair_unit_factor = 100.0_SP
       CASE ('PA','PASCAL')
         Pair_unit_factor = 1.0_SP
       CASE DEFAULT
         CALL FATAL_ERROR("UNKNOWN UNIT of AIR PRESSURE: " &
                 & //TRIM(ATT%CHR(1)))
       END SELECT
       ELSE
         ! There is no unit attribute for the air_pressure variable
         ! We assume the unit is 'Pa'.
         Pair_unit_factor = 1.0_SP
       END IF
       !<--- Siqi Li, 2021-01-27

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       AIR_PRESSURE_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN AIR PRESSURE")
       CALL NC_CONNECT_PVAR(AIR_PRESSURE_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       AIR_PRESSURE_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN AIR PRESSURE")
       CALL NC_CONNECT_PVAR(AIR_PRESSURE_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

!==============================================================
    CASE DEFAULT
!==============================================================
       CALL FATAL_ERROR("CAN NOT RECOGNIZE AIR PRESSURE FILE TYPE!")
!==============================================================
    END SELECT
!==============================================================

    ! ---------- new: 2016 , april, Karsten Lettmann after Hint by Qi and ayumi.fujisaki@noaa.gov------
    ! Initialize some variables 
    ! afm 20150914
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3. 
    ! PRECIPITATION
    AIR_PRESSURE_P%curr_stkcnt = 0 ; AIR_PRESSURE_N%curr_stkcnt = 0
    ! --------- end new ----------------------------------------------------------------



    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END SURFACE_AIRPRESSURE"
  END SUBROUTINE SURFACE_AIRPRESSURE

!================================================================
!================================================================

  SUBROUTINE SURFACE_PRECIPITATION
    IMPLICIT NONE
    ! SOME NC POINTERS
    TYPE(NCATT), POINTER :: ATT, ATT_DATE
    TYPE(NCDIM), POINTER :: DIM
    TYPE(NCVAR), POINTER :: VAR
    LOGICAL :: FOUND
    REAL(SP), POINTER :: STORAGE_ARR(:,:), storage_vec(:)
    TYPE(TIME) :: TIMETEST
    INTEGER :: LATS, LONS, I, Ntimes
    INTEGER :: STATUS
    CHARACTER(len=60) :: evpstr, prcstr

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "START SURFACE_PRECIPITATION"

    IF (.NOT. PRECIPITATION_ON ) THEN
       IF(DBG_SET(DBG_LOG)) write(IPT,*) "! SURFACE PRECIPITATION FORCING IS OFF!"
       ALLOCATE(PRECIP_FORCING_COMMENTS(1))
       PRECIP_FORCING_COMMENTS(1) = "SURFACE PRECIPITATION FORCING IS OFF" 
       RETURN
    END IF

    NULLIFY(ATT,DIM,VAR,STORAGE_ARR,STORAGE_VEC)

! DETERMINE HOW TO LOAD THE DATA
    SELECT CASE(PRECIPITATION_KIND)
    CASE (CNSTNT)

       write(evpstr,'(f8.4)') PRECIPITATION_EVP
       write(prcstr,'(f8.4)') PRECIPITATION_PRC

       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*)"! SETTING UP CONSTANT PRECIPITATION FORCING: "
          WRITE(IPT,*)"    EVAPORATION: "//trim(evpstr)
          WRITE(IPT,*)"  PRECIPITATION: "//trim(prcstr)
       END IF
       
       ALLOCATE(PRECIP_FORCING_COMMENTS(3))
       PRECIP_FORCING_COMMENTS(1) = "Using constant precipitation from run file:"
       PRECIP_FORCING_COMMENTS(2) = "Precipitation:"//trim(prcstr)
       PRECIP_FORCING_COMMENTS(3) = "Evaporation:"//trim(evpstr)
       RETURN

   CASE(STTC)

       CALL FATAL_ERROR("STATIC PRECIP Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT PRECIP Not Set Up Yet")

    CASE(PRDC)
    

       PRECIP_FILE => FIND_FILE(FILEHEAD,trim(PRECIPITATION_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE WIND BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE))

       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE   
       ATT => FIND_ATT(PRECIP_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(PRECIP_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
       

       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          PRECIP_FORCING_TYPE = PRECIP_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          PRECIP_FORCING_TYPE = PRECIP_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          PRECIP_FORCING_TYPE = PRECIP_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          PRECIP_FORCING_TYPE = PRECIP_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(PRECIP_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE PRECIP FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

        ALLOCATE(PRECIP_FORCING_COMMENTS(4))
       PRECIP_FORCING_COMMENTS(1) = "FVCOM periodic surface precip forcing:"
       PRECIP_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(PRECIPITATION_FILE)

       PRECIP_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(PRECIP_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          PRECIP_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          PRECIP_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF

       ! GET THE FILES LENGTH OF TIME AND SAVE FOR PERIODIC FORCING
       
       ! LOOK FOR THE DIMENSIONS
       DIM => FIND_UNLIMITED(PRECIP_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(WIND_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       PRECIP_PERIOD = get_file_time(PRECIP_FILE,ntimes)

       PRECIP_PERIOD = PRECIP_PERIOD - get_file_time(PRECIP_FILE,1)

       IF (PRECIP_PERIOD /= get_file_time(PRECIP_FILE,ntimes)) THEN

          CALL PRINT_REAL_TIME(get_file_time(PRECIP_FILE,1),IPT,"FIRST FILE TIME",TIMEZONE)
          CALL PRINT_REAL_TIME(get_file_time(PRECIP_FILE,ntimes),IPT,"LAST FILE TIME",TIMEZONE)

          CALL FATAL_ERROR&
               &("TO USE PERIODIC FORCING THE FILE TIME MUST COUNT FROM ZERO",&
               & "THE DIFFERENCE BETWEEN THE CURRENT MODEL TIME AND THE START TIME,",&
               & "MODULO THE FORCING PERIOD, DETERMINES THE CURRENT FORCING")
       END IF


       IF(DBG_SET(DBG_LOG)) THEN
          WRITE(IPT,*) "! USING PERIODIC PRECIP FORCING:"
          CALL PRINT_TIME(PRECIP_PERIOD,IPT,"PERIOD")
       END IF


    CASE(VRBL)

       PRECIP_FILE => FIND_FILE(FILEHEAD,trim(PRECIPITATION_FILE),FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("COULD NOT FIND SURFACE WIND BOUNDARY CONDINTION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE))
       
       ! DETERMINE GRID TYPE BASED ON SOURCE ATTRIBUTE   
       ATT => FIND_ATT(PRECIP_FILE,"source",FOUND)
       IF(.not. FOUND) ATT => FIND_ATT(PRECIP_FILE,"Source",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            &"COULD NOT FIND GLOBAL ATTRIBURE: 'source'")
       
       IF (ATT%CHR(1)(1:len_trim(wrf2fvcom_source)) ==&
            & TRIM(wrf2fvcom_source)) THEN
          PRECIP_FORCING_TYPE = PRECIP_IS_WRFGRID 

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_grid_source)) ==&
            & TRIM(fvcom_grid_source)) THEN
          PRECIP_FORCING_TYPE = PRECIP_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(fvcom_cap_grid_source)) ==&
            & TRIM(fvcom_cap_grid_source)) THEN
          PRECIP_FORCING_TYPE = PRECIP_IS_FVCOMGRID

       ELSE IF (ATT%CHR(1)(1:len_trim(wrf_grid_source)) ==&
            & TRIM(wrf_grid_source)) THEN
          PRECIP_FORCING_TYPE = PRECIP_IS_WRFGRID

       ELSE
          CALL PRINT_FILE(PRECIP_FILE)
          CALL FATAL_ERROR("CAN NOT RECOGNIZE PRECIP FILE!",&
               & "UNKNOWN SOURCE STRING:",TRIM(ATT%CHR(1)))
       END IF
       ! GOT GRID TYPE

       ALLOCATE(PRECIP_FORCING_COMMENTS(4))
       PRECIP_FORCING_COMMENTS(1) = "FVCOM periodic surface precip forcing:"
       PRECIP_FORCING_COMMENTS(2) = "FILE NAME:"//TRIM(PRECIPITATION_FILE)

       PRECIP_FORCING_COMMENTS(3) = "SOURCE:"//TRIM(ATT%CHR(1))

       ATT_DATE => FIND_ATT(PRECIP_FILE,"START_DATE",FOUND)
       IF (FOUND) THEN
          PRECIP_FORCING_COMMENTS(4) ="MET DATA START DATE:"//TRIM(ATT_DATE%CHR(1))
       ELSE
          PRECIP_FORCING_COMMENTS(4) = "Unknown start date meta data format"
       END IF

       DIM => FIND_UNLIMITED(PRECIP_FILE,FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            &"COULD NOT FIND THE UNLIMITED DIMENSION")

       NTIMES = DIM%DIM

       ! CHECK THE FILE TIME AND COMPARE WITH MODEL RUN TIME
       TIMETEST = get_file_time(PRECIP_FILE,1)
       IF(TIMETEST > STARTTIME) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            &"THE MODEL RUN STARTS BEFORE THE FORCING TIME SERIES")

       TIMETEST = get_file_time(PRECIP_FILE,ntimes)
       IF(TIMETEST < ENDTIME) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            &"THE MODEL RUN ENDS AFTER THE FORCING TIME SERIES")


    CASE DEFAULT
       CALL FATAL_ERROR("SURFACE_PRECIP: UNKNOWN WIND KIND?")
    END SELECT


! DEAL WITH DATA SET UP
!=====================================================================
    SELECT CASE(PRECIP_FORCING_TYPE)
!=====================================================================
    CASE(PRECIP_IS_WRFGRID)
!=====================================================================

       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP WIND STRESS FORCING FROM A 'wrf grid' FILE"

       ! LOOK FOR THE DIMENSIONS


       DIM => FIND_DIM(PRECIP_FILE,'south_north',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            & "COULD NOT FIND DIMENSION 'south_north'")

       LATS = DIM%DIM


       DIM => FIND_DIM(PRECIP_FILE,'west_east',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            & "COULD NOT FIND DIMENSION 'west_east'")
       LONS = DIM%DIM



       CALL SET_FILE_INTERP_bilinear(PRECIP_FILE,PRECIP_INTP_N,PRECIP_INTP_C)

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       ! WIND STRESS IN THE X or EAST-WEST DIRECTION
       VAR => FIND_VAR(PRECIP_FILE,"Precipitation",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            & "COULD NOT FIND VARIABLE 'PRECIPITATION'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       PRECIP_PRE_N => reference_var(var)
       CALL NC_CONNECT_PVAR(PRECIP_PRE_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       CALL NC_CONNECT_PVAR(PRECIP_PRE_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       PRECIP_PRE_P => reference_var(var)
       CALL NC_CONNECT_PVAR(PRECIP_PRE_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       CALL NC_CONNECT_PVAR(PRECIP_PRE_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! WIND STRESS IN THE X or EAST-WEST DIRECTION
       VAR => FIND_VAR(PRECIP_FILE,"Evaporation",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN SURFACE WIND BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            & "COULD NOT FIND VARIABLE 'Evaporation'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       PRECIP_EVP_N => reference_var(var)
       CALL NC_CONNECT_PVAR(PRECIP_EVP_N,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       CALL NC_CONNECT_PVAR(PRECIP_EVP_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       ALLOCATE(STORAGE_ARR(lons,lats), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       PRECIP_EVP_P => reference_var(var)
       CALL NC_CONNECT_PVAR(PRECIP_EVP_P,STORAGE_ARR)
       NULLIFY(STORAGE_ARR)

       ! MAKE SPACE FOR THE INTERPOLATED DATA
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       CALL NC_CONNECT_PVAR(PRECIP_EVP_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

!=====================================================================
    CASE(PRECIP_IS_FVCOMGRID)
!=====================================================================


       IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) &
            & "! SETTING UP PRECIPITATION FORCING FROM A 'FVCOM grid' FILE"

       ! LOOK FOR THE DIMENSIONS


       DIM => FIND_DIM(PRECIP_FILE,'node',FOUND)  
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            & "COULD NOT FIND DIMENSION 'node'")

      if (mgl /= dim%dim) CALL FATAL_ERROR&
            &("Surface PRECIP: the number of nodes in the file does not match the fvcom grid?")


       DIM => FIND_DIM(PRECIP_FILE,'nele',FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            & "COULD NOT FIND DIMENSION 'nele'")

       if (ngl /= dim%dim) CALL FATAL_ERROR&
            &("Surface PRECIP: the number of elements in the file does not match the fvcom grid?")

     !<----------ishid debug 20230104
     IF (PRECIPITATION_BULK_ON) THEN
          !precipitation
          VAR => FIND_VAR(PRECIP_FILE,"Precipitation",FOUND)
          IF(.not. FOUND) CALL FATAL_ERROR &
               & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
               & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
               & "COULD NOT FIND VARIABLE 'precip'")

          ! MAKE SPACE FOR THE DATA FROM THE FILE
          PRECIP_PRE_N => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
          CALL NC_CONNECT_PVAR(PRECIP_PRE_N,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)

                 ! MAKE SPACE FOR THE DATA FROM THE FILE
          PRECIP_PRE_P => reference_var(var)
          ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
          CALL NC_CONNECT_PVAR(PRECIP_PRE_P,STORAGE_VEC)
          NULLIFY(STORAGE_VEC)

         ! !wind speed
         ! !uwind
         ! VAR => FIND_VAR(PRECIP_FILE,"uwind_speed",FOUND)
         ! IF(.not. FOUND) CALL FATAL_ERROR &
         !      & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
         !      & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
         !      & "COULD NOT FIND VARIABLE 'uwind_speed'")
         ! 

         ! ! MAKE SPACE FOR THE DATA FROM THE FILE
         ! PRECIP_uwnd_N => reference_var(var)
         ! ALLOCATE(STORAGE_VEC(0:MT), stat = status)
         ! IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
         ! CALL NC_CONNECT_PVAR(PRECIP_PRE_N,STORAGE_VEC)
         ! NULLIFY(STORAGE_VEC)

                 ! MAKE SPACE FOR THE DATA FROM THE FILE
         ! PRECIP_uwnd_P => reference_var(var)
         ! ALLOCATE(STORAGE_VEC(0:MT), stat = status)
         ! IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
         ! CALL NC_CONNECT_PVAR(PRECIP_PRE_P,STORAGE_VEC)
         ! NULLIFY(STORAGE_VEC)

          !vwind
         ! VAR => FIND_VAR(PRECIP_FILE,"vwind_speed",FOUND)
         ! IF(.not. FOUND) CALL FATAL_ERROR &
         !      & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
         !      & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
         !      & "COULD NOT FIND VARIABLE 'vwind_speed'")

          ! MAKE SPACE FOR THE DATA FROM THE FILE
         ! PRECIP_vwnd_N => reference_var(var)
         ! ALLOCATE(STORAGE_VEC(0:MT), stat = status)
         ! IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
         ! CALL NC_CONNECT_PVAR(PRECIP_PRE_N,STORAGE_VEC)
         ! NULLIFY(STORAGE_VEC)

                 ! MAKE SPACE FOR THE DATA FROM THE FILE
         ! PRECIP_vwnd_P => reference_var(var)
         ! ALLOCATE(STORAGE_VEC(0:MT), stat = status)
         ! IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
         ! CALL NC_CONNECT_PVAR(PRECIP_PRE_P,STORAGE_VEC)
         ! NULLIFY(STORAGE_VEC)

          !air pressure
         ! VAR => FIND_VAR(PRECIP_FILE,"air_pressure",FOUND)
         ! IF(.not. FOUND) CALL FATAL_ERROR &
         !      & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
         !      & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
         !      & "COULD NOT FIND VARIABLE 'air_pressure'")

          ! MAKE SPACE FOR THE DATA FROM THE FILE
          !PRECIP_PAIR_N => reference_var(var)
          !ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          !IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
          !CALL NC_CONNECT_PVAR(PRECIP_PRE_N,STORAGE_VEC)
          !NULLIFY(STORAGE_VEC)

          !       ! MAKE SPACE FOR THE DATA FROM THE FILE
          !PRECIP_PAIR_P => reference_var(var)
          !ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          !IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
          !CALL NC_CONNECT_PVAR(PRECIP_PRE_P,STORAGE_VEC)
          !NULLIFY(STORAGE_VEC)

          !!relative humidity
          !VAR => FIND_VAR(PRECIP_FILE,"relative_humidity",FOUND)
          !IF(.not. FOUND) CALL FATAL_ERROR &
          !     & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
          !     & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
          !     & "COULD NOT FIND VARIABLE 'relative_humidity'")

          !! MAKE SPACE FOR THE DATA FROM THE FILE
          !PRECIP_RH_N => reference_var(var)
          !ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          !IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
          !CALL NC_CONNECT_PVAR(PRECIP_PRE_N,STORAGE_VEC)
          !NULLIFY(STORAGE_VEC)

           !      ! MAKE SPACE FOR THE DATA FROM THE FILE
          !PRECIP_RH_P => reference_var(var)
          !ALLOCATE(STORAGE_VEC(0:MT), stat = status)
          !IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
          !CALL NC_CONNECT_PVAR(PRECIP_PRE_P,STORAGE_VEC)
          !NULLIFY(STORAGE_VEC)

   
     ELSE
     !>---------ishid debug 20230104

       ! SETUP THE ACTUAL VARIABLES USED TO LOAD DATA!

       ! WIND STRESS IN THE X or EAST-WEST DIRECTION
       VAR => FIND_VAR(PRECIP_FILE,"precip",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            & "COULD NOT FIND VARIABLE 'precip'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       PRECIP_PRE_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       CALL NC_CONNECT_PVAR(PRECIP_PRE_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       PRECIP_PRE_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       CALL NC_CONNECT_PVAR(PRECIP_PRE_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)

       ! EVAP
       VAR => FIND_VAR(PRECIP_FILE,"evap",FOUND)
       IF(.not. FOUND) CALL FATAL_ERROR &
            & ("IN PRECIPITATION BOUNDARY CONDITION FILE OBJECT",&
            & "FILE NAME: "//TRIM(PRECIPITATION_FILE),&
            & "COULD NOT FIND VARIABLE 'evap'")

       ! MAKE SPACE FOR THE DATA FROM THE FILE
       PRECIP_EVP_N => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       CALL NC_CONNECT_PVAR(PRECIP_EVP_N,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)


       ! MAKE SPACE FOR THE DATA FROM THE FILE
       PRECIP_EVP_P => reference_var(var)
       ALLOCATE(STORAGE_VEC(0:MT), stat = status)
       IF(STATUS /= 0) CALL FATAL_ERROR("ALLOCATION ERROR IN SURFACE PRECIPITATION")
       CALL NC_CONNECT_PVAR(PRECIP_EVP_P,STORAGE_VEC)
       NULLIFY(STORAGE_VEC)
     !<----ishid debug
     ENDIF
     !>----ishid debug
       
!=====================================================================
    CASE DEFAULT
!=====================================================================
       CALL FATAL_ERROR("CAN NOT RECOGNIZE PRECIPITATION FILE TYPE")
!=====================================================================
    END SELECT
!=====================================================================


    ! ---------- new: 2016 , april, Karsten Lettmann after Hint by Qi and ayumi.fujisaki@noaa.gov------
    ! Initialize some variables 
    ! afm 20150914
    ! Need initialization. Otherwise, random values are asigned
    ! and cause a hanging problem of MPI job in UPDATE_VAR_BRACKET 
    ! This problem reported with Intel15.0.3. 
    ! PRECIPITATION
    PRECIP_PRE_P%curr_stkcnt = 0
    PRECIP_PRE_N%curr_stkcnt = 0
    write(IPT,*) "now we here"
    
    !<-------ishid debug
    IF (PRECIPITATION_BULK_ON) THEN
        ! UWIND
    !PRECIP_UWND_P%curr_stkcnt = 0
    !PRECIP_UWND_N%curr_stkcnt = 0
       ! VWND
    !PRECIP_VWND_P%curr_stkcnt = 0
    !PRECIP_VWND_N%curr_stkcnt = 0
       ! AIR PRESSEURE
    !PRECIP_PAIR_P%curr_stkcnt = 0
    !PRECIP_PAIR_N%curr_stkcnt = 0
       ! RELATIVE HUMIDITY
    !PRECIP_RH_P%curr_stkcnt = 0
    !PRECIP_RH_N%curr_stkcnt = 0
    !ELSE
    ! EVAPORATION
    ELSE
    PRECIP_EVP_N%curr_stkcnt = 0
    PRECIP_EVP_P%curr_stkcnt = 0
    ENDIF

    !<--------ishid dbg
    ! --------- end new ----------------------------------------------------------------

    IF(DBG_SET(DBG_SBR)) write(IPT,*) "END SURFACE_PRECIPITATION"
  END SUBROUTINE SURFACE_PRECIPITATION
  !==============================================================================|
  SUBROUTINE UPDATE_RIVERS(NOW,FLUX,TEMP,SALT,WQM,SED,BIO)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP), ALLOCATABLE :: FLUX(:)
    REAL(SP), ALLOCATABLE, OPTIONAL :: TEMP(:)
    REAL(SP), ALLOCATABLE, OPTIONAL :: SALT(:)
    REAL(SP), ALLOCATABLE, OPTIONAL :: WQM(:,:)
    REAL(SP), ALLOCATABLE, OPTIONAL :: SED(:,:)
    REAL(SP), ALLOCATABLE, OPTIONAL :: BIO(:,:)

    REAL(SP), POINTER :: VNP(:), VPP(:)

    REAL(SP), ALLOCATABLE :: CURRENT(:)
    TYPE(TIME) :: RIVTIME

    TYPE(NCFILE), POINTER :: NCF
    TYPE(NCVAR), POINTER :: VAR_N
    TYPE(NCVAR), POINTER :: VAR_P
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS, I, J, NRSF,IND,NS

    IF(.NOT. ALLOCATED(FLUX)) CALL FATAL_ERROR &
         &("THE RIVER FLUX VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

    IF(PRESENT(TEMP)) THEN
       IF(.NOT. ALLOCATED(TEMP)) CALL FATAL_ERROR &
            &("THE RIVER TEMP VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    END IF

    IF(PRESENT(SALT)) THEN
       IF(.NOT. ALLOCATED(SALT)) CALL FATAL_ERROR &
            &("THE RIVER SALT VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    END IF


    DO I = 1, SIZE(RIVER_FORCING) ! (NUMBER OF FILES)
       
       SELECT CASE (RIVER_KIND)
       CASE(PRDC)
          
          ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!          RIVTIME= NOW - RUNFILE_StartTime

          ! TO USE ZERO AS THE PHASE OF THE FORCING
          RIVTIME= NOW 
          
          
          RIVTIME = MOD(RIVTIME,RIVER_FORCING(I)%RIVER_PERIOD)
          
       CASE(VRBL)
          
          
          RIVTIME = NOW
       END SELECT
       

       NCF => RIVER_FORCING(I)%NCF
       FTM => NCF%FTIME

       NRSF = RIVER_FORCING(I)%RIVERS_IN_FILE

       ! RIVER FLUX
       VAR_N => RIVER_FORCING(I)%FLUX_N
       VAR_P => RIVER_FORCING(I)%FLUX_P
       CALL UPDATE_VAR_BRACKET(NCF,VAR_P,VAR_N,RIVTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE RIVER_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       ALLOCATE(CURRENT(NRSF))

       CALL NC_POINT_VAR(VAR_N,VNP)
       CALL NC_POINT_VAR(VAR_P,VPP)

       !====================================================
       ! Linear interpolation between time points
       !CURRENT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       !
       ! OR
       !
       ! Nearest time sets the value
       CURRENT = VNP
       if (FTM%PREV_WGHT .gt. 0.5_sp) CURRENT = VPP
       !====================================================

       DO J =1,NRSF
          IND = RIVER_FORCING(I)%RIV_FILE2LOC(J)
          IF(IND /= 0) FLUX(IND) = CURRENT(J)
       END DO

       DEALLOCATE(CURRENT)

       IF(PRESENT(SALT)) THEN

          ! RIVER SALT
          VAR_N => RIVER_FORCING(I)%SALT_N
          VAR_P => RIVER_FORCING(I)%SALT_P
          CALL UPDATE_VAR_BRACKET(NCF,VAR_P,VAR_N,RIVTIME,STATUS)
          IF (STATUS /= 0) THEN
             CALL FATAL_ERROR("COULD NOT UPATE RIVER_FILE TIME BRACKET: BOUNDS EXCEEDED?")
          end if

          ALLOCATE(CURRENT(NRSF))

          CALL NC_POINT_VAR(VAR_N,VNP)
          CALL NC_POINT_VAR(VAR_P,VPP)        

          !====================================================
          ! Linear interpolation between time points
          !CURRENT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
          !
          ! OR
          !
          ! Nearest time sets the value
          CURRENT = VNP
          if (FTM%PREV_WGHT .gt. 0.5_sp) current = VPP
          !====================================================


          DO J =1,NRSF
             IND = RIVER_FORCING(I)%RIV_FILE2LOC(J)
             IF(IND /= 0) SALT(IND) = CURRENT(J)
          END DO

          DEALLOCATE(CURRENT)
       END IF

       IF(PRESENT(TEMP)) THEN

          ! RIVER TEMP
          VAR_N => RIVER_FORCING(I)%TEMP_N
          VAR_P => RIVER_FORCING(I)%TEMP_P
          CALL UPDATE_VAR_BRACKET(NCF,VAR_P,VAR_N,RIVTIME,STATUS)
          IF (STATUS /= 0) THEN
             CALL FATAL_ERROR("COULD NOT UPATE RIVER_FILE TIME BRACKET: BOUNDS EXCEEDED?")
          end if

          ALLOCATE(CURRENT(NRSF))

          CALL NC_POINT_VAR(VAR_N,VNP)
          CALL NC_POINT_VAR(VAR_P,VPP)        

          !====================================================
          ! Linear interpolation between time points
          !CURRENT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
          !
          ! OR
          !
          ! Nearest time sets the value
          CURRENT = VNP
          if (FTM%PREV_WGHT .gt. 0.5_sp) current = VPP
          !====================================================

          DO J =1,NRSF
             IND = RIVER_FORCING(I)%RIV_FILE2LOC(J)
             IF(IND /= 0) TEMP(IND) = CURRENT(J)
          END DO

          DEALLOCATE(CURRENT)
       END IF




    END DO ! FOR EACH FILE

  END SUBROUTINE UPDATE_RIVERS
  !==============================================================================|
  SUBROUTINE UPDATE_GROUNDWATER(NOW,GW_FLUX,GW_TEMP,GW_SALT)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    TYPE(TIME)             :: GWTIME
    REAL(SP), ALLOCATABLE :: GW_FLUX(:)
    REAL(SP), ALLOCATABLE, OPTIONAL :: GW_SALT(:)
    REAL(SP), ALLOCATABLE, OPTIONAL :: GW_TEMP(:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS
    REAL(SP), POINTER :: VNP(:), VPP(:)

    IF(.NOT. ALLOCATED(GW_FLUX)) CALL FATAL_ERROR &
         &("THE GROUNDWATER FLUX VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")


!===================================================
    SELECT CASE(GROUNDWATER_KIND)
!===================================================
    CASE (CNSTNT)
       
       ! CONSTANT GROUND WATER FORCING IS ALWAYS A FLOW RATE (M/S)...
       ! CONVERT TO A FLUX
       GW_FLUX(1:MT) = GROUNDWATER_FLOW*ART1(1:MT)

       IF(GROUNDWATER_TEMP_ON .and. PRESENT(GW_TEMP)) THEN
          IF(.NOT. ALLOCATED(GW_TEMP)) CALL FATAL_ERROR &
               &("THE GROUNDWATER TEMPERATURE VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
          GW_TEMP(1:MT) = groundwater_temp
       END IF

       IF(GROUNDWATER_SALT_ON .and. PRESENT(GW_SALT)) THEN
          IF(.NOT. ALLOCATED(GW_SALT)) CALL FATAL_ERROR &
               &("THE GROUNDWATER SALINITY VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
          GW_SALT(1:MT) = groundwater_salt
       END IF

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC HEATING Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT HEATING Not Set Up Yet")

    CASE(PRDC)

       ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!      GWTIME= NOW - RUNFILE_StartTime
       
       ! TO USE ZERO AS THE PHASE OF THE FORCING
       GWTIME= NOW
       

       GWTIME = MOD(GWTIME,GWATER_PERIOD)

    CASE(VRBL)


       GWTIME = NOW
    END SELECT
!===================================================
!===================================================

!===================================================
    SELECT CASE(GWATER_FORCING_TYPE)
!===================================================
    CASE(GWATER_IS_FVCOMGRID)

       FTM => GWATER_FILE%FTIME

       ! GROUND WATER FLUX
       CALL UPDATE_VAR_BRACKET(GWATER_FILE,GWATER_FLUX_P,GWATER_FLUX_N,GWTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE GROUNDWATER_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(GWATER_FLUX_N,VNP)
       CALL NC_POINT_VAR(GWATER_FLUX_P,VPP)        

       GW_FLUX = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! IF THE GROUND WATER IS A FLOW RATE CONVERT TO A FLUX
       IF(GWATER_UNITS == GWATER_MS_1 ) GW_FLUX = GW_FLUX *ART1


       ! GROUND WATER TEMP
       IF(GROUNDWATER_TEMP_ON .and. PRESENT(GW_TEMP)) THEN

          IF(.NOT. ALLOCATED(GW_TEMP)) CALL FATAL_ERROR &
               &("THE GROUNDWATER TEMPERATURE VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

          CALL UPDATE_VAR_BRACKET(GWATER_FILE,GWATER_TEMP_P,GWATER_TEMP_N,GWTIME,STATUS)
          IF (STATUS /= 0) THEN
             CALL FATAL_ERROR("COULD NOT UPATE GROUNDWATER_FILE TIME BRACKET: BOUNDS EXCEEDED?")
          end if
          
          CALL NC_POINT_VAR(GWATER_TEMP_N,VNP)
          CALL NC_POINT_VAR(GWATER_TEMP_P,VPP)        
          
          GW_TEMP = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       END IF

       ! GROUND WATER SALT
       IF(GROUNDWATER_SALT_ON .and. PRESENT(GW_SALT)) THEN

          IF(.NOT. ALLOCATED(GW_SALT)) CALL FATAL_ERROR &
               &("THE GROUNDWATER SALINITY VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

          CALL UPDATE_VAR_BRACKET(GWATER_FILE,GWATER_SALT_P,GWATER_SALT_N,GWTIME,STATUS)
          IF (STATUS /= 0) THEN
             CALL FATAL_ERROR("COULD NOT UPATE GROUNDWATER_FILE TIME BRACKET: BOUNDS EXCEEDED?")
          end if
          
          CALL NC_POINT_VAR(GWATER_SALT_N,VNP)
          CALL NC_POINT_VAR(GWATER_SALT_P,VPP)        
          
          GW_SALT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       END IF



    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN GROUNDWATER_FORCING_TYPE IN UPDATE GROUNDWATER")
    END SELECT


  END SUBROUTINE UPDATE_GROUNDWATER
  !==============================================================================|
  !==============================================================================|
  !==============================================================================|
  SUBROUTINE UPDATE_HEAT_CALCULATED(NOW,HEAT_SWV,HEAT_NET,HEAT_SENSIBLE,HEAT_LATENT,HEAT_RLN,&
           Cwave,HSC1,TAU_WIND,U_star,Cd) ! Siqi Li, 2021-01-27
!  SUBROUTINE UPDATE_HEAT_CALCULATED(NOW,HEAT_SWV,HEAT_NET,HEAT_SENSIBLE,HEAT_LATENT,HEAT_RLN)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    TYPE(TIME)             :: HTIME
    REAL(SP), ALLOCATABLE :: HEAT_SWV(:), HEAT_NET(:)
    REAL(SP), ALLOCATABLE :: HEAT_SENSIBLE(:), HEAT_LATENT(:), HEAT_RLN(:)
    REAL(SP), ALLOCATABLE :: TAU_WIND(:), U_star(:), Cd(:), Cwave(:), HSC1(:) ! Siqi Li, 2021-01-27
    REAL(SP), DIMENSION(0:MT) :: TAU_node, USR_node, Cd_node  ! Siqi Li, 2021-01-27
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS,I,J
    REAL(SP), POINTER :: VNP(:), VPP(:)
    REAL(SP) :: WDS
    REAL(SP) :: HSB,HLB,TAU,USR,DTER  
    REAL(SP), DIMENSION(0:NT) :: WDSN
    REAL(SP) :: t1k,ulw_airf          ! ejw 8/16/2006

    IF(.NOT. ALLOCATED(HEAT_SWV)) CALL FATAL_ERROR &
         &("THE HEAT SHORTWAVE VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

    IF(.NOT. ALLOCATED(HEAT_NET)) CALL FATAL_ERROR &
         &("THE NET HEAT VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

    IF(.NOT. ALLOCATED(HEAT_SENSIBLE)) CALL FATAL_ERROR &
         &("THE SENSIBLE HEAT VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

    IF(.NOT. ALLOCATED(HEAT_LATENT)) CALL FATAL_ERROR &
         &("THE LATENT HEAT VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

    IF(.NOT. ALLOCATED(HEAT_RLN)) CALL FATAL_ERROR &
         &("THE LONGWAVE HEAT VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

    IF(WIND_TYPE /= 'speed')CALL FATAL_ERROR("WIND_TYPE must be 'speed' for heat flux calculating")
       
!===================================================
    SELECT CASE(HEATING_CALCULATE_KIND)
!===================================================
    CASE (CNSTNT)
    
       T_AIR(1:MT)   = AIR_TEMPERATURE
       RH_AIR(1:MT)  = RELATIVE_HUMIDITY
       PA_AIR(1:MT)  = SURFACE_PRESSURE
       DLW_AIR(1:MT) = LONGWAVE_RADIATION
       DSW_AIR(1:MT) = SHORTWAVE_RADIATION
       
       WDSN(1:NT) = SQRT(UUWIND(1:NT)*UUWIND(1:NT)+VVWIND(1:NT)*VVWIND(1:NT))

       DO I=1,MT  
         WDS = 0.0_SP
         DO J=1,NTVE(I)
           WDS = WDS + WDSN(NBVE(I,J))
         END DO
         WDS = WDS/FLOAT(NTVE(I))
  
         IF(TRIM(COARE_VERSION) == 'COARE26Z')THEN
!---> Siqi Li, 2021-01-27 : add Cd
           CALL COARE26Z(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1), &
                DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),USR,DTER,Cd_node(I))
!           CALL COARE26Z(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1),  &
!                DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),USR,DTER) !MDR include USR, DTER
           TAU_node(I) = TAU
           USR_node(I) = USR
!<--- Siqi Li, 2021-01-27

	   IF(.NOT. HEATING_FRESHWATER)THEN		
	     t1k=t1(i,1)+273.16_sp
             ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = DLW_AIR(I)+ulw_airf	     
	   ELSE
	     USRCOARE(I) = USR !MDR 3-19-2013 for wind wave mix
!             t1k=t1(i,1)+273.16_sp !MDR 4-24-2013
             t1k=t1(i,1)-DTER+273.16_sp !MDR 4-24-2013, include cool skin adjustment to surface temp
!MDR         ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
!MDR         HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             ulw_airf=-StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DSW_AIR(I)+HSB+HLB+emmiss*(DLW_AIR(I)+ulw_airf) !MDR, apply emmiss to NET long wave, as in COARE     
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = emmiss*(DLW_AIR(I)+ulw_airf)	     
	   END IF  
	 ELSE IF(TRIM(COARE_VERSION) == 'COARE40VN')THEN
           RAIN(I) = QPREC(I)*3600.0_SP*1000.0_SP  ! change unit to mm/hr  from m/s
	   ZI40 = fmiss
           CP40 = Cwave(I)
           SIGH = HSC1(I)
!---> Siqi Li, 2021-01-27
           CALL COARE40VN(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1),  &
                DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),ZI40(I),RAIN(I),CP40(I),SIGH(I),fmiss,USR,DTER,Cd_node(I)) 
           TAU_node(I) = TAU
           USR_node(I) = USR
!<--- Siqi Li, 2021-01-27

	   IF(.NOT. HEATING_FRESHWATER)THEN		
             t1k=t1(i,1)+273.16_sp
             ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = DLW_AIR(I)+ulw_airf	     
	   ELSE
             USRCOARE(I) = USR !MDR 3-19-2013 for wind wave mix
!             t1k=t1(i,1)+273.16_sp !MDR 4-24-2013
             t1k=t1(i,1)-DTER+273.16_sp !MDR 4-24-2013, include cool skin adjustment to surface temp
!MDR         ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
!MDR         HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             ulw_airf=-StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DSW_AIR(I)+HSB+HLB+emmiss*(DLW_AIR(I)+ulw_airf) !MDR, apply emmiss to NET long wave, as in COARE     
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = emmiss*(DLW_AIR(I)+ulw_airf)	     
	   END IF  
	 ELSE
	   CALL FATAL_ERROR("The value of COARE_VERSION should be 'COARE26Z' or 'COARE40VN'")
	 END IF  
       END DO

!---> Siqi Li, 2021-01-27
       DO I = 1, NT
         TAU_WIND(I) = (TAU_node(NV(I,1))+TAU_node(NV(I,2))+TAU_node(NV(I,3))) / 3.0_SP
         U_star(I) = (USR_node(NV(I,1))+USR_node(NV(I,2))+USR_node(NV(I,3))) / 3.0_SP
         Cd(I) = (Cd_node(NV(I,1))+Cd_node(NV(I,2))+Cd_node(NV(I,3))) / 3.0_SP
       END DO
!<--- Siqi Li, 2021-01-27

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC HEATING Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT HEATING Not Set Up Yet")

    CASE(PRDC)

       ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!      HTIME= NOW - RUNFILE_StartTime
       
       ! TO USE ZERO AS THE PHASE OF THE FORCING
       HTIME= NOW
       

       HTIME = MOD(HTIME,HEAT_PERIOD)

    CASE(VRBL)


       HTIME = NOW
    END SELECT
!===================================================
!===================================================


!===================================================
    SELECT CASE(HEAT_FORCING_TYPE)
!===================================================
    CASE(HEAT_IS_WRFGRID)

       FTM => HEAT_FILE%FTIME

!
       ! BULK AIR TEMPERATURE
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,T_AIR_P,T_AIR_N,HTIME,STATUS,HEAT_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(T_AIR_N,VNP)
       CALL NC_POINT_VAR(T_AIR_P,VPP)        

       T_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!
       ! RELATIVE HUMIDITY
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,RH_AIR_P,RH_AIR_N,HTIME,STATUS,HEAT_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(RH_AIR_N,VNP)
       CALL NC_POINT_VAR(RH_AIR_P,VPP)        

       RH_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!
       ! SURFACE PRESSURE
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,PA_AIR_P,PA_AIR_N,HTIME,STATUS,HEAT_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(PA_AIR_N,VNP)
       CALL NC_POINT_VAR(PA_AIR_P,VPP)        

       PA_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       PA_AIR = PA_AIR * Pair_unit_factor  ! Siqi Li, 2021-01-27
!
       ! DOWNWARD LONGWAVE RADIATION
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,DLW_AIR_P,DLW_AIR_N,HTIME,STATUS,HEAT_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(DLW_AIR_N,VNP)
       CALL NC_POINT_VAR(DLW_AIR_P,VPP)        

       DLW_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!
       ! DOWNWARD SHORTWAVE RADIATION
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,DSW_AIR_P,DSW_AIR_N,HTIME,STATUS,HEAT_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(DSW_AIR_N,VNP)
       CALL NC_POINT_VAR(DSW_AIR_P,VPP)        

       DSW_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!

       WDSN(1:NT) = SQRT(UUWIND(1:NT)*UUWIND(1:NT)+VVWIND(1:NT)*VVWIND(1:NT))

       DO I=1,M  
         WDS = 0.0_SP
         DO J=1,NTVE(I)
           WDS = WDS + WDSN(NBVE(I,J))
         END DO
         WDS = WDS/FLOAT(NTVE(I))
  
         IF(TRIM(COARE_VERSION) == 'COARE26Z')THEN
!---> Siqi Li, 2021-01-27 : add Cd
           CALL COARE26Z(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1), &
                DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),USR,DTER,Cd_node(I))
!           CALL COARE26Z(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1), &
!           &                DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),USR,DTER) !MDR include USR, DTER
           TAU_node(I) = TAU
           USR_node(I) = USR
!<--- Siqi Li, 2021-01-27

           IF(.NOT. HEATING_FRESHWATER)THEN
             t1k=t1(i,1)+273.16_sp
             ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = DLW_AIR(I)+ulw_airf	     
	   ELSE
             USRCOARE(I) = USR !MDR 3-19-2013 for wind wave mix
!             t1k=t1(i,1)+273.16_sp !MDR 4-24-2013
             t1k=t1(i,1)-DTER+273.16_sp !MDR 4-24-2013, include cool skin adjustment to surface temp
!MDR         ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
!MDR         HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             ulw_airf=-StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DSW_AIR(I)+HSB+HLB+emmiss*(DLW_AIR(I)+ulw_airf) !MDR, apply emmiss to NET long wave, as in COARE     
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = emmiss*(DLW_AIR(I)+ulw_airf)	     
	   END IF  
	 ELSE IF(TRIM(COARE_VERSION) == 'COARE40VN')THEN
           RAIN(I) = QPREC(I)*3600.0_SP*1000.0_SP  ! change unit to mm/hr  from m/s
	   ZI40 = fmiss
           CP40 = Cwave(I)
           SIGH = HSC1(I)
!---> Siqi Li, 2021-01-27
           CALL COARE40VN(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1),&
               DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),ZI40(I),RAIN(I),CP40(I),SIGH(I),fmiss,USR,DTER,Cd_node(I))
!           CALL COARE40VN(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1), &
!              DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),ZI40(I),RAIN(I),CP40(I),SIGH(I),fmiss,USR,DTER) !MDR include USR, DTER
           TAU_node(I) = TAU
           USR_node(I) = USR
!<--- Siqi Li, 2021-01-27

           IF(.NOT. HEATING_FRESHWATER)THEN
             t1k=t1(i,1)+273.16_sp
             ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = DLW_AIR(I)+ulw_airf	     
	   ELSE
	     USRCOARE(I) = USR !MDR 3-19-2013 for wind wave mix
!             t1k=t1(i,1)+273.16_sp !MDR 4-24-2013
             t1k=t1(i,1)-DTER+273.16_sp !MDR 4-24-2013, include cool skin adjustment to surface temp
!MDR         ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
!MDR         HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             ulw_airf=-StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DSW_AIR(I)+HSB+HLB+emmiss*(DLW_AIR(I)+ulw_airf) !MDR, apply emmiss to NET long wave, as in COARE     
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = emmiss*(DLW_AIR(I)+ulw_airf)	     
	   END IF  
	 ELSE
	   CALL FATAL_ERROR("The value of COARE_VERSION should be 'COARE26Z' or 'COARE40VN'")
	 END IF  
       END DO

    CASE(HEAT_IS_FVCOMGRID)

       FTM => HEAT_FILE%FTIME

!
       ! BULK AIR TEMPERATURE
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,T_AIR_P,T_AIR_N,HTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(T_AIR_N,VNP)
       CALL NC_POINT_VAR(T_AIR_P,VPP)        

       T_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!
       ! RELATIVE HUMIDITY
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,RH_AIR_P,RH_AIR_N,HTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(RH_AIR_N,VNP)
       CALL NC_POINT_VAR(RH_AIR_P,VPP)        

       RH_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!
       ! SURFACE PRESSURE
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,PA_AIR_P,PA_AIR_N,HTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(PA_AIR_N,VNP)
       CALL NC_POINT_VAR(PA_AIR_P,VPP)        

       PA_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       PA_AIR = PA_AIR * Pair_unit_factor  ! Siqi Li, 2021-01-27
!
       ! DOWNWARD LONGWAVE RADIATION
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,DLW_AIR_P,DLW_AIR_N,HTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(DLW_AIR_N,VNP)
       CALL NC_POINT_VAR(DLW_AIR_P,VPP)        

       DLW_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!
       ! DOWNWARD SHORTWAVE RADIATION
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,DSW_AIR_P,DSW_AIR_N,HTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE HEAT_FILE TIME BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(DSW_AIR_N,VNP)
       CALL NC_POINT_VAR(DSW_AIR_P,VPP)        

       DSW_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!

       WDSN(1:NT) = SQRT(UUWIND(1:NT)*UUWIND(1:NT)+VVWIND(1:NT)*VVWIND(1:NT))

       DO I=1,M  
         WDS = 0.0_SP
         DO J=1,NTVE(I)
           WDS = WDS + WDSN(NBVE(I,J))
         END DO
         WDS = WDS/FLOAT(NTVE(I))
  
         IF(TRIM(COARE_VERSION) == 'COARE26Z')THEN
!---> Siqi Li, 2021-01-27 : add Cd
           CALL COARE26Z(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1), &
                DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),USR,DTER,Cd_node(I))
!           CALL COARE26Z(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1),&
!           &       DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),USR,DTER) !MDR include USR, DTER
           TAU_node(I) = TAU
           USR_node(I) = USR
!<--- Siqi Li, 2021-01-27

           IF(.NOT. HEATING_FRESHWATER)THEN	      
             t1k=t1(i,1)+273.16_sp
             ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             HEAT_SWV(I) = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = DLW_AIR(I)+ulw_airf	     
	   ELSE
             USRCOARE(I) = USR !MDR 3-19-2013 for wind wave mix
!             t1k=t1(i,1)+273.16_sp !MDR 4-24-2013
             t1k=t1(i,1)-DTER+273.16_sp !MDR 4-24-2013, include cool skin adjustment to surface temp
!MDR         ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
!MDR         HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             ulw_airf=-StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DSW_AIR(I)+HSB+HLB+emmiss*(DLW_AIR(I)+ulw_airf) !MDR, apply emmiss to NET long wave, as in COARE     
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  =DLW_AIR(I)
             HEAT_RLN(I)  =emmiss*(DLW_AIR(I)+ulw_airf)	     
	   END IF  
	 ELSE IF(TRIM(COARE_VERSION) == 'COARE40VN')THEN
           RAIN(I) = QPREC(I)*3600.0_SP*1000.0_SP  ! change unit to mm/hr  from m/s
	   ZI40 = fmiss
           CP40 = Cwave(I)
           SIGH = HSC1(I)
!---> Siqi Li, 2021-01-27
           CALL COARE40VN(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1),&
                DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),ZI40(I),RAIN(I),CP40(I),SIGH(I),fmiss,USR,DTER,Cd_node(I))
!           CALL COARE40VN(WDS,ZUU,T_AIR(I),ZTT,RH_AIR(I),ZQQ,PA_AIR(I),T1(I,1), &
!                DLW_AIR(I),DSW_AIR(I),TAU,HSB,HLB,CORR(I),ZI40(I),RAIN(I),CP40(I),SIGH(I),fmiss,USR,DTER) !MDR include USR, DTER
           TAU_node(I) = TAU
           USR_node(I) = USR
!<--- Siqi Li, 2021-01-27


           IF(.NOT. HEATING_FRESHWATER)THEN	      
             t1k=t1(i,1)+273.16_sp
             ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             HEAT_SWV(I) = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  = DLW_AIR(I)
             HEAT_RLN(I)  = DLW_AIR(I)+ulw_airf	     
	   ELSE
             USRCOARE(I) = USR !MDR 3-19-2013 for wind wave mix
!             t1k=t1(i,1)+273.16_sp !MDR 4-24-2013
             t1k=t1(i,1)-DTER+273.16_sp !MDR 4-24-2013, include cool skin adjustment to surface temp
!MDR         ulw_airf=-emmiss*StefBo*t1k*t1k*t1k*t1k
!MDR         HEAT_NET(I) = DLW_AIR(I)+DSW_AIR(I)+HSB+HLB+ulw_airf        ! ejw 8/16/2006 based on ROMS
             ulw_airf=-StefBo*t1k*t1k*t1k*t1k
             HEAT_NET(I) = DSW_AIR(I)+HSB+HLB+emmiss*(DLW_AIR(I)+ulw_airf) !MDR, apply emmiss to NET long wave, as in COARE     
             HEAT_SWV(I)  = DSW_AIR(I)
	     !EJA edit - output sensible and latent heat flux
	     HEAT_SENSIBLE(I) = HSB
	     HEAT_LATENT(I) = HLB
!g.l. 10/10/2017         HEAT_RLN(I)  =DLW_AIR(I)
             HEAT_RLN(I)  =emmiss*(DLW_AIR(I)+ulw_airf)	     
	   END IF  
	 ELSE
	   CALL FATAL_ERROR("The value of COARE_VERSION should be 'COARE26Z' or 'COARE40VN'")
	 END IF  
       END DO

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN HEAT_FORCING_TYPE IN UPDATE HEAT")
    END SELECT

!---> Siqi Li, 2021-01-27
!#      if defined (MULTIPROCESSOR)
!       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,TAU_node)
!       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,USR_node)
!#      endif
    DO I = 1, N
       TAU_WIND(I) = (TAU_node(NV(I,1))+TAU_node(NV(I,2))+TAU_node(NV(I,3))) / 3.0_SP
       U_star(I) = (USR_node(NV(I,1))+USR_node(NV(I,2))+USR_node(NV(I,3))) / 3.0_SP
       Cd(I) = (Cd_node(NV(I,1))+Cd_node(NV(I,2))+Cd_node(NV(I,3))) / 3.0_SP
    END DO
!<--- Siqi Li, 2021-01-27


    RETURN
  END SUBROUTINE UPDATE_HEAT_CALCULATED
  !==============================================================================|
  !==============================================================================|
  !==============================================================================|

  SUBROUTINE UPDATE_WIND(NOW,wstrx,wstry)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    TYPE(TIME)             :: WTIME
    REAL(SP), ALLOCATABLE :: wstrx(:),wstry(:)
    REAL(SP), POINTER :: VNP(:), VPP(:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS


    IF(.NOT. ALLOCATED(wstrx)) CALL FATAL_ERROR &
         &("THE WIND VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(wstry)) CALL FATAL_ERROR &
         &("THE WIND VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")


!===================================================
    SELECT CASE(WIND_KIND)
!===================================================
    CASE (CNSTNT)

          wstrx(1:NT) = WIND_X
          wstry(1:NT) = WIND_Y

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC WIND Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT WIND Not Set Up Yet")

    CASE(PRDC)

       ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!      WTIME= NOW - RUNFILE_StartTime
       
       ! TO USE ZERO AS THE PHASE OF THE FORCING
       WTIME= NOW
       

       WTIME = MOD(WTIME,WINDS_PERIOD)

    CASE(VRBL)

       WTIME = NOW
    END SELECT
!===================================================
!===================================================


!===================================================
    SELECT CASE(WINDS_FORCING_TYPE)
!===================================================
    CASE(WINDS_ARE_WRFGRID)

       FTM => WINDS_FILE%FTIME

       ! THE X DIRECTION WIND STRESS
       CALL UPDATE_VAR_BRACKET(WINDS_FILE,WINDS_STRX_P,WINDS_STRX_N,WTIME,STATUS,WINDS_INTP_C)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND  X BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(WINDS_STRX_N,VNP)
       CALL NC_POINT_VAR(WINDS_STRX_P,VPP)   
       WSTRX = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! THE Y DIRECTION WIND STRESS
       CALL UPDATE_VAR_BRACKET(WINDS_FILE,WINDS_STRY_P,WINDS_STRY_N,WTIME,STATUS,WINDS_INTP_C)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND Y BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(WINDS_STRY_N,VNP)
       CALL NC_POINT_VAR(WINDS_STRY_P,VPP)   
       WSTRY = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!===================================================
    CASE(WINDS_ARE_FVCOMGRID)
!===================================================
       FTM => WINDS_FILE%FTIME

       ! THE X DIRECTION WIND STRESS
       CALL UPDATE_VAR_BRACKET(WINDS_FILE,WINDS_STRX_P,WINDS_STRX_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND X BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(WINDS_STRX_N,VNP)
       CALL NC_POINT_VAR(WINDS_STRX_P,VPP)   
       WSTRX = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! THE Y DIRECTION WIND STRESS
       CALL UPDATE_VAR_BRACKET(WINDS_FILE,WINDS_STRY_P,WINDS_STRY_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND Y BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(WINDS_STRY_N,VNP)
       CALL NC_POINT_VAR(WINDS_STRY_P,VPP)   
       WSTRY = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

!===================================================
    CASE(WINDS_ARE_PT_SOURCE)
!===================================================
       FTM => WINDS_FILE%FTIME

       ! THE X DIRECTION WIND STRESS
       CALL UPDATE_VAR_BRACKET(WINDS_FILE,WINDS_STRX_P,WINDS_STRX_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND X BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(WINDS_STRX_N,VNP)
       CALL NC_POINT_VAR(WINDS_STRX_P,VPP)
       WSTRX(1:NT) =   FTM%NEXT_WGHT * VNP(1) + FTM%PREV_WGHT * VPP(1)

       ! THE Y DIRECTION WIND STRESS
       CALL UPDATE_VAR_BRACKET(WINDS_FILE,WINDS_STRY_P,WINDS_STRY_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND Y BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(WINDS_STRY_N,VNP)
       CALL NC_POINT_VAR(WINDS_STRY_P,VPP)
       WSTRY(1:NT) = FTM%NEXT_WGHT * VNP(1) + FTM%PREV_WGHT * VPP(1)

!===================================================
    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN WINDS_FORCING_TYPE IN UPDATE WIND")
    END SELECT
!===================================================

  END SUBROUTINE UPDATE_WIND








!==============================================================================|
  SUBROUTINE UPDATE_PRECIPITATION(NOW,Qprec,Qevap)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    TYPE(TIME)             :: PTIME
    REAL(SP), ALLOCATABLE :: QEVAP(:),QPREC(:)
    !<-------ishid debug 20230104
    REAL(SP), ALLOCATABLE :: PAIR(:),RH(:)!,QS(:),ES(:)!,U(:),UWND(:),VWND(:),
    REAL(SP) :: Ce
    REAL(SP) :: WDS,ES,QS
    REAL(SP), DIMENSION(0:NT) :: WDSN
    INTEGER :: I,J
    !>-------ishid debug 20230104
    REAL(SP), POINTER :: VNP(:), VPP(:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS

    IF(.NOT. ALLOCATED(Qprec)) CALL FATAL_ERROR &
         &("THE PRECIPITATION VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(Qevap)) CALL FATAL_ERROR &
         &("THE EVAPORATION VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

!===================================================
    SELECT CASE(PRECIPITATION_KIND)
!===================================================
    CASE (CNSTNT)
       
       Qprec(1:MT) = PRECIPITATION_PRC
       Qevap(1:MT) = PRECIPITATION_EVP

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC PRECIP Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT PRECIP Not Set Up Yet")

    CASE(PRDC)

       ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!      PTIME= NOW - RUNFILE_StartTime
       
       ! TO USE ZERO AS THE PHASE OF THE FORCING
       PTIME= NOW
       

       PTIME = MOD(PTIME,PRECIP_PERIOD)

    CASE(VRBL)

       PTIME = NOW
    END SELECT
!===================================================
!===================================================




!===================================================
    SELECT CASE(PRECIP_FORCING_TYPE)
!===================================================
    CASE(PRECIP_IS_WRFGRID)

       FTM => PRECIP_FILE%FTIME
       
       ! PRECIPITATION
       CALL UPDATE_VAR_BRACKET(PRECIP_FILE,PRECIP_PRE_P,PRECIP_PRE_N,PTIME,STATUS,PRECIP_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE PRECIP BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(PRECIP_PRE_N,VNP)
       CALL NC_POINT_VAR(PRECIP_PRE_P,VPP)   
       Qprec = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       
       ! EVAPORATION
       CALL UPDATE_VAR_BRACKET(PRECIP_FILE,PRECIP_EVP_P,PRECIP_EVP_N,PTIME,STATUS,PRECIP_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE EVAP BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(PRECIP_EVP_N,VNP)
       CALL NC_POINT_VAR(PRECIP_EVP_P,VPP)   
       Qevap = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!===================================================
    CASE(PRECIP_IS_FVCOMGRID)
!===================================================

       FTM => PRECIP_FILE%FTIME
       
       ! PRECIPITATION
       CALL UPDATE_VAR_BRACKET(PRECIP_FILE,PRECIP_PRE_P,PRECIP_PRE_N,PTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE PRECIP BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(PRECIP_PRE_N,VNP)
       CALL NC_POINT_VAR(PRECIP_PRE_P,VPP)   
       Qprec = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       !<------ishid debug 20230104
       IF (PRECIPITATION_BULK_ON) THEN
          !uwnd
          !CALL UPDATE_VAR_BRACKET(PRECIP_FILE,PRECIP_UWND_P,PRECIP_UWND_N,PTIME,STATUS)
          !IF (STATUS /= 0) THEN
          !     CALL FATAL_ERROR("COULD NOT UPATE PRECIP BRACKET: BOUNDS EXCEEDED?")
          !end if
          
          !CALL NC_POINT_VAR(PRECIP_UWND_N,VNP)
          !CALL NC_POINT_VAR(PRECIP_UWND_P,VPP)   
          !uwnd = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

          !vwnd
          !CALL UPDATE_VAR_BRACKET(PRECIP_FILE,PRECIP_VWND_P,PRECIP_VWND_N,PTIME,STATUS)
          !IF (STATUS /= 0) THEN
          !     CALL FATAL_ERROR("COULD NOT UPATE PRECIP BRACKET: BOUNDS EXCEEDED?")
          !end if
          
          !CALL NC_POINT_VAR(PRECIP_VWND_N,VNP)
          !CALL NC_POINT_VAR(PRECIP_VWND_P,VPP)   
          !Vwnd = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

          !air pressure
          !CALL UPDATE_VAR_BRACKET(PRECIP_FILE,PRECIP_PAIR_P,PRECIP_PAIR_N,PTIME,STATUS)
          !IF (STATUS /= 0) THEN
          !     CALL FATAL_ERROR("COULD NOT UPATE PRECIP BRACKET: BOUNDS EXCEEDED?")
          !end if
          
          !CALL NC_POINT_VAR(PRECIP_PAIR_N,VNP)
          !CALL NC_POINT_VAR(PRECIP_PAIR_P,VPP)   
          !PAIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

          !relative humidity
          !CALL UPDATE_VAR_BRACKET(PRECIP_FILE,PRECIP_RH_P,PRECIP_RH_N,PTIME,STATUS)
          !IF (STATUS /= 0) THEN
          !     CALL FATAL_ERROR("COULD NOT UPATE PRECIP BRACKET: BOUNDS EXCEEDED?")
          !end if
          
          !CALL NC_POINT_VAR(PRECIP_RH_N,VNP)
          !CALL NC_POINT_VAR(PRECIP_RH_P,VPP)   
          !RH = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP


          !bulk 

          WDSN(1:NT) = SQRT(UUWIND(1:NT)*UUWIND(1:NT)+VVWIND(1:NT)*VVWIND(1:NT))

          DO I=1,M  
            WDS = 0.0_SP
            ES = 0.0_SP
            QS = 0.0_SP
            DO J=1,NTVE(I)
              WDS = WDS + WDSN(NBVE(I,J))
            END DO
            WDS = WDS/FLOAT(NTVE(I)) !WDS:wind speed at node from coare26Z 周囲の格子の風速を平均している。と言っても一様なので関係ない。
            ES = 6.1078_SP * 10.0_SP **(7.5_SP*(T1(I,1))/237.3_SP+T1(I,1))
            QS = 0.98_SP * 0.622_SP * (ES/PA_AIR(I)) / (1._SP - 0.378_SP * (ES/PA_AIR(I)))
            !Qevap(I) =10._SP**-10._SP** 1.293_SP * 1.2_SP*10._SP**-3._SP * WDS* (QS - RH_AIR(I))
            Qevap(I) = 0
          ENDDO

          
          

       ELSE
       !>------ishid debug 20230104
       ! EVAPORATION
       CALL UPDATE_VAR_BRACKET(PRECIP_FILE,PRECIP_EVP_P,PRECIP_EVP_N,PTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE EVAP BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(PRECIP_EVP_N,VNP)
       CALL NC_POINT_VAR(PRECIP_EVP_P,VPP)   
       Qevap = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       !<--------ishid debug20230104
       ENDIF
       !>--------ishid debug20230104
    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN WINDS_FORCING_TYPE IN UPDATE PRECIPITATION")
    END SELECT


  END SUBROUTINE UPDATE_PRECIPITATION
  !==============================================================================|


!==============================================================================|
  SUBROUTINE UPDATE_WAVE(NOW,WHS,WDIR,WPER,WLENGTH,WPER_BOT,WUB_BOT)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    TYPE(TIME)             :: PTIME
    REAL(SP), ALLOCATABLE  :: WHS(:),WDIR(:),WPER(:),WLENGTH(:),WPER_BOT(:),WUB_BOT(:)
    REAL(SP), POINTER :: VNP(:), VPP(:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS

    REAL    :: X1,X2,Y1,Y2,X0,Y0,ANGLE
    INTEGER :: I

    real:: a2, K2,H2,T2,Ub2,w2
    
    IF(.NOT. ALLOCATED(WHS)) CALL FATAL_ERROR &
         &("THE WAVE HEIGHT VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(WDIR)) CALL FATAL_ERROR &
         &("THE WAVE DIRECTION VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(WPER)) CALL FATAL_ERROR &
         &("THE WAVE PERIOD VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(WLENGTH)) CALL FATAL_ERROR &
         &("THE WAVE LENGTH VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(WPER_BOT)) CALL FATAL_ERROR &
         &("THE BOTTOM WAVE PERIOD VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(WUB_BOT)) CALL FATAL_ERROR &
         &("THE BOTTOM WAVE ORBITAL VELOCITY VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

!===================================================
    SELECT CASE(WAVE_KIND)
!===================================================
    CASE (CNSTNT)
       
       WHS(1:MT)      = WAVE_HEIGHT
       WDIR(1:MT)     = WAVE_DIRECTION
       WPER(1:MT)     = WAVE_PERIOD
       WLENGTH(1:MT)  = WAVE_LENGTH
       WPER_BOT(1:MT) = WAVE_PER_BOT
       WUB_BOT(1:MT)  = WAVE_UB_BOT

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC PRECIP Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT PRECIP Not Set Up Yet")

    CASE(PRDC)

       ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!      PTIME= NOW - RUNFILE_StartTime
       
       ! TO USE ZERO AS THE PHASE OF THE FORCING
       PTIME= NOW
       

       PTIME = MOD(PTIME,PRECIP_PERIOD)

    CASE(VRBL)

       PTIME = NOW
    END SELECT
!===================================================
!===================================================




!===================================================
    SELECT CASE(WAVES_FORCING_TYPE)
!===================================================
    CASE(WAVES_ARE_WRFGRID)

!===================================================
    CASE(WAVES_ARE_FVCOMGRID)
!===================================================

       FTM => WAVES_FILE%FTIME
       
       ! WAVE HEIGHT
       CALL UPDATE_VAR_BRACKET(WAVES_FILE,WAVES_HEIGHT_P,WAVES_HEIGHT_N,PTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WAVE HEIGHT BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(WAVES_HEIGHT_N,VNP)
       CALL NC_POINT_VAR(WAVES_HEIGHT_P,VPP)   
       WHS = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! WAVE DIRECTIION
       CALL UPDATE_VAR_BRACKET(WAVES_FILE,WAVES_DIRECTION_P,WAVES_DIRECTION_N,PTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WAVE DIRECTION BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(WAVES_DIRECTION_N,VNP)
       CALL NC_POINT_VAR(WAVES_DIRECTION_P,VPP)   
       WDIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       
       DO I=1,MT
          X1 = COS(VNP(I)*3.1415926/180.0)
          X2 = COS(VPP(I)*3.1415926/180.0)
          Y1 = SIN(VNP(I)*3.1415926/180.0)
          Y2 = SIN(VPP(I)*3.1415926/180.0)
          X0 = FTM%NEXT_WGHT * X1 + FTM%PREV_WGHT * X2
          Y0 = FTM%NEXT_WGHT * Y1 + FTM%PREV_WGHT * Y2
          ANGLE = ATAN2(Y0,X0)
          IF(ANGLE<0)ANGLE = ANGLE + 3.1415926*2.0
          WDIR(I) = ANGLE*180.0/3.1415926
       END DO
       

       ! WAVE LENGTH
       CALL UPDATE_VAR_BRACKET(WAVES_FILE,WAVES_LENGTH_P,WAVES_LENGTH_N,PTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WAVE LENGTH BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(WAVES_LENGTH_N,VNP)
       CALL NC_POINT_VAR(WAVES_LENGTH_P,VPP)   
       WLENGTH = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP


       ! WAVE PERIOD
       CALL UPDATE_VAR_BRACKET(WAVES_FILE,WAVES_PERIOD_P,WAVES_PERIOD_N,PTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WAVE PERIOD BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(WAVES_PERIOD_N,VNP)
       CALL NC_POINT_VAR(WAVES_PERIOD_P,VPP)   
       WPER = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       

       ! BOTTOM WAVE PERIOD 
       CALL UPDATE_VAR_BRACKET(WAVES_FILE,WAVES_PER_BOT_P,WAVES_PER_BOT_N,PTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE BOTTOM WAVE PERIOD BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(WAVES_PER_BOT_N,VNP)
       CALL NC_POINT_VAR(WAVES_PER_BOT_P,VPP)   
       WPER_BOT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP


       ! BOTTOM ORBITAL VELOCITY
       CALL UPDATE_VAR_BRACKET(WAVES_FILE,WAVES_UB_BOT_P,WAVES_UB_BOT_N,PTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE BOTTOM ORBITAL VELOCITY BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(WAVES_UB_BOT_N,VNP)
       CALL NC_POINT_VAR(WAVES_UB_BOT_P,VPP)   
       WUB_BOT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
!       print*,whs(36),wper(36),wlength(36)
!---------------------------------------------------------------------!
!    Test bottom orbital velocity 
!---------------------------------------------------------------------!
!      a2=0.5*whs(159)
!      K2=2.0*3.1415926/wlength(159)
!      H2=10.0
!      T2=wper(159)
!      w2=sqrt(9.8*K2*SINH(K2*H2)/COSH(K2*H2))
!      Ub2=a2*W2/sinh(H2*K2)
!     print*,whs(159),ub2
       

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN WAVES_FORCING_TYPE IN UPDATE WAVE")
    END SELECT


  END SUBROUTINE UPDATE_WAVE
  !==============================================================================|




!==============================================================================|
  SUBROUTINE UPDATE_AIRPRESSURE(NOW,PA_AIR)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    TYPE(TIME)             :: ATIME
    REAL(SP), ALLOCATABLE :: PA_AIR(:)
    REAL(SP), POINTER :: VNP(:), VPP(:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS

    IF(.NOT. ALLOCATED(PA_AIR)) CALL FATAL_ERROR &
         &("THE AIR PRESSURE VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

!===================================================
    SELECT CASE(AIRPRESSURE_KIND)
!===================================================
    CASE (CNSTNT)
       
       PA_AIR(1:MT) = AIRPRESSURE_VALUE

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC AIR PRESSURE Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT AIR PRESSURE Not Set Up Yet")

    CASE(PRDC)

       ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!      ATIME= NOW - RUNFILE_StartTime
       
       ! TO USE ZERO AS THE PHASE OF THE FORCING
       ATIME= NOW
       

       ATIME = MOD(ATIME,AIRPRESSURE_PERIOD)

    CASE(VRBL)

       ATIME = NOW
    END SELECT
!===================================================
!===================================================




!===================================================
    SELECT CASE(AIRPRESSURE_FORCING_TYPE)
!===================================================
    CASE(AIRPRESSURE_IS_WRFGRID)

       FTM => AIRPRESSURE_P_FILE%FTIME
       
       ! AIR PRESSURE
       CALL UPDATE_VAR_BRACKET(AIRPRESSURE_P_FILE,AIR_PRESSURE_P,AIR_PRESSURE_N,ATIME,STATUS,AIRPRESSURE_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE AIR PRESSURE BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(AIR_PRESSURE_N,VNP)
       CALL NC_POINT_VAR(AIR_PRESSURE_P,VPP)   
       PA_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       PA_AIR = PA_AIR * Pair_unit_factor  ! Siqi Li, 2021-01-27
       
!===================================================
    CASE(AIRPRESSURE_IS_FVCOMGRID)
!===================================================

       FTM => AIRPRESSURE_P_FILE%FTIME
       
       ! AIR PRESSURE
       CALL UPDATE_VAR_BRACKET(AIRPRESSURE_P_FILE,AIR_PRESSURE_P,AIR_PRESSURE_N,ATIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE AIR PRESSURE BRACKET: BOUNDS EXCEEDED?")
       end if
       
       CALL NC_POINT_VAR(AIR_PRESSURE_N,VNP)
       CALL NC_POINT_VAR(AIR_PRESSURE_P,VPP)   
       PA_AIR = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP
       PA_AIR = PA_AIR * Pair_unit_factor  ! Siqi Li, 2021-01-27 
    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN AIRPRESSURE_FORCING_TYPE IN UPDATE AIR PRESSURE")
    END SELECT


  END SUBROUTINE UPDATE_AIRPRESSURE
  !==============================================================================|
  SUBROUTINE UPDATE_TIDE(NOW,BND_ELV)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP), ALLOCATABLE :: BND_ELV(:)
    REAL(SP), POINTER :: VNP(:), VPP(:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS

    IF(.NOT. ALLOCATED(BND_ELV)) CALL FATAL_ERROR &
         &("THE BOUNDARY ELEVATION VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")


    SELECT CASE(TIDE_FORCING_TYPE)
    CASE(TIDE_FORCING_TIMESERIES)

       FTM => TIDE_FILE%FTIME

       ! PRECIPITATION
       CALL UPDATE_VAR_BRACKET(TIDE_FILE,TIDE_ELV_P,TIDE_ELV_N,NOW,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE TIDE ELVATION BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(TIDE_ELV_N,VNP)
       CALL NC_POINT_VAR(TIDE_ELV_P,VPP)   
       BND_ELV = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN TIDAL FORCING FILE TYPE IN UPDATE_TIDE")
    END SELECT


  END SUBROUTINE UPDATE_TIDE
  !==============================================================================|
  SUBROUTINE UPDATE_OBC_SALT(NOW,SALT)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP), ALLOCATABLE :: SALT(:,:)
    REAL(SP), POINTER :: VNP(:,:), VPP(:,:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS

    IF(.NOT. ALLOCATED(SALT)) CALL FATAL_ERROR &
         &("THE BOUNDARY SALINITY VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

    SELECT CASE(OBC_S_TYPE)
    CASE(OBC_S_SIGMA)

       FTM => OBC_S_FILE%FTIME

       ! OBC_SALT
       CALL UPDATE_VAR_BRACKET(OBC_S_FILE,OBC_S_P,OBC_S_N,NOW,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE OBC SALINITY BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(OBC_S_N,VNP)
       CALL NC_POINT_VAR(OBC_S_P,VPP)   
       SALT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN OBC SALINITY FILE TYPE IN UPDATE_OBC_SALT")
    END SELECT


  END SUBROUTINE UPDATE_OBC_SALT
  !==============================================================================|
  SUBROUTINE UPDATE_OBC_TEMP(NOW,TEMP)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    REAL(SP), ALLOCATABLE :: TEMP(:,:)
    REAL(SP), POINTER :: VNP(:,:), VPP(:,:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS

    IF(.NOT. ALLOCATED(TEMP)) CALL FATAL_ERROR &
         &("THE BOUNDARY TEMPERATURE VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")

    SELECT CASE(OBC_T_TYPE)
    CASE(OBC_T_SIGMA)

       FTM => OBC_T_FILE%FTIME

       ! PRECIPITATION
       CALL UPDATE_VAR_BRACKET(OBC_T_FILE,OBC_T_P,OBC_T_N,NOW,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE OBC TEMPERATURE BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(OBC_T_N,VNP)
       CALL NC_POINT_VAR(OBC_T_P,VPP)   
       TEMP = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN OBC TEMPERATURE FILE TYPE IN UPDATE_OBC_TEMP")
    END SELECT


  END SUBROUTINE UPDATE_OBC_TEMP
  !==============================================================================|
  !==============================================================================|
  !==============================================================================|
  !==============================================================================|
  SUBROUTINE UPDATE_ICE(NOW,SAT,SWV,SPQ,CLD)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    TYPE(TIME)             :: WTIME
    REAL(SP), ALLOCATABLE :: SAT(:)
    REAL(SP), ALLOCATABLE :: SWV(:)
    REAL(SP), ALLOCATABLE :: SLP(:)
    REAL(SP), ALLOCATABLE :: SPQ(:)
    REAL(SP), ALLOCATABLE :: CLD(:)
    REAL(SP), POINTER :: VNP(:), VPP(:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS

    IF(.NOT. ALLOCATED(SAT)) CALL FATAL_ERROR &
         &("THE Sea Surface Air Temperature VARIABLE PASSED TO UPDATE ICE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(SWV)) CALL FATAL_ERROR &
         &("THE SHORTWAVE RADIATION VARIABLE PASSED TO UPDATE ICE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(SPQ)) CALL FATAL_ERROR &
         &("THE SPECIFIC HUMIDIY VARIABLE PASSED TO UPDATE ICE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(CLD)) CALL FATAL_ERROR &
         &("THE CLOUD COVER VARIABLE PASSED TO UPDATE ICE IS NOT ALLOCATED")
   
!===================================================
    SELECT CASE(ICE_FORCING_KIND)
!===================================================
    CASE (CNSTNT)
       
       SAT(1:MT) = ICE_AIR_TEMP
       SPQ(1:MT) = ICE_SPEC_HUMIDITY
       CLD(1:MT) = ICE_CLOUD_COVER
       SWV(1:MT) = ICE_SHORTWAVE

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC ICE FORCING Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT ICE FORCING Not Set Up Yet")

    CASE(PRDC)

       ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!      WTIME= NOW - RUNFILE_StartTime
       
       ! TO USE ZERO AS THE PHASE OF THE FORCING
       WTIME= NOW
       

       WTIME = MOD(WTIME,ICE_PERIOD)

    CASE(VRBL)

       WTIME = NOW
    END SELECT

!===================================================
    SELECT CASE(ICE_FORCING_TYPE)
!===================================================
    CASE(ICE_IS_WRFGRID)
!===================================================

       FTM => ICE_FILE%FTIME

       ! THE SEA SURFACE AIR TEMP
       CALL UPDATE_VAR_BRACKET(ICE_FILE,ICE_SAT_P,ICE_SAT_N,WTIME,STATUS,ICE_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE ICE Surface Air Temp BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICE_SAT_N,VNP)
       CALL NC_POINT_VAR(ICE_SAT_P,VPP)   
       SAT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! SHORT WAVE
! afm 20151112 & EJA 20160921 - commented out for solar
! With SOLAR, use SOLAR-derived shortwave
! Without SOLAR, use shortwave from forcing data
!       CALL UPDATE_VAR_BRACKET(ICE_FILE,ICE_SWV_P,ICE_SWV_N,WTIME,STATUS,ICE_INTP_N)
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,ICE_SWV_P,ICE_SWV_N,WTIME,STATUS,ICE_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPDATE ICE SHORTWAVE BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICE_SWV_N,VNP)
       CALL NC_POINT_VAR(ICE_SWV_P,VPP)   
       SWV = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! THE SPECIFIC HUMIDITY
       CALL UPDATE_VAR_BRACKET(ICE_FILE,ICE_SPQ_P,ICE_SPQ_N,WTIME,STATUS,ICE_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE ICE SPECIFIC HUMIDITY BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICE_SPQ_N,VNP)
       CALL NC_POINT_VAR(ICE_SPQ_P,VPP) 
       SPQ = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP


       ! THE CLOUD COVER
       CALL UPDATE_VAR_BRACKET(ICE_FILE,ICE_CLD_P,ICE_CLD_N,WTIME,STATUS,ICE_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE ICE CLOUD COVER BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICE_CLD_N,VNP)
       CALL NC_POINT_VAR(ICE_CLD_P,VPP) 
       CLD = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP



!===================================================
!JQI    CASE(ICING_IS_FVCOMGRID)
    CASE(ICE_IS_FVCOMGRID)
!===================================================
!JQI       FTM => ICING_FILE%FTIME
       FTM => ICE_FILE%FTIME


       ! THE SEA SURFACE AIR TEMP
       CALL UPDATE_VAR_BRACKET(ICE_FILE,ICE_SAT_P,ICE_SAT_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE ICE Surface Air Temp BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICE_SAT_N,VNP)
       CALL NC_POINT_VAR(ICE_SAT_P,VPP)   
       SAT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! SHORT WAVE
! afm 20151001 commented out for solar
! With SOLAR, use SOLAR-derived shortwave
! Without SOLAR, use shortwave from forcing data
!       CALL UPDATE_VAR_BRACKET(ICE_FILE,ICE_SWV_P,ICE_SWV_N,WTIME,STATUS)
       CALL UPDATE_VAR_BRACKET(HEAT_FILE,ICE_SWV_P,ICE_SWV_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPDATE ICE SHORTWAVE BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICE_SWV_N,VNP)
       CALL NC_POINT_VAR(ICE_SWV_P,VPP)   
       SWV = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP


       ! THE SPECIFIC HUMIDITY
       CALL UPDATE_VAR_BRACKET(ICE_FILE,ICE_SPQ_P,ICE_SPQ_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE ICE SPECIFIC HUMIDITY BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICE_SPQ_N,VNP)
       CALL NC_POINT_VAR(ICE_SPQ_P,VPP) 
       SPQ = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP


       ! THE CLOUD COVER
       CALL UPDATE_VAR_BRACKET(ICE_FILE,ICE_CLD_P,ICE_CLD_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE ICE CLOUD COVER BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICE_CLD_N,VNP)
       CALL NC_POINT_VAR(ICE_CLD_P,VPP) 
       CLD = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

!===================================================
    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN ICING_FORCING_TYPE IN UPDATE ICING")
    END SELECT
!===================================================
  END SUBROUTINE UPDATE_ICE

!==============================================================================|
!==============================================================================|
  SUBROUTINE UPDATE_ICING(NOW,SAT,WSPDX,WSPDY)
    IMPLICIT NONE
    TYPE(TIME), INTENT(IN) :: NOW
    TYPE(TIME)             :: WTIME
    REAL(SP), ALLOCATABLE :: SAT(:)
    REAL(SP), ALLOCATABLE :: WSPDX(:)
    REAL(SP), ALLOCATABLE :: WSPDY(:)
    REAL(SP), POINTER :: VNP(:), VPP(:)
    TYPE(NCFTIME), POINTER :: FTM
    INTEGER :: STATUS
    REAL(SP), PARAMETER :: K2C     = 273.15_SP

    IF(.NOT. ALLOCATED(SAT)) CALL FATAL_ERROR &
         &("THE Sea Surface Air Temperature VARIABLE PASSED TO UPDATE IS NOT ALLOCATED")
    IF(.NOT. ALLOCATED(WSPDX) .or. .NOT.ALLOCATED(WSPDY)) CALL FATAL_ERROR &
         &("THE WIND SPEED VARIABLES PASSED TO UPDATE ARE NOT ALLOCATED")

!===================================================
    SELECT CASE(ICING_FORCING_KIND)
!===================================================
    CASE (CNSTNT)
       
       WSPDX(1:MT) = ICING_WSPD
       WSPDY=0.0_SP
       ! WEATHER DATA NEEDS TO HAVE WIND VELOCITY, MUST USE RECORD
       ! VECTOR BUT THE MODEL ONLY NEEDS A MAGNITUDE.

       SAT(1:MT) = ICING_AIR_TEMP

       RETURN

    CASE(STTC)

       CALL FATAL_ERROR("STATIC ICING Not Set Up Yet")

    CASE(TMDPNDNT)

       CALL FATAL_ERROR("TIME DEPENDANT ICING Not Set Up Yet")

    CASE(PRDC)

       ! TO SET ZERO TIME PHASE USING RUNFILE START TIME
!      WTIME= NOW - RUNFILE_StartTime
       
       ! TO USE ZERO AS THE PHASE OF THE FORCING
       WTIME= NOW
       

       WTIME = MOD(WTIME,ICING_PERIOD)

    CASE(VRBL)

       WTIME = NOW
    END SELECT


!===================================================
    SELECT CASE(ICING_FORCING_TYPE)
!===================================================
    CASE(ICING_IS_WRFGRID)
!===================================================

       FTM => ICING_FILE%FTIME

       ! THE X DIRECTION WIND SPEED
       CALL UPDATE_VAR_BRACKET(ICING_FILE,ICING_WSPX_P,ICING_WSPX_N,WTIME,STATUS,ICING_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND SPEED X BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICING_WSPX_N,VNP)
       CALL NC_POINT_VAR(ICING_WSPX_P,VPP)   
!       ALLOCATE(WSPDX(0:MT))
       WSPDX = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! THE Y DIRECTION WIND SPEED
       CALL UPDATE_VAR_BRACKET(ICING_FILE,ICING_WSPY_P,ICING_WSPY_N,WTIME,STATUS,ICING_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND SPEED Y BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICING_WSPY_N,VNP)
       CALL NC_POINT_VAR(ICING_WSPY_P,VPP) 
!       ALLOCATE(WSPDY(0:MT))
       WSPDY = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

!       WSPD = sqrt(wspdy**2 + wspdx**2)
!       wspd(0) = 0.0_sp
!       deallocate(wspdy,wspdx)

       ! THE SEA SURFACE AIR TEMP
       CALL UPDATE_VAR_BRACKET(ICING_FILE,ICING_SAT_P,ICING_SAT_N,WTIME,STATUS,ICING_INTP_N)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE Surface Air Temp BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICING_SAT_N,VNP)
       CALL NC_POINT_VAR(ICING_SAT_P,VPP)   
       SAT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP - K2C

!===================================================
    CASE(ICING_IS_FVCOMGRID)
!===================================================
       FTM => ICING_FILE%FTIME

       ! THE X DIRECTION WIND SPEED
       CALL UPDATE_VAR_BRACKET(ICING_FILE,ICING_WSPX_P,ICING_WSPX_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND SPEED X BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICING_WSPX_N,VNP)
       CALL NC_POINT_VAR(ICING_WSPX_P,VPP)   
!       ALLOCATE(WSPDX(0:MT))
       WSPDX = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

       ! THE Y DIRECTION WIND SPEED
       CALL UPDATE_VAR_BRACKET(ICING_FILE,ICING_WSPY_P,ICING_WSPY_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE WIND SPEED Y BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICING_WSPY_N,VNP)
       CALL NC_POINT_VAR(ICING_WSPY_P,VPP) 
!       ALLOCATE(WSPDY(0:MT))
       WSPDY = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP

!       WSPD = sqrt(wspdy**2 + wspdx**2)
!       wspd(0) = 0.0_sp
!       deallocate(wspdy,wspdx)

       ! THE SEA SURFACE AIR TEMP
       CALL UPDATE_VAR_BRACKET(ICING_FILE,ICING_SAT_P,ICING_SAT_N,WTIME,STATUS)
       IF (STATUS /= 0) THEN
          CALL FATAL_ERROR("COULD NOT UPATE Surface Air Temp BRACKET: BOUNDS EXCEEDED?")
       end if

       CALL NC_POINT_VAR(ICING_SAT_N,VNP)
       CALL NC_POINT_VAR(ICING_SAT_P,VPP)   
       SAT = FTM%NEXT_WGHT * VNP + FTM%PREV_WGHT * VPP -K2C

!===================================================
    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN ICING_FORCING_TYPE IN UPDATE ICING")
    END SELECT
!===================================================


  END SUBROUTINE UPDATE_ICING

!=====================================================================
!
!=====================================================================


!==============================================================================|
  !========================================================================

  SUBROUTINE GDAY2(IDD,IMM,IYY,ICC,KD)
!
!  given day,month,year and century(each 2 digits), gday returns
!  the day#, kd based on the gregorian calendar.
!  the gregorian calendar, currently 'universally' in use was
!  initiated in europe in the sixteenth century. note that gday
!  is valid only for gregorian calendar dates.
!
!  kd=1 corresponds to january 1, 0000
!
!  note that the gregorian reform of the julian calendar 
!  omitted 10 days in 1582 in order to restore the date
!  of the vernal equinox to march 21 (the day after
!  oct 4, 1582 became oct 15, 1582), and revised the leap 
!  year rule so that centurial years not divisible by 400
!  were not leap years.
!
!  this routine was written by eugene neufeld, at ios, in june 1990.
!
    integer idd, imm, iyy, icc, kd
    integer ndp(13)
    integer ndm(12)
    data ndp/0,31,59,90,120,151,181,212,243,273,304,334,365/
    data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
!
!  test for invalid input:
    if(icc.lt.0)then
!     write(11,5000)icc
      call pstop
    endif
    if(iyy.lt.0.or.iyy.gt.99)then
!     write(11,5010)iyy
      call pstop
    endif
    if(imm.le.0.or.imm.gt.12)then
!     write(11,5020)imm
      call pstop
      endif
    if(idd.le.0)then
!     write(11,5030)idd
      call pstop
    endif
    if(imm.ne.2.and.idd.gt.ndm(imm))then
!     write(11,5030)idd
      call pstop
    endif
    if(imm.eq.2.and.idd.gt.29)then
!     write(11,5030)idd
      call pstop
    endif
    if(imm.eq.2.and.idd.gt.28.and.((iyy/4)*4-iyy.ne.0.or.(iyy.eq.0.and.(icc/4)*4-icc.ne.0)))then
!     write(11,5030)idd
      call pstop
    endif
5000  format(' input error. icc = ',i7)
5010  format(' input error. iyy = ',i7)
5020  format(' input error. imm = ',i7)
5030  format(' input error. idd = ',i7)
!
!  calculate day# of last day of last century:
    kd = icc*36524 + (icc+3)/4
!
!  calculate day# of last day of last year:
    kd = kd + iyy*365 + (iyy+3)/4
!
!  adjust for century rule:
!  (viz. no leap-years on centurys except when the 2-digit
!  century is divisible by 4.)
    if(iyy.gt.0.and.(icc-(icc/4)*4).ne.0) kd=kd-1
!  kd now truly represents the day# of the last day of last year.
!
!  calculate day# of last day of last month:
    kd = kd + ndp(imm)
!
!  adjust for leap years:
    if(imm.gt.2.and.((iyy/4)*4-iyy).eq.0.and.((iyy.ne.0).or.(((icc/4)*4-icc).eq.0))) kd=kd+1
!  kd now truly represents the day# of the last day of the last
!  month.
!
!  calculate the current day#:
    kd = kd + idd

  RETURN
  END SUBROUTINE GDAY2 




END MODULE MOD_FORCE
