!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mercury_mod.F90
!
! !DESCRIPTION: Module MERCURY\_MOD contains variables and routines for the
!  GEOS-CHEM mercury simulation.  Many choices of reaction mechanism and
!  model processes can be selected with logical switches located in
!  INIT\_MERCURY.
!\\
!\\
! !INTERFACE:
!
MODULE MERCURY_MOD
!
! !USES:
!
  USE DEPO_MERCURY_MOD,  ONLY : ADD_HG2_SNOWPACK
  USE DEPO_MERCURY_MOD,  ONLY : LHGSNOW
  USE OCEAN_MERCURY_MOD, ONLY : LDYNSEASALT
  USE OCEAN_MERCURY_MOD, ONLY : LPOLARBR
  USE OCEAN_MERCURY_MOD, ONLY : LVEGEMIS
  USE OCEAN_MERCURY_MOD, ONLY : STRAT_BR_FACTOR
  USE OCEAN_MERCURY_MOD, ONLY : LAnthroHgOnly
  USE OCEAN_MERCURY_MOD, ONLY : LnoUSAemis
  USE OCEAN_MERCURY_MOD, ONLY : LOCEANCOEF
  USE PhysConstants           ! Physical constants
  USE PRECISION_MOD           ! For GEOS-Chem Precision (fp)


  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CHEMMERCURY
  PUBLIC  :: CLEANUP_MERCURY
  PUBLIC  :: INIT_MERCURY
  PUBLIC  :: EMISSMERCURY
  PUBLIC  :: Reset_Hg_Diags

!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: EMITHG
  PRIVATE :: SRCHg0
  PRIVATE :: SRCHg2
  PRIVATE :: SRCHgP
  PRIVATE :: MERCURY_READYR
  PRIVATE :: CALC_HG2_SEASALT_LOSSRATE
  PRIVATE :: OHNO3TIME
  PRIVATE :: Set_OxidPointers
  PRIVATE :: DiurnalHOx
  PRIVATE :: PolarBrOx
  PRIVATE :: PartNOx
  PRIVATE :: PartXOx
  PUBLIC  :: Set_HgOxidConc
!
! !REMARKS:

!
! !REVISION HISTORY:
!  06 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Parameters
  REAL(fp),  PARAMETER  :: SMALLNUM = 1e-20_fp

  ! Arrays
  INTEGER,  ALLOCATABLE :: AN_Hg0(:,:)    ! Index array for anth Hg0 regions
  INTEGER,  ALLOCATABLE :: AN_Hg2(:,:)    ! Index array for anth Hg2 regions
  INTEGER,  ALLOCATABLE :: AN_HgP(:,:)    ! Index array for anth HgP regions

  ! Arrays for tagged Hg simulation
  REAL(fp), ALLOCATABLE :: EHg0_an(:,:)   ! Anth Hg0 emis [kg/s] - Total
  REAL(fp), ALLOCATABLE :: EHg0_can(:,:)  ! Anth Hg0 emis [kg/s] - Canada
  REAL(fp), ALLOCATABLE :: EHg0_usa(:,:)  ! Anth Hg0 emis [kg/s] - USA
  REAL(fp), ALLOCATABLE :: EHg0_cam(:,:)  ! Anth Hg0 emis [kg/s] - C America
  REAL(fp), ALLOCATABLE :: EHg0_sam(:,:)  ! Anth Hg0 emis [kg/s] - S America
  REAL(fp), ALLOCATABLE :: EHg0_waf(:,:)  ! Anth Hg0 emis [kg/s] - W Africa
  REAL(fp), ALLOCATABLE :: EHg0_eaf(:,:)  ! Anth Hg0 emis [kg/s] - E Africa
  REAL(fp), ALLOCATABLE :: EHg0_saf(:,:)  ! Anth Hg0 emis [kg/s] - S Africa
  REAL(fp), ALLOCATABLE :: EHg0_naf(:,:)  ! Anth Hg0 emis [kg/s] - N Africa
  REAL(fp), ALLOCATABLE :: EHg0_eur(:,:)  ! Anth Hg0 emis [kg/s] - Europe
  REAL(fp), ALLOCATABLE :: EHg0_eeu(:,:)  ! Anth Hg0 emis [kg/s] - E Europe
  REAL(fp), ALLOCATABLE :: EHg0_mde(:,:)  ! Anth Hg0 emis [kg/s] - Mid. East
  REAL(fp), ALLOCATABLE :: EHg0_sov(:,:)  ! Anth Hg0 emis [kg/s] - Fmr USSR
  REAL(fp), ALLOCATABLE :: EHg0_sas(:,:)  ! Anth Hg0 emis [kg/s] - S Asia
  REAL(fp), ALLOCATABLE :: EHg0_eas(:,:)  ! Anth Hg0 emis [kg/s] - E Asia
  REAL(fp), ALLOCATABLE :: EHg0_sea(:,:)  ! Anth Hg0 emis [kg/s] - SE Asia
  REAL(fp), ALLOCATABLE :: EHg0_jpn(:,:)  ! Anth Hg0 emis [kg/s] - Japan
  REAL(fp), ALLOCATABLE :: EHg0_oce(:,:)  ! Anth Hg0 emis [kg/s] - Oceania
  REAL(fp), ALLOCATABLE :: EHg2_an(:,:)   ! Anth Hg2 emis [kg/s] - Total
  REAL(fp), ALLOCATABLE :: EHg2_can(:,:)  ! Anth Hg2 emis [kg/s] - Canada
  REAL(fp), ALLOCATABLE :: EHg2_usa(:,:)  ! Anth Hg2 emis [kg/s] - USA
  REAL(fp), ALLOCATABLE :: EHg2_cam(:,:)  ! Anth Hg2 emis [kg/s] - C America
  REAL(fp), ALLOCATABLE :: EHg2_sam(:,:)  ! Anth Hg2 emis [kg/s] - S America
  REAL(fp), ALLOCATABLE :: EHg2_waf(:,:)  ! Anth Hg2 emis [kg/s] - W Africa
  REAL(fp), ALLOCATABLE :: EHg2_eaf(:,:)  ! Anth Hg2 emis [kg/s] - E Africa
  REAL(fp), ALLOCATABLE :: EHg2_saf(:,:)  ! Anth Hg2 emis [kg/s] - S Africa
  REAL(fp), ALLOCATABLE :: EHg2_naf(:,:)  ! Anth Hg2 emis [kg/s] - N Africa
  REAL(fp), ALLOCATABLE :: EHg2_eur(:,:)  ! Anth Hg2 emis [kg/s] - Europe
  REAL(fp), ALLOCATABLE :: EHg2_eeu(:,:)  ! Anth Hg2 emis [kg/s] - E Europe
  REAL(fp), ALLOCATABLE :: EHg2_mde(:,:)  ! Anth Hg2 emis [kg/s] - Mid. East
  REAL(fp), ALLOCATABLE :: EHg2_sov(:,:)  ! Anth Hg2 emis [kg/s] - Fmr USSR
  REAL(fp), ALLOCATABLE :: EHg2_sas(:,:)  ! Anth Hg2 emis [kg/s] - S Asia
  REAL(fp), ALLOCATABLE :: EHg2_eas(:,:)  ! Anth Hg2 emis [kg/s] - E Asia
  REAL(fp), ALLOCATABLE :: EHg2_sea(:,:)  ! Anth Hg2 emis [kg/s] - SE Asia
  REAL(fp), ALLOCATABLE :: EHg2_jpn(:,:)  ! Anth Hg2 emis [kg/s] - Japan
  REAL(fp), ALLOCATABLE :: EHg2_oce(:,:)  ! Anth Hg2 emis [kg/s] - Oceania
  REAL(fp), ALLOCATABLE :: EHgP_an(:,:)   ! Anth HgP emis [kg/s] - Total
  REAL(fp), ALLOCATABLE :: EHgP_can(:,:)  ! Anth HgP emis [kg/s] - Canada
  REAL(fp), ALLOCATABLE :: EHgP_usa(:,:)  ! Anth HgP emis [kg/s] - USA
  REAL(fp), ALLOCATABLE :: EHgP_cam(:,:)  ! Anth HgP emis [kg/s] - C America
  REAL(fp), ALLOCATABLE :: EHgP_sam(:,:)  ! Anth HgP emis [kg/s] - S America
  REAL(fp), ALLOCATABLE :: EHgP_waf(:,:)  ! Anth HgP emis [kg/s] - W Africa
  REAL(fp), ALLOCATABLE :: EHgP_eaf(:,:)  ! Anth HgP emis [kg/s] - E Africa
  REAL(fp), ALLOCATABLE :: EHgP_saf(:,:)  ! Anth HgP emis [kg/s] - S Africa
  REAL(fp), ALLOCATABLE :: EHgP_naf(:,:)  ! Anth HgP emis [kg/s] - N Africa
  REAL(fp), ALLOCATABLE :: EHgP_eur(:,:)  ! Anth HgP emis [kg/s] - Europe
  REAL(fp), ALLOCATABLE :: EHgP_eeu(:,:)  ! Anth HgP emis [kg/s] - E Europe
  REAL(fp), ALLOCATABLE :: EHgP_mde(:,:)  ! Anth HgP emis [kg/s] - Mid. East
  REAL(fp), ALLOCATABLE :: EHgP_sov(:,:)  ! Anth HgP emis [kg/s] - Fmr USSR
  REAL(fp), ALLOCATABLE :: EHgP_sas(:,:)  ! Anth HgP emis [kg/s] - S Asia
  REAL(fp), ALLOCATABLE :: EHgP_eas(:,:)  ! Anth HgP emis [kg/s] - E Asia
  REAL(fp), ALLOCATABLE :: EHgP_sea(:,:)  ! Anth HgP emis [kg/s] - SE Asia
  REAL(fp), ALLOCATABLE :: EHgP_jpn(:,:)  ! Anth HgP emis [kg/s] - Japan
  REAL(fp), ALLOCATABLE :: EHgP_oce(:,:)  ! Anth HgP emis [kg/s] - Oceania

  ! Emissions from various sources
  REAL(fp), ALLOCATABLE :: EHg0_am(:,:)   ! Artisinal mining Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_oc(:,:,:) ! Ocean Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_ln(:,:,:) ! Hg reemission from land [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_dist(:,:) ! Spatial dist of terrestrial Hg0
                                          !  sources [unitless]
  REAL(fp), ALLOCATABLE :: EHg0_geo(:,:)  ! Geogenic Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_bb(:,:)   ! Biomass burning Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_vg(:,:)   ! Vegetation Hg0 emis [kg/s
  REAL(fp), ALLOCATABLE :: EHg0_so(:,:)   ! Soil Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_gtm(:,:)  ! GTMM Hg0 emis [kg/s]
  REAL(fp), ALLOCATABLE :: EHg0_snow(:,:,:) !Snow Hg0 emis [kg/s]

  ! Arrays for diurnal variation of HOx
  REAL(fp), ALLOCATABLE :: COSZM(:,:)     ! Max daily solar zenith angle
  REAL(fp), ALLOCATABLE :: TCOSZ(:,:)     ! Sum of solar zenith angle
  REAL(fp), ALLOCATABLE :: TTDAY(:,:)     ! Total daylight time at I,J [min]

  ! Deposition arrays
  REAL(fp), ALLOCATABLE :: ZERO_DVEL(:,:) ! Zero drydep velocity [cm/s]
  REAL(fp), ALLOCATABLE :: HG2_SEASALT_LOSSRATE(:,:)

  ! Arrays for photolysis rates
  REAL(fp), ALLOCATABLE :: JNO2_INST(:,:,:) ! Instantaneous JNO2
  REAL(fp), ALLOCATABLE :: JBrO_INST(:,:,:) ! Instantaneous JBrO
  REAL(fp), ALLOCATABLE :: JClO_INST(:,:,:) ! Instantaneous JClO

  ! Arrays for het rates
  REAL(fp), ALLOCATABLE :: AerHetRates(:,:,:,:)
  REAL(fp), ALLOCATABLE :: OrgHetRates(:,:,:,:)

  ! For now, we need an emission array for the HG simulation
  ! that can be passed to vdiff_mod.F90 since Trac_Tend does
  ! not exist anymore (ckeller, 10/21/2014).
  REAL(fp), ALLOCATABLE, PUBLIC :: HG_EMIS(:,:,:)

  ! Offline AOD arrays that will passed to Fast-JX
  REAL(fp), ALLOCATABLE, PUBLIC :: HEM_ODAer (:,:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: HEM_ODDust(:,:,:,:)

  ! Pointers to fields in the HEMCO data structure.
  ! These need to be declared REAL(f4), aka REAL*4.
  ! (NOTE: We can set them to NULL here because
  ! they are globally SAVEd variables (bmy, 4/29/16)

  REAL(f4), POINTER :: OCEAN_CONC(:,:,:)    => NULL()



  ! Hg species IDs
  INTEGER           :: N_Hg_CATS
  INTEGER           :: id_Hg0,     id_Hg2,     id_HgP
  INTEGER           :: id_O3,      id_OH,      id_HO2
  INTEGER           :: id_ClO,     id_Cl,      id_HOCl
  INTEGER           :: id_NO2,     id_Br,      id_BrO
  INTEGER           :: id_NO,      id_CO,      id_CH4
  INTEGER           :: id_HGBRNO2, id_HGBRHO2, id_HGBROH
  INTEGER           :: id_HGBRBRO, id_HGBRCLO, id_HGBR2
  INTEGER           :: id_HGCLNO2, id_HGCLHO2, id_HGCLOH
  INTEGER           :: id_HGCLBRO, id_HGCLCLO, id_HGCLBR
  INTEGER           :: id_HGCL2
  INTEGER           :: id_HGIIADS
  INTEGER           :: nHg2gasSpc

  INTEGER           :: Map_Hg2gas(25)

  INTEGER           :: ID_Hg_tot,  ID_Hg_can,  ID_Hg_usa
  INTEGER           :: ID_Hg_cam,  ID_Hg_sam,  ID_Hg_waf
  INTEGER           :: ID_Hg_eaf,  ID_Hg_saf,  ID_Hg_naf
  INTEGER           :: ID_Hg_eur,  ID_Hg_eeu,  ID_Hg_sov
  INTEGER           :: ID_Hg_mde,  ID_Hg_sas,  ID_Hg_eas
  INTEGER           :: ID_Hg_sea,  ID_Hg_jpn,  ID_Hg_oce
  INTEGER           :: ID_Hg_so,   ID_Hg_bb,   ID_Hg_geo
  INTEGER           :: ID_Hg_atl,  ID_Hg_nat,  ID_Hg_sat
  INTEGER           :: ID_Hg_npa,  ID_Hg_arc,  ID_Hg_ant
  INTEGER           :: ID_Hg_ocn,  ID_Hg_str

  ! Pointers for Hg indexing
  ! (NOTE: We can set them to NULL here because
  ! they are globally SAVEd variables (bmy, 4/29/16)
  INTEGER, POINTER  :: Hg0_Id_List(:) => NULL()
  INTEGER, POINTER  :: Hg2_Id_List(:) => NULL()
  INTEGER, POINTER  :: HgP_Id_List(:) => NULL()

  ! Species index of Hg oxidants in GEOS-Chem species database
  INTEGER               :: N_Hg_Oxid
  INTEGER, ALLOCATABLE  :: GC_HgOxid_TrID(:)

  ! Species index of aerosol species read from HEMCO
  INTEGER                        :: N_Aer, N_Dust
  CHARACTER(LEN=8), ALLOCATABLE  :: AerSpcNames(:)


  ! KPP Arrays
  INTEGER,  ALLOCATABLE :: PL_Kpp_ID (:)

  ! OxidPtr is a derived type to hold pointers to the oxidant fields
  TYPE :: HCOPtrObj
     REAL(f4), POINTER  :: Data(:,:,:) => NULL()
  END TYPE HCOPtrObj

  TYPE :: AeroPtrObj
     REAL(f4), POINTER  :: AOD(:,:,:)  => NULL()
     REAL(f4), POINTER  :: Area(:,:,:) => NULL()
     REAL(f4), POINTER  :: Radi(:,:,:) => NULL()
     REAL(f4), POINTER  :: Conc(:,:,:) => NULL()
  END TYPE AeroPtrObj

  ! Vectors holding the oxidant concentrations,
  ! which will be read and interpolated by HEMCO. The
  ! corresponding HEMCO fields must be specified in the HEMCO
  ! configuration file. The field names are assumed to be
  ! 'GLOBAL_XY', where XY is the species name.
  ! It is also assumed that the input data is in molec cm-3.
  TYPE(HCOPtrObj), POINTER :: HgOxidPtr(:)

  ! Vectors holding the AOD fields.
  TYPE(AeroPtrObj), POINTER :: AeroPtr(:)


CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemmercury
!
! !DESCRIPTION: Subroutine CHEMMERCURY is the driver routine for mercury
!  chemistry in the GEOS-CHEM module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEMMERCURY( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
!
! !USES:
!
    USE DEPO_MERCURY_MOD,   ONLY : ADD_HG2_DD
    USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_DD
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD
    USE DIAG03_MOD,         ONLY : AD03_Hg2_Hg0, AD03_Hg2_O3
    USE DIAG03_MOD,         ONLY : AD03_Hg2_Br,  AD03_Hg2_OH
    USE DIAG03_MOD,         ONLY : AD03_Hg2_BRY,  AD03_Hg2_CLY
    USE DIAG03_MOD,         ONLY : AD03_Hg2_Br_Y
    USE DIAG03_MOD,         ONLY : AD03_Hg2_SS,  LD03
    USE DIAG03_MOD,         ONLY : AD03_Hg2_SSR, ND03
    USE DIAG03_MOD,         ONLY : AD03_Br
#endif
    USE FAST_JX_MOD,          ONLY : FAST_JX
    USE CMN_FJX_MOD
    USE GCKPP_Monitor,        ONLY : SPC_NAMES, FAM_NAMES
    USE GCKPP_Parameters
    USE GCKPP_Integrator,     ONLY : INTEGRATE, NHnew
    USE GCKPP_Function
    USE GCKPP_Model
    USE GCKPP_Global
    USE GCKPP_Rates,          ONLY : UPDATE_RCONST, RCONST
    USE GCKPP_Initialize,     ONLY : Init_KPP => Initialize
    USE Timers_Mod
    USE PhysConstants,        ONLY : AVO
    USE PRESSURE_MOD
    USE Species_Mod,          ONLY : Species
    USE TIME_MOD,             ONLY : GET_TS_CHEM
    USE TIME_MOD,             ONLY : Get_Day
    USE TIME_MOD,             ONLY : Get_Month
    USE TIME_MOD,             ONLY : Get_Year
    USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
    USE TIME_MOD,           ONLY : ITS_TIME_FOR_A3
    USE UnitConv_Mod,         ONLY : Convert_Spc_Units
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE ERROR_MOD,          ONLY : SAFE_DIV
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:

!
! !REVISION HISTORY:
!  01 Oct 1995 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! For testing purposes
    LOGICAL            :: DO_HETCHEM
    LOGICAL            :: DO_PHOTCHEM

    ! Scalars
    LOGICAL            :: LDYNOCEAN
    LOGICAL            :: LGTMM
    LOGICAL            :: LNLPBL
    LOGICAL            :: prtDebug
    INTEGER            :: nAdvect, nSpecies
    INTEGER            :: I, J, L, K, N, NN, CN, Hg_Cat
    INTEGER                :: NA,        F,        SpcID,     KppID
    INTEGER                :: P,         MONTH,    YEAR,      WAVELENGTH
    INTEGER                :: TotSteps,  TotFuncs, TotJacob,  TotAccep
    INTEGER                :: TotRejec,  TotNumLU, HCRC,      IERR
    INTEGER                :: Day

    REAL(fp)           :: K_DRYD0, K_DRYD2G, K_DRYD2P, K_SALT
    REAL(fp)           :: K_OX, K_RED, K_DEP0, K_DEP2G, K_DEP2P
    REAL(fp)           :: DEP_Hg0, DEP_HG2G, DEP_HG2P, K_RED2P
    REAL(fp)           :: DEP_HG2G_DRY, DEP_HG2G_SALT, DEP_HG2P_DRY
    REAL(fp)           :: DEP_HG2_DRY
    REAL(fp)           :: K0,  Hg2CR,  R_DRY
    REAL(fp)           :: AREA_CM2, DEP_DRY_FLX
    REAL(fp)           :: F_HG0_DEP
    REAL(fp)               :: Start,     Finish,   rtim,      itim
    REAL(fp)               :: YLAT,      T,         TIN
    REAL(fp)               :: TOUT

    ! Strings
    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg,   ThisLoc
    CHARACTER(LEN=16)      :: ThisName

    ! Arrays
    INTEGER                :: ICNTRL     (                  20               )
    INTEGER                :: ISTATUS    (                  20               )
    REAL(dp)               :: RCNTRL     (                  20               )
    REAL(dp)               :: RSTATE     (                  20               )

    REAL(dp)               :: Vloc(NVAR), Aout(NREACT)

    ! Pointers
    REAL(fp), POINTER  :: Spc(:,:,:,:)
    REAL(fp), POINTER  :: TK(:,:,:   )

    ! Objects
    TYPE(Species), POINTER :: SpcInfo
!
! !DEFINED PARAMETERS:
!

    !================================================================
    ! CHEMMERCURY begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at ChemMercury (in GeosCore/mercury_mod.F90)'

    ! Copy values from Input_Opt
    LDYNOCEAN = Input_Opt%LDYNOCEAN
    LGTMM     = Input_Opt%LGTMM
    LNLPBL    = Input_Opt%LNLPBL
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Values from State_Chem
    nSpecies  = State_Chm%nSpecies

    ! Initialize KPP variables
    itim      =  0.0_fp
    rtim      =  0.0_fp
    totsteps  =  0
    totfuncs  =  0
    totjacob  =  0
    totaccep  =  0
    totrejec  =  0
    totnumLU  =  0
    Day       =  Get_Day()    ! Current day
    Month     =  Get_Month()  ! Current month
    Year      =  Get_Year()   ! Current year


    ! Initialize pointers
    Spc      => State_Chm%Species   ! Chemical species array [kg]
    TK       => State_Met%T         ! Temperature [K]
    SpcInfo  => NULL()

    ! Get info about this species from the database
    SpcInfo => State_Chm%SpcData(id_Hg2)%Info

    ! Get Henry's Law constant for Hg2 [mol/L/atm] (ewl, 1/5/16)
    K0    = SpcInfo%Henry_K0
    Hg2CR = SpcInfo%Henry_CR

    ! Calculate dry air gas constant in [L atm/K/mol] (ewl, 1/5/16)
    R_DRY = Rd * AIRMW / ATM


    !================================================================
    ! Set chemistry options and pointers to chemical inputs from HEMCO
    !=================================================================

    ! Turn heterogeneous chemistry and photolysis on/off for testing
    DO_HETCHEM  = .TRUE.
    DO_PHOTCHEM = .TRUE.
    IF ( FIRST ) THEN
       IF ( .not. DO_HETCHEM ) THEN
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
          WRITE( 6, '(a)' )  ' # Do_FlexChem: Heterogeneous chemistry' // &
                             ' is turned off for testing purposes.'
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
       ENDIF
       IF ( .not. DO_PHOTCHEM ) THEN
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
          WRITE( 6, '(a)' )  ' # Do_FlexChem: Photolysis chemistry' // &
                             ' is turned off for testing purposes.'
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
       ENDIF
    ENDIF

    ! Override the default value of LGCBROMINE and LRED_JNO2 depending
    ! on the settings in the HEMCO configuration file.  NOTE: This cannot
    ! be done in INIT_MERCURY because at that time the HEMCO configuration
    ! file has not been read yet. (bmy, 3/12/15)
    IF ( FIRST ) THEN
!       CALL SET_OPTIONS_FROM_HEMCO( Input_Opt, State_Grid, RC )

        ! Get pointers to oxidant fields via HEMCO
       CALL Set_OxidPointers ( Input_Opt, State_Chm, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "Set_OxidPointers"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
       ENDIF

       FIRST = .FALSE.
    ENDIF

!    State_Chm%Species(:,:,:,id_Hg2) = 0.0
!    DO CN = 2, State_Chm%N_Hg2_CATS
!       N = State_Chm%Hg2_Id_List(CN)
!       State_Chm%Species(:,:,:,id_Hg2) = State_Chm%Species(:,:,:,N) &
!            + State_Chm%Species(:,:,:,id_Hg2)
!    ENDDO

    IF ( ITS_A_NEW_MONTH() ) THEN

        ! Set AOD fields to pass to FastJX
!$OMP PARALLEL DO                                              &
!$OMP DEFAULT( SHARED )                                        &
!$OMP PRIVATE( I,       J,      L,       N      )
        DO L=1, State_Grid%NZ
        DO J=1, State_Grid%NY
        DO I=1, State_Grid%NX

            DO N=1, N_Dust
                HEM_ODDust(I,J,L,N)  = AeroPtr(N)%AOD(I,J,L)
            ENDDO

            DO N=1, N_Aer
                HEM_ODAer(I,J,L,N)  = AeroPtr(N_Dust+N)%AOD(I,J,L)
            ENDDO

        ENDDO
        ENDDO
        ENDDO
!$OMP END PARALLEL DO

    ENDIF

    !======================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'molec/cm3', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'mercury_mod.F90')
       RETURN
    ENDIF

    !=======================================================================
    ! Call photolysis routine to compute J-Values
    !=======================================================================
    IF ( DO_PHOTCHEM ) THEN
        !Do Photolysis
        CALL FAST_JX( WAVELENGTH, Input_Opt,  State_Chm, &
                      State_Diag, State_Grid, State_Met, RC )

        ! Trap potential errors
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "FAST_JX"!'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        !### Debug
        IF ( prtDebug ) THEN
           CALL DEBUG_MSG( '### ChemMercury: after FAST_JX' )
        ENDIF
    ENDIF

    !======================================================================
    ! Save J values for NO2, BrO, and ClO to partition species
    !======================================================================
!$OMP PARALLEL DO                                              &
!$OMP DEFAULT( SHARED )                                        &
!$OMP PRIVATE( I, J, L, N, P, SpcId,  ThisName      )
    DO L=1, State_Grid%NZ
    DO J=1, State_Grid%NY
    DO I=1, State_Grid%NX

        IF ( State_Met%SUNCOSmid(I,J) <= -0.1391731e+0_fp ) CYCLE

        ! Loop over the FAST-JX photolysis species
        DO N = 1, NRATJ

            ! GC photolysis species index
            P = GC_Photo_Id(N)

            ! Proceed only if species is in index
            IF ( P <= 0 ) CYCLE

            ! GEOS-Chem species ID
            SpcID = State_Chm%Map_Photol(P)

            ! Species name
            ThisName = State_Chm%SpcData(SpcID)%Info%Name

            ! Look for the relevant species
            SELECT CASE( TRIM( ThisName ) )
                CASE( 'NO2' )
                    JNO2_INST(I,J,L) = ZPJ(L,N,I,J)
                CASE( 'BrO' )
                    JBrO_INST(I,J,L) = ZPJ(L,N,I,J)
                CASE( 'ClO' )
                    JClO_INST(I,J,L) = ZPJ(L,N,I,J)
                CASE DEFAULT
                  ! Nothing
            END SELECT
        ENDDO

    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    !======================================================================
    ! Set instantaneous oxidant concentrations (molec cm-3)
    !======================================================================
    CALL Set_HgOxidConc( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Set_HgOxidConc"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### ChemMercury: after Set_HgOxidConc' )
    ENDIF

    !======================================================================
    ! Set heterogeneous uptake rates (s-1)
    !======================================================================
    CALL Set_HetRates( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Set_HetRates"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### ChemMercury: after Set_HetRates' )
    ENDIF

    !=======================================================================
    ! Set up integration convergence conditions and timesteps
    ! (cf. M. J. Evans)
    !=======================================================================

    !%%%%% TIMESTEPS %%%%%
    DT        = GET_TS_CHEM() ! [s]
    T         = 0d0
    TIN       = T
    TOUT      = T + DT

    !%%%%% CONVERGENCE CRITERIA %%%%%

    ! Absolute tolerance
    ATOL      = 1e-2_dp

    ! Relative tolerance
    RTOL      = 1e-2_dp

    !%%%%% SOLVER OPTIONS %%%%%

    ! Zero all slots of ICNTRL
    ICNTRL    = 0

    ! 0 - non-autonomous, 1 - autonomous
    ICNTRL(1) = 1

    ! 0 - vector tolerances, 1 - scalars
    ICNTRL(2) = 0

    ! Select Integrator
    ICNTRL(3) = 4 ! Rodas3

    ! 0 - adjoint, 1 - no adjoint
    ICNTRL(7) = 1

    !=======================================================================
    ! %%%%% SOLVE CHEMISTRY -- This is the main KPP solver loop %%%%%
    !=======================================================================
100 format('No. of function calls:', i6, /,                                 &
           'No. of jacobian calls:', i6, /,                                 &
           'No. of steps:         ', i6, /,                                 &
           'No. of accepted steps:', i6, /,                                 &
           'No. of rejected steps ', i6, /,                                 &
           '       (except at very beginning)',          /,                 &
           'No. of LU decompositions:             ', i6, /,                 &
           'No. of forward/backward substitutions:', i6, /,                 &
           'No. of singular matrix decompositions:', i6, /,                 &
            /,                                                              &
           'Texit, the time corresponding to the      ',        /,          &
           '       computed Y upon return:            ', f11.4, /,          &
           'Hexit, last accepted step before exit:    ', f11.4, /,          &
           'Hnew, last predicted step (not yet taken):', f11.4 )

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT  ( SHARED                                                 )&
    !$OMP PRIVATE  ( I,        J,        L,       N,      YLAT              )&
    !$OMP PRIVATE  ( IERR,     RCNTRL,   START,   FINISH, ISTATUS           )&
    !$OMP PRIVATE  ( RSTATE,   SpcID,    KppID,   F,      P                 )&
    !$OMP PRIVATE  ( Vloc,     Aout,    nn                                  )&
    !$OMP REDUCTION( +:ITIM                                                 )&
    !$OMP REDUCTION( +:RTIM                                                 )&
    !$OMP REDUCTION( +:TOTSTEPS                                             )&
    !$OMP REDUCTION( +:TOTFUNCS                                             )&
    !$OMP REDUCTION( +:TOTJACOB                                             )&
    !$OMP REDUCTION( +:TOTACCEP                                             )&
    !$OMP REDUCTION( +:TOTREJEC                                             )&
    !$OMP REDUCTION( +:TOTNUMLU                                             )&
    !$OMP SCHEDULE ( DYNAMIC,  1                                            )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !====================================================================
       ! For safety sake, initialize certain variables for each grid
       ! box (I,J,L), whether or not chemistry will be done there.
       !====================================================================
       HET       = 0.0_dp    ! Het chem array
       IERR      = 0         ! Success or failure flag
       ISTATUS   = 0.0_dp    ! Rosenbrock output
       PHOTOL    = 0.0_dp    ! Photolysis array
       RCNTRL    = 0.0_fp    ! Rosenbrock input
       RSTATE    = 0.0_dp    ! Rosenbrock output
       P         = 0         ! GEOS-Chem photolyis species ID

       ! Grid-box latitude [degrees]
       YLAT      = State_Grid%YMid(I,J)

       ! Temperature [K]
       TEMP      = State_Met%T(I,J,L)

       ! Pressure [hPa]
       PRESS     = GET_PCENTER( I, J, L )

       ! mje Calculate NUMDEN based on ideal gas law (# cm-3)
       NUMDEN    = State_Met%AIRNUMDEN(I,J,L)

       !====================================================================
       ! Get photolysis rates (daytime only)
       !====================================================================
       IF ( State_Met%SUNCOSmid(I,J) > -0.1391731e+0_fp ) THEN

          ! Loop over the FAST-JX photolysis species
          DO N = 1, JVN_

             IF ( DO_PHOTCHEM ) THEN
                ! Copy photolysis rate from FAST_JX into KPP PHOTOL array
                PHOTOL(N) = ZPJ(L,N,I,J)
             ENDIF

             !--------------------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Instantaneous photolysis rates [s-1] (aka J-values)
             ! and noontime photolysis rates [s-1]
             !
             !    NOTE: Attach diagnostics here instead of in module
             !    fast_jx_mod.F90 so that we can get the adjusted photolysis
             !    rates (output from routne PHOTRATE_ADJ above).
             !
             ! The mapping between the GEOS-Chem photolysis species and
             ! the FAST-JX photolysis species is contained in the lookup
             ! table in input file FJX_j2j.dat.
             !
             !    NOTE: Depending on the simulation, some GEOS-Chem
             !    species might not map to a of the FAST-JX species
             !    (e.g. SOA species will not be present in a tropchem run).
             !
             ! Some GEOS-Chem photolysis species may have multiple
             ! branches for photolysis reactions.  These will be
             ! represented by multiple entries in the FJX_j2j.dat
             ! lookup table.
             !
             ! To match the legacy bpch diagnostic, we archive the sum of
             ! photolysis rates for a given GEOS-Chem species over all of
             ! the reaction branches.
             !--------------------------------------------------------------

             ! GC photolysis species index
             P = GC_Photo_Id(N)

             ! If this FAST_JX photolysis species maps to a valid
             ! GEOS-Chem photolysis species (for this simulation)...
             IF ( P > 0 ) THEN

!                  IF ( I == 1 .AND. J == 25 .AND. L ==10 .AND. P == 1) THEN
!                    WRITE(*,'(A10,1X,ES10.2)') 'J_HgBrNO2',PHOTOL(N)
!                  ENDIF

                ! Archive the instantaneous photolysis rate
                ! (summing over all reaction branches)
                IF ( State_Diag%Archive_JVal ) THEN
                   State_Diag%JVal(I,J,L,P) = State_Diag%JVal(I,J,L,P)       &
                                            + PHOTOL(N)
                ENDIF

             ENDIF
          ENDDO
       ENDIF

       !====================================================================
       ! Test if we need to do the chemistry for box (I,J,L),
       ! otherwise move onto the next box.
       !====================================================================

       ! If we are not in the troposphere don't do the chemistry!
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       ! Skipping buffer zone (lzh, 08/10/2014)
       IF ( State_Grid%NestedGrid ) THEN
          IF ( J <=                 State_Grid%SouthBuffer ) CYCLE
          IF ( J >  State_Grid%NY - State_Grid%NorthBuffer ) CYCLE
          IF ( I <=                 State_Grid%EastBuffer  ) CYCLE
          IF ( I >  State_Grid%NX - State_Grid%WestBuffer  ) CYCLE
       ENDIF

       !====================================================================
       ! Intialize KPP solver arrays: CFACTOR, VAR, FIX, etc.
       !====================================================================
       CALL Init_KPP( )

       !====================================================================
       ! Get rates for heterogeneous chemistry
       !====================================================================
        IF ( DO_HETCHEM ) THEN
           DO N = 1, nSpecies
               HET(N,1) = AerHetRates( I, J, L, N )
               HET(N,2) = OrgHetRates( I, J, L, N )
           ENDDO
        ENDIF

       !====================================================================
       ! Initialize species concentrations
       !====================================================================

       ! Loop over KPP Species
       DO N = 1, NSPEC

          ! GEOS-Chem species ID
          SpcID = State_Chm%Map_KppSpc(N)

          ! Species info
          SpcInfo => State_Chm%SpcData(SpcID)%Info

          ! Initialize KPP species concentration array
          IF ( SpcID .eq. 0) THEN
             C(N) = 0.0_dp
          ELSE
             C(N) = State_Chm%Species(I,J,L,SpcID)
          ENDIF
          IF ( I == 1 .AND. J == 25 .AND. L ==20 ) THEN
            WRITE(*,'(A10,1X,A10,1X,ES10.2)') 'Before KPP', SpcInfo%Name, C(N)
          ENDIF
       ENDDO

       ! Zero out dummy species index in KPP
       DO F = 1, NFAM
         KppID = PL_Kpp_Id(F)
         IF ( KppID > 0 ) C(KppID) = 0.0_dp
       ENDDO

       !==================================================================
       ! Update KPP rates
       !==================================================================

       ! VAR and FIX are chunks of array C (mps, 2/24/16)
       VAR(1:NVAR) = C(1:NVAR)
       FIX         = C(NVAR+1:NSPEC)

       ! Update the array of rate constants
       CALL Update_RCONST( )

       !=================================================================
       ! Set options for the KPP Integrator (M. J. Evans)
       !=================================================================

       ! Zero all slots of RCNTRL
       RCNTRL    = 0.0_fp

       ! Starting value for integration time step
       RCNTRL(3) = State_Chm%KPPHvalue(I,J,L)

       !=================================================================
       ! Integrate the box forwards
       !=================================================================

       ! Call the KPP integrator
       CALL Integrate( TIN,    TOUT,    ICNTRL,      &
                       RCNTRL, ISTATUS, RSTATE, IERR )

       ! Print grid box indices to screen if integrate failed
       IF ( IERR < 0 ) THEN
          WRITE(6,*) '### INTEGRATE RETURNED ERROR AT: ', I, J, L
       ENDIF

       !------------------------------------------------------------------
       ! Try another time if it failed
       !------------------------------------------------------------------
       IF ( IERR < 0 ) THEN

          ! Reset first time step and start concentrations
          ! Retry the integration with non-optimized
          ! settings
          RCNTRL(3)  = 0e+0_fp
          CALL Init_KPP( )
          VAR = C(1:NVAR)
          FIX = C(NVAR+1:NSPEC)
          CALL Update_RCONST( )
          CALL Integrate( TIN,    TOUT,    ICNTRL,                           &
                          RCNTRL, ISTATUS, RSTATE, IERR                     )

          !------------------------------------------------------------------
          ! Exit upon the second failure
          !------------------------------------------------------------------
          IF ( IERR < 0 ) THEN
             WRITE(6,*) '## INTEGRATE FAILED TWICE !!! '
             WRITE(ERRMSG,'(a,i3)') 'Integrator error code :',IERR
             CALL ERROR_STOP(ERRMSG, 'INTEGRATE_KPP')
          ENDIF

       ENDIF

       !--------------------------------------------------------------------
       ! Continue upon successful return...
       !--------------------------------------------------------------------

       ! Copy VAR and FIX back into C (mps, 2/24/16)
       C(1:NVAR)       = VAR(:)
       C(NVAR+1:NSPEC) = FIX(:)

       ! Save for next integration time step
       State_Chm%KPPHvalue(I,J,L) = RSTATE(Nhnew)

       !====================================================================
       ! Check we have no negative values and copy the concentrations
       ! calculated from the C array back into State_Chm%Species
       !====================================================================
       ! Loop over KPP species
       DO N = 1, NSPEC

          ! GEOS-Chem species ID
          SpcID = State_Chm%Map_KppSpc(N)

          ! Species info
          SpcInfo => State_Chm%SpcData(SpcID)%Info

          ! Skip if this is not a GEOS-Chem species
          IF ( SpcID .eq. 0 ) CYCLE

          ! Set negative concentrations to zero
          C(N) = MAX( C(N), 0.0E0_dp )

          ! Copy concentrations back into State_Chm%Species
          State_Chm%Species(I,J,L,SpcID) = REAL( C(N), kind=fp )

          IF ( I == 1 .AND. J == 25 .AND. L ==20 ) THEN
            WRITE(*,'(A10,1X,A10,1X,ES10.2)') 'After KPP', SpcInfo%Name, C(N)
          ENDIF
       ENDDO

       !====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Prod and loss of families or species [molec/cm3/s]
       !
       ! NOTE: KppId is the KPP ID # for each of the prod and loss
       ! diagnostic species.  This is the value used to index the
       ! KPP "VAR" array (in module gckpp_Global.F90).
       !====================================================================

       ! Chemical loss of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Loss ) THEN
          DO F = 1, State_Chm%nLoss
             KppID                    = State_Chm%Map_Loss(F)
             State_Diag%Loss(I,J,L,F) = VAR(KppID) / DT
          ENDDO
       ENDIF

       ! Chemical production of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Prod ) THEN
          DO F = 1, State_Chm%nProd
             KppID                    = State_Chm%Map_Prod(F)
             State_Diag%Prod(I,J,L,F) = VAR(KppID) / DT
          ENDDO
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=======================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm,  State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90' )
       RETURN
    ENDIF

    !=======================================================================
    ! Partition HgII between gas and aerosol phase
    !=======================================================================
    CALL PARTITIONHG2( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ! Do deposition

    ! Record diagnostics


    ! Free pointer memory
    Spc     => NULL()
    TK      => NULL()
    SpcInfo => NULL()

  END SUBROUTINE CHEMMERCURY

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissmercury
!
! !DESCRIPTION: Subroutine EMISSMERCURY is the driver routine for mercury
!  emissions.
!\\
!\\
! NOTE/TODO: The mercury simulation is the only GEOS-Chem emission code that
! is not yet fully compatible with HEMCO. So far, only the anthropogenic
! emissions are included in HEMCO. For all other sources, the original
! mercury code is used.
!\\
!\\
! For the non-local PBL mixing, all emissions are written into module array
! HG\_EMIS (in kg m-2 s-1). These values are then used by vdiff\_mod.F90.
! This is just a workaround to ensure backwards compatibility of the mercury
! code. Once we have added all mercury emissions to HEMCO, HG\_EMIS is not
! used any more (ckeller, 10/21/2014).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSMERCURY( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
    USE DEPO_MERCURY_MOD,   ONLY : RESET_HG_DEP_ARRAYS
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE LAND_MERCURY_MOD,   ONLY : LAND_MERCURY_FLUX, VEGEMIS
    USE LAND_MERCURY_MOD,   ONLY : SOILEMIS, BIOMASSHG
    USE LAND_MERCURY_MOD,   ONLY : SNOWPACK_MERCURY_FLUX
    USE OCEAN_MERCURY_MOD,  ONLY : OCEAN_MERCURY_FLUX
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH, ITS_A_NEW_MONTH
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
    
    ! Added for GTMM (ccc, 11/19/09)
    !USE LAND_MERCURY_MOD,   ONLY : GTMM_DR
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  03 Jun 2013 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE     :: FIRST = .TRUE.
    INTEGER           :: THISMONTH, I, J
    CHARACTER(LEN=63) :: OrigUnit

    ! For fields from Input_Opt
    LOGICAL           :: LDYNOCEAN
    LOGICAL           :: LGTMM
    LOGICAL           :: LNLPBL
    LOGICAL           :: LPREINDHG
    LOGICAL           :: LEMIS
    LOGICAL           :: prtDebug

    ! Strings
    CHARACTER(LEN=255)   :: ErrMsg
    CHARACTER(LEN=255)   :: ThisLoc

    !=================================================================
    ! EMISSMERCURY begins here!
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at EMISSMERCURY (in module GeosCore/mercury_mod.F90)'

    ! Copy fields from Input_Opt
    LDYNOCEAN = Input_Opt%LDYNOCEAN
    LGTMM     = Input_Opt%LGTMM
    LNLPBL    = Input_Opt%LNLPBL
    LPREINDHG = Input_Opt%LPREINDHG
    LEMIS     = Input_Opt%LEMIS
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Convert species units to [kg] for EMISSMERCURY (ewl, 8/12/15)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg', RC, OrigUnit=OrigUnit )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Spc_Units" #1!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! First-time initialization
    IF ( FIRST ) THEN

       ! Check that emissions are turned on. Print error message
       ! and stop GEOS-Chem if emissions are turned off (ewl, 9/1/15)
       IF ( .not. LEMIS ) THEN
          ErrMsg = 'ERROR: Hg emissions are need for simulation ' // &
                   'but emissions are turned off in input.geos'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Read anthro, ocean, land emissions of Hg from disk
       CALL MERCURY_READYR( Input_Opt, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered call to "MERCURY_READYR"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !=================================================================
    ! Call emission routines for Hg(0), Hg(II), and Hg(P)
    !=================================================================

    ! Ocean flux of Hg(0)
    IF ( LDYNOCEAN ) THEN

       ! Set to zero to clear emissions from previous time step
       ! (cdh, 4/30/09)
       EHg0_oc = 0e+0_fp

       CALL OCEAN_MERCURY_FLUX( Input_Opt,  State_Chm, State_Diag, &
                                State_Grid, State_Met, EHg0_oc,    RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered call to "OCEAN_MERCURY_FLUX"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a OCEAN_FLUX' )

    ELSE
       EHg0_oc = 0e+0_fp
       CALL OFFLINEOCEAN_READMO( State_Chm, State_Diag, State_Grid, &
                                 State_Met, EHg0_oc, RC )
    ENDIF

    !==========================================================================
    ! Disable GTMM until it is brought up-to-date (mps, 3/10/19)
    !IF ( LGTMM ) THEN
    !   !--------------------------------------------------------------
    !   ! Here we are using the Global Terrstrial Mercury Model...
    !   !--------------------------------------------------------------
    !   IF ( ITS_A_NEW_MONTH() ) THEN
    !
    !      ! OLD CODE: CALL GTMM_DR( EHg0_gtm(:,:), State_Met )
    !      CALL GTMM_DR( Input_Opt, State_Chm, State_Grid, State_Met, &
    !                    EHg0_gtm,  RC )
    !
    !      IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a GTMM' )
    !   ENDIF
    !
    !ELSE
    !==========================================================================

    !--------------------------------------------------------------
    ! Here we are NOT using the Global Terrstrial Mercury Model...
    !--------------------------------------------------------------
    CALL LAND_MERCURY_FLUX ( EHg0_ln, LHGSNOW, State_Grid, State_Met )
    IF ( prtDebug )  CALL DEBUG_MSG( '### EMISSMERCURY: a LAND_FLUX' )

    CALL VEGEMIS( Input_Opt, State_Met, LVEGEMIS,  EHg0_dist, EHg0_vg,   RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a VEGEMIS' )

    CALL SOILEMIS( EHg0_dist, EHg0_so, State_Grid, State_Met )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SOILEMIS' )

    !==========================================================================
    !ENDIF
    !==========================================================================

    CALL SNOWPACK_MERCURY_FLUX( EHg0_snow,  LHGSNOW,  State_Chm, &
                                State_Grid, State_Met )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SNOW_FLUX' )

    CALL BIOMASSHG( Input_Opt, EHg0_bb, RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a BIOMASS' )

    IF ( LnoUSAemis ) THEN

       ! No Anthropogenic emissions over USA (cdh, 5/6/2009)
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          IF ( State_Grid%XMid(I,J) > -140 .AND. &
               State_Grid%XMid(I,J) < -40  .AND. &
               State_Grid%YMid(I,J) >  10  .AND. &
               State_Grid%YMid(I,J) <  70 ) THEN

             EHg0_An(I,J) = 0e+0_fp
             EHg2_An(I,J) = 0e+0_fp
             EHgP_An(I,J) = 0e+0_fp
             EHg0_bb(I,J) = 0e+0_fp
             EHg0_am(I,J) = 0e+0_fp

          ENDIF

       ENDDO
       ENDDO

    ENDIF

    CALL RESET_HG_DEP_ARRAYS
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: ' // &
         'a RESET_HG_DEP_ARRAYS' )

    ! If we are using the non-local PBL mixing,
    ! we need to initialize the EMIS_SAVE array (cdh, 08/27/09)
    ! EMIS_SAVE is now HG_EMIS (ckeller, 10/21/2014)
    IF ( LNLPBL ) HG_EMIS = 0.0e+0_fp

    ! Add Hg(0) source into State_Chm%Species [kg]
    CALL SRCHg0( Input_Opt,  State_Chm, State_Diag, State_Grid, State_Met, RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg0' )

    ! Add Hg(II) source into State_Chm%Species [kg]
    CALL SRCHg2( Input_Opt,  State_Chm, State_Diag, State_Grid, State_Met, RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg2' )

    ! Add HgP source into State_Chm%Species [kg]
    CALL SRCHgP( Input_Opt, State_Chm, State_Grid, State_Met, RC )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHgP' )

    ! Convert species units back to original unit
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'Routine EMISSMERCURY in mercury_mod.F90')
       RETURN
    ENDIF

  END SUBROUTINE EMISSMERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emithg
!
! !DESCRIPTION: Subroutine EMITHG directs emission either to the chemical
!  species array (State\_Chm%Species) directly or to Hg\_EMIS for use by the
!  non-local PBL mixing. This is a programming convenience.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMITHG( I, J, L, ID, E_HG, Input_Opt, State_Chm, State_Grid )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L, ID ! Grid boxes + species #
    REAL(fp),       INTENT(IN)    :: E_Hg        ! Hg emissions
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  27 Aug 2009 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(fp)          :: AM2, TS

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    !=================================================================
    ! EMITHG begins here!
    !=================================================================

    IF ( Input_Opt%LNLPBL ) THEN

       !--------------------------------------------------------------
       ! We are using FULL PBL MIXING (routine TURBDAY)
       !
       ! Save emissions for non-local PBL mixing or emit directly.
       ! Make sure that emitted mass is non-negative
       ! This is hear only for consistency with old code which warned
       ! of underflow error (cdh, 08/27/09)
       ! EMIS_SAVE is now HG_EMIS array. Convert kg to kg/m2/s
       ! (ckeller, 09/23/2014)
       !--------------------------------------------------------------

       ! Surface area [m2]
       AM2             = State_Grid%Area_M2(I,J)

       ! Emission timestep
       TS              = GET_TS_EMIS()

       ! Save emissions as [kg/m2/s].  These will be added
       ! to the chemical species array in routine DO_TEND
       ! (in mixing_mod.F90).
       HG_EMIS(I,J,ID) = HG_EMIS(I,J,ID) + ( MAX(E_HG,0e+0_fp) / AM2 / TS )

    ELSE

       !--------------------------------------------------------------
       ! We are using FULL PBL MIXING (routine TURBDAY)
       ! so add directly to the State_Chm%Species array
       !--------------------------------------------------------------

       ! Point to the chemical spcies array [kg]
       Spc             => State_Chm%Species

       ! Add emissions
       Spc(I,J,L,ID)   =  Spc(I,J,L,ID) + MAX( E_HG, 0e+0_fp )

       ! Free pointer
       Spc             => NULL()

    ENDIF

  END SUBROUTINE EMITHG
!EOP
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcHg0
!
! !DESCRIPTION: Subroutine SRCHg0 is the subroutine for Hg(0) emissions.
!  Emissions of Hg(0) will be distributed throughout the boundary layer.
!  (eck, cdh, bmy, 1/21/05, 4/6/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SRCHg0( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03, ND03, AD03_nat
#endif
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jan 2005 - N. (Eckley) Selin, C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      J,    L,        N,    NN,     PBL_MAX
    REAL(fp) :: DTSRCE, E_Hg, F_OF_PBL, T_Hg, T_Hg_An

    ! For fields from Input_Opt
    LOGICAL  :: LSPLIT, LPREINDHG, LGTMM

    !=================================================================
    ! SRCHg0 begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Copy values from Input_Opt
    LSPLIT    = Input_Opt%LSPLIT
    LPREINDHG = Input_Opt%LPREINDHG
    LGTMM     = Input_Opt%LGTMM

    ! Emission timestep [s]
    DTSRCE    = GET_TS_EMIS()

    ! Maximum extent of the PBL [model levels]
    PBL_MAX   = State_Met%PBL_MAX_L

    ! Loop over grid boxes
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, N, T_Hg_An, T_Hg, F_OF_PBL, E_Hg, NN)
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       IF ( LPREINDHG ) THEN  !eds

          ! Anthropogenic emissions
          T_Hg_An = 0e+0_fp

          ! No biomass burning emissions
          EHg0_bb(I,J) = 0e+0_fp

       ELSE

          ! Compute total anthropogenic Hg(0) emissions
          T_Hg_An = EHg0_an(I,J)

          IF ( LAnthroHgOnly ) THEN
             ! No other emissions
             EHg0_bb(I,J)      = 0e+0_fp
             EHg0_oc(I,J,:)    = 0e+0_fp
             EHg0_snow(I,J,:)  = 0e+0_fp
             IF ( LGTMM ) THEN
                EHg0_gtm(I,J) = 0e+0_fp
             ELSE
                EHg0_ln(I,J,:) = 0e+0_fp
                EHg0_vg(I,J)   = 0e+0_fp
                EHg0_so(I,J)   = 0e+0_fp
             ENDIF
          ENDIF

       ENDIF

       ! Compute total Hg(0) emissions (anthro+oceans+land+natural)
       IF ( LGTMM ) THEN
          T_Hg = T_Hg_An +                &
                 EHg0_bb(I,J) +           &
                 EHg0_oc(I,J,ID_Hg_tot) + &
                 EHg0_geo(I,J) +          &
                 EHg0_gtm(I,J) +          &
                 EHg0_snow(I,J,ID_Hg_tot)
       ELSE
          T_Hg = T_Hg_An +                &
                 EHg0_bb(I,J) +           &
                 EHg0_oc(I,J,ID_Hg_tot) + &
                 EHg0_ln(I,J,ID_Hg_tot) + &
                 EHg0_geo(I,J) +          &
                 EHg0_vg(I,J) +           &
                 EHg0_so(I,J) +           &
                 EHg0_snow(I,J,ID_Hg_tot)
       ENDIF

       !==============================================================
       ! Partition Hg0 throughout PBL; store into State_Chm%Species [kg]
       ! Now make sure State_Chm%Species does not underflow (cdh, bmy, 4/6/06)
       !==============================================================

       ! Loop up to max PBL level
       DO L = 1, PBL_MAX

          ! Fraction of box (I,J,L) w/in the PBL [unitless]
          F_OF_PBL = State_Met%F_OF_PBL(I,J,L)

          !-----------------
          ! Total Hg species
          !-----------------
          N    = Hg0_Id_List(ID_Hg_tot)
          E_Hg = F_OF_PBL * T_Hg * DTSRCE
          CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

          !-----------------
          ! Tagged species
          !-----------------
          IF ( LSPLIT ) THEN

             !--------------------
             ! Primary emissions
             !--------------------

             ! Anthro Canada Hg0
             N    = Hg0_Id_List(ID_Hg_can)
             E_Hg = F_OF_PBL * EHg0_can(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro United States Hg0
             N    = Hg0_Id_List(ID_Hg_usa)
             E_Hg = F_OF_PBL * EHg0_usa(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Central America Hg0
             N    = Hg0_Id_List(ID_Hg_cam)
             E_Hg = F_OF_PBL * EHg0_cam(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South America Hg0
             N    = Hg0_Id_List(ID_Hg_sam)
             E_Hg = F_OF_PBL * EHg0_sam(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro West Africa Hg0
             N    = Hg0_Id_List(ID_Hg_waf)
             E_Hg = F_OF_PBL * EHg0_waf(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Africa Hg0
             N    = Hg0_Id_List(ID_Hg_eaf)
             E_Hg = F_OF_PBL * EHg0_eaf(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Africa Hg0
             N    = Hg0_Id_List(ID_Hg_saf)
             E_Hg = F_OF_PBL * EHg0_saf(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro North Africa Hg0
             N    = Hg0_Id_List(ID_Hg_naf)
             E_Hg = F_OF_PBL * EHg0_naf(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro OECD Europe Hg0
             N    = Hg0_Id_List(ID_Hg_eur)
             E_Hg = F_OF_PBL * EHg0_eur(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Eastern Europe Hg0
             N    = Hg0_Id_List(ID_Hg_eeu)
             E_Hg = F_OF_PBL * EHg0_eeu(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Middle East Hg0
             N    = Hg0_Id_List(ID_Hg_mde)
             E_Hg = F_OF_PBL * EHg0_mde(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Former USSR Hg0
             N    = Hg0_Id_List(ID_Hg_sov)
             E_Hg = F_OF_PBL * EHg0_sov(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Asia Hg0
             N    = Hg0_Id_List(ID_Hg_sas)
             E_Hg = F_OF_PBL * EHg0_sas(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Asia Hg0
             N    = Hg0_Id_List(ID_Hg_eas)
             E_Hg = F_OF_PBL * EHg0_eas(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Southeast Asia Hg0
             N    = Hg0_Id_List(ID_Hg_sea)
             E_Hg = F_OF_PBL * EHg0_sea(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Japan Hg0
             N    = Hg0_Id_List(ID_Hg_jpn)
             E_Hg = F_OF_PBL * EHg0_jpn(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Oceania Hg0
             N    = Hg0_Id_List(ID_Hg_oce)
             E_Hg = F_OF_PBL * EHg0_oce(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Soil emissions
             N    = Hg0_Id_List(ID_Hg_so)
             E_Hg = F_OF_PBL * EHg0_so(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Biomass burning emissions
             N    = Hg0_Id_List(ID_Hg_bb)
             E_Hg = F_OF_PBL * EHg0_bb(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Geogenic emissions
             N    = Hg0_Id_List(ID_Hg_geo)
             E_Hg = F_OF_PBL * EHg0_geo(I,J) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             !--------------------
             ! Ocean re-emissions
             !--------------------

             ! Anthro Canada re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_can)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_can) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro United States re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_usa)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_usa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Central America re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_cam)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_cam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South America re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sam)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro West Africa re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_waf)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_waf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Africa re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_eaf)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eaf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Africa re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_saf)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_saf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro North Africa re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_naf)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_naf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro OECD Europe re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_eur)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eur) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Eastern Europe re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_eeu)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eeu) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Middle East re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_mde)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_mde) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Former Soviet Union re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sov)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sov) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Asia re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sas)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Asia re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_eas)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_eas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Southeast Asia re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sea)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sea) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Japan re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_jpn)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_jpn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Oceania re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_oce)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_oce) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Soil re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_so)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_so) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Biomass burning re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_bb)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_bb) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Geogenic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_geo)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_geo) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Mid-Atlantic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_atl)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_atl) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Atlantic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_nat)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_nat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer South Atlantic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_sat)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_sat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Pacific re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_npa)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_npa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Arctic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_arc)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_arc) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt,  State_Chm, State_Grid )

             ! Below Mixed Layer Antarctic re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_ant)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_ant) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer All Other Ocean re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_ocn)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_ocn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Residual Stratosphere from Spin-up re-emission from ocean
             N    = Hg0_Id_List(ID_Hg_str)
             E_Hg = F_OF_PBL * EHg0_oc(I,J,ID_Hg_str) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             !--------------------
             ! Land re-emissions
             !--------------------

             ! Anthro Canada re-emission from land
             N    = Hg0_Id_List(ID_Hg_can)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_can) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro United States re-emission from land
             N    = Hg0_Id_List(ID_Hg_usa)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_usa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Central America re-emission from land
             N    = Hg0_Id_List(ID_Hg_cam)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_cam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South America re-emission from land
             N    = Hg0_Id_List(ID_Hg_sam)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro West Africa re-emission from land
             N    = Hg0_Id_List(ID_Hg_waf)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_waf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Africa re-emission from land
             N    = Hg0_Id_List(ID_Hg_eaf)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eaf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Africa re-emission from land
             N    = Hg0_Id_List(ID_Hg_saf)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_saf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro North Africa re-emission from land
             N    = Hg0_Id_List(ID_Hg_naf)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_naf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro OECD Europe re-emission from land
             N    = Hg0_Id_List(ID_Hg_eur)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eur) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Eastern Europe re-emission from land
             N    = Hg0_Id_List(ID_Hg_eeu)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eeu) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Middle East re-emission from land
             N    = Hg0_Id_List(ID_Hg_mde)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_mde) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Former Soviet Union re-emission from land
             N    = Hg0_Id_List(ID_Hg_sov)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sov) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Asia re-emission from land
             N    = Hg0_Id_List(ID_Hg_sas)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Asia re-emission from land
             N    = Hg0_Id_List(ID_Hg_eas)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_eas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Southeast Asia re-emission from land
             N    = Hg0_Id_List(ID_Hg_sea)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sea) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Japan re-emission from land
             N    = Hg0_Id_List(ID_Hg_jpn)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_jpn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Oceania re-emission from land
             N    = Hg0_Id_List(ID_Hg_oce)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_oce) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Soil re-emission from land
             N    = Hg0_Id_List(ID_Hg_so)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_so) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Biomass burning re-emission from land
             N    = Hg0_Id_List(ID_Hg_bb)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_bb) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Geogenic re-emission from land
             N    = Hg0_Id_List(ID_Hg_geo)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_geo) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Mid-Atlantic re-emission from land
             N    = Hg0_Id_List(ID_Hg_atl)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_atl) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Atlantic re-emission from land
             N    = Hg0_Id_List(ID_Hg_nat)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_nat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer South Atlantic re-emission from land
             N    = Hg0_Id_List(ID_Hg_sat)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_sat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Pacific re-emission from land
             N    = Hg0_Id_List(ID_Hg_npa)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_npa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Arctic re-emission from land
             N    = Hg0_Id_List(ID_Hg_arc)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_arc) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Antarctic re-emission from land
             N    = Hg0_Id_List(ID_Hg_ant)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_ant) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer All Other Ocean re-emission from land
             N    = Hg0_Id_List(ID_Hg_ocn)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_ocn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Residual Stratosphere from Spin-up re-emission from land
             N    = Hg0_Id_List(ID_Hg_str)
             E_Hg = F_OF_PBL * EHg0_ln(I,J,ID_Hg_str) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             !--------------------
             ! Snow re-emissions
             !--------------------

             ! Anthro Canada re-emission from snow
             N    = Hg0_Id_List(ID_Hg_can)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_can) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro United States re-emission from snow
             N    = Hg0_Id_List(ID_Hg_usa)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_usa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Central America re-emission from snow
             N    = Hg0_Id_List(ID_Hg_cam)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_cam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South America re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sam)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sam) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro West Africa re-emission from snow
             N    = Hg0_Id_List(ID_Hg_waf)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_waf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Africa re-emission from snow
             N    = Hg0_Id_List(ID_Hg_eaf)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eaf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Africa re-emission from snow
             N    = Hg0_Id_List(ID_Hg_saf)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_saf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro North Africa re-emission from snow
             N    = Hg0_Id_List(ID_Hg_naf)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_naf) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro OECD Europe re-emission from snow
             N    = Hg0_Id_List(ID_Hg_eur)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eur) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Eastern Europe re-emission from snow
             N    = Hg0_Id_List(ID_Hg_eeu)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eeu) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Middle East re-emission from snow
             N    = Hg0_Id_List(ID_Hg_mde)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_mde) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Former Soviet Union re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sov)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sov) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro South Asia re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sas)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro East Asia re-emission from snow
             N    = Hg0_Id_List(ID_Hg_eas)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_eas) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Southeast Asia re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sea)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sea) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Japan re-emission from snow
             N    = Hg0_Id_List(ID_Hg_jpn)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_jpn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Anthro Oceania re-emission from snow
             N    = Hg0_Id_List(ID_Hg_oce)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_oce) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Soil re-emission from snow
             N    = Hg0_Id_List(ID_Hg_so)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_so) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Biomass burning re-emission from snow
             N    = Hg0_Id_List(ID_Hg_bb)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_bb) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Geogenic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_geo)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_geo) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Mid-Atlantic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_atl)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_atl) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Atlantic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_nat)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_nat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer South Atlantic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_sat)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_sat) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer North Pacific re-emission from snow
             N    = Hg0_Id_List(ID_Hg_npa)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_npa) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Arctic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_arc)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_arc) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer Antarctic re-emission from snow
             N    = Hg0_Id_List(ID_Hg_ant)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_ant) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Below Mixed Layer All Other Ocean re-emission from snow
             N    = Hg0_Id_List(ID_Hg_ocn)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_ocn) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             ! Residual Stratosphere from Spin-up re-emission from snow
             N    = Hg0_Id_List(ID_Hg_str)
             E_Hg = F_OF_PBL * EHg0_snow(I,J,ID_Hg_str) * DTSRCE
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

          ENDIF
       ENDDO

#ifdef BPCH_DIAG
       !==============================================================
       ! %%%%% ND03 (bpch) DIAGNOSTIC %%%%%
       !
       ! Total Hg(0) emissions [kg]
       ! 1=anthro; 3=from ocean; 4=land re-emission; 5=natural src
       !==============================================================

       IF ( ND03 > 0 ) THEN

          IF ( LGTMM ) THEN

             !eds 9/9/10
             DO NN = 1, N_HG_CATS
                AD03(I,J, 1,NN)=AD03(I,J, 1,NN)+(T_Hg_An          *DTSRCE)
                AD03(I,J, 3,NN)=AD03(I,J, 3,NN)+(EHg0_oc(I,J,NN)  *DTSRCE)
                AD03(I,J, 4,NN)=AD03(I,J, 4,NN)+(EHg0_gtm(I,J)    *DTSRCE)
                AD03(I,J, 5,NN)=AD03(I,J, 5,NN)+(EHg0_geo(I,J)    *DTSRCE)
                AD03(I,J,13,NN)=AD03(I,J,13,NN)+(EHg0_bb(I,J)     *DTSRCE)
                AD03(I,J,18,NN)=AD03(I,J,18,NN)+(EHg0_snow(I,J,NN)*DTSRCE)
             ENDDO

          ELSE

             !eds 9/9/10
             DO NN = 1, N_HG_CATS
                AD03(I,J, 1,NN)=AD03(I,J, 1,NN)+(T_Hg_An          *DTSRCE)
                AD03(I,J, 3,NN)=AD03(I,J, 3,NN)+(EHg0_oc(I,J,NN)  *DTSRCE)
                AD03(I,J, 4,NN)=AD03(I,J, 4,NN)+(EHg0_ln(I,J,NN)  *DTSRCE)
                AD03(I,J, 5,NN)=AD03(I,J, 5,NN)+(EHg0_geo(I,J)    *DTSRCE)
                AD03(I,J,13,NN)=AD03(I,J,13,NN)+(EHg0_bb(I,J)     *DTSRCE)
                AD03(I,J,14,NN)=AD03(I,J,14,NN)+(EHg0_vg(I,J)     *DTSRCE)
                AD03(I,J,15,NN)=AD03(I,J,15,NN)+(EHg0_so(I,J)     *DTSRCE)
                AD03(I,J,18,NN)=AD03(I,J,18,NN)+(EHg0_snow(I,J,NN)*DTSRCE)
             ENDDO

          ENDIF

          ! for preindustrial simulation, archive only soil,
          ! CDH- WHY ONLY SOIL??
          ! for present day archive soil, geogenic, biomass burning,
          ! vegetation, and rapid recycing
          IF ( LPREINDHG ) THEN !eds
             AD03_nat(I,J,:)=MAX(EHg0_so(I,J)*DTSRCE, SMALLNUM)
          ELSE
             IF ( LGTMM ) THEN
                DO N = 1, N_HG_CATS
                   AD03_nat(I,J,N) = DTSRCE * ( EHg0_geo(I,J)     + &
                                                EHg0_bb(I,J)      + &
                                                EHg0_gtm(I,J)     + &
                                                EHg0_snow(I,J,N)  )
                ENDDO
             ELSE
                DO N = 1, N_HG_CATS
                   AD03_nat(I,J,N) = DTSRCE * ( EHg0_ln(I,J,N)   + &
                                                EHg0_geo(I,J)    + &
                                                EHg0_bb(I,J)     + &
                                                EHg0_vg(I,J)     + &
                                                EHg0_so(I,J)     + &
                                                EHg0_snow(I,J,N) )
                ENDDO
             ENDIF
          ENDIF

       ENDIF
#endif

       !==============================================================
       ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
       !
       ! Save the various Hg0 emissions to fields of State_Diag
       ! in units of [kg/s].
       !
       ! NOTES (bmy, 10/24/18):
       ! (1) Normally these would go to the HEMCO_diagnostics*.nc
       !     file, but several of these diagnostics were defined as
       !     HEMCO manual diagnostics.  Therefore to preserve
       !     backwards compatibility with the existing Hg simulation,
       !     we save the emission fields into State_Diag.
       !
       ! (2) For now, save emissions only for total species
       !     and not the tagged species.  Colin Thackray says that
       !     the tagged Hg simulation is in some disrepair and is
       !     not widely used.
       !
       ! (3) The total Hg0 emissions diagnostic is not defined;
       !     users can sum up the individual categories in
       !     post-processing.
       !
       ! (4) In the current Hg simulation, HgP anthro emissions are
       !     added into the Hg2 species.  Therefore we only have
       !     a single emissions field of State_Diag for the
       !     combined Hg2 + HgP anthro emissions.
       !==============================================================

       ! Anthropogenic Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0anthro ) THEN
          State_Diag%EmisHg0anthro(I,J) = EHg0_an(I,J)
       ENDIF

       ! Biomass Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0biomass ) THEN
          State_Diag%EmisHg0biomass(I,J) = EHg0_bb(I,J)
       ENDIF

       ! Geogenic Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0geogenic ) THEN
          State_Diag%EmisHg0geogenic(I,J) = EHg0_geo(I,J)
       ENDIF

       ! Land Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0ocean ) THEN
          State_Diag%EmisHg0land(I,J) = EHg0_ln(I,J,1)
       ENDIF

       ! Oceanic Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0ocean ) THEN
          State_Diag%EmisHg0ocean(I,J) = EHg0_oc(I,J,1)
       ENDIF

       ! Snow Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0snow ) THEN
          State_Diag%EmisHg0snow(I,J) = EHg0_snow(I,J,1)
       ENDIF

       ! Soil Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0soil ) THEN
          State_Diag%EmisHg0soil(I,J) = EHg0_so(I,J)
       ENDIF

       ! Vegetation Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0vegetation ) THEN
          State_Diag%EmisHg0vegetation(I,J) = EHg0_vg(I,J)
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE SRCHg0
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcHg2
!
! !DESCRIPTION: Subroutine SRCHg2 is the subroutine for Hg(II) emissions.
!  Emissions of Hg(II) will be distributed throughout the boundary layer.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SRCHg2( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03, ND03
#endif
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  07 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      J,        L,    N,   PBL_MAX
    REAL(fp) :: DTSRCE, F_OF_PBL, E_Hg

    ! For values from Input_Opt
    LOGICAL  :: LSPLIT, LPREINDHG

    !=================================================================
    ! SRCHg2 begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Copy values from Input_Opt
    LSPLIT    = Input_Opt%LSPLIT
    LPREINDHG = Input_Opt%LPREINDHG

    ! Emission timestep [s]
    DTSRCE    = GET_TS_EMIS()

    ! Maximum extent of the PBL [model levels]
    PBL_MAX   = State_Met%PBL_MAX_L

    IF (.NOT. LPREINDHG ) THEN

       ! Loop over grid boxes
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, F_OF_PBL, E_Hg, N )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Loop up to the max PBL layer
          DO L = 1, PBL_MAX

             ! Fraction of box (I,J,L) w/in the PBL [unitless]
             F_OF_PBL = State_Met%F_OF_PBL(I,J,L)

             ! Partition total Hg2 into box (I,J,L) [kg]
             E_Hg     = F_OF_PBL * EHg2_an(I,J) * DTSRCE

             !---------------------------
             ! Total anthro Hg(II) [kg]
             !---------------------------
             N    = Hg2_Id_List(ID_Hg_tot)
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, State_Grid )

             !---------------------------
             ! Tagged anthro Hg(II) [kg]
             !---------------------------
             IF ( LSPLIT ) THEN

                ! Anthro Canada Hg2
                N    = Hg2_Id_List(ID_Hg_can)
                E_Hg = F_OF_PBL * EHg2_can(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro United States Hg2
                N    = Hg2_Id_List(ID_Hg_usa)
                E_Hg = F_OF_PBL * EHg2_usa(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Central America Hg2
                N    = Hg2_Id_List(ID_Hg_cam)
                E_Hg = F_OF_PBL * EHg2_cam(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro South America Hg2
                N    = Hg2_Id_List(ID_Hg_sam)
                E_Hg = F_OF_PBL * EHg2_sam(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro West Africa Hg2
                N    = Hg2_Id_List(ID_Hg_waf)
                E_Hg = F_OF_PBL * EHg2_waf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro East Africa Hg2
                N    = Hg2_Id_List(ID_Hg_eaf)
                E_Hg = F_OF_PBL * EHg2_eaf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro South Africa Hg2
                N    = Hg2_Id_List(ID_Hg_saf)
                E_Hg = F_OF_PBL * EHg2_saf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro North Africa Hg2
                N    = Hg2_Id_List(ID_Hg_naf)
                E_Hg = F_OF_PBL * EHg2_naf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro OECD Europe Hg2
                N    = Hg2_Id_List(ID_Hg_eur)
                E_Hg = F_OF_PBL * EHg2_eur(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Eastern Europe Hg2
                N    = Hg2_Id_List(ID_Hg_eeu)
                E_Hg = F_OF_PBL * EHg2_eeu(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Middle East Hg2
                N    = Hg2_Id_List(ID_Hg_mde)
                E_Hg = F_OF_PBL * EHg2_mde(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Former USSR Hg2
                N    = Hg2_Id_List(ID_Hg_sov)
                E_Hg = F_OF_PBL * EHg2_sov(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro South Asia Hg2
                N    = Hg2_Id_List(ID_Hg_sas)
                E_Hg = F_OF_PBL * EHg2_sas(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro East Asia Hg2
                N    = Hg2_Id_List(ID_Hg_eas)
                E_Hg = F_OF_PBL * EHg2_eas(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Southeast Asia Hg2
                N    = Hg2_Id_List(ID_Hg_sea)
                E_Hg = F_OF_PBL * EHg2_sea(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Japan Hg2
                N    = Hg2_Id_List(ID_Hg_jpn)
                E_Hg = F_OF_PBL * EHg2_jpn(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro Oceania Hg2
                N    = Hg2_Id_List(ID_Hg_oce)
                E_Hg = F_OF_PBL * EHg2_oce(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

             ENDIF
          ENDDO

#ifdef BPCH_DIAG
          !==============================================================
          ! %%%%% ND03 (bpch) DIAGNOSTICS %%%%%
          !
          ! Anthro Hg(II) [kg]
          ! NOTE: HgP is emitted into Hg2
          !==============================================================
          IF ( ND03 > 0 ) THEN
             AD03(I,J,6,1) = AD03(I,J,6,1) + ( EHg2_an(I,J) * DTSRCE )
          ENDIF
#endif

          !==============================================================
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Save the combined Hg2 + HgP emissions to State_Diag
          ! in units of [kg/s].
          !
          ! NOTES (bmy, 10/24/18):
          !
          ! (1) In the current Hg simulation, HgP anthro emissions
          !     are added into the Hg2 species.  Therefore we only
          !     have defined a single emissions field of State_Diag
          !     for the combined Hg2 + HgP anthro emissions.
          !==============================================================
          IF ( State_Diag%Archive_EmisHg2HgPanthro ) THEN
             State_Diag%EmisHg2HgPanthro(I,J) = EHg2_an(I,J)
          ENDIF

       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

  END SUBROUTINE SRCHg2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcHgp
!
! !DESCRIPTION: Subroutine SRCHgP is the subroutine for HgP emissions.
!  Emissions of HgP will be distributed throughout the boundary layer.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SRCHgP( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE DIAG03_MOD,         ONLY : AD03, ND03
#endif
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  07 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      J,        L,    N,   PBL_MAX
    REAL(fp) :: DTSRCE, F_OF_PBL, E_Hg

    ! For values from Input_Opt
    LOGICAL  :: LSPLIT, LPREINDHG

    !=================================================================
    ! SRCHgP begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Copy values from Input_Opt
    LSPLIT    = Input_Opt%LSPLIT
    LPREINDHG = Input_Opt%LPREINDHG

    ! Chemistry timestep [s]
    DTSRCE    = GET_TS_EMIS()

    ! Maximum extent of the PBL [model levels]
    PBL_MAX   = State_Met%PBL_MAX_L

    IF (.NOT. LPREINDHG) THEN

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, F_OF_PBL, E_Hg, N )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Loop up to PBL top layer
          DO L = 1, PBL_MAX

             ! Fraction of box (I,J,L) w/in the PBL [unitless]
             F_OF_PBL           = State_Met%F_OF_PBL(I,J,L)

             ! Partition HgP into box (I,J,L) [kg]
             E_Hg               = F_OF_PBL * EHgP_an(I,J) * DTSRCE

             !------------------------
             ! Total anthro HgP [kg]
             !------------------------
             N               = HgP_Id_List(ID_Hg_tot)
             CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                          State_Grid )

             !------------------------
             ! Tagged anthro HgP [kg]
             !------------------------
             IF ( LSPLIT ) THEN

                ! Anthro Canada HgP
                N            = HgP_Id_List(ID_Hg_can)
                E_Hg         = F_OF_PBL * EHgP_can(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

                ! Anthro United States HgP
                N            = HgP_Id_List(ID_Hg_usa)
                E_Hg         = F_OF_PBL * EHgP_usa(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Central America HgP
                N            = HgP_Id_List(ID_Hg_cam)
                E_Hg         = F_OF_PBL * EHgP_cam(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro South America HgP
                N            = HgP_Id_List(ID_Hg_sam)
                E_Hg         = F_OF_PBL * EHgP_sam(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro West Africa HgP
                N            = HgP_Id_List(ID_Hg_waf)
                E_Hg         = F_OF_PBL * EHgP_waf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro East Africa HgP
                N            = HgP_Id_List(ID_Hg_eaf)
                E_Hg         = F_OF_PBL * EHgP_eaf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro South Africa HgP
                N            = HgP_Id_List(ID_Hg_saf)
                E_Hg         = F_OF_PBL * EHgP_saf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro North Africa HgP
                N            = HgP_Id_List(ID_Hg_naf)
                E_Hg         = F_OF_PBL * EHgP_naf(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro OECD Europe HgP
                N            = HgP_Id_List(ID_Hg_eur)
                E_Hg         = F_OF_PBL * EHgP_eur(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Eastern Europe HgP
                N            = HgP_Id_List(ID_Hg_eeu)
                E_Hg         = F_OF_PBL * EHgP_eeu(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Middle East HgP
                N            = HgP_Id_List(ID_Hg_mde)
                E_Hg         = F_OF_PBL * EHgP_mde(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Former USSR HgP
                N            = HgP_Id_List(ID_Hg_sov)
                E_Hg         = F_OF_PBL * EHgP_sov(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro South Asia HgP
                N            = HgP_Id_List(ID_Hg_sas)
                E_Hg         = F_OF_PBL * EHgP_sas(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro East Asia HgP
                N            = HgP_Id_List(ID_Hg_eas)
                E_Hg         = F_OF_PBL * EHgP_eas(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Southeast Asia HgP
                N            = HgP_Id_List(ID_Hg_sea)
                E_Hg         = F_OF_PBL * EHgP_sea(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Japan HgP
                N            = HgP_Id_List(ID_Hg_jpn)
                E_Hg         = F_OF_PBL * EHgP_jpn(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )
                
                ! Anthro Oceania HgP
                N            = HgP_Id_List(ID_Hg_oce)
                E_Hg         = F_OF_PBL * EHgP_oce(I,J) * DTSRCE
                CALL EMITHG( I, J, L, N, E_Hg, Input_Opt, State_Chm, &
                             State_Grid )

             ENDIF
          ENDDO

#ifdef BPCH_DIAG
          !-----------------------------------------------------------
          ! %%%%% ND03 (bpch) diagnostics
          !
          ! Anthro HgP [kg]
          ! NOTE: This is zero, because in the current Hg simulation
          ! HgP is emitted into the Hg2 species (bmy, 10/26/18)
          !-----------------------------------------------------------
          IF ( ND03 > 0 ) THEN
             AD03(I,J,9,1) = AD03(I,J,9,1) + (EHgP_an(I,J) * DTSRCE)
          ENDIF
#endif

       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

  END SUBROUTINE SRCHgP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mercury_readyr
!
! !DESCRIPTION: Subroutine MERCURY\_READYR reads the year-invariant emissions
!  for Mercury from anthropogenic, ocean, and land sources.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MERCURY_READYR( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ERROR_MOD
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
    USE HCO_INTERFACE_MOD,  ONLY : GetHcoDiagn
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  06 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS
!
    REAL(fp), PARAMETER      :: SEC_PER_YR = 365.25e+0_fp*86400e+0_fp
!
! !LOCAL VARIABLES:
!
    ! For values from Input_Opt
    LOGICAL                  :: LDYNOCEAN
    LOGICAL                  :: LPREINDHG
    LOGICAL                  :: LSPLIT

    ! HEMCO update
    CHARACTER(LEN=63)        :: DgnName
    REAL(f4),        POINTER :: Ptr2D(:,:)

    ! Strings
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=512)       :: ErrMsg

    !=================================================================
    ! MERCURY_READYR begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = '-> at MERCURY_READYR (in GeosCore/mercury_mod.F90)'

    ! Initialize pointers
    Ptr2D     => NULL()

    ! Copy values from Input_Opt
    LDYNOCEAN = Input_Opt%LDYNOCEAN
    LPREINDHG = Input_Opt%LPREINDHG
    LSPLIT    = Input_Opt%LSPLIT

    !=================================================================
    ! Anthropogenic Emissions
    !=================================================================

    !=================================================================
    ! Now get the emission fields through HEMCO. The emissions
    ! settings and source data are now set in the HEMCO configuration
    ! file (specified in input.geos) and emissions are calculated
    ! based on the content of the configuration file.
    ! Here, we just import the HEMCO diagnostics fields (defined
    ! in hcoi_gc_diagn_mod.F90) and pass the emission fields to
    ! the corresponding arrays (EHg0_an, etc.). Those internal
    ! arrays are in kg/s, so we need to convert the diagnostics
    ! from kg/m2/s to kg/s.
    !
    ! NOTE: the current implementation is a hybrid version between
    ! HEMCO and the old mercury code, primarily to make sure that
    ! all the diagnostics based on EHg* still work. Because of
    ! that, we convert all emissions from kg/m2/s to kg/s and then
    ! back to kg/m2/s when we pass them to HG_EMIS (via SRCHg0 and
    ! SRCHg2). HG_EMIS is used for the non-local PBL mixing scheme.
    ! (ckeller, 09/24/2014)
    !=================================================================

    ! Get HEMCO state object
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       ErrMsg = 'HcoState not associated!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! No anthropogenic emissions in preindustrial simulation
    IF ( LPREINDHG ) THEN

       EHg0_an = 0e+0_fp  !eds
       EHg2_an = 0e+0_fp
       EHgP_an = 0e+0_fp

    ELSE

       ! ---------------
       ! Hg0 emissions
       ! ---------------

       ! Anthropogenic emissions
       DgnName = 'HG0_ANTHRO'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not get HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       EHg0_an =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
       Ptr2D   => NULL()

       ! Artisanal emissions
       DgnName = 'HG0_ARTISANAL'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not get HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       EHg0_am =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
       Ptr2D   => NULL()

       ! ---------------------------------------------
       ! Hg2 emissions
       ! Note: HgP emissions are added to Hg2 by HEMCO
       ! ---------------------------------------------

       ! Anthropogenic emissions
       ! DgnName = 'HG2_ANTHRO'
       ! CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       ! IF ( RC /= HCO_SUCCESS ) THEN
       !    ErrMsg = 'Could not get HEMCO field ' // TRIM( DgnName )
       !    CALL GC_Error( ErrMsg, RC, ThisLoc )
       !    RETURN
       ! ENDIF
       ! EHg2_an =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
       ! Ptr2D   => NULL()

    ENDIF

    ! Natural emissions
    DgnName = 'HG0_NATURAL'
    CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Could not get HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    EHg0_geo =  Ptr2D(:,:) * HcoState%Grid%AREA_M2%Val(:,:)
    Ptr2D    => NULL()

    ! Soil distribution
    CALL HCO_GetPtr( HcoState, 'HG0_SOILDIST', Ptr2D, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Could not get pointer to HEMCO field HG0_SOILDIST!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    EHg0_dist =  Ptr2D(:,:)
    Ptr2D     => NULL()

    !=================================================================
    ! Print totals to the screen in [Gg/yr]
    !=================================================================
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, 210   )
    WRITE( 6, '(a)' )
    WRITE( 6, 211   ) SUM( EHg0_an ) * SEC_PER_YR * 1e-6_fp
    WRITE( 6, 213   ) SUM( EHg0_ln ) * SEC_PER_YR * 1e-6_fp
    WRITE( 6, 214   ) SUM( EHg0_geo ) * SEC_PER_YR * 1e-6_fp

    ! Only write ocean total if we are doing offline ocean
    IF ( .not. LDYNOCEAN ) THEN
       WRITE( 6, 217   ) SUM( EHg0_oc ) * SEC_PER_YR * 1e-6_fp
    ENDIF

    WRITE( 6, 215   ) SUM( EHg2_an ) * SEC_PER_YR * 1e-6_fp
    WRITE( 6, 216   ) SUM( EHgP_an ) * SEC_PER_YR * 1e-6_fp
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )

    ! FORMAT strings
210 FORMAT( 'M E R C U R Y   E M I S S I O N S' )
211 FORMAT( 'Total Anthro     Hg(0)  : ', f7.3, ' [Gg/yr]' )
213 FORMAT( 'Total Re-Emitted Hg(0)  : ', f7.3, ' [Gg/yr]' )
214 FORMAT( 'Total Natural    Hg(0)  : ', f7.3, ' [Gg/yr]' )
215 FORMAT( 'Total Anthro     Hg(II) : ', f7.3, ' [Gg/yr]' )
216 FORMAT( 'Total Anthro     HgP    : ', f7.3, ' [Gg/yr]' )
217 FORMAT( 'Total Ocean      Hg(0)  : ', f7.3, ' [Gg/yr]' )

  END SUBROUTINE MERCURY_READYR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_OxidPointers
!
! !DESCRIPTION: Subroutine Set_OxidPointers gets the Hg oxidant concentration data
! read by HEMCO. The pointers only need to be established once. Target data
! is automatically updated through HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_OxidPointers( Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorological State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (based on set_brypointers)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=16)   :: ThisName
    CHARACTER(LEN=255)  :: PREFIX, FIELDNAME, ThisLoc
    CHARACTER(LEN=1024) :: ErrMsg
    INTEGER             :: N

    !=================================================================
    ! Set_OxidPointers begins here
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Set_OxidPointers (in module GeosCore/mercury_mod.F90)'

    ! Do for each Hg oxidant species
    DO N = 1, N_Hg_Oxid

       ! Skip if species is not defined
       IF ( GC_HgOxid_TrID(N) <= 0 ) CYCLE

       ! Get oxidant name
       ThisName = State_Chm%SpcData(GC_HgOxid_TrID(N))%Info%Name

       ! Construct field name using species name
       FIELDNAME = 'GLOBAL_'//TRIM(ThisName)

       ! Get pointer to this field. These are the concentrations (molec cm-3).
       CALL HCO_GetPtr( HcoState, FIELDNAME, HgOxidPtr(N)%Data, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO! Oxidant data '//&
                   'is expected to be listed in the HEMCO configuration '  //&
                   'file. This error occured when trying to get field '     //&
                  TRIM( FieldName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDDO !N

    ! Do for each aerosol species
    DO N = 1, N_Aer + N_Dust

       !------------------------------
       ! AOD
       !------------------------------

       ! Get aerosol species name
       FIELDNAME = 'AOD_' // TRIM( AerSpcNames(N) )

       ! Get pointer to this field. These are AODs.
       CALL HCO_GetPtr( HcoState, FIELDNAME, AeroPtr(N)%AOD, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO for ' //&
                  TRIM( FieldName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !------------------------------
       ! Area (cm2 cm-3)
       !------------------------------

       ! Get aerosol species name
       FIELDNAME = 'Area_' // TRIM( AerSpcNames(N) )

       ! Get pointer to this field. These are AODs.
       CALL HCO_GetPtr( HcoState, FIELDNAME, AeroPtr(N)%Area, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO for ' //&
                  TRIM( FieldName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !------------------------------
       ! Radi (cm)
       !------------------------------

       ! Get aerosol species name
       FIELDNAME = 'Radi_' // TRIM( AerSpcNames(N) )

       ! Get pointer to this field. These are AODs.
       CALL HCO_GetPtr( HcoState, FIELDNAME, AeroPtr(N)%Radi, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO for ' //&
                  TRIM( FieldName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !------------------------------
       ! Mass concentration (ug m-3)
       !------------------------------

       ! Get aerosol species name
       FIELDNAME = 'Conc_' // TRIM( AerSpcNames(N) )

       ! Mass conc of some species are not read
       IF ( FIELDNAME .eq. 'Conc_DST5' .or. FIELDNAME .eq. 'Conc_DST6' .or. &
            FIELDNAME .eq. 'Conc_DST7' .or. FIELDNAME .eq. 'Conc_BGSULF' .or. &
            FIELDNAME .eq. 'Conc_ICEI' ) CYCLE

       ! Get pointer to this field.
       CALL HCO_GetPtr( HcoState, FIELDNAME, AeroPtr(N)%Conc, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO for ' //&
                  TRIM( FieldName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDDO !N

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Set_OxidPointers
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_HgOxidConc
!
! !DESCRIPTION: Subroutine Set_HgOxidConc transfers oxidant concentration fields
!               to State_Chm after applying diurnal variation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_HgOxidConc( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!

    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH


!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    ! Scalars
    LOGICAL            :: prtDebug
    INTEGER            :: I, J, L, N, NN     ! lon, lat, lev, indexes
    CHARACTER(LEN=60)  :: Prefix             ! utility string
    CHARACTER(LEN=255) :: ThisLoc            ! routine location
    CHARACTER(LEN=255) :: ErrMsg             ! message

    ! Pointers
    REAL(fp),  POINTER        :: Spc(:,:,:,:)

    ! Objects
    TYPE(Species),    POINTER :: SpcInfo

    ! Initialize pointers
    SpcInfo     => NULL()
    Spc         => NULL()

    !================================================================
    ! Get_HgOxConc begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Get_HgOxConc (in GeosCore/mercury_mod.F90)'

    ! Copy values from Input_Opt
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Point to the chemical spcies array
    Spc             => State_Chm%Species

    !=================================================================
    ! Set instantaneous species concentrations
    !=================================================================

    ! Set species concentration to monthly mean value
    DO N=1, N_Hg_Oxid
        ! Get speciesID
        NN = GC_HgOxid_TrID(N)

        ! Proceed only if species is in the database
        IF ( NN <= 0 ) CYCLE

        ! Get value from pointer to monthly mean field
        Spc(:,:,:,NN)  = HgOxidPtr(N)%Data(:,:,:)
    ENDDO

    ! Impose diurnal cycle
    ! Compute sum of cosine of the solar zenith angle over a 24 hour day
     CALL OHNO3TIME( State_Grid )
     IF ( prtDebug ) CALL DEBUG_MSG( 'Set_HgOxConc: a OHNO3TIME' )

    ! Calculate instantaneous HOx
    CALL DiurnalHOx( State_Chm, State_Grid, State_Met )

    ! Partition NOx based on O3 and JNO2_Inst
    CALL PartNOx( State_Chm, State_Grid, State_Met )

    ! Apply dirunal cycle and partition XOx based on O3, NO and J_XO
    CALL PartXOx( State_Chm, State_Grid, State_Met )

    ! Add BrOx from springtime polar bromine explosion events
    IF ( LPOLARBR ) CALL PolarBrOx( State_Chm, State_Grid, State_Met )

  END SUBROUTINE Set_HgOxidConc
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
  SUBROUTINE DiurnalHOx( State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE TIME_MOD,       ONLY : GET_TS_CHEM
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState

!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
!
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (modified from GET_OH)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I, J, L
    REAL(fp) :: DiurnalFac, C_OH, C_HO2

    ! Pointer to Species array
    REAL(fp),  POINTER    :: Spc(:,:,:,:)

    ! Initialize pointers
    Spc         => NULL()

    !=================================================================
    ! DirunalHOx begins here!
    !=================================================================

    ! Point to the chemical spcies array
    Spc             => State_Chm%Species

    ! Loop over gridcells
!$OMP PARALLEL DO                                              &
!$OMP DEFAULT( SHARED )                                        &
!$OMP PRIVATE( I,       J,       L                           ) &
!$OMP PRIVATE( DiurnalFac,       C_OH,      C_HO2            )
    DO L=1, State_Grid%NZ
    DO J=1, State_Grid%NY
    DO I=1, State_Grid%NX

        ! Get HOx concentrations
        C_OH   = Spc( I,J,L,id_OH  )
        C_HO2  = Spc( I,J,L,id_HO2 )

        ! Test for sunlight...
        IF ( State_Met%SUNCOS(I,J) > 0e+0_fp .and.  TCOSZ(I,J) > 0e+0_fp ) THEN

           DiurnalFac = ( State_Met%SUNCOS(I,J)  / TCOSZ(I,J) )* &
                        ( 86400e0_fp             / GET_TS_CHEM() )

           ! Make sure factor is not negative
           DiurnalFac = MAX( DiurnalFac, 0e+0_fp )

        ELSE

           ! At night, OH goes to zero
           DiurnalFac = 0e+0_fp

        ENDIF

        Spc( I,J,L,id_OH  ) = C_OH  * DiurnalFac
        Spc( I,J,L,id_HO2 ) = C_HO2 * DiurnalFac

    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE DiurnalHOx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
! FUNCTION GetDiurnalXOx( I, J, State_Grid, State_Met ) RESULT( SCALEF )
!
! !USES:
!
!   USE ERROR_MOD,          ONLY : SAFE_DIV
!   USE State_Grid_Mod,     ONLY : GrdState
!   USE State_Met_Mod,      ONLY : MetState
!   USE TIME_MOD,           ONLY : GET_TS_CHEM
!   USE TIME_MOD,           ONLY : GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
!   INTEGER,        INTENT(IN)  :: I              ! Longitude index
!   INTEGER,        INTENT(IN)  :: J              ! Latitude index
!   TYPE(GrdState), INTENT(IN)  :: State_Grid   ! Grid State object
!   TYPE(MetState), INTENT(IN)  :: State_Met      ! Meteorology State object
!
! !RETURN VALUE:
!
!   REAL(fp)                    :: SCALEF   ! Scale mean to inst.
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
!
!   REAL(fp)              :: LOCALTIME, HOUR
!   ! Test for sunlight...
!   IF ( (State_Met%SUNCOS(I,J) > 0e+0_fp) .and. (TTDAY(I,J) > 0e+0_fp) ) THEN
!      ! Use a constant function if daylight is < 2 hours or > 22hours
!      IF ( ( TTDAY(I,J) < 120e+0_fp  ) .OR. ( TTDAY(I,J) > 1320e+0_fp ) ) THEN
!         SCALEF = SAFE_DIV( 1440e+0_fp, TTDAY(I,J), 0e+0_fp )
!      ELSE
!         ! Local time: 0-23.999
!         LOCALTIME = GET_LOCALTIME( I, J, 1, State_Grid )
!         ! Interpolate the real hour to lie within an equinoctal day
!         ! i.e. between 6-18
!         HOUR = 12e+0_fp + ( LOCALTIME - 12e+0_fp ) * 720e+0_fp / ( TTDAY(I,J))
!         ! Multiplicative factor for the diurnal cycle
!         ! This formula comes from the Holmes et al (2009) MBL box model
!         ! which fit RGM observations for Okinawa and a Pacific cruise
!         SCALEF = (1e+6_fp - 1e+2_fp * (HOUR - 6e+0_fp) ** 4e+0_fp + &
!                  ( 18e+0_fp - HOUR ) * 1e+6_fp / 12e+0_fp ) / 2e+6_fp
!         ! Normalize the multiplicative factor to have a 24-h mean of 1
!         SCALEF = SAFE_DIV( SCALEF, (4e-4_fp * TTDAY(I,J)), 0e+0_fp )
!      ENDIF
!      ! Make sure that diurnal scaling factor is non-negative
!      SCALEF = MAX( SCALEF, 0e+0_fp )
!   ELSE
!      ! At night, OH goes to zero
!      SCALEF = 0e+0_fp
!   ENDIF
! END FUNCTION GetDiurnalXOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PolarBrOx
!
! !DESCRIPTION: Subroutine PolarBr calculates BrOx during bromine explosion events
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PolarBrOx( State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE TIME_MOD,           ONLY : GET_MONTH
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
!
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (modified from cdh's and jaf's
!                            GET_BR)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

!
! !LOCAL VARIABLES:
!    ! Scalars
    INTEGER              :: I, J, L            ! lon, lat, lev indexes

    REAL(fp)              :: FPBL
    ! Make polar BrO variable rather than fixed, varying as a factor of
    ! solar radiation and temperatures. (jaf, 11/29/11)
    REAL(fp)              :: BRO_POLAR_PPTV, BR_POLAR_PPTV, O3_POLAR, SWRAD
    REAL(fp)              :: BRO_POLAR_CONC, BR_POLAR_CONC, BRO_CONC, BR_CONC

    LOGICAL               :: IS_MOSTLY_ICE
!
! !DEFINED PARAMETERS:
!

    ! Assume 5 ppb O3 when BrO present. Pohler et al. (2010) shows 5 ppb O3
    ! can exist for 3 ppt < [BrO] < 40 ppt. (jaf, 11/29/11)
    REAL(fp), PARAMETER  :: O3_POLAR_PPBV = 5e+0_fp
    ! Parameters for calculating Br/BrO photostationary state
    ! BrO J value, /s
    ! Rate coefficient BrO + NO -> Br + NO2, cm3/molec/s
    REAL(fp), PARAMETER   :: K_BRO_NO = 2.1e-11_fp
    ! Rate coefficient Br + O3 -> BrO + O2, cm3/molec/s
    REAL(fp), PARAMETER   :: K_BR_O3  = 1.2e-12_fp
    ! Concentration of NO, based on 10pptv, molec/cm3
    REAL(fp), PARAMETER   :: C_NO     = 2.5e+8_fp

    ! Pointer to Species array
    REAL(fp),  POINTER    :: Spc(:,:,:,:)

    ! Initialize pointers
    Spc         => NULL()

    !=================================================================
    ! PolarBrOx begins here!
    !=================================================================
    ! Point to the chemical spcies array
    Spc             => State_Chm%Species

    ! Loop over gridcells and add BrOx in polar regions
!$OMP PARALLEL DO                                              &
!$OMP DEFAULT( SHARED )                                        &
!$OMP PRIVATE( I,       J,      L     )
    DO L=1, State_Grid%NZ
    DO J=1, State_Grid%NY
    DO I=1, State_Grid%NX

        IF ( ((State_Grid%YMid(I,J) > 50e+0_fp) .AND.           &
            (GET_MONTH() >= 2) .AND. (GET_MONTH() <= 6)) .OR.  &
            ((State_Grid%YMid(I,J) < -50e+0_fp) .AND.            &
            (GET_MONTH() >= 8) .AND. (GET_MONTH() <= 12)) )  THEN

            !----------------------------------------------------------------
            ! Add Br in the polar PBL
            !----------------------------------------------------------------
            FPBL = State_Met%F_UNDER_PBLTOP(I,J,L)

            !
            ! Bromine in the polar boundary layer requires the following
            ! criteria be met:
            ! - In the PBL
            ! - Downward shortwave radiation > 100 W/m2 (Pohler et al. 2010)
            ! - Sea ice exists (>50% native boxes have >10% ice cover)
            ! - Breaks in sea ice exist (<100% native boxes have >90% ice cover)
            ! - Month is between Feb & June (Arctic) or Aug & Dec (Antarctic)
            !   based on http://bro.aeronomie.be/level3_monthly.php?cmd=map
            ! - Temperature is less than 0C
            !
            ! If these criteria are met, BrO is a function of ambient temp.
            ! with [BrO] based on findings from Pohler et al. (2010) and
            ! Prados-Roman et al. (2011). O3 used to convert BrO to Br is 5
            ! ppb, based on data from Pohler et al. (2010).
            !----------------------------------------------------------------
            SWRAD         = State_Met%SWGDN(I,J)

            IS_MOSTLY_ICE = ( State_Met%SEAICE00(I,J) <= 0.5e+0_fp .AND. &
                              State_Met%SEAICE90(I,J) < 1e+0_fp )

            IF ( (FPBL > 0e+0_fp) .AND. (IS_MOSTLY_ICE) .AND. &
                 (SWRAD > 1e+2_fp) .AND. (State_Met%TS(I,J) <= 273e+0_fp) ) THEN

              ! Get BrOx concentration from species data
              BRO_CONC = Spc(I,J,L,id_BrO)
              BR_CONC  = Spc(I,J,L,id_Br)

              ! [BrO] is a linear function of temperature derived based on
              ! results from Pohler et al. (2010), Prados-Roman et al. (2011)
              ! and ability to match Hg0 seasonal cycle at Alert. (jaf,
              ! 12/24/11)
              IF ( State_Met%TS(I,J) <= 253e+0_fp ) THEN
                 BRO_POLAR_PPTV = 20e+0_fp
              ELSE
                 BRO_POLAR_PPTV = -1e+0_fp * ( State_Met%TS(I,J) - 253e+0_fp ) + 20e+0_fp
              ENDIF

              ! Convert O3 to molec/cm3
              O3_POLAR  = O3_POLAR_PPBV * 1.0e-9_fp * State_Met%AIRNUMDEN(I,J,L)

              ! Compute polar Br, BrO concentrations in pptv
              BrO_POLAR_PPTV = BrO_POLAR_PPTV * FPBL
              Br_POLAR_PPTV  = BRO_POLAR_PPTV * ( JBrO_INST(I,J,L) + K_BRO_NO * C_NO ) / &
                                                ( K_BR_O3 * O3_POLAR )

              ! Convert Br and BrO to molec/cm3
              BRO_POLAR_CONC = BRO_POLAR_PPTV * 1.0e-12_fp * State_Met%AIRNUMDEN(I,J,L)
              BR_POLAR_CONC  = BR_POLAR_PPTV  * 1.0e-12_fp * State_Met%AIRNUMDEN(I,J,L)

              ! Add Br to State_Chem array
              Spc(I,J,L,id_BrO) = BRO_CONC + BRO_POLAR_CONC
              Spc(I,J,L,id_Br ) = BR_CONC  + BR_POLAR_CONC

            ENDIF! Polar BrO criteria

        ENDIF ! Month

    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE PolarBrOx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PartXOx
!
! !DESCRIPTION: Subroutine PartXOx partitions XOx based on PSS with ozone
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PartXOx( State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : SAFE_DIV
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState

!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
!
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (modified from hmh's GET_NO2)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I, J, L
    REAL(fp) :: C_Br,    C_BrO,    C_O3,     C_BrOx,   C_NO
    REAL(fp) :: C_Cl,    C_ClO,    C_ClOx
    REAL(fp) :: k_Br_O3, k_BrO_NO, F_Br_BrO, J_BrO
    REAL(fp) :: k_Cl_O3, k_ClO_NO, F_Cl_ClO, J_ClO
    REAL(fp) :: DiurnalFac
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: A_BrO_NO     = 8.8e-12_fp
    REAL(fp), PARAMETER :: EdivR_BrO_NO = -260e+0_fp
    REAL(fp), PARAMETER :: A_ClO_NO     = 6.4e-12_fp
    REAL(fp), PARAMETER :: EdivR_ClO_NO = -290e+0_fp
    REAL(fp), PARAMETER :: A_Br_O3      = 1.6e-11_fp
    REAL(fp), PARAMETER :: EdivR_Br_O3  = 780e+0_fp
    REAL(fp), PARAMETER :: A_Cl_O3      = 2.3e-11_fp
    REAL(fp), PARAMETER :: EdivR_Cl_O3  = 200e+0_fp

    !  =======================================================================
    !  NOTES:
    !
    !  (R1) XO + hv -> X + O3,  j
    !  (R2) XO + NO -> X + NO2, k_XO_NO
    !  (R2) O3 + X ->  XO + O2, k_X_O3
    !
    !  NOx steady state:
    !  [X]/[XO] = (j+k_XO_NO * C_NO) / (k_X_O3 * C_O3)
    !
    !***************************************************************************

    ! Pointer to Species array
    REAL(fp),  POINTER    :: Spc(:,:,:,:)

    ! Initialize pointers
    Spc         => NULL()

    !=================================================================
    ! PartXOx begins here!
    !=================================================================

    ! Point to the chemical spcies array
    Spc             => State_Chm%Species

    ! Loop over gridcells
!$OMP PARALLEL DO                                               &
!$OMP DEFAULT( SHARED )                                         &
!$OMP PRIVATE( I,       J,        L                           ) &
!$OMP PRIVATE( C_Br,    C_BrO,    C_O3,     C_BrOx,   C_NO    ) &
!$OMP PRIVATE( C_Cl,    C_ClO,    C_ClOx                      ) &
!$OMP PRIVATE( k_Br_O3, k_BrO_NO, F_Br_BrO, J_BrO             ) &
!$OMP PRIVATE( k_Cl_O3, k_ClO_NO, F_Cl_ClO, J_ClO             ) &
!$OMP PRIVATE( DiurnalFac                                     )
    DO L=1, State_Grid%NZ
    DO J=1, State_Grid%NY
    DO I=1, State_Grid%NX

        ! Get species concentrations
        C_Br   = Spc( I,J,L,id_Br  )
        C_Cl   = Spc( I,J,L,id_Cl  )
        C_BrO  = Spc( I,J,L,id_BrO )
        C_ClO  = Spc( I,J,L,id_ClO )
        C_NO   = Spc( I,J,L,id_NO  )
        C_O3   = Spc( I,J,L,id_O3  )

        C_BrOx = C_Br + C_BrO
        C_ClOx = C_Cl + C_ClO

        ! Test for sunlight...
        IF ( (State_Met%SUNCOS(I,J) > 0e+0_fp) .and. (TTDAY(I,J) > 0e+0_fp) ) THEN

            ! Use a constant function for XOx
            DiurnalFac = SAFE_DIV( 1440e+0_fp, TTDAY(I,J), 0e+0_fp )

            ! Apply dirunal scale to XOx
            C_BrOx = C_BrOx * DiurnalFac
            C_ClOx = C_ClOx * DiurnalFac

            ! Calculate temperature dependent reaction rates for partitioning
            k_BrO_NO = A_BrO_NO * exp( -EdivR_BrO_NO / State_Met%T(I,J,L))
            k_ClO_NO = A_ClO_NO * exp( -EdivR_ClO_NO / State_Met%T(I,J,L))
            k_Br_O3  = A_Br_O3  * exp( -EdivR_Br_O3  / State_Met%T(I,J,L))
            k_Cl_O3  = A_Cl_O3  * exp( -EdivR_Cl_O3  / State_Met%T(I,J,L))

            ! Instantaneous J
            J_BrO    = JBrO_Inst( I, J, L )
            J_ClO    = JClO_Inst( I, J, L )

            ! Fraction of [X]/[XO]
            F_Br_BrO = SAFE_DIV( J_BrO+k_BrO_NO*C_NO, k_Br_O3*C_O3, 0e+0_fp )
            F_Cl_ClO = SAFE_DIV( J_ClO+k_ClO_NO*C_NO, k_Cl_O3*C_O3, 0e+0_fp )

            ! Species concentrations
            C_Br     =  C_BrOx * F_Br_BrO/( 1e+0_fp+F_Br_BrO )
            C_BrO    =  C_BrOx - C_Br
            C_Cl     =  C_ClOx * F_Cl_ClO/( 1e+0_fp+F_Cl_ClO )
            C_ClO    =  C_ClOx - C_Cl
        ELSE
            ! At night, XOx goes to zero
            C_Br     =  0e+0_fp
            C_BrO    =  0e+0_fp
            C_Cl     =  0e+0_fp
            C_ClO    =  0e+0_fp
        ENDIF

        Spc(I,J,L,id_Br ) = C_Br
        Spc(I,J,L,id_Cl ) = C_Cl
        Spc(I,J,L,id_BrO) = C_BrO
        Spc(I,J,L,id_ClO) = C_ClO

    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE PartXOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PartNOx
!
! !DESCRIPTION: Subroutine PartNOx partitions NOx based on PSS with ozone
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PartNOx( State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : SAFE_DIV
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState

!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
!
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (modified from hmh's GET_NO2)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I, J, L
    REAL(fp) :: C_O3, C_NOx
    REAL(fp) :: J_NO2, k3,   F_NO2
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: A = 3.0e-12_fp
    REAL(fp), PARAMETER :: EdivR = 1.5e+3_fp

    !  =======================================================================
    !  NOTES:
    !
    !  (R1) NO2 + hv -> NO + O,  j
    !  (R2) O +O2 -> O3,         k2
    !  (R3) O3 + NO -> NO2 + O2, k3
    !
    !  [NOx] = [NO] + [NO2]
    !
    !  NOx steady state:
    !  j[NO2] = k3[O3][NO]
    !  j[NO2] = k3[O3]([NOx]-[NO2])
    !  [NO2](j+k3[O3] = k3[O3][NOx]
    !  [NO2]/[NOx] = k3[O3]/(j+k3[O3])
    !
    !   k3 = A exp(-E / RT) Arrhenius Equation
    !   A = 3.0e-12 Seinfeld & Pandis
    !   E/R = 1500
    !
    !***************************************************************************

    ! Pointer to Species array
    REAL(fp),  POINTER    :: Spc(:,:,:,:)

    ! Initialize pointers
    Spc         => NULL()
    !=================================================================
    ! PartNOx begins here!
    !=================================================================

    ! Point to the chemical spcies array
    Spc             => State_Chm%Species

    ! Loop over gridcells
!$OMP PARALLEL DO                                              &
!$OMP DEFAULT( SHARED )                                        &
!$OMP PRIVATE( I,       J,      L       )                      &
!$OMP PRIVATE( J_NO2,  k3,      F_NO2,   C_O3,  C_NOx   )
    DO L=1, State_Grid%NZ
    DO J=1, State_Grid%NY
    DO I=1, State_Grid%NX
        ! Test for sunlight...
        IF ( (State_Met%SUNCOS(I,J) > 0e+0_fp) .and. (TTDAY(I,J) > 0e+0_fp) ) THEN

            k3 = A*exp(-EdivR/State_Met%T(I,J,L))

            ! Get NOx and O3 concentrations (molec cm-3)
            C_NOx = Spc(I,J,L,id_NO) + Spc(I,J,L,id_NO2)
            C_O3  = Spc(I,J,L,id_O3)

            ! Instantaneous JNO2
            J_NO2 = JNO2_INST( I, J, L )

            ! Fraction of NO2/NOx
            F_NO2 = SAFE_DIV( k3*C_O3, J_NO2+k3*C_O3, 0e+0_fp )

        ELSE
            ! NO goes to zero
            F_NO2 = 1e+0_fp
        ENDIF

        Spc(I,J,L,id_NO2) = F_NO2 * C_NOx
        Spc(I,J,L,id_NO)  = ( 1e+0_fp - F_NO2 ) * C_NOx

    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE PartNOx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the
!  solar zenith angle over a 24 hour day, as well as the total length of
!  daylight.  This is needed to scale the offline OH and NO3 concentrations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OHNO3TIME( State_Grid )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : GET_NHMSb,   GET_ELAPSED_SEC
    USE TIME_MOD,       ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Park & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE :: FIRST = .TRUE.
    INTEGER       :: I, J, L, N, NT, NDYSTEP
    REAL(fp)      :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)      :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)      :: SUNTMP(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! OHNO3TIME begins here!
    !=================================================================

    !  Solar declination angle (low precision formula, good enough for us):
    A0 = 0.006918
    A1 = 0.399912
    A2 = 0.006758
    A3 = 0.002697
    B1 = 0.070257
    B2 = 0.000907
    B3 = 0.000148
    R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

    DEC = A0 - A1*cos(  R) + B1*sin(  R) &
             - A2*cos(2*R) + B2*sin(2*R) &
             - A3*cos(3*R) + B3*sin(3*R)

    LHR0 = int(float( GET_NHMSb() )/10000.)

    ! Only do the following at the start of a new day
    IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN

       ! Zero arrays
       TTDAY(:,:) = 0e+0_fp
       TCOSZ(:,:) = 0e+0_fp
       COSZM(:,:) = 0e+0_fp

       ! NDYSTEP is # of chemistry time steps in this day
       NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 3600 / GET_TS_CHEM()

       ! NT is the elapsed time [s] since the beginning of the run
       NT = GET_ELAPSED_SEC()

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Zero SUNTMP array
          SUNTMP = 0e+0_fp

          ! Loop over surface grid boxes
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R( I, J )

             TIMLOC = real(LHR0) + real(NT)/3600.0 + State_Grid%XMid(I,J)/15.0

             DO WHILE (TIMLOC < 0)
                TIMLOC = TIMLOC + 24.0
             ENDDO

             DO WHILE (TIMLOC > 24.0)
                TIMLOC = TIMLOC - 24.0
             ENDDO

             AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP(I,J) = sin(YMID_R) * sin(DEC) + &
                           cos(YMID_R) * cos(DEC) * cos(AHR)

             ! TCOSZ is the sum of SUNTMP at location (I,J)
             ! Do not include negative values of SUNTMP
             TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(I,J), 0e+0_fp )

             ! COSZM is the peak value of SUMTMP during a day at (I,J)
             ! (rjp, bmy, 3/30/04)
             COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(I,J) )

             ! TTDAY is the total daylight time at location (I,J)
             IF ( SUNTMP(I,J) > 0e+0_fp ) THEN
                TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() ) / 60e+0_fp
             ENDIF
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

          ! Increment elapsed time [sec]
          NT = NT + GET_TS_CHEM()
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

  END SUBROUTINE OHNO3TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_hg2_seasalt_lossrate
!
! !DESCRIPTION: Subroutine CALC\_HG2\_SEASALT\_LOSSRATE calculates the loss
!  rate of RGM (/s) by uptake of RGM into sea salt aerosol for each model
!  grid. Return value is a loss frequency (/s)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_HG2_SEASALT_LOSSRATE( State_Grid, State_Met )
!
! !USES:
!
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  The formula used here is a least-squares fit to the full-physics model of
!  sea-salt aerosol emissions, hydroscopic growth, mass-transport limited
!  uptake of Hg(II), and aerosol deposition presented by Holmes et al. (2009)
!  See Holmes et al. 2010 for evaluation of this parameterization.
!  (cdh, 11/25/09)
!
! !REVISION HISTORY:
!  25 Nov 2009 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)      :: U10M, S
    REAL(fp)      :: LOSS_FREQ
    INTEGER       :: I, J
    REAL(fp),SAVE :: TABLE_S(21), TABLE_U10(20)
    REAL(fp),SAVE :: TABLE_UPTAKE(21,20)
    REAL(fp)      :: SFCWINDSQR

    ! Flag for first call
    LOGICAL, SAVE :: FIRST=.TRUE.

    !=================================================================
    ! HG2_SEASALT_LOSSRATE begins here!
    !=================================================================

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, U10M, S, LOSS_FREQ, SFCWINDSQR )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Only calculate deposition via sea salt over water
       IF ( State_Met%IsWater(I,J) ) THEN

          ! Wind speed at 10m altitude [m/s]
          SFCWINDSQR = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2
          U10M       = SQRT( SFCWINDSQR )

          ! Don't allow wind >20 m/s which is the limit of this parameterization
          U10M       = MAX( MIN( U10M, 20e+0_fp ), 1e+0_fp )

          ! Relative humidity as a saturation ratio
          ! Use the relative humidity of the lowest layer, although this is
          ! lower than the higher layers of the MBL
          !
          ! Don't allow supersaturation, as [Cl-] is undefined for RH>=1
          ! Cap RH at 99%, Don't allow RH < 75% as happens in coastal areas
          S = MAX( MIN( State_Met%RH(I,J,1), 99e+0_fp ), 75e+0_fp ) * 1e-2_fp

          LOSS_FREQ = 1e-10_fp * ( 1e+0_fp - EXP( -57.758e+0_fp * &
                      (1e+0_fp-S) ) ) * &
                      EXP( -1.9351e+0_fp  * U10M + &
                            9.0047e+0_fp  * SQRT( U10M ) + &
                            0.14788e+0_fp * U10M**1.5e+0_fp )

          ! Loss frequency must be positive
          LOSS_FREQ = MAX( LOSS_FREQ, 1e-10_fp )

       ELSE

          ! No loss over land
          LOSS_FREQ = 0e+0_fp

       ENDIF

       HG2_SEASALT_LOSSRATE(I,J) = LOSS_FREQ

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE CALC_HG2_SEASALT_LOSSRATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: Set_Hetrates
!
! !DESCRIPTION: Main heterogenous chemistry driver routine.  Sets up the
!  vector of heterogeneous chemistry rates for the KPP chemistry solver.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Set_HetRates( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!

    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState


!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure

! !REMARKS:
!
! !REVISION HISTORY:

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! Scalars
      LOGICAL            :: prtDebug
      CHARACTER(LEN=60)  :: Prefix             ! utility string
      CHARACTER(LEN=255) :: ThisLoc            ! routine location
      CHARACTER(LEN=255) :: ErrMsg             ! message


      INTEGER  :: I, J, L, N
      INTEGER  :: NAEROTYPE
      REAL(fp) :: QICE,     QLIQ

      LOGICAL, SAVE :: FIRST = .TRUE.
      LOGICAL       :: STRATBOX

      ! Cloud parameters
      Real(fp)         :: rLiq, ALiq, VLiq, CLDFr
      Real(fp)         :: rIce, AIce, VIce
      Real(fp)         :: TEMPK, XTEMP, XDENA, SUNCOS

      ! Volume of air (cm3)
      Real(fp)         :: VAir

      ! Hg(II) sticking coefficient
      Real(fp), Parameter :: ALPHA_Hg2 = 0.1e+0_fp

        ! Molecular weights [g/mole]
      Real(fp), Parameter :: MW_HgBrNO2 = 327e+0_fp
      Real(fp), Parameter :: MW_HgBrHO2 = 313e+0_fp
      Real(fp), Parameter :: MW_HgBrBrO = 376.4e+0_fp
      Real(fp), Parameter :: MW_HgBrClO = 332e+0_fp
      Real(fp), Parameter :: MW_HgBrOH  = 297.5e+0_fp
      Real(fp), Parameter :: MW_HgBr2   = 360.4e+0_fp
      Real(fp), Parameter :: MW_HgClNO2 = 282e+0_fp
      Real(fp), Parameter :: MW_HgClHO2 = 269e+0_fp
      Real(fp), Parameter :: MW_HgClBrO = 332e+0_fp
      Real(fp), Parameter :: MW_HgClClO = 287.5e+0_fp
      Real(fp), Parameter :: MW_HgClOH  = 253e+0_fp
      Real(fp), Parameter :: MW_HgClBr  = 316e+0_fp



      ! Arrays
      Real(fp)           :: XArea(15), XRadi(15), XVol(15)

      ! Pointer to Species array
      REAL(fp),  POINTER    :: Spc(:,:,:,:)

      ! Initialize pointers
      Spc         => NULL()

      !====================================================================
      ! Set_Hetrates begins here!
      !====================================================================

      ! Assume success
      RC        = GC_SUCCESS
      ErrMsg    = ''
      ThisLoc   = ' -> at Set_HetRates (in GeosCore/mercury_mod.F90)'

      ! Copy values from Input_Opt
      prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

      ! Zero scalars and arrays
      NAEROTYPE     = N_Dust + N_Aer

      ! Initialize logicals
      STRATBOX      = .FALSE.

      ! Zero array
      AerHetRates   = 0e+0_fp
      OrgHetRates   = 0e+0_fp

      ! Point to the chemical spcies array
      Spc             => State_Chm%Species

      ! Loop over gridcells
!$OMP PARALLEL DO                                               &
!$OMP DEFAULT( SHARED )                                         &
!$OMP PRIVATE( I,       J,      L,     N       )                &
!$OMP PRIVATE( XArea,   XRadi,  XVol, rLiq, ALiq, VLiq, CLDFr ) &
!$OMP PRIVATE( rIce, AIce, VIce, TEMPK, XTEMP, XDENA, SUNCOS  ) &
!$OMP PRIVATE( QICE, QLIQ, Vair )
      DO L=1, State_Grid%NZ
      DO J=1, State_Grid%NY
      DO I=1, State_Grid%NX

          ! If we are not in the troposphere don't proceed
          IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

          !--------------------------------------------------------------------
          ! Get species concentrations [molec/cm3]
          !--------------------------------------------------------------------
          ! !EEM: Add for ClNO2 yield calculation
          ! IND = Ind_( 'SALA' )
          ! IF (IND .le. 0) THEN
          !    SPC_SALA    = 0.0e+0_fp
          ! ELSE
          !    SPC_SALA    = spcVec(IND)
          ! ENDIF

          !--------------------------------------------------------------------
          ! Aerosol Physical Properties
          !--------------------------------------------------------------------

          DO N=1, nAeroType

              ! Aerosol specific surface area, cm2(aerosol)/cm3(air)
              XAREA(N) = AeroPtr(N)%Area(I,J,L)

              ! Aerosol effective radius, cm
              XRADI(N) = AeroPtr(N)%Radi(I,J,L)

              ! Aerosol specific volume, cm3(aerosol)/cm3(air)
              XVOL(N)  = XAREA(N) * XRADI(N) / 3e+0_fp

          ENDDO

          !--------------------------------------------------------------------
          ! Get fields from State_Met, State_Chm, and Input_Opt
          !--------------------------------------------------------------------

          TEMPK  = State_Met%T(I,J,L)              ! Temperature [K]
          XTEMP  = sqrt(State_Met%T(I,J,L))        ! Square root of temperature
          XDENA  = State_Met%AIRNUMDEN(I,J,L)      ! Dry air density [molec/cm3]
          SUNCOS = State_Met%SUNCOSmid(I,J)        ! COS(SZA),midpt of chem timestep
          VAir   = State_Met%AIRVOL(I,J,L)*1.0e6_fp! Volume of air (cm3)
          QICE   = State_Met%QI(I,J,L)             ! Ice   mix ratio [kg/kg dry air]
          QLIQ   = State_Met%QL(I,J,L)             ! Water mix ratio [kg/kg dry air]

          !--------------------------------------------------------------------
          !  Calculate parameters for cloud halogen chemistry
          !  under the new scheme (SDE 2016-12-21)
          !--------------------------------------------------------------------

          ! Get cloud physical parameters
          !      CALL Cld_Params( I, J, L, XDenA, VAir, TempK, QLiq, QIce, State_Met, &
          !                       rLiq,  ALiq,  VLiq, rIce,  AIce,  VIce, CLDFr )


          !--------------------------------------------------------------------
          ! Calculate heterogeous uptake in clouds
          !--------------------------------------------------------------------

          !       HET(id_HGBRNO2, 1) = CloudHet( 'HG2', CldFr, Aliq, &
          !                       Aice, rLiq, rIce, TempK, XDenA )

          !--------------------------------------------------------------------
          ! Calculate heterogeneous uptake on aqueous aerosols
          !--------------------------------------------------------------------

          IF ( id_HGBRNO2 > 0 ) THEN
              AerHetRates( I,J,L,id_HGBRNO2 ) = AerHet( MW_HgBrNO2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGBRNO2 ) = AerHet( MW_HgBrNO2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGBRHO2 > 0 ) THEN
              AerHetRates( I,J,L,id_HGBRHO2 ) = AerHet( MW_HgBrHO2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGBRHO2 ) = AerHet( MW_HgBrHO2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGBROH > 0 ) THEN
              AerHetRates( I,J,L,id_HGBROH  ) = AerHet( MW_HgBrOH, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGBROH  ) = AerHet( MW_HgBrOH, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGBRBRO > 0 ) THEN
              AerHetRates( I,J,L,id_HGBRBRO ) = AerHet( MW_HgBrBrO, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGBRBRO ) = AerHet( MW_HgBrBrO, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGBRCLO > 0 ) THEN
              AerHetRates( I,J,L,id_HGBRCLO ) = AerHet( MW_HgBrClO, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGBRCLO ) = AerHet( MW_HgBrClO, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGBR2 > 0 ) THEN
              AerHetRates( I,J,L,id_HGBR2   ) = AerHet( MW_HgBr2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGBR2   ) = AerHet( MW_HgBr2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGCLNO2 > 0 ) THEN
              AerHetRates( I,J,L,id_HGCLNO2 ) = AerHet( MW_HGCLNO2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGCLNO2 ) = AerHet( MW_HgClNO2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGCLHO2 > 0 ) THEN
              AerHetRates( I,J,L,id_HGCLHO2 ) = AerHet( MW_HGCLHO2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGCLHO2 ) = AerHet( MW_HgClHO2, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGCLOH > 0 ) THEN
              AerHetRates( I,J,L,id_HGCLOH  ) = AerHet( MW_HGCLOH, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGCLOH ) = AerHet( MW_HgClOH, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGCLBRO > 0 ) THEN
              AerHetRates( I,J,L,id_HGCLBRO ) = AerHet( MW_HGCLBrO, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGCLBRO ) = AerHet( MW_HgClBrO, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGCLCLO > 0 ) THEN
              AerHetRates( I,J,L,id_HGCLCLO ) = AerHet( MW_HGClClO, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGCLCLO ) = AerHet( MW_HgClClO, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
          IF ( id_HGCLBR > 0 ) THEN
              AerHetRates( I,J,L,id_HGCLBR  ) = AerHet( MW_HgClBr, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 1 )
              OrgHetRates( I,J,L,id_HGCLBR ) = AerHet( MW_HgClBr, ALPHA_Hg2, &
                                                        XAREA, XRADI, XDENA, XTEMP, 2 )
          ENDIF
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      ! Free pointer memory
      Spc => NULL()

  END SUBROUTINE Set_HetRates
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AerHet
!
! !DESCRIPTION: Set the heterogenous chemistry rate for NO3.
!\\
!\\
! !INTERFACE:
!
FUNCTION AerHet( MW, XSTKCF, XAREA, XRADI, XDENA, XTEMP, TYPE ) &
    RESULT( HetRate )
  !
  ! !Uses
  USE GCKPP_HetRates,       ONLY : ARSL1K
  !
  ! !INPUT PARAMETERS:
  !
  ! Rate coefficients
  REAL(fp), INTENT(IN) :: MW, XSTKCF, XDENA, XTEMP
  REAL(fp), INTENT(IN) :: XAREA(:), XRADI(:)
  INTEGER              :: TYPE !1:Sea-salt & sulfate ; 2:OA
  !
  ! !RETURN VALUE:
  !
  REAL(fp)             :: HetRate
  !
  ! !REMARKS:
  !
  ! !REVISION HISTORY:
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  ! !LOCAL VARIABLES:
  !
  INTEGER  :: N
  REAL(fp) :: ADJUSTEDRATE

  ! Initialize
  HetRate      = 0.0_fp
  ADJUSTEDRATE = 0.0_fp


  IF ( TYPE .EQ. 1 ) THEN
      ! Loop over aerosol types
      DO N = 1, N_Dust + N_Aer

          ! Uptake happens only on sulfate and sea-salt aerosols
          IF ( N.eq.8 .or. N.eq.11 .or. N.eq.12 ) THEN
              ! Reaction rate for surface of aerosol
              ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                  (MW**0.5_FP))
             !PRINT *, N, ADJUSTEDRATE
          ELSE
              ADJUSTEDRATE = 0e+0_fp
          ENDIF

          ! Add to overall reaction rate
          HetRate = HetRate + ADJUSTEDRATE

      ENDDO
  ENDIF

  IF ( TYPE .EQ. 2 ) THEN
      ! Organic aerosol only
      HetRate =ARSL1K(XAREA(10),XRADI(10),XDENA,XSTKCF,XTEMP, &
                  (MW**0.5_FP))
  ENDIF
END FUNCTION AerHet
!EOC
  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Transport Model                  !
  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: partitionhg2
  !
  ! !DESCRIPTION: Subroutine PARTITIONHG2 splits Hg(II) into gas and aerosol
  !  portions  according to the thermodynamic equilibrium determined by
  !  temperature and aerosol surface area. This is done only for dust, BC and OCPO.
  !
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE PARTITIONHG2( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!

    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState


!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  12-Jul-2010 - H. Amos     - Add option to partition Hg2 according to Fg/Fp
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
! Scalars
  INTEGER           :: I, J, L, N
  REAL(fp)          :: Hg2gasTOT, Hg2TOT, DryPM, k_ad, RatPartGas, Fgas

  REAL(fp)          :: Hg2SpcConc(20)

  CHARACTER(len=255):: ErrMsg, ThisLoc

  ! Pointers
  REAL(fp), POINTER :: Spc(:,:,:,:)

  !=================================================================
  ! PARTITIONHG2 begins here!
  !=================================================================

      ! Assume success
      RC  =  GC_SUCCESS

      ErrMsg    = ''
      ThisLoc   = ' -> at PartitionHg2 (in GeosCore/mercury_mod.F90)'


      ! Point to the chemical species array [kg]
      Spc => State_Chm%Species

!$OMP PARALLEL DO       &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( I, J, L, N, Hg2gasTOT, Hg2TOT, Hg2SpcConc ) &
!$OMP PRIVATE( DryPM, k_ad, RatPartGas, Fgas )
      DO L=1, State_Grid%NZ
      DO J=1, State_Grid%NY
      DO I=1, State_Grid%NX

          !Initialize
          Hg2gasTOT = 0e+0_fp
          Hg2TOT  = 0e+0_fp
          DryPM   = 0e+0_fp
          Hg2SpcConc = 0e+0_fp

          ! Dry PM concentration
          DO N=1, n_Dust+n_Aer

            ! Add dust and BC concentration
            ! BC includes OCPO
            IF ( N .le. 4 .or. N .eq. 9 ) &

            ! Otherwise add to PM conc (ug/m3)
            DryPM = DryPM + AeroPtr(N)%Conc(I,J,L)

          ENDDO

          ! Calculate adsorption coefficient (m-3/ug)
          ! This is from Amos et al 2012
          k_ad = 10e+0_fp**( ( 2.5e+3_fp / State_Met%T(I,J,L)) - 10e+0_fp )

          ! HgII particle-gas ratio (Rat = HgII_ads/HgII_gas)
          RatPartGas = k_ad * DryPM


          ! Loop over all HgII gas spcies
          DO N=1, nHg2gasSpc

              ! Total Hg(II) (gas +aerosol)
              Hg2SpcConc(N) = Spc(I,J,L,Map_Hg2gas(N))

          ENDDO

          ! Total HgII gas
          Hg2gasTOT = sum(Hg2SpcConc)

          ! Caluclate total Hg2 (in gas-phase + adsorbed on dust,BC & OCPO)
          ! intitialize with only aerosol phase HgP
          HG2TOT = Spc(I,J,L,id_HgIIads) + Hg2gasTOT

          ! Gas fraction
          Fgas = 1e+0_fp / (1e+0_fp + RatPartGas)

          ! Aerosol portion
          Spc(I,J,L,id_HgIIads) = Hg2TOT * (1e+0_fp - Fgas)

          ! Gas portion
          ! Loop over all HgII gas spcies
          DO N=1, nHg2gasSpc

              ! Total Hg(II) (gas +aerosol)
              Spc(I,J,L,Map_Hg2gas(N)) = Hg2TOT * Fgas * Hg2SpcConc(N)/Hg2gasTOT

          ENDDO

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer memory
      Spc => NULL()

END SUBROUTINE PARTITIONHG2
!EOC
  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Transport Model                  !
  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: offlineocean_readmo
  !
  ! !DESCRIPTION: Subroutine OFFLINEOCEAN\_READMO reads the monthly varying
  !     offline ocean evasion emissions if LDYNOCEAN is FALSE. Will not actually
  !     need mixed layer depth when i get stuff from Yanxu
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE OFFLINEOCEAN_READMO( State_Chm, State_Diag, State_Grid, &
      State_Met, FLUX, RC )
      !
      ! !USES:
      !
      USE ErrCode_Mod
      USE TIME_MOD,          ONLY : EXPAND_DATE,GET_YEAR, GET_TS_EMIS
      USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH, GET_MONTH
#ifdef BPCH_DIAG
    USE DIAG03_MOD,        ONLY : AD03, ND03
#endif
      USE State_Chm_Mod,     ONLY : ChmState
      USE State_Diag_Mod,    ONLY : DgnState
      USE State_Grid_Mod,    ONLY : GrdState
      USE State_Met_Mod,     ONLY : MetState
      USE HCO_INTERFACE_MOD, ONLY : HcoState
      USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
      !
      ! !INPUT PARAMETERS:
      !
      TYPE(ChmState), INTENT(IN)    :: State_Chm    ! Chemistry State object
      TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
      TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
      !
      ! !INPUT/OUTPUT PARAMETERS:
      !
      TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State Object
      !
      ! !OUTPUT PARAMETERS:
      !
      REAL*8,         INTENT(OUT)   :: FLUX(State_Grid%NX,State_Grid%NY,1)
      INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
      !
      ! !REVISION HISTORY:
      !  12 Aug 2015 - H. Horowitz - Initial version
      !  See https://github.com/geoschem/geos-chem for complete history
      !------------------------------------------------------------------------------
      !BOC
      !
      ! !LOCAL VARIABLES:
      !
      INTEGER                 :: I, J, M, DN(12) ! M is mon., DN is days in mon
      INTEGER                 :: NNN
      INTEGER                 :: THISYEAR
      INTEGER                 :: THISMONTH

      LOGICAL, SAVE           :: FIRST = .TRUE.
      LOGICAL                 :: IS_OCEAN_BOX
      INTEGER                 :: NN, N
      REAL(fp)                :: A_M2,     DTSRCE
      REAL(fp)                :: CHg0aq,   CHg0,     vi, Hg0aqtemp
      REAL(fp)                :: TC,       TK,       Kw
      REAL(fp)                :: Sc,       ScCO2,    USQ
      REAL(fp)                :: FRAC_L,   FRAC_O,   H, D
      REAL(fp)                :: FUP(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
      REAL(fp)                :: FDOWN(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
      REAL(fp)                :: Hg0aq(State_Grid%NX,State_Grid%NY,N_Hg_CATS)
      REAL(fp)                :: MHg0_air

      ! Conversion factor from [cm/h * ng/L] --> [kg/m2/s]
      REAL(fp),  PARAMETER    :: TO_KGM2S = 1.0e-11_fp / 3600.0e+0_fp

      ! Small numbers to avoid dividing by zero
      REAL(fp),  PARAMETER    :: SMALLNUM = 1.0e-32_fp

      REAL(fp)                :: SFCWINDSQR

      ! Pointers
      ! We need to define local arrays to hold corresponding values
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER       :: STT(:,:,:,:)

      ! Characters
      CHARACTER(LEN=255)      :: ThisLoc
      CHARACTER(LEN=512)      :: ErrMsg

      !=================================================================
      ! OFFLINEOCEAN_READMO begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at OFFLINEOCEAN_READMO (in GeosCore/mercury_mod.F90)'

      ! Get month
      THISMONTH = GET_MONTH()
      M         = THISMONTH

      ! Days in each month (will use later) 9/16/15 hmh
      DN =  (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

      !-----------------------------------------------------------------
      ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
      !
      ! Zero flux arrays of State_Diag to prevent leftover values
      ! from the last timestep from being included in the averaging
      !-----------------------------------------------------------------
      IF ( State_Diag%Archive_FluxHg0fromOceanToAir ) THEN
          State_Diag%FluxHg0fromOceanToAir = 0.0_f4
      ENDIF

      IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
          State_Diag%FluxHg0fromAirToOcean = 0.0_f4
      ENDIF

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! read monthly ocean evasion  !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF ( ITS_A_NEW_MONTH() ) THEN

          CALL HCO_GetPtr( HcoState, 'GLOBAL_OCEAN', OCEAN_CONC, RC )
          IF ( RC /= GC_SUCCESS ) THEN
              ErrMsg = 'Cannot get pointer to HEMCO field GLOBAL_OCEAN!'
              CALL GC_Error( ErrMsg, RC, ThisLoc )
              RETURN
          ENDIF

      ENDIF

      ! Only doing Hg0 overall, should add trap for LSPLIT (cpt - 2017)
      DO NNN = 1, N_Hg_CATS
          IF ( NNN .eq. 1 ) THEN
              Hg0aq(:,:,NNN) = OCEAN_CONC(:,:,1)
          ELSE
              Hg0aq(:,:,NNN) = 0.0e+0_fp
          ENDIF
      ENDDO

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS()

      STT => State_Chm%Species

      ! Loop over surface boxes
      !$OMP PARALLEL DO       &
      !$OMP DEFAULT( SHARED ) &
      !$OMP PRIVATE( I,   A_M2,    vi,     ScCO2 )               &
      !$OMP PRIVATE( J,   NN,      TC,     TK )                  &
      !$OMP PRIVATE( N,   CHg0,    FRAC_L, FRAC_O  )             &
      !$OMP PRIVATE( H,   Kw,      CHg0aq, Hg0aqtemp, MHg0_air ) &
      !$OMP PRIVATE( IS_OCEAN_BOX, Sc,     Usq, D )              &
      !$OMP SCHEDULE( DYNAMIC )
      DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

              ! Grid box surface area [m2]
              A_M2       = State_Grid%Area_M2( I, J )

              ! Initialize values
              Kw         = 0e0_fp
              TK         = 0e0_fp
              TC         = 0e0_fp

              ! Get fractions of land and ocean in the grid box [unitless]
              ! Use fractional land type information in MERRA. Also make sure
              ! we do not use boxes that are mostly sea ice for consistency
              ! FROCEAN is a constant, so to get correct ocean fraction we
              ! need to subtract the sea ice fraction. Don't let the fraction
              ! be less than zero (jaf, 4/26/11)
              FRAC_L       = State_Met%FRLAND(I,J)
              FRAC_O       = MAX( State_Met%FROCEAN(I,J) - &
                  State_Met%FRSEAICE(I,J), 0e0_fp )
              IS_OCEAN_BOX = ( ( FRAC_O > 0e0_fp ) .and. &
                  ( State_Met%SEAICE00(I,J)  > 0.5e0_fp ) )

              IF ( (IS_OCEAN_BOX) ) THEN

                  !--------------------------------------------------------------
                  ! Sea surface temperature in both [K] and [C]
                  !--------------------------------------------------------------
                  ! where TSKIN is the temperature (K) at the ground/sea surface
                  ! (Use as surrogate for SST, cap at freezing point)
                  TK     = MAX( State_Met%TSKIN(I,J), 273.15e0_fp )
                  TC     = TK - 273.15e0_fp

                  !==============================================================
                  ! Volatilisation of Hg0
                  !==============================================================

                  ! Henry's law constant (gas->liquid) [unitless] [L water/L air]
                  ! (ref: Andersson et al. 2008)
                  H      = EXP( ( -2404.3e0_fp / TK ) + 6.92e0_fp )

                  ! Viscosity as a function of changing temperatures
                  ! (ref: Loux 2001)
                  ! The paper says the viscosity is given in cP but us really P
                  ! and we therefor multiply with 100 to get cP.
                  vi    = ( 10**( ( 1301.0e0_fp / ( 998.333e0_fp + 8.1855e0_fp &
                      * ( TC - 20.0e0_fp )+ 0.00585e0_fp &
                      * ( TC - 20.0e0_fp )**2 ) ) - 3.30233e0_fp ) ) * 100.0e0_fp

                  ! Schmidt # for Hg [unitless]
                  ! Sc = v/D = kinematic viscosity/diffusivity
                  ! (ref: Poissant et al 2000; Wilke and Chang 1995)
                  ! to correct for seawater D0 is decreased by 6% as suggested
                  ! by Wanninkhof (1992)
                  D = 7.4e-8_fp * sqrt( 2.26e0_fp * 18.0e0_fp ) * TK / &
                      ( ( 14.8e0_fp**0.6e0_fp ) *vi )
                  Sc   = ( 0.017e0_fp * EXP( -0.025e0_fp * TC ) ) / D

                  ! Schmidt # of CO2 [unitless] for CO2 in seawater at 20 degrees C
                  ! The value is set to a constant based on other ocean studies
                  ! (Gardfeld et al. 2003, Rolfhus & Fitzgerald2004,Mason et al.2001)
                  !
                  ! Correction of the Schmidt # with temperature based on Poissant
                  ! et al. (2000) (for freshwatersystems).
                  ScCO2  = 644.7e0_fp + TC * ( -6.16e0_fp + TC * ( 0.11e0_fp))

                  ! Square of surface (actually 10m) wind speed [m2/s2]
                  Usq    = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2

                  !------------------------------------------------------
                  ! Parameterizations for calculating water side mass trasfer
                  ! coefficient
                  !------------------------------------------------------
                  ! Mass transfer coefficient [cm/h], from Nightingale et al. 2000
                  Kw     = ( 0.25e0_fp * Usq ) / SQRT( Sc / ScCO2 )

                  ! Loop over all Hg categories
                  DO NN = 1, N_Hg_CATS

                      ! Hg0 tracer number (for STT)
                      N = Hg0_Id_List(NN)

                      !--------------------------------------------------------
                      ! Calculate oceanic and gas-phase concentration of Hg(0)
                      !--------------------------------------------------------

                      ! Concentration of Hg(0) in the ocean [ng/L]
                      ! now converting from Hg0aq in mol/m3 to ng/L
                      CHg0aq = Hg0aq(I,J,NN) *200.59e0_fp * 1.0e9_fp / 1.0e3_fp

                      ! Gas phase Hg(0) concentration: convert [kg] -> [ng/L]
                      MHg0_air = STT(I,J,1,N)
                      CHg0     = MHg0_air *  1.0e9_fp /State_Met%AIRVOL(I,J,1)

                      !--------------------------------------------------------
                      ! Compute flux of Hg(0) from the ocean to the air
                      !--------------------------------------------------------

                      ! Compute ocean flux of Hg0 [cm/h*ng/L]
                      FLUX(I,J,NN)  = Kw * ( CHg0aq - ( CHg0 / H ) )

                      !Extra diagnostic: compute flux up and flux down
                      FUP(I,J,NN)   = ( Kw * CHg0aq )
                      FDOWN(I,J,NN) = ( Kw * CHg0 / H )

                      !--------------------------------------------------
                      ! Convert [cm/h*ng/L] --> [kg/m2/s] --> [kg/s]
                      ! Also account for ocean fraction of grid box
                      FLUX(I,J,NN)  = FLUX(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O

                      ! hmh 5/11/16 reverting to old version and uncommenting here
                      FUP(I,J,NN)  = FUP(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O
                      FDOWN(I,J,NN)  = FDOWN(I,J,NN) * TO_KGM2S * A_M2 * FRAC_O
                      !--------------------------------------------------

                      !--------------------------------------------------------
                      ! Flux limited by ocean and atm Hg(0)
                      !--------------------------------------------------------

                      ! Cap the flux w/ the available Hg(0) ocean mass
                      ! THIS IS THE PROBLEM!
                      Hg0aqtemp = CHg0aq * A_M2 * FRAC_O *1.0e-8_fp

                      IF ( FLUX(I,J,NN) * DTSRCE > Hg0aqtemp ) THEN
                          FLUX(I,J,NN) = Hg0aqtemp / DTSRCE
                          FUP(I,J,NN)  = FLUX(I,J,NN)-FDOWN(I,J,NN)
                      ENDIF

                      ! Cap the neg flux w/ the available Hg(0) atm mass
                      IF ( (-FLUX(I,J,NN) * DTSRCE ) > MHg0_air ) THEN
                          FLUX(I,J,NN) = -MHg0_air / DTSRCE
                      ENDIF

                      ! make sure Fup and Fdown do not underflow either
                      ! debug 2x2.5 diagnostic?
                      FUP(I,J,NN) = MAX(FUP(I,J,NN), SMALLNUM )
                      FDOWN(I,J,NN) = MAX(FDOWN(I,J,NN), SMALLNUM )

#ifdef BPCH_DIAG
                      !--------------------------------------------------------
                      ! %%%%% ND03 (bpch) DIAGNOSTICS %%%%%
                      !
                      ! Fluxes of Hg0 from air to ocean and ocean to air [kg]
                      !--------------------------------------------------------
             IF ( ND03 > 0 ) THEN
                AD03(I,J,16,NN) = AD03(I,J,16,NN) + FUP(I,J,NN) * DTSRCE
                AD03(I,J,17,NN) = AD03(I,J,17,NN) + FDOWN(I,J,NN) *DTSRCE
             ENDIF
#endif

                      !--------------------------------------------------------
                      ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
                      !
                      ! Fluxes of Hg0 from air to ocean and ocean to air [kg/s]
                      ! NOTE: Implement for total Hg species at ths time
                      !--------------------------------------------------------
                      IF ( NN == 1 ) THEN

                          ! Flux of Hg0 from ocean to air [kg/s]
                          IF ( State_Diag%Archive_FluxHg0fromOceanToAir ) THEN
                              State_Diag%FluxHg0fromOceanToAir(I,J) = FUP(I,J,NN)
                          ENDIF

                          IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
                              State_Diag%FluxHg0fromAirToOcean(I,J) = FDOWN(I,J,NN)
                          ENDIF

                      ENDIF

                  ENDDO

              ELSE

                  DO NN = 1, N_Hg_CATS
                      FLUX(I,J,NN)  = 0e0_fp
                      FUP(I,J,NN)   = 0e0_fp
                      FDOWN(I,J,NN) = 0e0_fp
                  ENDDO

              ENDIF
          ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

  END SUBROUTINE OFFLINEOCEAN_READMO
  !EOC
  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Transport Model                  !
  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: Reset_Hg_Diags
  !
  ! !DESCRIPTION: Zeroes the relevant diagnostic fields of State_Diag for
  !  the Hg and tagged Hg simulations.  Some of these need to be done
  !  for example, at the start of each timestep.
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE Reset_Hg_Diags( Input_Opt, State_Diag, RC )
      !
      ! !USES:
      !
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Diag_Mod, ONLY : DgnState
      USE Time_Mod,       ONLY : Its_Time_For_Chem
      !
      ! !INPUT PARAMETERS:
      !
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      !
      ! !INPUT/OUTPUT PARAMETERS:
      !
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
      !
      ! !REMARKS:
      !  NOTE: Not all netCDF diagnostic fields need to be zeroed.  Some fields
      !  of State_Diag are always
      !
      ! !REVISION HISTORY:
      !  06 Jan 2015 - R. Yantosca - Initial version
      !  See https://github.com/geoschem/geos-chem for complete history
      !EOP
      !------------------------------------------------------------------------------
      !BOC
      !
      ! !LOCAL VARIABLES:
      !
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Assume success
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at Reset_Hg_Diags (in module GeosCore/mercury_mod.F90'

      !=================================================================
      ! Reset_Hg_Diags begins here!
      !=================================================================

      ! Exit if it's not a Hg sim
      IF ( .not. Input_Opt%ITS_A_MERCURY_SIM ) THEN
          ErrMsg = 'Routine "Reset_Hg_Diags" was called, but this ' // &
              'routine is only valid for Hg or tagHg simulations!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
      ENDIF

      !--------------------------------------------------------------
      ! Zero diagnostics in depo_mercury_mod.F90
      !--------------------------------------------------------------
      IF ( State_Diag%Archive_FluxHg2HgPfromAirToSnow ) THEN
          State_Diag%FluxHg2HgPfromAirToSnow = 0.0_fp
      ENDIF

      ! These are called once per chemistry or emissions timestep
      IF ( Its_Time_For_Chem() ) THEN

          !--------------------------------------------------------------
          ! Zero diagnostics in mercury_mod.F90
          !--------------------------------------------------------------
          IF ( State_Diag%Archive_DryDepChm .or. State_Diag%Archive_DryDep ) THEN
              State_Diag%DryDepChm = 0.0_f4
          ENDIF

          !-------------------------------------------------------------
          ! Zero diagnostics in ocean_mercury_mod.F90
          !-------------------------------------------------------------
          IF ( State_Diag%Archive_EmisHg2rivers ) THEN
              State_Diag%EmisHg2rivers = 0.0_f4
          ENDIF

          IF ( State_Diag%Archive_EmisHg2snowToOcean ) THEN
              State_Diag%EmisHg2snowToOcean = 0.0_f4
          ENDIF

          IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
              State_Diag%FluxHg0fromAirToOcean = 0.0_f4
          ENDIF

          IF ( State_Diag%Archive_FluxHg0fromOceantoAir ) THEN
              State_Diag%FluxHg0fromOceanToAir = 0.0_f4
          ENDIF

          IF ( State_Diag%Archive_FluxHg2HgPfromAirToOcean ) THEN
              State_Diag%FluxHg2HgPfromAirToOcean = 0.0_f4
          ENDIF

          IF ( State_Diag%Archive_FluxOCtoDeepOcean ) THEN
              State_Diag%FluxOCtoDeepOcean = 0.0_f4
          ENDIF

      ENDIF

  END SUBROUTINE Reset_Hg_Diags
  !EOC
  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Transport Model                  !
  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: init_mercury
  !
  ! !DESCRIPTION: Subroutine INIT\_MERCURY allocates and zeroes all
  !  module arrays.
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE INIT_MERCURY( Input_Opt, State_Chm, State_Grid, RC )
      !
      ! !USES:
      !
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Grid_Mod,     ONLY : GrdState
      USE Gckpp_Monitor,      ONLY : Eqn_Names, Fam_Names
      USE Gckpp_Parameters,   ONLY : nFam, nReact
      USE FAST_JX_MOD,        ONLY : INIT_FJX

      !
      ! !INPUT PARAMETERS:
      !
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
      TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
      !
      ! !OUTPUT PARAMETERS:
      !
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
      !
      ! !REVISION HISTORY:
      !  02 Dec 2004 - N. (Eckley) Selin - Initial version
      !  See https://github.com/geoschem/geos-chem for complete history
      !EOP
      !------------------------------------------------------------------------------
      !BOC
      !
      ! !LOCAL VARIABLES:
      !
      ! Scalars
      LOGICAL, SAVE          :: IS_INIT = .FALSE.
      LOGICAL                :: LSPLIT
      LOGICAL                :: LDRYD
      LOGICAL                :: LNLPBL
      LOGICAL                :: LGTMM
      LOGICAL                :: LHALOGENCHEM
      INTEGER                :: nAdvect, nSpecies
      INTEGER                :: AS,    N
      INTEGER                :: KppId, I

      ! Strings
      CHARACTER(LEN=255)     :: ThisLoc
      CHARACTER(LEN=512)     :: ErrMsg

      ! Pointers
      TYPE(Species), POINTER :: SpcInfo

      !=================================================================
      ! INIT_MERCURY begins here!
      !=================================================================

      ! Assume success
      RC =  GC_SUCCESS

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Initialize
      SpcInfo  => NULL()
      LDRYD    = Input_Opt%LDRYD            ! Use drydep?
      LGTMM    = Input_Opt%LGTMM            ! Use GTMM model?
      LNLPBL   = Input_Opt%LNLPBL           ! Use non-local PBL?
      LSPLIT   = Input_Opt%LSPLIT           ! Tagged simulation?
      nAdvect  = State_Chm%nAdvect          ! # of Hg advected species
      nSpecies = State_Chm%nSpecies         ! # of species

      ! Location string for error messages
      ErrMsg   = ''
      ThisLoc  = '-> at DEFINE_TAGGED_Hg (in GeosCore/mercury_mod.F90)'

      ! Store the # of tagged Hg categories in a module variable
      N_Hg_CATS   = State_Chm%N_Hg_CATS

      ! Hg species index corresponding to a given Hg category number
      Hg0_Id_List => State_Chm%Hg0_Id_List
      Hg2_Id_List => State_Chm%Hg2_Id_List
      HgP_Id_List => State_Chm%HgP_Id_List

      ! Locate various species flags
      DO N = 1, State_Chm%nSpecies

          ! Point to the species database entry for species # N
          SpcInfo => State_Chm%SpcData(N)%Info

          SELECT CASE( TRIM( SpcInfo%Name ) )
              CASE( 'Hg0'     )
                  id_Hg0    = SpcInfo%ModelId
                  ID_Hg_tot = SpcInfo%Hg_Cat
              CASE( 'Hg2'     )
                  id_Hg2    = SpcInfo%ModelId
              CASE( 'HgP'     )
                  id_HgP    = SpcInfo%ModelId
              CASE( 'Hg0_can' )
                  ID_Hg_can = SpcInfo%Hg_Cat
              CASE( 'Hg0_usa' )
                  ID_Hg_usa = SpcInfo%Hg_Cat
              CASE( 'Hg0_cam' )
                  ID_Hg_cam = SpcInfo%Hg_Cat
              CASE( 'Hg0_sam' )
                  ID_Hg_sam = SpcInfo%Hg_Cat
              CASE( 'Hg0_waf' )
                  ID_Hg_waf = SpcInfo%Hg_Cat
              CASE( 'Hg0_eaf' )
                  ID_Hg_eaf = SpcInfo%Hg_Cat
              CASE( 'Hg0_saf' )
                  ID_Hg_saf = SpcInfo%Hg_Cat
              CASE( 'Hg0_naf' )
                  ID_Hg_naf = SpcInfo%Hg_Cat
              CASE( 'Hg0_eur' )
                  ID_Hg_eur = SpcInfo%Hg_Cat
              CASE( 'Hg0_eeu' )
                  ID_Hg_eeu = SpcInfo%Hg_Cat
              CASE( 'Hg0_sov' )
                  ID_Hg_sov = SpcInfo%Hg_Cat
              CASE( 'Hg0_mde' )
                  ID_Hg_mde = SpcInfo%Hg_Cat
              CASE( 'Hg0_sas' )
                  ID_Hg_sas = SpcInfo%Hg_Cat
              CASE( 'Hg0_eas' )
                  ID_Hg_eas = SpcInfo%Hg_Cat
              CASE( 'Hg0_sea' )
                  ID_Hg_sea = SpcInfo%Hg_Cat
              CASE( 'Hg0_jpn' )
                  ID_Hg_jpn = SpcInfo%Hg_Cat
              CASE( 'Hg0_oce' )
                  ID_Hg_oce = SpcInfo%Hg_Cat
              CASE( 'Hg0_so'  )
                  ID_Hg_so  = SpcInfo%Hg_Cat
              CASE( 'Hg0_bb'  )
                  ID_Hg_bb  = SpcInfo%Hg_Cat
              CASE( 'Hg0_geo' )
                  ID_Hg_geo = SpcInfo%Hg_Cat
              CASE( 'Hg0_atl' )
                  ID_Hg_atl = SpcInfo%Hg_Cat
              CASE( 'Hg0_nat' )
                  ID_Hg_nat = SpcInfo%Hg_Cat
              CASE( 'Hg0_sat' )
                  ID_Hg_sat = SpcInfo%Hg_Cat
              CASE( 'Hg0_npa' )
                  ID_Hg_npa = SpcInfo%Hg_Cat
              CASE( 'Hg0_arc' )
                  ID_Hg_arc = SpcInfo%Hg_Cat
              CASE( 'Hg0_ant' )
                  ID_Hg_ant = SpcInfo%Hg_Cat
              CASE( 'Hg0_ocn' )
                  ID_Hg_ocn = SpcInfo%Hg_Cat
              CASE( 'Hg0_str' )
                  ID_Hg_str = SpcInfo%Hg_cat
              CASE DEFAULT
             ! Do nothing
          END SELECT

          ! Free pointer
          SpcInfo => NULL()
      ENDDO

      WRITE( 6 ,'(a)' ) ' KPP Reaction Reference '
      DO N = 1, NREACT
          WRITE( 6, '(i8,a3,a85)' ) N,' | ',EQN_NAMES(N)
      END DO

      !--------------------------------------------------------------------
      ! Pre-store the KPP indices for each KPP prod/loss species or family
      !--------------------------------------------------------------------

      IF ( nFam > 0 ) THEN

          ! Allocate mapping array for KPP Ids for ND65 bpch diagnostic
          ALLOCATE( PL_Kpp_Id( nFam ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:PL_Kpp_Id', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN

          ! Loop over all KPP prod/loss species
          DO N = 1, nFam
              ! NOTE: KppId is the KPP ID # for each of the prod and loss
              ! diagnostic species.  This is the value used to index the
              ! KPP "VAR" array (in module gckpp_Global.F90).
              KppID = Ind_( TRIM ( Fam_Names(N) ), 'K' )

              ! If the species ID is OK, save in ND65_Kpp_Id
              PL_Kpp_Id(N) = KppId
          ENDDO

      ENDIF

      ! Set oxidant species IDs
      id_O3     = Ind_( 'O3'   )
      id_CO     = Ind_( 'CO'   )
      id_OH     = Ind_( 'OH'   )
      id_HO2    = Ind_( 'HO2'  )
      id_ClO    = Ind_( 'ClO'  )
      id_Cl     = Ind_( 'Cl'   )
      id_HOCl   = Ind_( 'HOCl' )
      id_NO2    = Ind_( 'NO2'  )
      id_NO     = Ind_( 'NO'   )
      id_Br     = Ind_( 'Br'   )
      id_BrO    = Ind_( 'BrO'  )
      id_CH4    = Ind_( 'CH4'  )

      ! Locate Hg(II) gas species
      id_HGBRNO2  = Ind_( 'HGBRNO2' )
      id_HGBRHO2  = Ind_( 'HGBRHO2' )
      id_HGBROH   = Ind_( 'HGBROH ' )
      id_HGBRBRO  = Ind_( 'HGBRBRO' )
      id_HGBRCLO  = Ind_( 'HGBRCLO' )
      id_HGBR2    = Ind_( 'HGBR2  ' )
      id_HGCLNO2  = Ind_( 'HGCLNO2' )
      id_HGCLHO2  = Ind_( 'HGCLHO2' )
      id_HGCLOH   = Ind_( 'HGCLOH ' )
      id_HGCLBRO  = Ind_( 'HGCLBRO' )
      id_HGCLCLO  = Ind_( 'HGCLCLO' )
      id_HGCLBR   = Ind_( 'HGCLBR'  )
      id_HGCL2    = Ind_( 'HGCL2'   )
      id_HgIIads  = Ind_( 'HGIIADS' )

      ! Initialize variables
      nHg2gasSpc = 0
      Map_Hg2gas = 0

      IF (  id_HGBRNO2 > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGBRNO2
      ENDIF
      IF (  id_HGBRHO2 > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGBRHO2
      ENDIF
      IF (  id_HGBROH  > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGBROH
      ENDIF
      IF (  id_HGBRBRO > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGBRBRO
      ENDIF
      IF (  id_HGBRCLO > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGBRCLO
      ENDIF
      IF (  id_HGBR2 > 0   ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGBR2
      ENDIF
      IF (  id_HGCLNO2 > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGCLNO2
      ENDIF
      IF (  id_HGCLHO2 > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGCLHO2
      ENDIF
      IF (  id_HGCLOH  > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HgCLOH
      ENDIF
      IF (  id_HGCLBRO > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGCLBRO
      ENDIF
      IF (  id_HGCLCLO > 0 ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGCLCLO
      ENDIF
      IF (  id_HGCLBR > 0  ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGCLBR
      ENDIF
      IF (  id_HGCL2  > 0  ) THEN
        nHg2gasSpc           = nHg2gasSpc + 1
        Map_Hg2gas(nHg2gasSpc) = id_HGCL2
      ENDIF

      !=================================================================
      ! Allocate arrays
      !=================================================================
      WRITE( 6, 101 )
101   FORMAT( '     - INIT_MERCURY: Allocating arrays for CHEMISTRY' )

      ALLOCATE( COSZM( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:COSZM', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      COSZM = 0e+0_fp

      ALLOCATE( EHg0_an( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHg0_an', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHg0_an = 0e+0_fp

      ALLOCATE( EHg0_am( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHg0_am', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHg0_am = 0e+0_fp

      ALLOCATE( EHg2_an( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHg2_an', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHg2_an = 0e+0_fp

      ALLOCATE( EHgP_an( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHgP_an', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHgP_an = 0e+0_fp

      ALLOCATE( EHg0_oc( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHg0_oc', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHg0_oc = 0e+0_fp

      ALLOCATE( EHg0_dist( State_Grid%NX, State_Grid%NY), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHg0_dist', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHg0_dist = 0e+0_fp

      ALLOCATE( EHg0_geo( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHg0_Geo', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHg0_geo = 0e+0_fp

      ALLOCATE( EHg0_bb( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHg0_bb', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHg0_bb = 0e+0_fp

      ALLOCATE( EHg0_snow( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:EHg0_snow', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      EHg0_snow = 0e+0_fp

      ! Allocate the following if tagged Hg simulation
      IF ( LSPLIT ) THEN

          ! Tagged Hg0 arrays (eds 8/31/10)
          ALLOCATE( EHg0_can( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_can', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_can = 0e+0_fp

          ALLOCATE( EHg0_usa( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_usa', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_usa = 0e+0_fp

          ALLOCATE( EHg0_cam( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_cam', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_cam = 0e+0_fp

          ALLOCATE( EHg0_sam( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_sam', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_sam = 0e+0_fp

          ALLOCATE( EHg0_waf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_waf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_waf = 0e+0_fp

          ALLOCATE( EHg0_eaf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_eaf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_eaf = 0e+0_fp

          ALLOCATE( EHg0_saf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_saf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_saf = 0e+0_fp

          ALLOCATE( EHg0_naf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_naf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_naf = 0e+0_fp

          ALLOCATE( EHg0_eur( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_eur', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_eur = 0e+0_fp

          ALLOCATE( EHg0_eeu( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_eeu', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_eeu = 0e+0_fp

          ALLOCATE( EHg0_mde( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_mde', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_mde = 0e+0_fp

          ALLOCATE( EHg0_sov( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_sov', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_sov = 0e+0_fp

          ALLOCATE( EHg0_sas( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_sas', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_sas = 0e+0_fp

          ALLOCATE( EHg0_eas( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_eas', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_eas = 0e+0_fp

          ALLOCATE( EHg0_sea( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_sea', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_sea = 0e+0_fp

          ALLOCATE( EHg0_jpn( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_jpn', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_jpn = 0e+0_fp

          ALLOCATE( EHg0_oce( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_oce', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_oce = 0e+0_fp

          ! Tagged Hg2 arrays (eds 8/31/10)
          ALLOCATE( EHg2_can( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_can', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_can = 0e+0_fp

          ALLOCATE( EHg2_usa( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_usa', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_usa = 0e+0_fp

          ALLOCATE( EHg2_cam( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_cam', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_cam = 0e+0_fp

          ALLOCATE( EHg2_sam( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_sam', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_sam = 0e+0_fp

          ALLOCATE( EHg2_waf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_waf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_waf = 0e+0_fp

          ALLOCATE( EHg2_eaf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_eaf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_eaf = 0e+0_fp

          ALLOCATE( EHg2_saf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_saf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_saf = 0e+0_fp

          ALLOCATE( EHg2_naf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_naf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_naf = 0e+0_fp

          ALLOCATE( EHg2_eur( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_eur', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_eur = 0e+0_fp

          ALLOCATE( EHg2_eeu( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_eeu', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_eeu = 0e+0_fp

          ALLOCATE( EHg2_mde( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_mde', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_mde = 0e+0_fp

          ALLOCATE( EHg2_sov( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_sov', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_sov = 0e+0_fp

          ALLOCATE( EHg2_sas( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_sas', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_sas = 0e+0_fp

          ALLOCATE( EHg2_eas( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:Hg2_eas', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_eas = 0e+0_fp

          ALLOCATE( EHg2_sea( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_sea', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_sea = 0e+0_fp

          ALLOCATE( EHg2_jpn( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_jpn', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_jpn = 0e+0_fp

          ALLOCATE( EHg2_oce( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg2_oce', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg2_oce = 0e+0_fp

          ! Tagged HgP arrays (eds 8/31/10)
          ALLOCATE( EHgP_can( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_can', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_can = 0e+0_fp

          ALLOCATE( EHgP_usa( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_usa', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_usa = 0e+0_fp

          ALLOCATE( EHgP_cam( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_cam', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_cam = 0e+0_fp

          ALLOCATE( EHgP_sam( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_sam', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_sam = 0e+0_fp

          ALLOCATE( EHgP_waf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_waf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_waf = 0e+0_fp

          ALLOCATE( EHgP_eaf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_eaf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_eaf = 0e+0_fp

          ALLOCATE( EHgP_saf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_saf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_saf = 0e+0_fp

          ALLOCATE( EHgP_naf( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_naf', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_naf = 0e+0_fp

          ALLOCATE( EHgP_eur( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_eur', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_eur = 0e+0_fp

          ALLOCATE( EHgP_eeu( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_eeu', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_eeu = 0e+0_fp

          ALLOCATE( EHgP_mde( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_mde', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_mde = 0e+0_fp

          ALLOCATE( EHgP_sov( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_sov', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_sov = 0e+0_fp

          ALLOCATE( EHgP_sas( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_sas', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_sas = 0e+0_fp

          ALLOCATE( EHgP_eas( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_eas', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_eas = 0e+0_fp

          ALLOCATE( EHgP_sea( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_sea', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_sea = 0e+0_fp

          ALLOCATE( EHgP_jpn( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_jpn', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_jpn = 0e+0_fp

          ALLOCATE( EHgP_oce( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHgP_oce', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHgP_oce = 0e+0_fp

      ENDIF

      IF ( LGTMM ) THEN

          ALLOCATE( EHg0_gtm( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_gtm', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_gtm = 0e+0_fp

      ELSE

          ALLOCATE( EHg0_ln( State_Grid%NX, State_Grid%NY, N_Hg_CATS ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_ln', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_ln = 0e+0_fp

          ALLOCATE( EHg0_vg( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_vg', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_vg = 0e+0_fp

          ALLOCATE( EHg0_so( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:EHg0_so', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          EHg0_so = 0e+0_fp

      ENDIF

      ALLOCATE( TCOSZ( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:TCOSZ', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      TCOSZ = 0e+0_fp

      ALLOCATE( TTDAY( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:TTDAY', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      TTDAY = 0e+0_fp

      ! Allocate ZERO_DVEL if we use non-local PBL mixing or
      ! if drydep is turned off
      IF ( LNLPBL .OR. (.not. LDRYD) ) THEN
          ALLOCATE( ZERO_DVEL( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:ZERO_DVEL', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          ZERO_DVEL = 0e+0_fp
      ENDIF

      ALLOCATE( HG2_SEASALT_LOSSRATE( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:HG2_SEASALT_LOSSRATE', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      HG2_SEASALT_LOSSRATE = 0e+0_fp

      ALLOCATE( JNO2_Inst( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:JNO2_Inst', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      JNO2_Inst = 0e+0_fp

      ALLOCATE( JBrO_Inst( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:JBrO_Inst', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      JBrO_Inst = 0e+0_fp

      ALLOCATE( JClO_Inst( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:JClO_Inst', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      JClO_Inst = 0e+0_fp

      !=================================================================
      ! Allocate & initialize arrays for tagged species indices
      !=================================================================
      IF ( LSPLIT ) THEN

          ! Species indices for tagged anthro regions
          ALLOCATE( AN_Hg0( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:AN_Hg0', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          AN_Hg0 = 0e+0_fp

          ALLOCATE( AN_Hg2( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:AN_Hg2', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          AN_Hg2 = 0e+0_fp

          ALLOCATE( AN_HgP( State_Grid%NX, State_Grid%NY ), STAT=RC )
          CALL GC_CheckVar( 'mercury_mod.F90:AN_HgP', 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          AN_HgP = 0e+0_fp

      ENDIF

      ! HG_EMIS is needed for non-local PBL mixing
      ALLOCATE ( HG_EMIS( State_Grid%NX, State_Grid%NY, nAdvect ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:HG_EMIS', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      HG_EMIS = 0e+0_fp



      !=================================================================
      ! Allocate and initialize oxidant concentration pointer
      !=================================================================
      ! Number of species in Hg oxidants list
      N_Hg_Oxid = 11   ! Br, BrO, Cl, ClO, OH, HO2, NO, NO2, O3, CO, & CH4


      ! Allocate and initialize oxidant species map to State_Chem
      ALLOCATE( GC_HgOxid_TrID ( N_Hg_oxid ), STAT=AS )
      CALL GC_CheckVar( 'mercury_mod.F90:GC_HgOxid_TrID', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      GC_HgOxid_TrID(1:N_Hg_Oxid) = (/id_Br, id_BrO, id_Cl, id_ClO, &
          id_OH, id_HO2, id_NO, id_NO2, &
          id_O3, id_CO,  id_CH4/)


      ALLOCATE( HgOxidPtr( N_Hg_Oxid ), STAT=AS )
      CALL GC_CheckVar( 'mercury_mod.F90:HgOxidPtr', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ! Number of aerosol & dust speices
      N_Aer  =  7
      N_Dust =  7

      ! Aerosol species name
      ALLOCATE( AerSpcNames ( N_Dust + N_Aer ), STAT=AS )
      CALL GC_CheckVar( 'mercury_mod.F90:AerSpcNames', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      AerSpcNames = (/'DST1','DST2','DST3','DST4','DST5','DST6',  'DST7', &
          'SO4', 'BC',  'OC','SSA', 'SSC', 'BGSULF','ICEI'/)

      ALLOCATE( AeroPtr( N_Dust + N_Aer ), STAT=AS )
      CALL GC_CheckVar( 'mercury_mod.F90:AeroPtr', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ALLOCATE( HEM_ODDust( State_Grid%NX, State_Grid%NY, State_Grid%NZ, N_Dust ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:HEM_ODDust', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      HEM_ODDust = 0e+0_fp

      ALLOCATE( HEM_ODAer( State_Grid%NX, State_Grid%NY, State_Grid%NZ, N_Aer ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:HEM_ODAer', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      HEM_ODAer = 0e+0_fp

      ALLOCATE( AerHetRates( State_Grid%NX, State_Grid%NY, State_Grid%NZ, nSpecies ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:AerHetRates', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      AerHetRates = 0e+0_fp

      ALLOCATE( OrgHetRates( State_Grid%NX, State_Grid%NY, State_Grid%NZ, nSpecies ), STAT=RC )
      CALL GC_CheckVar( 'mercury_mod.F90:OrgHetRates', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      OrgHetRates = 0e+0_fp

      !=================================================================
      ! Settings
      !=================================================================

      ! Switch uses ocean rate coefficients from parameter inversion,
      ! ref. Song et al. 2015 ACP
      LOCEANCOEF=.FALSE.

      ! Switch determines whether uptake of Hg2 by sea-salt aerosol
      ! is calculated dynamically (TRUE) or uses a constant rate (FALSE)
      LDYNSEASALT = .TRUE.

      ! no Hg emitted through transpiration (VEGEMIS off)
      LVEGEMIS=.FALSE.

      ! Switch adds bromine explosion in Northern springtime
      LPOLARBR=.TRUE.

      ! Switch turns on snowpack Hg storage until snowmelt
      LHGSNOW = .TRUE.

      ! Multiplicative factor for increasing stratospheric Br and BrO
      STRAT_BR_FACTOR = 1e+0_fp

      ! Switch turns off all emissions except direct anthropogenic
      LAnthroHgOnly = .FALSE.

      ! Switch turns off all anthropogenic emissions from contiguous USA
      LnoUSAemis = .FALSE.

      !=================================================================
      ! Done
      !=================================================================

      ! Reset IS_INIT, since we have already allocated arrays
      IS_INIT = .TRUE.

  END SUBROUTINE INIT_MERCURY
  !EOC

  !------------------------------------------------------------------------------
  !                  GEOS-Chem Global Chemical Transport Model                  !
  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: cleanup_mercury
  !
  ! !DESCRIPTION: Subroutine CLEANUP\_MERCURY deallocates all module arrays.
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE CLEANUP_MERCURY
      !
      ! !REVISION HISTORY:
      !  06 Dec 2004 - N. (Eckley) Selin - Initial version
      !  See https://github.com/geoschem/geos-chem for complete history
      !EOP
      !------------------------------------------------------------------------------
      !BOC
      !
      ! !LOCAL VARIABLES:
      !
      INTEGER                :: I

      IF ( ALLOCATED( AN_Hg0   ) ) DEALLOCATE( AN_Hg0   )
      IF ( ALLOCATED( AN_Hg2   ) ) DEALLOCATE( AN_Hg2   )
      IF ( ALLOCATED( AN_HgP   ) ) DEALLOCATE( AN_HgP   )
      IF ( ALLOCATED( COSZM    ) ) DEALLOCATE( COSZM    )
      IF ( ALLOCATED( EHg0_an  ) ) DEALLOCATE( EHg0_an  )
      IF ( ALLOCATED( EHg0_can ) ) DEALLOCATE( EHg0_can )
      IF ( ALLOCATED( EHg0_usa ) ) DEALLOCATE( EHg0_usa )
      IF ( ALLOCATED( EHg0_sam ) ) DEALLOCATE( EHg0_sam )
      IF ( ALLOCATED( EHg0_eaf ) ) DEALLOCATE( EHg0_eaf )
      IF ( ALLOCATED( EHg0_waf ) ) DEALLOCATE( EHg0_waf )
      IF ( ALLOCATED( EHg0_saf ) ) DEALLOCATE( EHg0_saf )
      IF ( ALLOCATED( EHg0_naf ) ) DEALLOCATE( EHg0_naf )
      IF ( ALLOCATED( EHg0_eur ) ) DEALLOCATE( EHg0_eur )
      IF ( ALLOCATED( EHg0_eeu ) ) DEALLOCATE( EHg0_eeu )
      IF ( ALLOCATED( EHg0_sov ) ) DEALLOCATE( EHg0_sov )
      IF ( ALLOCATED( EHg0_mde ) ) DEALLOCATE( EHg0_mde )
      IF ( ALLOCATED( EHg0_sas ) ) DEALLOCATE( EHg0_sas )
      IF ( ALLOCATED( EHg0_eas ) ) DEALLOCATE( EHg0_eas )
      IF ( ALLOCATED( EHg0_sea ) ) DEALLOCATE( EHg0_sea )
      IF ( ALLOCATED( EHg0_jpn ) ) DEALLOCATE( EHg0_jpn )
      IF ( ALLOCATED( EHg0_oce ) ) DEALLOCATE( EHg0_oce )
      IF ( ALLOCATED( EHg0_am  ) ) DEALLOCATE( EHg0_am  )
      IF ( ALLOCATED( EHg2_an  ) ) DEALLOCATE( EHg2_an  )
      IF ( ALLOCATED( EHg2_can ) ) DEALLOCATE( EHg2_can )
      IF ( ALLOCATED( EHg2_usa ) ) DEALLOCATE( EHg2_usa )
      IF ( ALLOCATED( EHg2_sam ) ) DEALLOCATE( EHg2_sam )
      IF ( ALLOCATED( EHg2_eaf ) ) DEALLOCATE( EHg2_eaf )
      IF ( ALLOCATED( EHg2_waf ) ) DEALLOCATE( EHg2_waf )
      IF ( ALLOCATED( EHg2_saf ) ) DEALLOCATE( EHg2_saf )
      IF ( ALLOCATED( EHg2_naf ) ) DEALLOCATE( EHg2_naf )
      IF ( ALLOCATED( EHg2_eur ) ) DEALLOCATE( EHg2_eur )
      IF ( ALLOCATED( EHg2_eeu ) ) DEALLOCATE( EHg2_eeu )
      IF ( ALLOCATED( EHg2_sov ) ) DEALLOCATE( EHg2_sov )
      IF ( ALLOCATED( EHg2_mde ) ) DEALLOCATE( EHg2_mde )
      IF ( ALLOCATED( EHg2_sas ) ) DEALLOCATE( EHg2_sas )
      IF ( ALLOCATED( EHg2_eas ) ) DEALLOCATE( EHg2_eas )
      IF ( ALLOCATED( EHg2_sea ) ) DEALLOCATE( EHg2_sea )
      IF ( ALLOCATED( EHg2_jpn ) ) DEALLOCATE( EHg2_jpn )
      IF ( ALLOCATED( EHg2_oce ) ) DEALLOCATE( EHg2_oce )
      IF ( ALLOCATED( EHgP_an  ) ) DEALLOCATE( EHgP_an  )
      IF ( ALLOCATED( EHgP_can ) ) DEALLOCATE( EHgP_can )
      IF ( ALLOCATED( EHgP_usa ) ) DEALLOCATE( EHgP_usa )
      IF ( ALLOCATED( EHgP_sam ) ) DEALLOCATE( EHgP_sam )
      IF ( ALLOCATED( EHgP_eaf ) ) DEALLOCATE( EHgP_eaf )
      IF ( ALLOCATED( EHgP_waf ) ) DEALLOCATE( EHgP_waf )
      IF ( ALLOCATED( EHgP_saf ) ) DEALLOCATE( EHgP_saf )
      IF ( ALLOCATED( EHgP_naf ) ) DEALLOCATE( EHgP_naf )
      IF ( ALLOCATED( EHgP_eur ) ) DEALLOCATE( EHgP_eur )
      IF ( ALLOCATED( EHgP_eeu ) ) DEALLOCATE( EHgP_eeu )
      IF ( ALLOCATED( EHgP_sov ) ) DEALLOCATE( EHgP_sov )
      IF ( ALLOCATED( EHgP_mde ) ) DEALLOCATE( EHgP_mde )
      IF ( ALLOCATED( EHgP_sas ) ) DEALLOCATE( EHgP_sas )
      IF ( ALLOCATED( EHgP_eas ) ) DEALLOCATE( EHgP_eas )
      IF ( ALLOCATED( EHgP_sea ) ) DEALLOCATE( EHgP_sea )
      IF ( ALLOCATED( EHgP_jpn ) ) DEALLOCATE( EHgP_jpn )
      IF ( ALLOCATED( EHgP_oce ) ) DEALLOCATE( EHgP_oce )
      IF ( ALLOCATED( EHg0_oc  ) ) DEALLOCATE( EHg0_oc  )
      IF ( ALLOCATED( EHg0_ln  ) ) DEALLOCATE( EHg0_ln  )
      IF ( ALLOCATED( EHg0_snow) ) DEALLOCATE( EHg0_snow)
      IF ( ALLOCATED( EHg0_geo ) ) DEALLOCATE( EHg0_geo )
      IF ( ALLOCATED( EHg0_bb  ) ) DEALLOCATE( EHg0_bb  )
      IF ( ALLOCATED( EHg0_gtm ) ) DEALLOCATE( EHg0_gtm )
      IF ( ALLOCATED( EHg0_vg  ) ) DEALLOCATE( EHg0_vg  )
      IF ( ALLOCATED( EHg0_so  ) ) DEALLOCATE( EHg0_so  )
      IF ( ALLOCATED( EHg0_dist) ) DEALLOCATE( EHg0_dist)
      IF ( ALLOCATED( TCOSZ    ) ) DEALLOCATE( TCOSZ    )
      IF ( ALLOCATED( TTDAY    ) ) DEALLOCATE( TTDAY    )
      IF ( ALLOCATED( ZERO_DVEL) ) DEALLOCATE( ZERO_DVEL)
      IF ( ALLOCATED( HG_EMIS  ) ) DEALLOCATE( HG_EMIS  )
      IF ( ALLOCATED( HG2_SEASALT_LOSSRATE ) ) DEALLOCATE( HG2_SEASALT_LOSSRATE )
      IF ( ALLOCATED( JNO2_INST ) ) DEALLOCATE( JNO2_INST )
      IF ( ALLOCATED( JBrO_INST ) ) DEALLOCATE( JBrO_INST )
      IF ( ALLOCATED( JClO_INST ) ) DEALLOCATE( JClO_INST )
      IF ( ALLOCATED( GC_HgOxid_TrID ) ) DEALLOCATE( GC_HgOxid_TrID )
      IF ( ALLOCATED( AerSpcNames ) ) DEALLOCATE( AerSpcNames )
      IF ( ALLOCATED( PL_Kpp_Id ) ) DEALLOCATE( PL_Kpp_Id  )
      IF ( ALLOCATED( HEM_ODAer ) ) DEALLOCATE( HEM_ODAer  )
      IF ( ALLOCATED( HEM_ODDust) ) DEALLOCATE( HEM_ODDust )
      IF ( ALLOCATED( AerHetRates ) ) DEALLOCATE( AerHetRates )
      IF ( ALLOCATED( OrgHetRates ) ) DEALLOCATE( OrgHetRates )

      ! Cleanup OxidPtr array of pointers
      IF ( ASSOCIATED( HgOxidPtr ) ) THEN
          DO I = 1, N_Hg_Oxid
              HgOxidPtr(I)%Data => NULL()
          ENDDO
          DEALLOCATE( HgOxidPtr )
      ENDIF

      ! Cleanup OxidPtr array of pointers
      IF ( ASSOCIATED( AeroPtr ) ) THEN
          DO I = 1, N_Aer + N_Dust
              AeroPtr(I)%AOD  => NULL()
              AeroPtr(I)%Area => NULL()
              AeroPtr(I)%Radi => NULL()
              AeroPtr(I)%Conc => NULL()
          ENDDO
          DEALLOCATE( AeroPtr )
      ENDIF

      ! Free Hg indexing pointers
      Hg0_Id_List => NULL()
      Hg2_Id_List => NULL()
      HgP_Id_List => NULL()

      ! Free other pointers
      OCEAN_CONC => NULL()

  END SUBROUTINE CLEANUP_MERCURY
  !EOC
  END MODULE MERCURY_MOD
