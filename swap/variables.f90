! File VersionID:
!   $Id: variables.f90 326 2017-04-22 10:52:23Z kroes006 $
! ----------------------------------------------------------------------
! File Name:    variables.for
! Content:      common variables of modular SWAP code
! Sections:     time & control, meteo, irrigation, crop, soilwater, macropore, surfacewater, heat, snow, solute
      module  variables
      implicit none
      save
      include 'arrays.fi'

! --- time & control variables
      integer   nprintday          ! Number of output times during one day
      logical   flprintdt          ! Flag indicating output every dt 
      logical   flprintshort       ! Flag indicating several output times during a day
      logical   floutputshort      ! Flag indicating time for output during a day is reached
      integer   nprintcount        ! Counter for output during a day
      integer   cntper             ! Day number of intermediate period
      integer   daycum             ! Day number from start of simulation
      integer   daynr              ! Day number of calendar year
      integer   imonth             ! Month number of calendar year
      integer   ioutdat            ! Counter of output date for water and solute balance
      integer   ioutdatint         ! Counter of intermediate output date
      integer   isteps             ! Number of time steps from the start of the day
      integer   iyear              ! Year number of calendar year
      integer   iyearm1            ! Year number of previous calendar year
      integer   logf               ! Internal number of logbook output file *.LOG
      integer   period             ! Length of prescribed output interval (T)
      integer   swheader           ! Switch for printing of header in output files at each balance period: 0 = no; 1 = yes
      integer   swodat             ! Switch for extra, specific output dates in the input file: 0 = no; 1 = yes
      integer   swres              ! Switch for counter of output interval: 0 = no reset; 1 = reset at start of calendar year
      integer   swscre             ! Switch of screen display: 0 = no display; 1 = summary water balance; 2 = daynumber
      real(8)   dt                 ! Time step (T)
      real(8)   dtmax              ! Maximum time step (T)
      real(8)   dtmin              ! Minimum time step (T)
      real(8)   outdat(maout)      ! Array with output dates for water and solute balances
      real(8)   outdatint(maout)   ! Array with intermediate output dates
      real(8)   outper             ! Length of actual output interval (T)
      real(8)   t                  ! Time since start of calendar year (T)
      real(8)   t1900              ! Time since 1900 (T)
      real(8)   tcum               ! Time since start of simulation (T)
      real(8)   tend               ! End date of simulation run
      real(8)   tstart             ! Start date of simulation run
      logical   flbaloutput        ! Flag indicating time for output of water and solute balance
      logical   fldayend           ! Flag indicating end of day
      logical   fldaystart         ! Flag indicating that this time step is the first one of a day
      logical   fldecdt            ! Flag indicating decrease of time step
      logical   fldecdtmin         ! Flag indicating that the time step should be reset to the minimum time step
      logical   fldtmin            ! Flag indicating that the time step is equal to the minimum time step
      logical   fldtreduce
      logical   flheader           ! Flag indicating that header should be printed in output file
      logical   floutput           ! Flag indicating time for ouput
      logical   flrunend           ! Flag indicating end of run
      logical   flSwapShared       ! Flag to indicate the shared simultaneous simulation with other applications
      logical   flyearstart        ! Flag indicating the beginning of a new year
      logical   flzerocumu         ! Flag indicating that cumulative fluxes should be reset to zero
      logical   flzerointr         ! Flag indicating that intermediate fluxes should be reset to zero
      character(len=11) date       ! Current date
      character(len=16) outfil     ! Name of output file
      character(len=80) pathwork   ! Path to work directory
      character(len=80) project    ! Name of project

! --- meteo variables
      integer   daymeteo           ! Calendar day number for which meteorological data should be read from input file
      integer   daynrfirst         ! First calendar day number for which meteorological data is available in current year
      integer   daynrlast          ! Last calendar day number for which meteorological data is available in current year
      integer   detrecord(nmetfile) ! Record number of meteo file with detailed meteo data (-)
      integer   irectotal          ! Total record number with detailed meteo input for new weather file (-)
      integer   nmetdetail         ! Number of detailed records for ET and rainfall per day (-)
      integer   nmrain             ! Number of rain event records (-)
      integer   nofd               ! number of days for running average Tmin (-)
      integer   rainrec            ! Actual rain record when detailed rainfall are used from separate rain file
      integer   swdivide           ! Switch on division ET into E and T: 0 = according to the SWAP traditional way; 1 = according to PMdirect
      integer   swetr              ! Switch: 0 = use daily meteorological basic data; 1 = use daily Etref values
      integer   swetsine           ! Switch: 0 = Tp and Ep uniform during a day; 1 = Tp and Ep are distributed as sine waves during a day
      integer   swinter            ! Switch for interception method: 0 = no interception; 1 = agricultural crops; 2 = trees and forests
      integer   swmetdetail        ! Switch: 0 = daily meteorological records; 1 = detailed records for both ET and rainfall
      integer   swmeteo            ! Switch: 1 = no detailed meteo data; 2 = detailed meteo data are required for crop growth simulation
      integer   swrain             ! Switch: 0 = use daily rain amounts; 1 = use daily amounts + mean intensity; 
                                   !         2 = use daily amounts + duration; 3 = use detailed rainfall data from separate file
      integer   wrecord            ! Actual number of weather record in case of detailed weather input (< 1 day)
      integer   yearmeteo          ! Year for which meteorological data should be read from input file
      real(8)   aetr(366)          ! Array with daily ETref input data (L/T)
      real(8)   ahum(366)          ! Array with daily humidity input data (M/L/T2)  
      real(8)   aintcdt            ! Interception flux of ONLY Rain during iteration timesteps (L/T) 
      real(8)   alt                ! Altitude of meteorological station (L)
      real(8)   altw               ! Height of wind speed measurement (L)
      real(8)   angstroma          ! first  angstrom coefficient [-]
      real(8)   angstromb          ! second angstrom coefficient [-]
      real(8)   arad(366)          ! Array with daily radiation input data (M/T2)
      real(8)   arai(366)          ! Array with daily precipitation sum input data (L/T)
      real(8)   atav(96)           ! In case of detailed weather input, air temperature of each weather record (L/T)
      real(8)   atmdem             ! Atmospheric demand = daily potential transpiration of a dry crop (L/T)
      real(8)   atmin7(7)          ! Array with minimum temperatures of last week (C)
      real(8)   atmn(366)          ! Array with daily minimum temperature input data (  )
      real(8)   atmx(366)          ! Array with daily maximum temperature input data (  )
      real(8)   awin(366)          ! Array with daily wind speed input data (L/T)     
      real(8)   caintc             ! Cumulative amount of rainfall interception (L)
      real(8)   cevap              ! Cumulative amount of actual soil evaporation (L)
      real(8)   cgrai              ! Cumulative amount of gross precipitation (L)
      real(8)   cnrai              ! Cumulative amount of net precipitation (L)
      real(8)   cpeva              ! Cumulative amount of potential soil evaporation (L)
      real(8)   cptra              ! Cumulative amount of potential transpiration (L)
      real(8)   daylp              ! Photoperiodic daylength in hours (T)
      real(8)   dethum(nmetfile)   ! Array with detailed humidity input data (M/L/T2) 
      real(8)   detrad(nmetfile)   ! Array with detailed radiation input data (M/T2
      real(8)   detrain(nmetfile)  ! Array with detailed precipitation sum input data (L/T)
      real(8)   dettav(nmetfile)   ! Array with detailed temperature input data (  )
      real(8)   dettime(nmetfile)  ! Array with dates of detailed meteo input data
      real(8)   detwind(nmetfile)  ! Array with detailed wind speed input data (L/T)   
      real(8)   dtEventRain        ! Time step length for next precipitation event (T)
      real(8)   empreva            ! Reduced soil evaporation flux according to empirical functions (L/T)
      real(8)   epot(96)           ! In case of detailed weather input, calculated Epot of each weather record (L/T)
      real(8)   cfevappond         ! Parameter equal to the ratio ponding layer evaporation / ETref (-)
      real(8)   finterception      ! Ratio net / gross rain flux in case of detailed rainfall data (-)
      real(8)   fprecnosnow        ! Ratio rain (excl. snow and rain on snow) / gross rain flux in case of detailed rainfall data (-)
      real(8)   grai               ! Daily gross rain flux (L/T), without rain on snow
      real(8)   graidt             ! Gross precipitation flux during iteration timesteps (L/T) 
      real(8)   grain(96)          ! In case of detailed weather input, gross rain flux of each weather record (L/T)
      real(8)   ievap              ! Intermediate amount of actual soil evaporation (L)
      real(8)   inrai              ! Intermediate amount of net precipitation (L)
      real(8)   ipeva              ! Intermediate amount of potential soil evaporation (L)
      real(8)   iptra              ! Intermediate amount of potential transpiration (L)
      real(8)   lat                ! Latitude of meteorological station (degrees)
      real(8)   metperiod          ! Length of weather record in case of detailed ET and rainfall input
      real(8)   nraida             ! Daily average net precipitation flux (L/T)
      real(8)   nraidt             ! Net precipitation flux during iteration timesteps (L/T)
      real(8)   nrain(96)          ! In case of detailed weather input, calculated netto rain of each weather record (L/T)
      real(8)   peva               ! Potential soil evaporation flux (L/T)
      real(8)   pevaday            ! Potential soil evaporation of one day (L)
      real(8)   ptra               ! Potential transpiration flux (L/T)
      real(8)   ptraday            ! Potential transpiration of one day (L)
      real(8)   rad                ! Global solar radiation (J/m2/d)
      real(8)   rainamount(mrain)  ! Array with short duration rainfall amounts (L)
      real(8)   raintab(60)        ! Array with mean rainfall intensity (L/T) as function of time (T)
      real(8)   rainfluxarray(mrain) ! Array with short duration rainfall intensities (L/T)
      real(8)   raintimearray(mrain) ! Array with times (T) at which rainfall intensity changes
      real(8)   rh                 ! Relative air humidity (-)
      real(8)   tav                ! Average air temperature of a day ( C)
      real(8)   tavd               ! Average air temperature during day time ( C)
      real(8)   timjan1            ! January first of meteo year as number of day since 1-1-1900(T) 
      real(8)   tmn                ! Minimum air temperature of current day ( C)
      real(8)   tmnr               ! Average of minimum air temperature during past 7 days ( C)
      real(8)   tmx                ! Maximum air temperature of current day ( C)
      real(8)   tpot(96)           ! In case of detailed weather input, calculated Tpot of each weather record (L/T)
      real(8)   tra                ! Actual transpiration flux (L/T)
      real(8)   wet(366)           ! Fraction of each day the crop is wet (L)    
      logical   fletsine           ! Flag indicating that Tp and Ep are distributed according to a sine wave 
      logical   flmeteodt          ! Flag indicating that calculations of meteo variables on time step basis is required    !!!!!!!!!!!  Robnew
      logical   flmetdetail        ! Flag indicating that detailed meteorological records for both ET and rainfall are used 
      logical   flrainintens       ! Flag indicating that rainfall intensity info in combination with daily ET records is used              !!!!!!!!!!!  Robnew
      logical   flupdmetdet        ! Flag indicating that update of detailed meteorological input within the day is required               !!!!!!!!!!!  Robnew
      character(len=200) metfil    ! Name of meteorological input file
      character(len=80) pathatm    ! Path to folder with meteorological input files
      character(len=200) rainfil   ! Name of input file with detailed rainfall intensities
!   - meteo output variables for PEARL
      real(4)   out_etr            ! Reference evapotranspiration  of current day (m/d)
      real(4)   out_hum            ! Air humidity  of current day (kPa)
      real(4)   out_rad            ! Global solar radiation (KJ/m2)
      real(4)   out_tmn            ! Minimum air temperature of current day ( C)
      real(4)   out_tmx            ! Maximum air temperature of current day ( C)
      real(4)   out_wet            ! Rainfall duration of current day (d)
      real(4)   out_win            ! Average wind speed of current day (m/s)

! --- irrigation variables
      integer   irg                ! Internal number of irrigation output file *.IRG
      integer   irrigevent         ! Switch: 0 = no irrigation; 1 = fixed irrigation event; 2 = scheduled irrigation event
      integer   irtype(mairg)      ! Type of fixed irrigation: 0 = sprinkling irrigation; 1 = surface irrigation
      integer   isua               ! Switch for type of irrigation: 0 = sprinkling irrigation, 1 = surface irrigation
      integer   isuas              ! Switch for type of scheduled irrigation: 0 = sprinkling irrigation, 1 = surface irrigation
      integer   nirri              ! Number of irrigation event
      integer   phormc             ! Swith for 5th irrigation criterion: 0 = use pressure head; 1 = use water content
      integer   schedule           ! Switch for simulation of irrigation scheduling: 0 = no, 1 = yes
      integer   swirfix            ! Switch for fixed irrigation: 0 = no applications prescribed; 1 = applications are prescribed
      integer   swcirrthres        ! Switch to allow over irrigation when a conc-threshold is exceeded: 0 = no; 1 = yes/allowed
      real(8)   cirrs              ! Solute concentration of irrigation water (M/L3)
      real(8)   cirrthres          ! Threshold value (M/L3) indicating the concentration that initiates over irrigation
      real(8)   dcrit              ! Depth (L) of sensor for soil water pressure head or water content
      real(8)   ditab(14)          ! Array with amount of under- or over-irrigation (L) as function of crop development stage
      real(8)   dwatab(14)         ! Array with maximum amounts of water depleted as function of crop development stage
      real(8)   fidtab(14)         ! Array with prescribed fixed irrigation depth (L) as function of crop development stage
      real(8)   gird               ! Gross irrigation depth (L)
      real(8)   hcritab(14)        ! Array with minimum soil water pressure heads (L) as function of crop development stage
      real(8)   igird              ! Intermediate depth of gross irrigation (L)
      real(8)   inird              ! Intermediate depth of net irrigation (L)
      real(8)   irconc(mairg)      ! Array with irrigation concentrations (M/L3) in case of fixed irrigation
      real(8)   irdate(mairg)      ! Array with fixed irrigation dates
      real(8)   irdepth(mairg)     ! Array with fixed irrigation depths (L)
      real(8)   nird               ! Net irrigation depth (L)
      real(8)   perirrsurp         ! percentage (-) of the scheduled irrigation depths that may be over irrigated
      real(8)   raithreshold       ! Threshold value (L) indicating the amount of rainfall which is substracted from scheduled irrigation depths
      real(8)   rawtab(14)         ! Array with minimum of readily available water as function of crop development stage
      real(8)   tawtab(14)         ! Array with minimum of totally available water as function of crop development stage
      real(8)   tcritab(14)        ! Array with minimum volumetric soil water contents as function of crop development stage
      real(8)   tstairrig          ! Date after which scheduled irrigation is allowed
      real(8)   tendirrig          ! Date after which scheduled irrigation is NOT allowed
      real(8)   treltab(14)        ! Array with minimum of ratio actual/potential transpiration as function of crop development stage
      logical   flheadirg          ! Flag indicating that header should be printed in irrigation output file
      logical   flirrigate         ! Flag indicating irrigation in simulation run (either fixed or scheduled)
      logical   flIrg1Start        ! Flag indicating irrigation output of 1st crop of the simulation period
      logical   FlIrrigationOutput ! Flag indication irrigation output

! --- crop variables
      integer   swCrop             ! Switch for simulating crop: 0 = No; 1 = Yes
      logical   flCropCalendar     ! Flag indicating that crop season is active (but currently might be bare or cropped)
      logical   flCropEmergence    ! Flag indicating period from crop emergence until harvest
      logical   flCropHarvest      ! Flag indicating period from crop harvest until the end of crop season
      logical   flCropReadFile     ! Flag indicating reading of input.crp
      logical   flCropOpenFile     ! Flag indicating to create output.crp
      logical   flcropoutput       ! Flag indicating data should be written to crop output file
      logical   fltsumttd          ! Flag indicating suppress grass growth until criteria for temperaature, time and depth are met

      integer   croptype(macrop)   ! Switch for type of crop model: 1 = simple; 2 = general detailed model; 3 = detailed for grass
      integer   crp                ! Internal number of crop output file *.CRP
      integer   daycrop            ! Number of days that a crop exists
      integer   icrop              ! Current crop number
      integer   daymow             ! Number of days that grass is growing after management event or emergence, actual run (d)
      integer   daymowpot          ! Number of days that grass is growing after management event or emergence, potential run (d)
      integer   idaysgraz          ! Day number of grass grazing, actual run (d)
      integer   idaysgrazpot       ! Day number of grass grazing, potential run (d)
      integer   idev               ! Switch for length of growth period in case of simple crop: 1 = fixed; 2 = depends on temperature sum
      integer   idregr             ! Number of days for regrowth of grassland, actual growth
      integer   idregrpot          ! Number of days for regrowth of grassland, potential growth
      integer   idsl               ! Switch for crop development before anthesis: 0 = depends on temperature; 
                                   !   1 = depends on temperature and day length; 2 = depends on temperature, day length and vernalisation factor
      integer   iharvest           ! Grass harvest number of grass crop when harvest dates are fixed
      integer   initcrp(macrop)    ! Switch for crop-specific type of initialization (1 = default=emergence; 2 = sowing); 
      integer   ilvold             ! Age of oldest leaf (d) of actual crop
      integer   ilvoldpot          ! Age of oldest leaf (d) of potential crop
      integer   iseqgm             ! Counter in sequence of grass grazing and mowing, actual crop
      integer   iseqgmpot          ! Counter in sequence of grass grazing and mowing, potential crop
      integer   noddrz             ! Compartment number at bottom root zone (-)
      integer   seqgrazmow(366)    ! Sequence of grass grazing and mowing, actual crop
      integer   seqgrazmowpot(366) ! Sequence of grass grazing and mowing, potential crop
      integer   swcf               ! Switch for simple crop: 1 = crop factor is input; 2 = crop height is input
      integer   swdmi2rd           ! rooting depth limitation by relative dry matter increase (dmi/dmipot)
      integer   swdrought          ! Switch for drought stress (1 = Feddes et al., 1978; 2 = De Jong van Lier et al., 2008)
      integer   swgc               ! Switch for simple crop: 1 = leaf area index is input; 2 = soil cover fraction is input
      integer   swrootradius       ! Switch for root radius to calculate Oxygen stress: 1 = calculate, 2 = input
      integer   swsalinity         ! Switch for salinity stress: 0 = no stress; 1 = Maas and Hoffman (1977); 2 = Osmotic head
      integer   tsumtime           ! time (nrs of sequential days) with temp above tsumtemp for grass growth [1..20 days, I]

      real(8)   atmtr              ! Daily atmospheric transmission (-)
      real(8)   agerm              ! Coefficient a  of germination
      real(8)   cgerm              ! Coefficient c  of germination
      real(8)   bgerm              ! Coefficient b  of germination
      real(8)   adcrh              ! Level of high atmospheric demand (L/T)
      real(8)   adcrl              ! Level of low atmospheric demand (L/T)
      real(8)   air_filled_root_por ! Air filled root porosity [0..1.0 -, R]
      real(8)   albedo             ! Crop reflection coefficient (-)
      real(8)   alphacrit          ! Critical stress index for compensation of root water uptake (-)
      real(8)   alpJvLier          ! Uniform drought reduction factor based on concept Jong van Lier (-)
      real(8)   amaxtb(30)         ! Maximum CO2 assimilation rate (kg/ha/hr) as function of development stage (-)
      real(8)   avevaptb(2*magrs)  ! Gash interception model: average evaporation intensity during shower (-) as function of time (T)
      real(8)   avprectb(2*magrs)  ! Gash interception model: average rainfall intensity (-) as function of time (T)
      real(8)   cf                 ! Crop factor (-)
      real(8)   cfeic              ! Crop factor wet (-)
      real(8)   cftb(2*magrs)      ! Array with either crop factors (-) or crop height (L) as function of development stage
      real(8)   cfeictb(2*magrs)   ! Array with crop factors, wet (-) 
      real(8)   ch                 ! Crop height (cm)
      real(8)   chtb(2*magrs)      ! Array with crop heights (cm) as function of development stage
      real(8)   c_mroot            ! Maintenance coefficient of root [0.0..1.0 kg O2/kg/d, R]
      real(8)   cofab              ! Interception coefficient Von Hoyningen-Hune and Braden (L)
      real(8)   cropend(macrop)    ! Array with crop end dates
      real(8)   cropstart(macrop)  ! Array with crop start dates
      real(8)   cumdens(202)       ! Cumulative root density as function of relative soil depth (-)
      real(8)   cuptgraz           ! Cumulative dry weight of grass consumed with animal grazing for actual run (kg/ha)
      real(8)   cuptgrazpot        ! Cumulative dry weight of grass consumed with animal grazing for potential run (kg/ha)
      real(8)   cvl                ! Efficiency of assimilate conversion into leaves (kg/kg)
      real(8)   cvo                ! Efficiency of assimilate conversion into storage organs (kg/kg)
      real(8)   cvr                ! Efficiency of assimilate conversion into roots (kg/kg)
      real(8)   cvs                ! Efficiency of assimilate conversion into stems (kg/kg)
      real(8)   cwdm               ! Dry weight of dead and living plant organs (kg/ha)
      real(8)   cwdmpot            ! Dry weight of dead and living plant organs for potential growth (kg/ha)
      real(8)   difpp              ! Diffuse irradiation perpendicular to direction of light (J/m2/s)
      real(8)   dlc                ! Shortest day length (T) for any crop development
      real(8)   dlo                ! Minimum day length (T) for optimal crop development
      real(8)   dry_mat_cont_roots ! Dry matter content of roots [0..1.0 -, R]
      real(8)   dsinbe             ! Daily total of effective solar height (s)
      real(8)   dtsmtb(30)         ! Increase of temperature sum ( C) as function of daily average temperature ( C)
      real(8)   dvs                ! Crop development stage (-)
      real(8)   dvsend             ! Crop development stage at harvest
      real(8)   dwlv               ! Dry weight of plant leafs of actual crop (kg/ha)
      real(8)   dwlvCrop           ! Dry weight of deceased leafs that remain attached to plant (kg/ha)
      real(8)   dwlvSoil           ! Dry weight of deceased leafs allocated to soil (kg/ha)
      real(8)   dwlvpot            ! Dry weight of plant leafs of potential crop (kg/ha)
      real(8)   dwrt               ! Dry weight of plant roots of actual crop (kg/ha)
      real(8)   dwrtpot            ! Dry weight of plant roots of potential crop (kg/ha)
      real(8)   dwso               ! Dry weight of plant storage organs of actual crop (kg/ha)
      real(8)   dwst               ! Dry weight of plant stem of actual crop (kg/ha)
      real(8)   dwstpot            ! Dry weight of plant stem of potential crop (kg/ha)
      real(8)   eff                ! Light use efficiency of a leaf (kg CO2 / J adsorbed)
      real(8)   f_senes            ! Reduction factor for senescence, used for maintenance respiration [0..1.0 -, R]
      real(8)   fltb(30)           ! Fraction of total dry matter increase partitioned to the leaves (-) as function of dvs
      real(8)   fotb(30)           ! Fraction of total dry matter increase partitioned to the storage organs (-) as function of dvs
      real(8)   frtb(30)           ! Fraction of total dry matter increase partitioned to the roots (-) as function of dvs
      real(8)   fstb(30)           ! Fraction of total dry matter increase partitioned to the stems (-) as function of dvs
      real(8)   gasst              ! Total gross assimilation for actual crop (kg/ha)
      real(8)   gasstpot           ! Total gross assimilation for potential crop (kg/ha)
      real(8)   gc                 ! Ground cover in case of a crop
      real(8)   HarLosOrm_tot      ! Harvest losses added to the soil (roots+fr_shoots_fr_stems+fr_stor.organs) at harvest (kg/ha DM) 
      real(8)   hdrygerm           ! Criterium Hdry of germination
      real(8)   hwetgerm           ! Criterium Hwet of germination
      real(8)   hlim1              ! Pressure head above which root water uptake stops (L)
      real(8)   hlim2l             ! Pressure head below which optimum water uptake starts for sub layer (L)
      real(8)   hlim2u             ! Pressure head below which optimum water uptake starts for top layer (L)
      real(8)   hlim3h             ! Pressure head below which water uptake reduction starts at high Tpot (L)
      real(8)   hlim3l             ! Pressure head below which water uptake reduction starts at low Tpot (L)
      real(8)   hlim4              ! Wilting point, no root water uptake at lower soil water pressure heads (L)
      real(8)   kdif               ! Extinction coefficient for diffuse visible light (-)
      real(8)   kdir               ! Extinction coefficient for direct visible light (-)
      real(8)   lai                ! Leaf area index
      real(8)   laiem              ! Leaf area index (-) at crop emergence
      real(8)   laiexp             ! Leaf area index (-) in exponential growth stage of actual crop
      real(8)   laiexppot          ! Leaf area index (-) in exponential growth stage of potential crop
      real(8)   laimax             ! Maximum leaf area index reached during growth of actual grass crop (-)
      real(8)   laipot             ! Leaf area index for potential run (-)
      real(8)   lv(366)            ! Array with leaf weight (kg/ha) as function of crop day number of actual crop
      real(8)   lvpot(366)         ! Array with leaf weight (kg/ha) as function of crop day number of potential crop
      real(8)   lvage(366)         ! Array with leaf age (d) as function of crop day number of actual crop
      real(8)   lvagepot(366)      ! Array with leaf age (d) as function of crop day number of potential crop
      real(8)   max_resp_factor    ! Ratio root total respiration / maintenance respiration [1..5.0 -, R]
      real(8)   mowrest            ! Dry weight of above ground grass (leaves + stems) after mowing(kg/ha)
      real(8)   mrest              ! Total maintenance respiration for actual crop (kg/ha)
      real(8)   mrestpot           ! Total maintenance respiration for potential crop (kg/ha)
      real(8)   mrftb(2*magrs)     ! Array with ratio root total respiration / maintenance respiration as function of DVS (kg/m3)
      real(8)   perdl              ! Maximum relative death rate of leaves due to water stress (/T)
      real(8)   pfreetb(2*magrs)   ! Gash interception model: free throughfall coefficient (-) as function of time (T)
      real(8)   pstemtb(2*magrs)   ! Gash interception model: stem flow coefficient (-) as function of time (T)
      real(8)   siccaptb(2*magrs)  ! NHI interception model: interception capacity as a function of time (T)
      real(8)   fimin              ! start-up saturation fraction for relative interception evaporation (-)
      real(8)   siccapact          ! interceprion storage capacity of canopy (cm)
      real(8)   sicact             ! amount of water stored on canopy (cm)
      real(8)   siccaplai          ! interception storage per unit of LAI (cm/LAI)
      real(8)   q10                ! Relative increase of respiration rate with temperature (/10 C)
      real(8)   q10_microbial      ! Relative increase in microbial respiration at temperature increase of 10 ºC [1.0..4.0 -, R]
      real(8)   q10_root           ! Relative increase in root respiration at temperature increase of 10 ºC [1.0..4.0 -, R]
      real(8)   hrz1               ! Soil water pressure head in root zone at depth z1 (L)
      real(8)   rd                 ! Rooting depth (L)
      real(8)   rdc                ! Maximum rooting depth of particular crop (L)
      real(8)   rdctb(22)          ! Array with relative root density (-) as function of relative root depth (-)
      real(8)   rdi                ! Initial rooting depth (L)
      real(8)   rdpot              ! Rooting depth for potential run (L)
      real(8)   rdrrtb(30)         ! Array with relative death rates of roots (kg/kg/d) as function of development stage (-)
      real(8)   rdrstb(30)         ! Array with relative death rates of stems (kg/kg/d) as function of development stage (-)
      real(8)   rdmax              ! Maximum rooting depth in soil profile (L)
      real(8)   reltr              ! relative transpiration factor that reduces crop growth (-)
      real(8)   rfsetb(30)         ! Reduction factor for senescence (-) as function of development stage (-)
      real(8)   rgrlai             ! Maximum relative increase in leaf area index (/T)
      real(8)   rid                ! Real day number of detailed grass crop (d)
      real(8)   rml                ! Relative maintenance respiration rate of leaves (kg CH2O)/kg/d)
      real(8)   rmo                ! Relative maintenance respiration rate of storage organs (kg CH2O)/kg/d)
      real(8)   rmr                ! Relative maintenance respiration rate of roots (kg CH2O)/kg/d)
      real(8)   rms                ! Relative maintenance respiration rate of stems (kg CH2O)/kg/d)
      real(8)   rootcoefa          ! Defines relative distance at which mean soil water content occurs between roots
      real(8)   rooteff            ! Root system efficiency factor [0..1.0 -, R]
      real(8)   rootradius         ! Root radius drought stress (cm)
      real(8)   root_radiusO2      ! Root radius oxygen stress (m)
      real(8)   rri                ! Maximum daily increase of rooting depth (L/T)
      real(8)   rsc                ! Minimum canopy resistance of dry crop (T/L)
      real(8)   rsw                ! Canopy resistance of intercepted water (T/L)
      real(8)   scanopytb(2*magrs) ! Gash interception model: storage capacity of canopy (-) as function of time (T)
      real(8)   shape_factor_rootr ! Shape factor for exponential decrease of root respiration rate with depth [0..1.0 -, R]
      real(8)   sla(366)           ! Array with specific leaf area (ha/kg) as function of crop day number of actual crop
      real(8)   slapot(366)        ! Array with specific leaf area (ha/kg) as function of crop day number of potential crop
      real(8)   slatb(30)          ! Array with specific leaf area (ha/kg) as function of development stage
      real(8)   spa                ! Specific pod area (ha/kg)
      real(8)   span               ! Life span of leaves at optimum conditions (T)
      real(8)   spec_weight_root_tissue ! Specific weight of non-airfilled root tissue [0.d0..1.d5 kg root/m3 root, R]
      real(8)   specific_resp_humus ! Respiration rate of humus at 25 ºC [0.0..1.0 kg O2/kg C/d, R] 
      real(8)   srl                ! Specific root length [0.d0..1.d10 m root/kg root, R]      
      real(8)   ssa                ! Specific stem area (ha/kg)
      real(8)   tadw               ! Dry weight of plant minus roots of actual growth (kg/ha)
      real(8)   tadwpot            ! Dry weight of plant minus roots of potential growth (kg/ha)
      real(8)   tagp               ! Dry weight of dead and living grass organs (kg/ha)
      real(8)   tagppot            ! Dry weight of dead and living grass organs for potential run (kg/ha)
      real(8)   tagpt              ! Dry weight of harvested grass (kg/ha)
      real(8)   tagptpot           ! Dry weight of harvested grass for potential run (kg/ha)
      real(8)   tbase              ! Lower threshold temperature for ageing of leaves ( C)
      real(8)   tdwi               ! Initial total crop dry weight (kg/ha)
      real(8)   tmnftb(30)         ! Reduction factor for maximum assimilation rate (-) as function of minimum day temperature ( C)
      real(8)   tmpftb(30)         ! Reduction factor for maximum assimilation rate (-) as function of average day temperature ( C)
      real(8)   tsum               ! Temperature sum from cropstart to cropend ( C)
      real(8)   tsumam             ! Temperature sum from anthesis to maturity ( C)
      real(8)   tsumea             ! Temperature sum from emergence to anthesis ( C)
      real(8)   tsum200            ! Temperature sum from 1st day of calendar year (C)
      real(8)   tsumemeopt         ! Temperature sum for crop emergence under optimal conditions
      real(8)   tsumgerm           ! Temperature sum during germination
      real(8)   tsumtemp           ! temperature limit to initiate grass growth  [0.0..20.0 grC, R]
      real(8)   tsumdepth          ! depth at which temp above tsumtemp for grass growth [0.0..100.0 cm below soil surface, R]
      real(8)   TBASEM             ! Lower threshold temp. for emergence (C)
      real(8)   TEFFMX             ! max. eff. temp. for emergence (C)

      real(8)   var_a              ! Variance of root radius [0.d0..1.d0 -, R]
      real(8)   w_root_ss          ! Dry weight of roots at soil surface [0.0..10.0 kg/m3, R]
      real(8)   wiltpoint          ! Minimum pressure head at interface soil-root (cm)

      real(8)   wlv                ! Dry weight of plant leaves (kg/ha)
      real(8)   wlvpot             ! Dry weight of plant leaves for potential growth (kg/ha)
      real(8)   wrtb(2*magrs)      ! Array with dry weight of root at soil surface as function of DVS (kg/m3)
      real(8)   wrt                ! Dry weight of plant root (kg/ha)
      real(8)   wrtpot             ! Dry weight of plant root for potential growth (kg/ha)
      real(8)   wso                ! Dry weight of storage organ (kg/ha)
      real(8)   wsopot             ! Dry weight of storage organ for potential growth (kg/ha)
      real(8)   wst                ! Dry weight of plant stem (kg/ha)
      real(8)   wstpot                ! Dry weight of plant stem for potential growth (kg/ha)
      logical   flanthesis         ! Flag indicating anthesis stage of a crop
      logical   flHarvest           ! Flag indicating that grass should be harvested, actual crop
      logical   flHarvestDay        ! Flag indicating that current day is harvest day
      logical   flHarvestpot        ! Flag indicating that grass should be harvested, potential crop
      logical   flGrazing           ! Flag indicating that cattle grazes the grass, actual run
      logical   flGrazingpot        ! Flag indicating that cattle grazes the grass, potential run
      character(len=40) cropfil(macrop)   ! Array with names of crop files
      character(len=40) cropname(macrop)  ! Array with crop names
      character(len=80) pathcrop          ! Path to folder with crop input files
      character(len=200) inifil           ! Name of file with output data of previous day which is used for initialization
!     CO2
      logical   flCO2              ! Flag indicating correction of CO2
      real(8)   fco2amax           ! factor to correct AMAX for CO2
      real(8)   fco2eff            ! factor to correct EFF for CO2
      real(8)   fco2tra            ! factor to correct TRA for CO2
      real(8)   co2amaxtb(30)      ! table with factors to correct AMAX for CO2
      real(8)   co2efftb(30)       ! table with factors to correct EFF for CO2
      real(8)   co2tratb(30)       ! table with factors to correct TRA for CO2
      integer   co2year(mayrs)     ! table with years for which CO2 concentrations are given
      real(8)   co2ppm(mayrs)      ! table with CO2 concentrations (ppm), for each year in co2year
!     vernalisation
      real(8)   verndvs            ! critical development stage after which the effect of vernalisation is halted [-]
      real(8)   vernsat            ! saturated vernalisation requirement [d]
      real(8)   vernbase           ! base vernalisation requirement [d]
      real(8)   verntb(30)         ! table with rate of vernalisation as function of tav [days/degrees]

! --- Nitrogen: crop and soil management
      logical   flCropNut          ! Flag indicating simulation of crop nutrient stress
      real(8)   amFERT             ! amount of applied Fertilizer (kg/ha/d N)

! --- soilwater variables
      integer   MaxIterTime        ! Maximum cputime (secs), introduced to be able to interrupt (near) endless iterations
      integer   MaxBackTr
      integer   MaxIt
      integer   Itnumb(100,2)
      real(8)   CritDevh1Cp        ! Convergence criterium for Richards equation: relative difference in pressure heads (-)
      real(8)   CritDevh2Cp        ! Convergence criterium for Richards equation: absolute difference in pressure heads (L)
      real(8)   CritDevPondDt
      logical   fldumpconvcrit     ! flag to generate additional output about convergence-warnings from subr Headcalc
      logical   flksatexm          ! flag to use values, examined in lab or field, of saturated hydraulic conductivity (L/T) for each soil layer
      logical   flMaxIterTime      ! flag to enable input of Maximum cputime
      logical   FlRunoff
      logical   swcaprise          ! flag to minimize cap.rise to rootzone (for experts only)
      logical   swcapriseoutput    ! flag to generate an output file with cap.rise to/form rootzone (for experts only)

      real(8)   h0max
      real(8)   k1max
      real(8)   q0

      integer   afo                ! Internal number of output file *.AFO with formatted hydrologic data for soil water quality models
      integer   aun                ! Internal number of output file *.AUN with unformatted hydrologic data for soil water quality models
      integer   bal                ! Internal number of output file *.BAL with overview of water balance
      integer   blc                ! Internal number of output file *.BLC with all water balance components in detail
      integer   bma                ! Internal number of output file *.BMA with detailed yearly water balance in case of macropores
      integer   bpegwl             ! Node at bottom of perched groundwater
      integer   botcom(maho)       ! Array with number of bottom compartments in each soil layer
      integer   dra                ! Internal number of drainage input file *.DRA
      integer   dramet             ! Switch for lateral drainage: 1 = table of flux - groundwater level; 2 = Hooghoudt or Ernst; 
                                   !                              3 = drainage/infiltration resistance
      integer   inc                ! Internal number of output file *.INC with incremental water balance data
      integer   indeks(macp)       ! Index denoting wetting or drying curve in case of hysteresis: 1 = wetting; -1 = drying
      integer   ipos               ! Switch for position of drain (see *.DRA input file for overview)
      integer   isoillay(maho)     ! Number of soil layer, starting with 1 at the soil surface
      integer   layer(macp)        ! Array with soil layer number for each compartment
      integer   msteps             ! Maximum number of iteration steps during a day to solve Richards equation
      integer   ncomp(macp)        ! Array with number of compartments in each sublayer
      integer   nhead              ! Number of initial soil water pressure heads as provided in the input
      integer   nodgwl             ! Node directly above groundwater level
      integer   nod1lay(maho)      ! node nr of first node of each soil layer (from top to bottom)
      integer   nodfrostbot        ! Node nr of deepest node with frost conditions
      integer   npegwl             ! Node directly above perched groundwater level
      integer   nrlevs             ! Number of drainage levels
      integer   nrstaring          ! Number of soil type [1..18] according to Staring series (Wosten et al., 2001)
      integer   nsublay            ! Number of sublayers in the soil profile
      integer   numbit             ! Iteration number for solving Richards equation
      integer   numlay             ! Number of (physical) soil layers
      integer   numnod             ! Number of nodes or compartments
      integer   numnodnew          ! Number of desired nodes for soil water quality models
      integer   numtab(macp)       ! Number of table entries of soil physical values for each model compartment
      integer   numtablay(maho)    ! Number of table entries of soil physical values for each soil layer
      integer   rot                ! Internal number of output file *.ROT with microscopic root water extraction data 
      integer   str                ! Internal number of output file *.STR with stress factors for transpiration
      integer   sw2                ! Switch for prescribed bottom flux: 1 = sine function; 2 = table
      integer   sw3                ! Switch for prescribed hydraulic head of deep aquifer: 1 = sine function; 2 = table
      integer   sw4                ! Switch for extra groundwater flux as function of time: 0 = no extra flux; 1 = include extra flux
      integer   swafo              ! Switch for extra output file with formatted data for water quality models: 
                                   !       0 = no output; 1 = output to file *.AFO; 2 = output to file *.BFO
      integer   swaun              ! Switch for extra output file with unformatted data for water quality models: 
                                   !       0 = no output; 1 = output to file *.AUN; 2 = output to file *.BUN
      integer   swblc              ! Switch for output file with detailed yearly water balance *.BLC: 0 = no; 1 = yes
      integer   swbotb             ! Switch for bottom boundary condition (see *.SWP input file for overview)
      integer   swbotb3Impl        ! Switch for implicit solution with lower boundary option 3 (Cauchy): 0 = explicit, 1 = implicit
      integer   SwBotb3ResVert     ! Switch to suppress addition of vertical resistance between bottom of model and groundwater level
      integer   swcfbs             ! Switch for use of coefficient CFBS to convert potential ET into potential E: 0 = no; 1 = yes
      integer   swdiscrvert        ! Switch to convert vertical discretization for soil water quality models: 0 = no; 1 = yes
      integer   swdivd             ! Switch to distribute drainage flux vertically according to transmissivity: 0 = no; 1 = yes
      integer   swdivdinf          ! Switch to distribute infiltration flux vertically according to transmissivity, separatly from drainage fluxes: 0 = no; 1 = yes
      integer   swdislay           ! Switch to distribute drainage flux vertically with a given position of the top of the model discharge layers: 0 = no; 1 = yes
      integer   swtopdislay(madr)  ! Switch, for each drainage level, to distribute drainage flux vertically with a given position of the top of the model discharge layers: 0 = no; 1 = yes
      integer   swdra              ! Switch for simulation of lateral drainage: 0 = no drainage; 1 = use basic drainage routine; 
                                   !                                            2 = simulate drainage and surface water
      integer   swfrost            ! Switch for reduction of hydraulic conductivity in case of frost: 0 = no; 1 = yes
      integer   swhyst             ! Switch for hysteresis of soil moisture retention function: 0 = no; 1 = yes
      integer   swinco             ! Switch for initial soil moisture condition: 1 = pressure heads; 2 = hydrostatic equilibrium; 
                                   !                                             3 = final pressure heads from previous simulation
      integer   swkmean            ! Switch for mean of hydraulic conductivity: 1 = unweighted arithmic mean, 2 = weighted arithmic mean
                                   !                                            3 = unweighted geometric mean,4 = weighted geometric mean
      integer   swkimpl            ! Switch for implicit solution with hydraulic conductivity: 0 = explicit, 1 = implicit
      integer   swliminf           ! Switch for limit of infiltration head to the waterdepth in the channel: 0 = nolimit, 1 = limitation
      integer   swoutputmodflow    ! Switch for extra output file with data for Modflow
      integer   swoxygen           ! Switch oxygen stress: 1 = concept Feddes et al. (1978); 2 = concept Bartholomeus et al. (2008)
      integer   swoxygentype       ! Switch for method oxygen stress calculation: 1 = physical processes; 2 = repro functions
      integer   swpondmx           ! Switch for time dependent maximum amount of ponding (L) on soil surface before runoff starts
      integer   swqhbot            ! Switch for flux-groundwater level relationship: 1 = exponential function; 2 = tabular function
      integer   swcofqhc           ! Switch for additional flux added to exponential flux-groundwater level relationship: 0 = no, 1 = yes
      integer   swredu             ! Switch for reduction of soil evaporation: 0 = no empirical function; 1 = use function of Black; 
                                   !                                           2 = use function of Boesten/Stroosnijder
      integer   swsophy            ! Switch for input of soil hydraulica properties as function parameters (0) or as table (1)
      integer   ientrytab(macp,0:matabentries)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment
      integer   ientrytablay(maho,0:matabentries) ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each soil layer

      integer   swtopsub           ! Switch for topsoil or subsoil: 1 = topsoil, 2 = subsoil
      integer   swvap              ! Switch for output file *.VAP with soil profile data (water content, pressure head, 
                                   !        solute concentration, temperature): 0 = no; 1 = yes
      integer   vap                ! Internal number of output file *.VAP with soil profile data 
      integer   wba                ! Internal number of output file *.WBA with cumulative water balance data
      real(8)   aqamp              ! Amplitude of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(8)   aqave              ! Average hydraulic head in deep aquifer (L)
      real(8)   aqper              ! Period of prescribed sine wave of hydraulic head in deep aquifer (T)
      real(8)   aqtmax             ! Time with maximum hydraulic head in deep aquifer (T)
      real(8)   atop(18,6)         ! Parameters of reprofunctions for oxygen stress
      real(8)   basegw             ! Depth of impervious layer (L) for drainage according to Hooghoudt or Ernst
      real(8)   bdens(maho)        ! Array with dry bulk density for each soil layer (M/L3)
      real(8)   c_top(macp)        ! Oxygen concentration at top of compartment(kg/m3)
      real(8)   cfbs               ! Coefficient (-) to convert potential evapotranspiration into potential evaporation
      real(8)   cgird              ! Cumulative amount of gross irrigation (L)
      real(8)   cinund             ! Cumulative amount of inundation (L)
      real(8)   cnird              ! Cumulative amount of net irrigation (L)
      real(8)   cofani(maho)       ! Anisotropy coefficient (horizontal / vertical saturated hydraulic conductivity) (-)
      real(8)   cofgen(12,macp)    ! Array (adjusted by hysteresis) soil hydraulic parameters according to Mualem - van Genuchten for each soil layer
      real(8)   cofqha             ! Coefficient A in exponential relationship between drainage flux and groundwater level (L/T)
      real(8)   cofqhb             ! Coefficient B in exponential relationship between drainage flux and groundwater level (/T)
      real(8)   cofqhc             ! Coefficient C (flux) in exponential relationship between drainage flux and groundwater level (L/T)
      real(8)   cofred             ! Soil evaporation coefficient of Black or Boesten/Stroosnijder
      real(8)   cqbot              ! Cumulative amount of water flow through bottom of simulated soil column (L)
      real(8)   cqbotdo            ! Cumulative amount of water (L) passed through the soil column bottom in downward direction
      real(8)   cqbotup            ! Cumulative amount of water (L) passed through the soil column bottom in upward direction
      real(8)   cqdra              ! Cumulative amount of lateral drainage (L)
      real(8)   cqdrain(Madr)      ! Cumulative amount of lateral drainage for each drainage level (L)
      real(8)   cqdrainin(Madr)    ! Cumulative infiltration flux (L) for each drainage level
      real(8)   cqdrainout(Madr)   ! Cumulative drainage flux (L) for each drainage level
      real(8)   cqprai             ! Cumulative amount of net rain (L)
      real(8)   cqrot              ! Cumulative amount of extracted water by roots (L)
      real(8)   cqtdo              ! Cumulative amount of water (L) passed through the soil surface in downward direction
      real(8)   cqtup              ! Cumulative amount of water (L) passed through the soil surface in upward direction
      real(8)   CritDevMasBal      ! Maximum error in water balance (L)
      real(8)   CriterHr           ! Maximum difference of Hroot between iterations; convergence criterium  (L)
      real(8)   crunoff            ! Cumulative runoff (L)
      real(8)   crunon             ! Cumulative amount of runon (L)
      real(8)   deepgw             ! hydraulic head in aquifer (L)
      real(8)   dimoca(macp)       ! Differential soil moisture capacity (/L)
      real(8)   disnod(macp+1)     ! Distance between actual node and upper node (L)
      real(8)   drainl(Madr)       ! Drainage level (maximum of depth of drain and surface water level) (L)
      real(8)   drares(Madr)       ! Array with drainage resistance (T) for each drainage level
      real(8)   dz(macp)           ! Compartment thickness (L)
      real(8)   dznew(macp)        ! Desired thickness of compartments for soil water quality models (L)
      real(8)   entres             ! Drain entry resistance (T)
      real(8)   es0                !  potential evaporation rate from a wet bare soil [mm/d]
      real(8)   et0                !  potential transpiration rate from a dry crop [mm/d]
      real(8)   ew0                !  potential transpiration rate from a wet crop [mm/d]
      real(8)   evp(macp)          ! Internal evaporation flux of top soil compartments (L/T)
      real(8)   FacDpthInf         ! Factor to reduce Depth of Infiltration layer as fraction of max depth Discharge layer (-)
      real(8)   ftopdislay(madr)   ! Array with factor for function to determine depth of top of model discharge layer for each drain level, see also swtopdislay (L)
      real(8)   geofac             ! Geometry factor (-) for analytical drainage formula of Ernst
      real(8)   gwl                ! Groundwater level (L)
      real(8)   gwlconv            ! Maximum difference of groundwater levels between iterations to solve Richards equation
      real(8)   gwli               ! Groundwater level (L) at start of simulation
      real(8)   gwlinp             ! Prescribed groundwater level (L) for current time
      real(8)   gwlm1              ! Groundwater level (L) at former time level
      real(8)   gwltab(mabbc*2)    ! Array with prescribed groundwater level (L) as function of time (T)
      real(8)   h(macp)            ! Soil water pressure head (L)
      real(8)   h_enpr(macp)       ! Soil water Entry Pressure head for Modified MualemVanGenuchten curve (L)
      real(8)   haqtab(mabbc*2)    ! Array with specified hydraulic head in deep aquifer (L) as function of time (T)
      real(8)   hatm               ! Pressure head of air (L) near the soil surface
      real(8)   hbot               ! Soil water pressure head (L) at bottom of soil column
      real(8)   hbotab(mabbc*2)    ! Array with specified pressure head of lowest compartment (L) as function of time (T)
      real(8)   hcomp(macp)        ! Array with prescribed height of numerical compartments (L) for each sublayer
      real(8)   hdrain             ! Mean drainage level (L) to derive regional average groundwater level for bottom boundary condition
      real(8)   hleaf              ! Pressure head inside leaves (cm)
      real(8)   hm1(macp)          ! Soil water pressure head (L) at former time level
      real(8)   hroot(macp)        ! Pressure head of a compartment at the root-soil interface (L)
      real(8)   hsublay(macp)      ! Array with prescribed height of sublayers (L)
      real(8)   hsurf              ! Soil water pressure head at the soil surface (cm)
      real(8)   Hxylem             ! Pressure head in root xylem (L)
      real(8)   igrai              ! Intermediate amount of gross rainfall (L)
      real(8)   ies0               ! Intermediate potential evaporation rate from a wet bare soil [cm/d]
      real(8)   iet0               ! Intermediate potential transpiration rate from a dry crop [cm/d]
      real(8)   iew0               ! Intermediate potential transpiration rate from a wet crop [cm/d]
      real(8)   iintc              ! Intermediate amount of rainfall interception by vegetation (L)
      real(8)   infres(Madr)       ! Array with infiltration resistance (T) for each drainage level
      real(8)   inpola(macp)       ! Weight for interpolation between current node and upper node
      real(8)   inpolb(macp)       ! Weight for interpolation between current node and lower node
      real(8)   inq(macp+1)        ! Array with intermediate amounts of water flow between current and upper compartment (L)
      real(8)   inqdra(Madr,macp)  ! Array with intermediate amounts of lateral drainage for each level and compartment (L)
      real(8)   inqrot(macp)       ! Array with intermediate amounts of extracted water by roots for each compartment (L)
      real(8)   ipondbeg           ! Ponding water layer (L) on soil surface at start of current intermediate period
      real(8)   iprec              ! Intermediate amount of gross precipitation + gross irrigation (L)
      real(8)   iqbot              ! Intermediate amount of water flow through bottom of simulated soil column (L)
      real(8)   iqdra              ! Intermediate amount of lateral drainage (L)
      real(8)   iqrot              ! Intermediate amount of extracted water by roots (L)
      real(8)   iqredwet           ! Intermediate amount of reduced root water extraction due to wet conditions (L)
      real(8)   iqreddry           ! Intermediate amount of reduced root water extraction due to dry conditions (L)
      real(8)   iqredsol           ! Intermediate amount of reduced root water extraction due to salt conditions (L)
      real(8)   iqredfrs           ! Intermediate amount of reduced root water extraction due to frost conditions (L)
      real(8)   iruno              ! Intermediate amount of runoff (L)
      real(8)   irunon             ! Intermediate amount of runon (L)
      real(8)   issnowbeg          ! Amount of snow in soil water equivalent (L) at start of current intermediate period
      real(8)   ithetabeg(macp)    ! Array with volumetric soil water contents (-) for each compartment at start of intermediate period
      real(8)   k(macp+1)          ! Array with soil hydraulic conductivity (L/T) for each numerical compartment
      real(8)   khbot              ! Horizontal hydraulic conductivity of bottom layer (L/T)
      real(8)   khtop              ! Horizontal hydraulic conductivity of top layer (L/T)
      real(8)   kmean(macp+1)      ! Array with mean soil hydraulic conductivity (L/T) at the interface of current and upper compartment
      real(8)   Kroot              ! Hydraulic, radial conductivity of root tissue (L/T)
      real(8)   ksatfit(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: fitted on VG based on lab data
      real(8)   ksatexm(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: examined in lab or field 
      real(8)   ksatthr(maho)      ! Array with saturated hydraulic conductivity (L/T) for each soil layer: to interpolate VG and Ksatexm
      real(8)   kstem              ! Conductance in the path from leaf to root xylem (/d)
      real(8)   kvbot              ! Vertical hydraulic conductivity of bottom layer (L/T)
      real(8)   kvtop              ! Vertical hydraulic conductivity of top layer (L/T)
      real(8)   L(Madr)            ! Array with spacing between drains for each drainage level (L)
      real(8)   ldwet              ! Length of dry period (L) as used in Black's model for reduction of soil evaporation
      real(8)   mfluxtable(maho,801)  ! Reference table with matric flux potential of each soil layer (L2/T)
      real(8)   mflux(macp)        ! Actual matric flux potential of each node (L2/T)
      real(8)   mroot(macp)        ! Matrix flux head of a compartment at the root-soil interface (L2/T)
      real(8)   OxygenIntercept(6) ! Parameters of reproduction function for oxygen stress according to Bartholomeus
      real(8)   OxygenSlope(6)     ! Parameters of reproduction function for oxygen stress according to Bartholomeus
      real(8)   paramvg(10,maho)   ! Array with input values of soil hydraulic parameters according to Mualem - van Genuchten for each soil layer
      real(8)   pegwl              ! Perched groundwater level (L)
      real(8)   pond               ! Height of ponding layer (L)
      real(8)   pondini            ! Ponding water layer (L) on soil surface at start of current water balance period
      real(8)   pondm1             ! Ponding water layer (L) on soil surface at former time level
      real(8)   pondmx             ! Maximum amount of ponding (L) on soil surface before runoff starts
      real(8)   pondmxtab(2*mairg) ! Table with time-dependent input (date,value) for maximum amount of ponding (L) on soil surface before runoff starts
      real(8)   q(macp+1)          ! Soil water flux between current compartment and upper compartment (L/T)
      real(8)   qbot               ! Water flux through bottom of simulated soil column (L/T)
      real(8)   qbotab(mabbc*2)    ! Array with specified bottom flux (L/T) as function of time (T)
      real(8)   qbot_nonfrozen     ! Water flux through bottom of non-frozen soil column (L/T)
      real(8)   qdra(Madr,macp)    ! Array with lateral drainage flux (L/T) for each drainage level and compartment
      real(8)   qdrain(Madr)       ! Total lateral drainage flux (L/T) for each drainage level
      real(8)   qdraincomp(macp)   ! Total lateral drainage flux (L/T) for each compartment
      real(8)   qdrtab(50)         ! Array with lateral drainage flux (L/T) as function of groundwater level (L)
      real(8)   qdrtot             ! Total lateral drainage flux (L/T)
      real(8)   qimmob(macp)       ! Soil water flux between mobile and immobile fraction in case of fingered flow (L/T)
      real(8)   qrosum             ! Total root water extraction flux (L/T)
      real(8)   qredwetsum         ! Total reduction of root water extraction due to wet conditions (L/T)
      real(8)   qreddrysum         ! Total reduction of root water extraction due to dry conditions (L/T)
      real(8)   qredsolsum         ! Total reduction of root water extraction due to salt conditions (L/T)
      real(8)   qredfrssum         ! Total reduction of root water extraction due to frost conditions (L/T)
      real(8)   qrot(macp)         ! Array with root water extraction flux for each compartment (L/T)
      real(8)   qtop               ! Water flux through soil surface (L/T)
      real(8)   relsatthr(maho)    ! Array with relative saturation (-) for each soil layer: to interpolate VG and Ksatexm
      real(8)   ResultsOxygenStress(19,macp) ! array with results for OxygenStress; for output only 
      real(8)   reva               ! Actual soil evaporation rate (L/T)
      real(8)   rfcp(macp)         ! Reduction factor for frozen conditions in each model compartment (-)
      real(8)   rimlay             ! Vertical resistance of aquitard (T)
      real(8)   rmax(macp)         ! Radius around roots in which water is extracted (L)
      real(8)   RootPhi(macp)      ! Factor Phi of a compartment used in drought reduction of De Jong van Lier et al. (T/L)
      real(8)   RootRho(macp)      ! Factor Rho of a compartment used in drought reduction of De Jong van Lier et al. (/L2)
      real(8)   rsigni             ! Minimum amount of rainfall (L) which resets the empirical soil evaporation reduction models
      real(8)   rsoil              ! Soil resistance of wet soil of PMdirect (T/L)
      real(8)   rsro               ! Drainage resistance for surface runoff (T)
      real(8)   rsroexp            ! Exponent to calculate surface runoff (T)
      real(8)   Rxylem             ! Mean radius of xylem tube inside roots (L)
      real(8)   runon              ! Water runon flux (L/T)
      real(8)   runonarr(maday)    ! Array with runon (L) data for each day
      real(8)   runots             ! Amount of runoff during a time step (L)
      real(8)   saev               ! Cumulative actual evaporation (L) as used in the model of Boesten/Stroosnijder for reduction of E
      real(8)   shape              ! Shape factor: ratio between the mean and the maximum groundwater level elevation above the drainage base (-)
      real(8)   sinamp             ! Amplitude of prescribed bottom flux (L/T) in case of sine function
      real(8)   sinave             ! Average value of prescribed bottom flux (L/T) in case of sine function
      real(8)   sinmax             ! Time of the year with maximum bottom flux in case of prescribed sine function
      real(8)   spev               ! Cumulative potential evaporation (L) as used in the model of Boesten/Stroosnijder for reduction of E
      real(8)   sptab(5,macp,matab)    ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each model compartment
      real(8)   sptablay(5,maho,matab) ! Soil Physical functions (h,theta,k,dthetadh,dkdtheta) tabulated for each soil layer
      real(8)   StepHr             ! Maximum difference of Hroot and Hxylem between iterations; convergence criterium  (L)
      real(8)   taccur             ! Maximum absolute difference between simulated and calculated potential transpiration rate (cm/d)
      real(8)   Tactual            ! Actual transpiration at former iteration in JongvanLier (cm/d)
      real(8)   tau                ! Minimum pressure head difference (L) to change from wetting to drying in case of hysteresis
      real(8)   theta(macp)        ! Volumic soil water content (-)
      real(8)   thetar(macp)       ! Residual volumic soil water content (-) for each numerical compartment
      real(8)   thetas(macp)       ! Saturated volumic soil water content (-) for each numerical compartment
      real(8)   thetm1(macp)       ! Volumic soil water content (-) at former time level
      real(8)   thetsl(maho)       ! Saturated volumic water content (-) for each soil layer
      real(8)   twilt(macp)        ! Pressure head of a compartment at wilting point (L)
      real(8)   volact             ! Water storage (L) of soil column at current time level
      real(8)   volini             ! Water storage (L) of soil column at start of simulation
      real(8)   volm1              ! Water storage (L) of soil column at former time level
      real(8)   vtair              ! Total depth of air in the soil column (L)
      real(8)   wbalance           ! Cumulative water balance error (L)
      real(8)   wetper(Madr)       ! Array with wet perimeter of drain (L) for each drainage level
      real(8)   z(macp)            ! Depth of a node (L)
      real(8)   zfrostbot          ! Depth of bottom of frost layer (L)
      real(8)   zfrosttop          ! Depth of top of frost layer (L)
      real(8)   zbotdr(Madr)       ! Array with depth of drain bottom for each drain level
      real(8)   zi(macp)           ! Array with soil depths (L) used to specify initial soil water pressure heads
      real(8)   zintf              ! Depth (L) at which fine top layer ends and coarse sub layer starts
      real(8)   ztopdislay(Madr)   ! Array with depth of top of model discharge layer for each drain level, see also swtopdislay (L)
      logical   fldrain            ! Flag indicating basic drainage
      logical   FlHydrLift         ! Flag indicating release of water from root to soil is allowed
      logical   fllowgwl           ! Flag indicating precribed groundwater level below bottom soil column
      logical   flrunon            ! Flag indicating the existance of runon
      logical   ftoph              ! Flag indicating that the pressure head is prescribed at the soil surface
      character(len=16) drfil      ! Name of drainage input file
      character(len=80) pathdrain  ! Path to folder with drainage input files

! --- heat variables
      integer   nheat              ! Number of initial soil temperatures as provided in the input
      integer   swbotbhea          ! Switch for bottom boundary condition: 1 = heat flux is zero; 2 = prescribed temperature
      integer   swtopbhea          ! Switch for top boundary condition: 1 = use air temperatures; 2 = read measured surface temperatures
      integer   swcalt             ! Switch for method of soil water heat flow simulation: 1 = analytical method; 2 = numerical method
      integer   swhea              ! Switch for simulation of soil heat flow: 0 = no; 1 = yes
      integer   swtem              ! Switch for output file *.TEM with soil temperatures: 0 = no; 1 = yes
      integer   tem                ! Internal number of output file *.TEM with soil temperatures
      real(8)   ddamp              ! Damping depth (L) of temperature wave in soil
      real(8)   fclay(macp)        ! Array with gravimatric content of clay (g/g mineral parts) of each numerical compartment
      real(8)   forg(macp)         ! Array with gravimatric organic matter content (g/g mineral parts) of each numerical compartment
      real(8)   fquartz(macp)      ! Array with gravimatric content of sand+silt (g/g mineral parts) of each numerical compartment
      real(8)   orgmat(maho)       ! Array with gravimetric organic matter content (g/g mineral parts) for each soil layer
      real(8)   pclay(maho)        ! Array with gravimetric clay content (g/g mineral parts) for each soil layer
      real(8)   psand(maho)        ! Array with gravimetric sand content (g/g mineral parts) for each soil layer
      real(8)   psilt(maho)        ! Array with gravimetric silt content (g/g mineral parts) for each soil layer
      real(8)   tampli             ! Amplitude of prescribed annual temperature wave (ºC) at soil surface
      real(8)   tebot              ! Temperatures (ºC) at bottom of soil profile
      real(8)   tembtab(mabbc*2)   ! Array with specified bottom temperature (ºC) as function of time (T)
      real(8)   temtoptab(mabbc*2) ! Array with specified soil surface temperature (ºC) as function of time (T)
      real(8)   tfroststa          ! Soil temperature (ºC) where reduction of water fluxes starts
      real(8)   tfrostend          ! Soil temperature (ºC) where reduction of water fluxes ends
      real(8)   timref             ! Time in the year (T) with top of prescribed sine temperature wave
      real(8)   tmean              ! Prescribed mean annual temperature (ºC) at soil surface
      real(8)   tsoil(macp)        ! Array with soil temperatures (ºC) for each compartment
      real(8)   tetop              ! Temperatures (ºC) at top of soil profile (under snow cover)
      real(8)   zh(macp)           ! Array with soil depths (L) used to specify initial soil temperatures
      logical   fltemperature      ! Flag indicating simulation of soil heat flow

! --- snow variables
      integer   snw                ! Internal number of output file *.SNW with snow pack data
      integer   swsnow             ! Switch for simulation of snow accumulation and melt: 0 = no; 1 = yes
      integer   swsublim           ! Switch for suppressing simulation of sublimation of snow: 1 = suppress ! Adaptation 3 for PEARL-MACRO
      real(8)   cgsnow             ! Cumulative amount of gross snow fall (L water)
      real(8)   cmelt              ! Cumulative amount of melted snow (L water)
      real(8)   csnrai             ! Cumulative amount of net snow fall (L water)
      real(8)   csubl              ! Cumulative amount of sublimated snow (L water)
      real(8)   gsnow              ! Gross snow rate (L/T)
      real(8)   igsnow             ! Incremental amount of gross snow fall (L water)
      real(8)   isnrai             ! Incremental amount of net snow fall (L water)
      real(8)   isubl              ! Incremental amount of sublimated snow (L water)
      real(8)   melt               ! Melting rate (L/T)
      real(8)   slw                ! liquid water stored in snow pack           ! RobPearlMacro
      real(8)   snowcoef           ! Snow melt factor (-)
      real(8)   snowinco           ! Amount of snow (L water) at start of balance period
      real(8)   snrai              ! Net rain rate on snow pack (L/T)
      real(8)   ssnow              ! Amount of snow (L water)
      real(8)   subl               ! Sublimation rate (L/T)
      real(8)   TePrRain           ! Temperature above which all precipitation is rain,[ 0.0...5.0 oC, R]
      real(8)   TePrSnow           ! Temperature below which all precipitation is snow,[-5.0...0.0 oC, R]
      logical   flsnow             ! Flag indicating simulation of snow accumulation and melt

! --- solute variables
      integer   nconc              ! Number of initial solute concentrations as provided in the input
      integer   sba                ! Internal number of output file *.SBA with cumulative solute balance components
      integer   swbr               ! Switch to consider mixed reservoir for solute breakthrough in the saturated zone: 0 = no; 1 = yes
      integer   swbotbc            ! Switch for bottom boundary condition of solute-concentration (see *.SWP input file for overview)
      integer   swsolu             ! Switch for simulation of solute transport: 0 = no; 1 = yes
      integer   swsp               ! Switch (in case of solute transport) for simulation of sorption  0 = no; 1 = yes
      real(8)   AgeGwl1m           ! Age (d) of groundwater in upper 1 meter of saturated zone
      real(8)   bexp               ! Exponent in decomposition reduction factor due to dryness (-)
      real(8)   cdrain             ! Mean solute concentration in aquifer or drainage system (M/L3 water)
      real(8)   cirr               ! Solute concentration (M/L3) in irrigation water
      real(8)   cml(macp)          ! Array with solute concentration (M/L3 water) in mobile region
      real(8)   cmsy(macp)         ! Array with dissolved + adsorbed solute concentration (M/L3 soil volume) in mobile region
      real(8)   cpond              ! Mean solute concentration (M/L3) in ponding layer on soil surface
      real(8)   cpre               ! Solute concentration (M/L3) in precipitation
      real(8)   cref               ! Reference solute concentration (M/L3) for Freundlich adsorption
      real(8)   cseep              ! Mean solute concentration in upward seepage water at bottom of profile (M/L3 water)
      real(8)   cseeptab(mabbc*2)  ! Array with Mean solute concentration in upward seepage water at bottom of profile (M/L3 water) as function of time (T)
      real(8)   csurf              ! Total amount of solutes (M/L2) in ponding layer on soil surface
      real(8)   daquif             ! Thickness of saturated aquifer (L) to calculate solute breakthrough to surface water
      real(8)   ddif               ! Molecular diffusion coefficient (L2/T)
      real(8)   decpot(maho)       ! Array with Potential decomposition rate (/T) for each soil layer
      real(8)   decsat             ! Decomposition rate in aquifer (/T)
      real(8)   dectot             ! Cumulative amount of solute decomposition (M/L2)
      real(8)   dtsolu             ! Maximum time step (T) for accurate numerical solution of solute transport equation
      real(8)   fdepth(maho)       ! Array with reduction factor for decomposition (-) for each soil layer
      real(8)   frexp              ! Array with Freundlich exponent (-) for solute adsorption
      real(8)   gampar             ! Reduction factor for decomposition due to low temperatures (/C)
      real(8)   icAgeBot           ! Incremental (over output interval) age (d) of groundwater leaving bottom comp.
      real(8)   icAgeDra(madr)     ! Incremental (over output interval) age (d) of groundwater leaving bottom comp.
      real(8)   icAgeRot           ! Incremental (over output interval) age (d) of groundwater leaving by root uptake
      real(8)   icAgeSur           ! Incremental (over output interval) age (d) of groundwater leaving by surface runoff
      real(8)   isqbot             ! Solute flux at the bottom of the soil column (M/L2/T) 
      real(8)   isqtop             ! Solute flux through the soil top surface (M/L2/T) 
      real(8)   kf(maho)           ! Array with Freundlich coefficient (L3/M) for solute adsorption for each soil layer
      real(8)   kfsat              ! Linear adsorption coefficient in aquifer (L3/M)
      real(8)   ldis(maho)         ! Array with Solute dispersion length (L) for each soil layer
      real(8)   poros              ! Porosity of aquifer (-) to calculate solute breakthrough
      real(8)   rottot             ! Cumulative amount of solutes (M/L2) extracted by plant roots
      real(8)   rtheta             ! Minimum volumetric water content (-) for potential decomposition
      real(8)   salthead           ! Conversion salt concentration (mg/cm3) into osmotic head (cm) [0..1000.0 cm/(mg/cm3), R]
      real(8)   saltmax            ! Threshold salt concentration in soil water  [0..100 mg/cm3, R]
      real(8)   saltslope          ! Decline of rootwater uptake above threshold [0..1.0 cm3/mg, R]
      real(8)   samcra             ! Total amount of solutes (M/L2) entrapped in cracks
      real(8)   samini             ! Total amount of solutes (M/L2) in soil profile at start of current balance period
      real(8)   sampro             ! Total amount of solutes (M/L2) in soil column
      real(8)   solbal             ! Cumulative solute balance (M/L2) for current balance period
      real(8)   sqbot              ! Cumulative amount of solutes (M/L2) passed through the soil column bottom
      real(8)   sqdra              ! Total amount of solutes (M/L2) transported to drainage canals
      real(8)   sqirrig            ! Cumulative amount of solutes (M/L2) in irrigation water
      real(8)   sqprec             ! Cumulative amount of solutes (M/L2) in precipitation
      real(8)   sqrap              ! Cumulative amount of solutes (M/L2) in rapid drainage
      real(8)   sqsur              ! Cumulative amount of solutes (M/L2) transported to surface water
      real(8)   tscf               ! Relative uptake of solutes by roots (-)
      real(8)   zc(macp)           ! Array with soil depths (L) used to specify initial solute concentrations
      logical   flsolute           ! Flag indicating simulation of solute transport
      logical   flAgeTracer        ! Flag indicating simulation of Ageing (groundwater age)

! --- macropore Input parameters
      integer NumSbDm              ! Number of Subdomains in IC domain (-)
      integer SwDarcy              ! Switch for eliminating Darcy flow unsaturated zone
      integer SwDrRap              ! Switch for kind of drainage function (-) TEMPORARY: TEST optio
      integer SwPowM               ! Switch for double convex/concave freq. distr. curve (-)
      integer SwShrInp(MaHo)       ! Switch for determining shrinkage curve (-)
      integer SwSoilShr(MaHo)      ! Switch for kind of soil for determining shrinkage curve (-)
      integer SwSorp(MaHo)         ! Switch for kind of sorptivity function (-) 
      real(8) DiPoMa               ! Maximal diameter soil polygones (deep) (L) 
      real(8) DiPoMi               ! Minimal diameter soil polygones (shallow) (L) 
      real(8) CritUndSatVol        ! Critical value for undersaturation volume (L)
      real(8) GeomFac(MaHo)        ! Geometry factor for (an)isotropic shrinkage (-)
      real(8) PndmxMp              ! Threshold value for ponding (L) on soil surface before overland flow into macropores starts    
      real(8) PowM                 ! Power M for frequency distribut. curve IC domain (-)   
      real(8) PpIcSs               ! Proportion of IC domain at Soil Surface (-)
      real(8) RapDraReaExp         ! Reaction coefficient for rapid drainage (-)
      real(8) RapDraResRef(Madr)   ! Reference rapid drainage resistance (T)  
      real(8) Rzah                 ! Fraction macropores ended at bottom A-horizon (-)
      real(8) ShapeFacMp           ! Shape factor for description of macropore water exchange with matrix (-)
      real(8) SorpAlfa(MaHo)       ! Fitting parameter for emperical sorptivity curve (-)
      real(8) SorpMax(MaHo)        ! Maximal sorptivity at theta residual (L/T^0.5)
      real(8) SorpFacParl(MaHo)    ! Factor for modifying Parlange function (-)
      real(8) ShrParA(MaHo)        ! Parameter 1 for describing shrinkage curve (depending on SWSoilShr and SwShrInp)
      real(8) ShrParB(MaHo)        ! Parameter 2 for describing shrinkage curve (depending on SWSoilShr and SwShrInp)
      real(8) ShrParC(MaHo)        ! Parameter 3 for describing shrinkage curve (depending on SWSoilShr and SwShrInp)
      real(8) ShrParD(MaHo)        ! Parameter 4 for describing shrinkage curve (depending on SWSoilShr and SwShrInp)
      real(8) ShrParE(MaHo)        ! Parameter 5 for describing shrinkage curve (depending on SWSoilShr and SwShrInp)
      real(8) Spoint               ! Symmetry Point for freq. distr. curve (-)
      real(8) ThetCrMp(MaHo)       ! Critical volumetric water content below which cracks are formed (L^3/L^3)
      real(8) VlMpStSs             ! Volume of Static Macropores at Soil Surface (L^3/L^3)
      real(8) Z_Ah                 ! Depth bottom A-horizon (L)
      real(8) Z_Ic                 ! Depth bottom Internal Catchment (IC) domain (L)
      real(8) Z_MB50               ! Depth where static volume of MB domain has decreased to half of original(L) ! Adaptation 4 for PEARL-MACRO  
      real(8) Z_St                 ! Depth bottom Static macropores (L)
      real(8) ZDiPoMa              ! Depth below which diameter of soil polygons is maximum (L)
      real(8) ZnCrAr               ! Depth at which crack area of soil surface is calculated (L)
! --- macropore variables
      integer ICpBtDmPot(MaDm)     ! Compartment number containing potential bottom depth of domain (-) 
      integer IDecMpRat            ! Counter for number of times macropore fluxes are decreased with factor 10 because convergence is not reached (-)
      integer NodGWlFlCpZo         ! Node directly above groundwater level of full capillary zone (-)
      integer NumDm                ! Number of macropore domains (-)     
      integer NumLevRapDra         ! Number of drainage level 
      integer SwBma                ! Switch for output file with detailed yearly macropore water balance *.BMA: 0 = no; 1 = yes
      integer SwMacro              ! Switch for simulation of macropore flow: 0 = no; 1 = yes
      real(8) ArMpSs               ! Area fraction of macropores at soil surface (-) 
      real(8) cQMpLatSs            ! Cumulative amount of macropore inflow at soil surface by lateral overland flow (L)
      real(8) cQMpInIntSatDm1      ! Cumulative amount of interflow out off perched groundwater into macropores of domain 1 (MB) (L)
      real(8) cQMpInIntSatDm2      ! Cumulative amount of interflow out off perched groundwater into macropores of domain 2 (IC) (L)
      real(8) cQMpInMtxSatDm1      ! Cumulative amount of exfiltration out off saturated matrix into macropores of domain 1 (MB) (L)
      real(8) cQMpInMtxSatDm2      ! Cumulative amount of exfiltration out off saturated matrix into macropores of domain 2 (IC) (L)
      real(8) cQMpInTopLatDm1      ! Cumulative amount of lateral overland flow into macropores of domain 1 (MB) (L)
      real(8) cQMpInTopLatDm2      ! Cumulative amount of lateral overland flow into macropores of domain 2 (IC) (L)
      real(8) cQMpInTopPreDm1      ! Cumulative amount of direct precipiation into macropores of domain 1 (MB) (L)
      real(8) cQMpInTopPreDm2      ! Cumulative amount of direct precipiation into macropores of domain 2 (IC) (L)
      real(8) cQMpOutDrRap         ! Cumulative amount of rapid drainage out off macropores of domain 1 (MB) (L)
      real(8) cQMpOutMtxSatDm1     ! Cumulative amount of infiltration into saturated matrix out off macropores of domain 1 (MB) (L) 
      real(8) cQMpOutMtxSatDm2     ! Cumulative amount of infiltration into saturated matrix out off macropores of domain 2 (IC) (L)
      real(8) cQMpOutMtxUnsDm1     ! Cumulative amount of infiltration into unsaturated matrix out off macropores of domain 1 (MB) (L) 
      real(8) cQMpOutMtxUnsDm2     ! Cumulative amount of infiltration into unsaturated matrix out off macropores of domain 2 (IC) (L) 
      real(8) DiPoCp(MaCp)         ! Diameter of soil matrix polygon per compartment (L)
      real(8) dFdhMp(MaCp)         ! Contribution of macropores to derivative of compartment (1/T)
      real(8) dtold                ! Length of previous Time step (T)
      real(8) QMpLatSs             ! Macropore inflow flux at soil surface by lateral overland flow (L/T)
      real(8) FrArMtrx(MaCp)       ! Fraction of horizontal area of soil matrix per compartment (-)
      real(8) GWlFlCpZo            ! Groundwater level of full capillary zone (L) (only unsaturated zones with less than CritUndSatVol air)
      real(8) IAvFrMpWlWtDm1(MaCp) ! Incremental sum of average wet macropore wall fraction weighted for time step, for domain 1 (MB) (-)
      real(8) IAvFrMpWlWtDm2(MaCp) ! Incremental sum of average wet macropore wall fraction weighted for time step, for domain 2 (IC) (-)
      real(8) iQExcMtxDm1Cp(MaCp)  ! Incremental amount of water exchange between matrix and macropores of domain 1 (MB) (L)
      real(8) iQExcMtxDm2Cp(MaCp)  ! Incremental amount of water exchange between matrix and macropores of domain 2 (IC) (L)
      real(8) iQInTopLatDm1        ! Incremental amount of lateral overland flow into macropores of domain 1 (MB) (L)
      real(8) iQInTopLatDm2        ! Incremental amount of lateral overland flow into macropores of domain 2 (IC) (L)
      real(8) iQInTopPreDm1        ! Incremental amount of direct precipiation into macropores of domain 1 (MB) (L)
      real(8) iQInTopPreDm2        ! Incremental amount of direct precipiation into macropores of domain 2 (IC) (L)
      real(8) iQMpOutDrRap         ! Incremental amount of rapid drainage out off macropores of domain 1 (MB) (L)
      real(8) iQOutDrRapCp(MaCp)   ! Incremental amount of rapid drainage out off macropores per compartment (L)
      real(8) IWaSrDm1Beg          ! Incremental amount of water storage in macropores of domain 1 (MB) at beginning of period (L)
      real(8) IWaSrDm2Beg          ! Incremental amount of water storage in macropores of domain 2 (IC) at beginning of period (L)
      real(8) KsMpSs               ! Vertical hydraulic conductivity of macropores at soil surface (L/T) 
      real(8) PpDmCp(MaDm,MaCp)    ! Volumetric proportion of macropore domain per compartment (-)
      real(8) QExcMpMtx(MaCp)      ! Water exchange flux between matrix and macropores per compartment (L/T) 
      real(8) QInTopLatDm1         ! Lateral overland flow flux into macropores of domain 1 (MB) (L/T)
      real(8) QInTopLatDm2         ! Lateral overland flow flux into macropores of domain 2 (IC) (L/T)
      real(8) QInTopPreDm1         ! Direct precipiation flux into macropores of domain 1 (MB) (L/T)
      real(8) QInTopPreDm2         ! Direct precipiation flux into macropores of domain 2 (IC) (L/T)
      real(8) QMaPo                ! Total exchange flux between matrix and macropores (L/T)
      real(8) QRapDra              ! Total rapid drainage flux (L/T)
      real(8) SubsidCp(MaCp)       ! Vertical subsidence of matrix per compartment (L)
      real(8) VlMp                 ! Total macropore volume (L)
      real(8) VlMpDm1              ! Total macropore volume of domain 1 (MB) (L)
      real(8) VlMpDm2              ! Total macropore volume of domain 2 (IC) (L)
      real(8) VlMpDyCp(MaCp)       ! Dynamic macropore volume per compartment (L)
      real(8) VlMpStCp(MaCp)       ! Static macropore volume per compartment (L)
      real(8) VlMpStDm1(MaCp)      ! Static macropore volume of domain 1 (MB) per compartment (L)
      real(8) VlMpStDm2(MaCp)      ! Static macropore volume of domain 2 (IC) per compartment (L)
      real(8) WaLevDm1             ! Water level in domain 1 (MB) (L)
      real(8) WaSrDm1              ! Water storage in domain 1 (MB) (L)
      real(8) WaSrDm1Ini           ! Initial water storage in domain 1 (MB) (L)
      real(8) WaSrDm2              ! Water storage in domain 2 (IC) (L)
      real(8) WaSrDm2Ini           ! Initial water storage in domain 2 (IC) (L)
      real(8) ZDraBas              ! Level of drainage basis (drain depth or surface water level) for rapid drainage calculations (L)
      logical FlDecMpRat           ! Flag indicating decrease of macropore fluxes when convergence is not reached
      logical flInitDraBas         ! Flag indicating initialization of basis for rapid drainage
      logical flmacropore          ! Flag indicating simulation of macropore flow

! --- surface water variables
      integer swswb,swdrf,swsrf,swallo(Madr),swdtyp(Madr),swnrsrf
      integer swqhr,swsec,nrpri,nrsec,nmper,swman(mamp),SwTopnrsrf
      integer nqh(mamp),drf,swb,nphase(mamp),nodhd(mamp),numadj
      integer intwl(mamp),imper
      real(8) widthr(Madr),taludr(Madr),rdrain(Madr),rsurfdeep
      real(8) rsurfshallow,rinfi(Madr),rentry(Madr),rexit(Madr)
      real(8) gwlinf(Madr),wlptab(2*mawlp)
      real(8) impend(mamp)
      real(8) wldip(mamp),wscap(mamp),hbweir(mamp)
      real(8) osswlm,wlstar,wlp,alphaw(mamp),betaw(mamp)
      real(8) dropr(mamp*mamte),hdepth(mamp*mamte)
      real(8) gwlcrit(mamp,mamte),hcrit(mamp,mamte),vcrit(mamp,mamte)
      real(8) hqhtab(mamp,mamte),qqhtab(mamp,mamte)
      real(8) wlsman(mamp,mamte)
      real(8) wlstab(2*mawls),sttab(22,2)
      real(8) swstini,swst,wlsbak(4)
      real(8) cofintfl,expintfl,cqdrd,cwsupp,cwout,wls
      real(8) owltab(Madr,2*maowl),hwlman,wlsold,qdrd
      logical flsurfacewater,overfl

      end module variables
