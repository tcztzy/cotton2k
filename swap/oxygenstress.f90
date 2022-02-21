! File VersionID:
!   $Id: oxygenstress.f90 298 2016-07-25 20:07:31Z kroes006 $
!
! ----------------------------------------------------------------------
      subroutine OxygenStress(node,rwu_factor,ResultsOxStr) 
! ----------------------------------------------------------------------
!     Last modified      : January 2014              
!     Purpose            : calculates oxygen stress according to Bartholomeus et al. (2008)
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local
      integer glit,lay,node,i,j
      real(8) resp_factor,rwu_factor,xi,accuracy
      real(8) c_min_micro, c_macro,depth
      real(8) theta0
      real(8) matric_potential, gas_filled_porosity
      real(8) sat_water_cont,res_water_cont,gen_n,alpha,gen_m
      real(8) perc_org_mat,soil_density,percentage_sand,w_root_z0
      real(8) soil_temp,air_temp,root_radius,w_root
      real(8) shape_factor_microbialr,r_microbial_z0
      real(8) d_o2inwater,d_root,d_gassfreeair,o2_atmosphere
      real(8) surface_tension_water,bunsencoeff
      real(8) waterfilm_thickness,SOLVE,pi
      real(8) soilphystab(5,matab)
      real(8) diff_water_cap_actual
      real(8) theta100,theta500,campbell_b,d_soil,gfp100
      integer numrec_tab
      real(8) watcon   
      real(8) rdepth,rdens,cumrdens,mdepth,afgen !RB20140110
      integer node_temp
      real(8) rdepth_top,rdens_top !RB20140114 
          
      real(8)   ResultsOxStr(19,macp)    ! result oxygen stress; tabulated for each model compartment
      
      parameter (pi = 3.1415926535d0)
      
!     save values of locals
      save   cumrdens,waterfilm_thickness,r_microbial_z0,d_soil

! --- Get max_resp_factor, i.e. the ratio between total respiration and maintenance respiration
      if (node.eq.1) call GET_MAX_RESP_FACTOR(max_resp_factor)
      !write(*,*) max_resp_factor

! --- initialize
      c_min_micro = 0.1d0
      c_macro = 0.2d0
      resp_factor = 1.d0

! --- soil layer in SWAP
      lay = layer(node)

! --- dry weight of root per unit length of root [kg/m]
      w_root = 1.0d0/SRL
! --- root radius [m]
      if (swrootradius .eq. 1) then
        root_radius = ((w_root/(pi*dry_mat_cont_roots*                  &
     &              (1-air_filled_root_por)*spec_weight_root_tissue))-  &
     &              (var_a))**0.5
      endif
      if (swrootradius .eq. 2) then
        root_radius = root_radiusO2
      endif

! --- RB20140117 get wofost parameters
      if ((croptype(icrop) .eq. 2).or.(croptype(icrop) .eq. 3)) then
          q10_root = q10
          c_mroot = rmr*32.0d0/30.0d0 !CH2O --> O2
      endif
      if (croptype(icrop) .eq. 2) then
          f_senes=afgen(rfsetb,30,dvs)    
      endif
      if (croptype(icrop) .eq. 3) then
          f_senes=afgen(rfsetb,30,rid)    
      endif
   
! --- extract a number of variables from Swap for local use in module OxygenStress

! --- set soil density [kg m-3]      
      soil_density = bdens(lay) 
! --- set parameter n of soil hydraulic functions           
      gen_n = cofgen(6,node) 
! --- set saturated water content [-]      
      sat_water_cont = cofgen(2,node) 
! --- set residual water content [-]
      res_water_cont = cofgen(1,node) 
! --- set parameter alpha [1/Pa] of soil hydraulic functions, so divide main swap alpha by 100      
      alpha = cofgen(4,node)/100.0d0 
! --- set percentage organic matter [%]
      perc_org_mat = orgmat(lay)*100.0d0 
! --- set percentage sand in % of total soil            
      percentage_sand = (psand(lay)*(1.0d0-orgmat(lay)))*100.0d0 
! --- get soil moisture content as defined in further calculations within this routine [-]      
      theta0 = theta(node) 
! --- gas filled porosity      
      gas_filled_porosity = sat_water_cont-theta0    

! --- thickness of the soil compartment [m]
      depth = dz(node)/100.0d0   
! --- microbial respiration calculated from organic matter content in actual soil compartment; keep this value fixed
      shape_factor_microbialr = 0.9d0 
!RB20140109 start shape_factor_rootr added
! --- microbial respiration calculated from organic matter content in actual soil compartment; keep this value fixed
      shape_factor_rootr = 0.9d0 
!RB20140109 start shape_factor_rootr added
! --- temperature in the soil compartment [K]
      soil_temp = tsoil(node)+273.d0  
! --- dry weight of roots at nodal depth
!RB20140109 start new calculation of w_root_z0
!previous:  w_root_z0 = w_root_ss * exp(0.01*z(node)/shape_factor_rootr)  
!new:  
! --- static crop. w_root_z0 relative to value of top layer              
      if (croptype(icrop) .eq. 1) then
            rdepth_top = (-z(1)-0.5*dz(1))/rd
            rdens_top = afgen(rdctb,22,rdepth_top)
            rdepth = (-z(node)-0.5*dz(node))/rd
            rdens = afgen(rdctb,22,rdepth) 
            w_root_z0 = w_root_ss * rdens/rdens_top !static crop
      endif
! --- dynamic crop. wrt [kg/ha] = 10-4 kg/m2; 
      if ((croptype(icrop) .eq. 2) .or. (croptype(icrop) .eq. 3)) then
! --- estimate sum of densities; cumrdens
            if (node.eq.1) then !.and.(flDayStart)
                cumrdens = 0.0d0        
                do node_temp = 1,noddrz-1
                    rdepth = (-z(node_temp)-0.5*dz(node_temp))/rd
                    rdens = afgen(rdctb,22,rdepth)
                    cumrdens = cumrdens+rdens
                enddo
! --- last node, partly filled with roots
                node_temp = noddrz
                mdepth = (z(node_temp)+0.5*dz(node_temp)-rd)/2.d0 !mean depth
                rdepth = mdepth/(-rd)
                rdens = afgen(rdctb,22,rdepth)
                cumrdens = cumrdens+rdens
            endif  
                rdepth = (-z(node)-0.5*dz(node))/rd
                rdens = afgen(rdctb,22,rdepth)
                w_root_z0 = (0.0001*wrt/(0.01*rd))*rdens/cumrdens !dynamic crop (wofost)
      endif
!RB20140109 end new calculation of w_root_z0      

! --- Calculate matric potential [Pa]
       matric_potential = -100.d0 * h(node)          
! --- if gas filled porosity = 0, then root water uptake = 0. Store results and go to end of routine.        
      if (gas_filled_porosity .lt. 1.0d-4) then !RB20131106 .eq. 0
         ! RB 20140106 start if statement added; if max_resp_factor = 1 then no stress so rwufactor = 1        
          if (max_resp_factor .gt. 1.0d0) then
            rwu_factor = 0.0d0
          else
            rwu_factor = 1.0d0
          endif
          c_macro = 0.d0 !RB20140825 added
! --- store value for output results
        if (node.eq.1) ResultsOxStr = -999.d0 !RB20140114
          ResultsOxStr(1,node) = node
          ResultsOxStr(2,node) = z(node)*0.01d0
          ResultsOxStr(3,node) = h(node)
          ResultsOxStr(4,node) = matric_potential
          ResultsOxStr(5,node) = gas_filled_porosity
          ResultsOxStr(6,node) = -999.d0
          ResultsOxStr(7,node) = -999.d0
          ResultsOxStr(8,node) = rwu_factor
! --- Rpotz0 kg/m3/d
          ResultsOxStr(9,node)= f_senes * (c_mroot * w_root_z0 *        &
     &      max_resp_factor) *                                          &
     &      ( q10_root ** ( (soil_temp - 298.0d0) / 10.0d0 ) )
! --- RpotLay kg/m2/d; integrate over function potresp_top_of_compartment*exp(x/shape_factor_rootR); analytical solution
          ResultsOxStr(10,node)= -1.0d0*((ResultsOxStr(9,node)*         &
     &      shape_factor_rootR*exp(-depth/shape_factor_rootr))-         &
     &      (ResultsOxStr(9,node)*shape_factor_rootr*                   &
     &      exp((0.0d0)/shape_factor_rootr)))
! --- Ractz0 kg/m3/d
          if (resp_factor.gt.max_resp_factor) then
            resp_factor = max_resp_factor
          endif
          if (resp_factor.lt.1.0d0) then
            resp_factor = 0.0d0
          endif      
          ResultsOxStr(11,node)= f_senes * (c_mroot * w_root_z0 *       &
     &      resp_factor) *                                              &
     &      ( q10_root ** ( (soil_temp - 298.0d0 ) / 10.d0 ) ) 
! --- RactLay kg/m2/d; integrate over function actresp_top_of_compartment*exp(x/shape_factor_rootR); analytical solution
        ResultsOxStr(12,node)= -1.d0*((ResultsOxStr(11,node)*           &
     &      shape_factor_rootR*exp(-depth/shape_factor_rootr))-         &
     &      (ResultsOxStr(11,node)*shape_factor_rootr*                  &
     &    exp((0.d0)/shape_factor_rootr)))
! --- RredLay kg/m2/d
          ResultsOxStr(13,node)= ResultsOxStr(10,node)-                 &
     &      ResultsOxStr(12,node)    
! --- RsoilActLay kg/m2/d     
          ResultsOxStr(14,node)= -999.d0
! --- dsoil
          ResultsOxStr(15,node) = -999.d0     
! --- c_macro
          ResultsOxStr(16,node)=-999.d0 !kg/m3
          ResultsOxStr(17,node)=-999.d0       
! --- Max resp factor RB20140115
          ResultsOxStr(18,node)=max_resp_factor 
! --- wrt from wofost
          ResultsOxStr(19,node)=-999.d0 !kg/ha          
      
      else   ! (if (gas_filled_porosity .lt. 1.0d-6)) 
        
! --- In case of tabular soil hydraulic functions
        if(swsophy.eq.1) then
! --- Get tabular soil hydraulic function for node
          do i = 1,5
            do j = 1, numtablay(lay)  !check this
              soilphystab(i,j) = sptab(i,node,j)
            end do
          end do
! --- Get differential water capacity at actual node, (/L --> /Pa)
          diff_water_cap_actual = 0.01d0*dimoca(node)
          numrec_tab = j-1
        endif
! --- get theta at h=-100 cm and at h=-500 cm; for diffusion coef
        theta100 = watcon(node,-100.d0,cofgen,swsophy,numtab,sptab,     &
     &                    ientrytab)
        theta500 = watcon(node,-500.d0,cofgen,swsophy,numtab,sptab,     &
     &                    ientrytab)    

! --- atmosphere oxygen concentration [kg/m3] according to general gas law
        if (node.eq.1) then
! ---   atmospheric temperature [K]      
          air_temp = tav + 273.0d0
          o2_atmosphere = (672.0d0) / (8.314472d0 * air_temp)
! ---   PLAATS O2_atmosphere IN DE VECTOR VOOR C_TOP 
          C_top(1) = o2_atmosphere
        endif    

! --- Calculate temperature dependent parameters
        call TEMP_DEPENDENT_PARAMETERS (d_o2inwater,                    &
     &          d_root,d_gassfreeair,                                   &
     &          surface_tension_water,bunsencoeff,soil_temp)

! --- Calculate diffusivity Dsoil
        gfp100 = sat_water_cont - theta100  
        campbell_b = (log10(500.0d0)-log10(100.0d0))/                   &
     &               (log10(theta100)-log10(theta500))
        d_soil = d_gassfreeair*(2.0d0*(gfp100**3)+0.04d0*gfp100)*       &
     &          ((gas_filled_porosity/gfp100)**(2.0d0+3.0d0/campbell_b))

! --- Calculate the thickness of the water film that surrounds the roots      
        call waterfilmthickness (waterfilm_thickness,                   &
     &   matric_potential,sat_water_cont,res_water_cont,                &
     &       alpha,gen_n,gen_m,surface_tension_water,glit,              &
     &          soilphystab,diff_water_cap_actual,numrec_tab)

! --- Calculate microbial respiration rate
        call microbial_resp (r_microbial_z0,soil_temp,perc_org_mat,     &
     &   soil_density,percentage_sand,matric_potential,                 &
     &   specific_resp_humus,q10_microbial)

! --- Calculate sink term variable
! --- Define input for the solving procedure
        xi = 0.5d0*max_resp_factor
        accuracy = 1.0d-4
! --- Calculate actual respiration factor from the solving procedure
        resp_factor =  SOLVE(xi,accuracy,                               &
     &                 c_mroot,w_root,w_root_z0,f_senes,q10_root,       &
     &                 soil_temp,                                       &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro,                                     &
     &                 c_macro,                                         &
     &                 c_top(node),                                     &
     &                 shape_factor_rootr,r_microbial_z0,               &
     &                 max_resp_factor,d_soil)

! --- PLAATS C_MACRO IN DE VECTOR VOOR C_TOP. BEREKENDE WAARDE IS INPUT VOOR VOLGENDE COMPARTIMENT
        C_top(node+1) = C_macro

! --- Calculate the sink term (Root Water Uptake) variable due to oxygen stress.
! --- The decrease in root water uptake is assumed proportional to the decrease
! --- in respiration (given by the maximum and actual respiration factor).

! RB 20140106 start if statement added; if max_resp_factor = 1 then no stress so rwufactor = 1        
          if (max_resp_factor .gt. 1.0d0) then
            rwu_factor = (1.0d0/(max_resp_factor-1.0d0)) * resp_factor  &
     &                 - (1.0d0/(max_resp_factor-1.0d0))
          else
            rwu_factor = 1.0d0
          endif
! RB 20140106 end if statement added

          if (rwu_factor .gt. 1.0d0) then
             rwu_factor = 1.0d0
          endif
          if (rwu_factor .lt. 0.d0) then
              rwu_factor = 0.d0
          endif
      
      endif !if (gas_filled_porosity .lt. 1.0d-6) !RB20131216 goto removed
      
! --- store value for output results
      if (node.eq.1) ResultsOxStr = -999.0d0 !RB20140114
      ResultsOxStr(1,node) = node
      ResultsOxStr(2,node) = z(node)*0.01d0
      ResultsOxStr(3,node) = h(node)
      ResultsOxStr(4,node) = matric_potential
      ResultsOxStr(5,node) = gas_filled_porosity
      ResultsOxStr(6,node) = w_root_z0
      ResultsOxStr(7,node) = waterfilm_thickness
      ResultsOxStr(8,node) = rwu_factor
! --- Rpotz0 kg/m3/d
      ResultsOxStr(9,node)= f_senes * (c_mroot * w_root_z0 *            &
     &    max_resp_factor) *                                            &
     &    ( q10_root ** ( (soil_temp - 298.0d0 ) / 10.0d0 ) )
! --- RpotLay kg/m2/d; integrate over function potresp_top_of_compartment*exp(x/shape_factor_rootR); analytical solution
      ResultsOxStr(10,node)= -1.0d0*((ResultsOxStr(9,node)*             &
     &  shape_factor_rootR*exp(-depth/shape_factor_rootr))-             &
     &  (ResultsOxStr(9,node)*shape_factor_rootr*                       &
     &  exp((0.0d0)/shape_factor_rootr)))
! --- Ractz0 kg/m3/d
      if (resp_factor.gt.max_resp_factor) then
        resp_factor = max_resp_factor
      endif
      if (resp_factor.lt.1.0d0) then
        resp_factor = 0.0d0
      endif      
      ResultsOxStr(11,node)= f_senes * (c_mroot * w_root_z0 *           &
     &  resp_factor) *                                                  &
     &  ( q10_root ** ( (soil_temp - 298.0d0 ) / 10.0d0 ) ) 
! --- RactLay kg/m2/d; integrate over function actresp_top_of_compartment*exp(x/shape_factor_rootR); analytical solution
      ResultsOxStr(12,node)= -1.d0*((ResultsOxStr(11,node)*             &
     &  shape_factor_rootR*exp(-depth/shape_factor_rootr))-             &
     &  (ResultsOxStr(11,node)*shape_factor_rootr*                      &
     &  exp((0.d0)/shape_factor_rootr)))
! --- RredLay kg/m2/d
      ResultsOxStr(13,node)= ResultsOxStr(10,node)-                     &
     & ResultsOxStr(12,node)    
! --- RsoilActLay kg/m2/d     
      ResultsOxStr(14,node)= -1.d0*((r_microbial_z0*                    &
     &   shape_factor_microbialr*exp(-depth/shape_factor_microbialr))-  &
     &     (r_microbial_z0*shape_factor_microbialr*                     &
     &      exp((0.d0)/shape_factor_microbialr)))
! --- dsoil
      ResultsOxStr(15,node) = d_soil     
! --- c_macro
      ResultsOxStr(16,node)=c_macro !kg/m3
      ResultsOxStr(17,node)=100.d0*                                     &
     & c_macro/((0.032*1d5)/(8.314472d0*soil_temp)) !% or kPa       
! --- Max resp factor RB20140115
      ResultsOxStr(18,node)=max_resp_factor    
! --- wrt from wofost RBf20140120
      if ((croptype(icrop) .eq. 2).or.(croptype(icrop) .eq. 3)) then
        ResultsOxStr(19,node)=wrt    
      endif
      return
      end

! --- End of main module OxygenStress ---------------------------------------------------------------------------------------
      
      subroutine GET_MAX_RESP_FACTOR (max_resp_factor_gmrf)
      use Variables
      implicit none
! --- Procedure to derive max_resp_factor, 
! --- i.e. the ratio between total respiration and maintenance respiration [-]
! --- This ratio is either given in the input file (for a static crop) 
! --- or calculated from a series of equations taken from WOFOST (for a dynamic crop)

      real(8) dayl,sinld,cosld
      real(8) amax_gmrf,afgen,dtga_gmrf,pgass_gmrf,gass_gmrf
      real(8) rmres_gmrf,teff_gmrf,mres_gmrf,asrc_gmrf
      real(8) fr_gmrf,fl_gmrf,fs_gmrf,fo_gmrf   
      real(8) cvf_gmrf
      real(8) Froots, Rg_roots,Rm_roots,Max_resp_factor_gmrf        
      
! --- static crop        
! --- static crop: max_resp_factor is given in the input file
      if (croptype(icrop) .eq. 1) then
        max_resp_factor_gmrf = max_resp_factor
      endif !if (croptype(icrop) .eq. 1)

! --- dynamic crop: max_resp_factor is calculated following the procedure
! --- for the calculation of root maintenance respiration and root growth respiration as 
! --- used in WOFOST.
! --- dynamic crop, not grass 
      if (croptype(icrop) .eq. 2) then
! --- phenological development rate  
        call astro(daynr,lat,rad,dayl,daylp,sinld,cosld,difpp,          &
     &             atmtr,dsinbe)
! --- gross assimilation
        amax_gmrf = afgen (amaxtb,30,dvs)
! --- correction for sub-optimum average daytemperature
        amax_gmrf = amax_gmrf * afgen (tmpftb,30,tavd)
        call totass(dayl,amax_gmrf,eff,lai,kdif,rad,                    &
     &              difpp,dsinbe,sinld,cosld,dtga_gmrf)
! --- correction for low minimum temperature
        dtga_gmrf = dtga_gmrf * afgen (tmnftb,30,tmnr)
! --- potential assimilation in kg ch2o per ha
        pgass_gmrf = dtga_gmrf * 30.0d0/44.0d0
! --- water stress reduction of pgass to gass
! --- not for calculation of max_resp_factor as max_resp_factor only occurs without stress
      !reltr = 1.0d0
        gass_gmrf = pgass_gmrf !* reltr
! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration, based on actual plant state variables
        rmres_gmrf = (rmr*wrt+rml*wlv+rms*wst+rmo*wso)*                 &
     &            afgen(rfsetb,30,dvs)
        !teff_gmrf = q10**((tsoil(10)-25.0d0)/10.0d0) !TEMPORARY!!!! ONLY TO CHECK EFFECT OF USING TSOIL INSTEAD OF TAV
        teff_gmrf = q10**((tav-25.0d0)/10.0d0)      
        mres_gmrf = min(gass_gmrf,rmres_gmrf*teff_gmrf)
        asrc_gmrf = gass_gmrf - mres_gmrf
! --- partitioning factors
        fr_gmrf = afgen(frtb,30,dvs) !rid for grass, dvs for wofost
        fl_gmrf = afgen(fltb,30,dvs)
        fs_gmrf = afgen(fstb,30,dvs)
        fo_gmrf = afgen(fotb,30,dvs)
! --- dry matter increase, only part in which cvf is calculated
        cvf_gmrf = 1.0d0/((fl_gmrf /cvl+fs_gmrf /cvs+fo_gmrf/cvo)*      &
     &  (1.0d0-fr_gmrf)+fr_gmrf/cvr)
! --- cvf: factor used in wofost to calculate the increase in biomass (dmi) from the
! ---  net assimilation of the whole plant (asrc); dmi = cvf*asrc. 
! ---  What is left is the growth respiration (i.e. asrc = dmi + growth respiration). 
! ---  Therefore, growth respiration of the whole plant = asrc*(1-cvf)
! --- Froots: contribution of the roots to cvf
        Froots = (fr_gmrf/cvr)*cvf_gmrf
! --- Rg_roots: growth respiration roots        
        Rg_roots = Froots*(1.0d0-cvf_gmrf)*asrc_gmrf
! --- Rm_roots: maintenance respiration roots        
        Rm_roots = min(Froots*(1.0d0-cvf_gmrf)*gass_gmrf,               &
     &      rmr*wrt*afgen(rfsetb,30,dvs)*teff_gmrf)
! --- Max_resp_factor: ratio total respiration / maintenance respiration        
        if (Rm_roots.gt.0.0d0) then
            Max_resp_factor_gmrf = (Rg_roots+Rm_roots)/Rm_roots
        else 
            Max_resp_factor_gmrf = 1.0d0  
        endif         
      endif !if (croptype(icrop) .eq. 2)
        
! --- dynamic crop, grass 
      if (croptype(icrop) .eq. 3) then        
        Max_resp_factor_gmrf = 1.0d0  !RB20140317
! --- skip in case of regrowth, equal to wofost detailed grass
! --- note: daycrop.ge.idregrpot (wofost) --> daycrop.gt.idregrpot, because idregrpot is result of wofost of previous day
        if (daycrop.eq.0 .or.daycrop.gt.idregr                          &
     &                .or. seqgrazmow(iseqgm).eq.1) then           
! --- equations from cropgrowth.for subroutine grass      
! --- gross assimilation
          amax_gmrf = afgen (amaxtb,30,rid)
! --- correction for sub-optimum average daytemperature
          amax_gmrf = amax_gmrf * afgen (tmpftb,30,tavd)
! --- gross assimilation  
          call astro(daynr,lat,rad,dayl,daylp,sinld,cosld,difpp,        &
     &               atmtr,dsinbe)
          call totass(dayl,amax_gmrf,eff,lai,kdif,rad,                  &
     &                difpp,dsinbe,sinld,cosld,dtga_gmrf)
! --- correction for low minimum temperature
          dtga_gmrf = dtga_gmrf * afgen (tmnftb,30,tmnr)
! --- potential assimilation in kg ch2o per ha
          pgass_gmrf = dtga_gmrf * 30.0d0/44.0d0
! --- water stress reduction of pgass to gass
! --- not for calculation of max_resp_factor as max_resp_factor only occurs without stress
      !reltr = 1.0d0
          gass_gmrf = pgass_gmrf !* reltr
! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration, based on actual plant state variables
          rmres_gmrf = (rmr*wrt+rml*wlv+rms*wst)*afgen(rfsetb,30,rid) 
!        teff_gmrf = q10**((tsoil(10)-25.0d0)/10.0d0) !TEMPORARY!!!! ONLY TO CHECK EFFECT OF USING TSOIL INSTEAD OF TAV
          teff_gmrf = q10**((tav-25.0d0)/10.0d0)
          mres_gmrf = min(gass_gmrf,rmres_gmrf*teff_gmrf)
          asrc_gmrf = gass_gmrf - mres_gmrf
! --- partitioning factors
          fr_gmrf = afgen(frtb,30,rid) !rid for grass, dvs for wofost
          fl_gmrf = afgen(fltb,30,rid)
          fs_gmrf = afgen(fstb,30,rid)
! --- dry matter increase, only part in which cvf is calculated
          cvf_gmrf = 1.0d0/((fl_gmrf /cvl+fs_gmrf /cvs)*                &
     &        (1.0d0-fr_gmrf)+fr_gmrf/cvr)
! --- cvf: factor used in wofost to calculate the increase in biomass (dmi) from the
! ---  net assimilation of the whole plant (asrc); dmi = cvf*asrc. 
! ---  What is left is the growth respiration (i.e. asrc = dmi + growth respiration). 
! ---  Therefore, growth respiration of the whole plant = asrc*(1-cvf)
! --- Froots: contribution of the roots to cvf
          Froots = (fr_gmrf/cvr)*cvf_gmrf
! --- Rg_roots: growth respiration roots            
          Rg_roots = Froots*(1.0d0-cvf_gmrf)*asrc_gmrf
! --- Rm_roots: maintenance respiration roots     
          Rm_roots = min(Froots*(1.0d0-cvf_gmrf)*gass_gmrf,             &
     &        rmr*wrt*afgen(rfsetb,30,rid)*teff_gmrf)
! --- Max_resp_factor: ratio total respiration / maintenance respiration        
          if (Rm_roots.gt.0.d0) then
              Max_resp_factor_gmrf = (Rg_roots+Rm_roots)/Rm_roots
          else 
              Max_resp_factor_gmrf = 1.0d0  
          endif                  
        endif !RB20140317 #skip in case of regrowth                   
      endif !if (croptype(icrop) .eq. 3)

      return
      end
      
      subroutine TEMP_DEPENDENT_PARAMETERS (d_o2inwater,                &
     &          d_root,d_gassfreeair,                                   &
     &          surface_tension_water,bunsencoeff,soil_temp)

      implicit none
      real(8) d_o2inwater,d_root,d_gassfreeair
      real(8) surface_tension_water,bunsencoeff,soil_temp
! --- Calculate parameters that depend on temperature
! --- diffusion coefficient for oxygen in water [m2/d] Lango et al 1996
      d_o2inwater = 24.0d0*3600.0d0*((1.0d-5*1.2d0/(100.0d0*100.0d0))*  &
     &              exp(0.026d0 * (soil_temp - 273.0d0) ) )
! --- diffusion coefficient for oxygen in root tissue [m2/d]
! --- scaled to d_o2inwater, according the value given by
! --- van noordwijk & de willigen 1987 at 293 k
      d_root = 0.4d0*d_o2inwater
! --- diffusion coefficient for oxygen in free air [m2/d]
! --- hirschfelder et al 1964 molecular theory of gases and liquids
      d_gassfreeair = 1.74528d0 * (soil_temp**3)/(293.0d0**3)
! --- surface tension of water [n/m] eotvos rule
      surface_tension_water = 0.07275d0 *                               &
     &                      ( 1.0d0 - 0.002d0 * (soil_temp - 291.0d0 ) )
! --- bunsen solubility coeff for oxygen Lango et al 1996
      bunsencoeff =                                                     &
     &         (1413.d0*(exp(-144.397d0+7775.18d0*(soil_temp**(-1.d0))  &
     &              +18.3977d0*log(soil_temp)+0.0094437d0*soil_temp)))  &
     &              *273.15d0/soil_temp 
      return
      end

      real(8) function FUNC(x,sat_water_cont,res_water_cont,alpha,       &
     &                   gen_n,gen_m,surface_tension_water )

      implicit none
      real(8) x
      real(8) sat_water_cont,res_water_cont,alpha
      real(8) gen_n,gen_m,surface_tension_water
      real(8) pi
      parameter (pi = 3.1415926535897932d0)
! --- function needed for calculation of length density of gas filled pores in
! ---   subroutine 'water film thickness'                                          
      func= -(( -(sat_water_cont - res_water_cont) * alpha *            &    
     &      ( (alpha * x)**(gen_n - 1.d0) ) *                           &
     &      ( ( 1.d0 + (alpha * x)**gen_n )**(gen_m - 1.d0) )*          &
     &      gen_m*gen_n ) /                                             &
     &      (((( (alpha * x)**gen_n ) + 1.d0 )**gen_m)**2) ) /          &
     &      (pi*(4.d0*((surface_tension_water**2)/(x**2))))
      return
      end

      real(8) function FUNCtab(x,diff_water_cap,surface_tension_water )

      implicit none
      real(8) x
      real(8) diff_water_cap,surface_tension_water
      real(8) pi
      parameter (pi = 3.1415926535897932d0)
! --- function needed for calculation of length density of gas filled pores in
! ---   subroutine 'water film thickness'                                          
      functab= diff_water_cap /                                         &
     &      (pi*(4.d0*((surface_tension_water**2)/(x**2))))
      return
      end     

      subroutine TRAPZD(a,b,s,n,sat_water_cont,res_water_cont,alpha,    &
     &                   gen_n,gen_m,surface_tension_water,glit)
     
      implicit none
      real(8) FUNC
      real(8) a,b,s
      integer n
! --- procedure from numerical recipes. calculate integral numerically
! --- Needed for procedure 'water film thickness'                                
! --- programs calling trapzd must provide a function
! --- func(x:real(8)):real(8)which is to be integrated. they must
! --- also define the variable
! --- var glit: integer;
! --- in the main routine. *)
!RB20140312 sum replaced by mysum
      integer j,glit
      real(8) x,tnm,mysum,del
      real(8) sat_water_cont,res_water_cont,alpha
      real(8) gen_n,gen_m,surface_tension_water
     
! ---  initialize j
      j=0
      if (n .eq. 1) then
         s = 0.5d0*(b-a)*(func(a,sat_water_cont,res_water_cont,alpha,   &
     &                   gen_n,gen_m,surface_tension_water)             &
     &              +func(b,sat_water_cont,res_water_cont,alpha,        &
     &                   gen_n,gen_m,surface_tension_water))
         glit = 1
      else
         tnm = dfloat(glit)
         del = (b-a)/tnm
         x = a+0.5d0*del
         mysum = 0.0d0
         do while (j .lt. glit) 
            mysum = mysum+func(x,sat_water_cont,res_water_cont,alpha,   &
     &                     gen_n,gen_m,surface_tension_water)
            x = x+del
            j = j + 1
         enddo
         s = 0.5d0*(s+(b-a)*mysum/tnm)
         glit = 2*glit
      endif
      return
      end

      subroutine TRAPZDtab(a,b,s,n,diff_water_cap,                      &
     &                   surface_tension_water,glit)
     
      implicit none
      real(8) functab
      real(8) a,b,s
      integer n
! --- procedure from numerical recipes. calculate integral numerically
! --- Needed for procedure 'water film thickness'                                
! --- programs calling trapzd must provide a function
! --- func(x:real(8)):real(8)which is to be integrated. they must
! --- also define the variable
! --- var glit: integer;
! --- in the main routine. *)
!RB20140312 sum replaced by mysum
      integer j,glit
      real(8) x,tnm,mysum,del
      real(8) diff_water_cap,surface_tension_water
! ---  initialize j
      j=0
      if (n .eq. 1) then
         s = 0.5d0*(b-a)*(functab(a,diff_water_cap,                     &
     &                   surface_tension_water)                         &
     &              +functab(b,diff_water_cap,                          &
     &                   surface_tension_water))
         glit = 1
      else
         tnm = dfloat(glit)
         del = (b-a)/tnm
         x = a+0.5d0*del
         mysum = 0.0d0
         do while (j .lt. glit) 
            mysum = mysum+functab(x,diff_water_cap,                     &
     &                     surface_tension_water)
            x = x+del
            j = j + 1
         enddo
         s = 0.5d0*(s+(b-a)*mysum/tnm)
         glit = 2*glit
      endif
      return
      end

      subroutine waterfilmthickness (waterfilm_thickness,               &
     & matric_potential,sat_water_cont,res_water_cont,                  &
     &       alpha,gen_n,gen_m,surface_tension_water,glit,              &
     &       soilphystab,diff_water_cap_actual,numrec_tab)
! --- calculate water film thickness. method according to simojoki 2000
      use Variables
      implicit none
      
      real(8) waterfilm_thickness, matric_potential,sat_water_cont
      real(8) res_water_cont, alpha,gen_n,gen_m,surface_tension_water
      real(8) soilphystab(5,matab)
      real(8) diff_water_cap_actual
      real(8) lowlim,upplim,diff_water_cap
      real(8) length_density_gas_pores_sub
      integer glit    
      integer i,ii,numrec_tab
      real(8) new,old,s
      real(8) length_density_gas_pores
      
      save s

      real(8) pi
      parameter (pi = 3.1415926535897932d0)
      
      if (swsophy.eq.0) then
! --- m from van genuchten parameters
        gen_m = 1-(1/gen_n)
! --- calculate length density air filled (gas) pores. [number per m2]
! --- subroutine trapdz is used to solve the integral defined in
! --- function 'func'.
        i = 1
        call TRAPZD(1.d-10,matric_potential,s,i,sat_water_cont,         &
     &    res_water_cont,alpha,gen_n,gen_m,surface_tension_water,glit)
        new = s
        old = new + new

        do while (abs(new-old) .gt. 0.00001*new) 
          i = i + 1
          call TRAPZD(1.d-10,matric_potential,s,i,sat_water_cont,       &
     &        res_water_cont,alpha,gen_n,gen_m,                         &
     &        surface_tension_water,glit)
          old = new
          new = s
        enddo
! --- result from trapzd and func:      
        length_density_gas_pores = new
! --- calculated water film thickness [m]
        waterfilm_thickness = 2.d0 *                                    &
     &   ( ( ( 1 / ( pi * length_density_gas_pores ) ) ** 0.5d0 )       &
     &   - 2.d0 * surface_tension_water / matric_potential )
      endif
      !!!!sptab(1,node,ii) = soilphystab(1,ii)
      if(swsophy.eq.1) then
          ii = numrec_tab !34  !run over points in soil hydraulic table
! --- lower value of matric potential for integration interval
         lowlim = 0.000001d0 !-100*soilphystab(1,ii) !initial value should be zero !RB20140725, very close to zero, otherwise/0 in functab
! --- upper value of matric potential for integration interval
         upplim = -100.d0*0.5d0*(soilphystab(1,ii)+soilphystab(1,ii-1)) !middle of points 1 and 2
         upplim = MIN(upplim, matric_potential) !upplim can be higher than matpot          
! --- initialize length density gas filled pores
         length_density_gas_pores = 0.0d0 !initial value
! --- start loop, i.e. run over full integration interval;
! --- in each loop, the length density of gas filled pores of that interval is calculated;
! --- The sum of all loops gives the total length density.
! --- The integration interval runs from 0 to the matric potential at the actual node
         do while ((upplim.lt.matric_potential) .and. (ii.gt.2)) !adjusted RB20150714 ii.gt.2 added
            diff_water_cap = 0.01d0*soilphystab(4,ii) !at matric potential centre of lowlim and upplim, except for first point
            call TRAPZDtab(lowlim,upplim,s,1,                           & !n=1 as integral is over a linear function, i.e. convergence is reached after one step
     &                     diff_water_cap,                              &
     &                     surface_tension_water,glit)
            length_density_gas_pores_sub = s
            length_density_gas_pores = length_density_gas_pores +       &
     &                                      length_density_gas_pores_sub
            ii = ii - 1
            lowlim = upplim
            upplim=-100.d0*0.5d0*(soilphystab(1,ii)+soilphystab(1,ii-1))
          enddo
! --- in the last step (upplim >= matric_potential), always use the differential water capacity corresponding to the actual matric_potential of the node
          diff_water_cap = diff_water_cap_actual 
          upplim = MIN(upplim, matric_potential)
          call TRAPZDtab(lowlim,upplim,s,1,                             &  !n=1 as integral is over a linear function, i.e. convergence is reached after one step
     &                   diff_water_cap,surface_tension_water,glit)
          length_density_gas_pores_sub = s
          length_density_gas_pores = length_density_gas_pores +         &
     &                                      length_density_gas_pores_sub
! --- calculated water film thickness [m]
          waterfilm_thickness = 2 *                                     &
     &     ( ( ( 1.d0 / ( pi * length_density_gas_pores ) ) ** 0.5d0 )  &
     &       - 2.d0 * surface_tension_water / matric_potential )
      endif
      !ResultsOxStr(7,node) = waterfilm_thickness !RB20140115
      return
      end

      subroutine microbial_resp (r_microbial_z0,soil_temp,perc_org_mat, &
     & soil_density,percentage_sand,matric_potential,                   &
     & specific_resp_humus,q10_microbial)
! --- Calculate microbial respiration rate, based on available amount of
! --- organic matter, soil moisture and temperature

      implicit none
      real(8) r_microbial_z0,soil_temp,perc_org_mat
      real(8) soil_density,percentage_sand,matric_potential
      real(8) specific_resp_humus,q10_microbial
      real(8) f_moisture_humus,t_humus,carbon_humuspools
      real(8) saturated_matric_potential
      
! --- humus temperature [K]
      t_humus = soil_temp
! --- available amount of organic carbon [kg/m3]
      carbon_humuspools = 0.48d0 * ( perc_org_mat/100.d0 ) *soil_density
! --- saturated matric potential [pa] calculated according to cosby et al. 1984
! --- * 100: cm -> pa
      saturated_matric_potential =(10.d0**(-0.0131d0*percentage_sand+   &
     &                            1.88d0))*100
! --- reduction function for soil moisture:
      if ( matric_potential .lt. saturated_matric_potential) then
         f_moisture_humus = 0.5d0
      endif
      if (( matric_potential .ge. saturated_matric_potential)           &
     &                      .and. (matric_potential .lt. 25000.d0)) then
         f_moisture_humus = 1.d0 - 0.5d0 *                              &
     &                      ((log10(25000.d0)-log10(matric_potential))/ &
     &                      (log10(25000.d0)-                           &
     &                      log10(saturated_matric_potential) ) )
      endif
      if ((matric_potential .ge. 25000d0)                               &
     &                    .and. (matric_potential .le. 762500d0)) then
         f_moisture_humus = 1.0d0
      endif
      if ((matric_potential .gt. 762500d0)                              &
     &                    .and. (matric_potential .le. 1500000d0)) then
         f_moisture_humus = 1.d0 -                                      &
     &                     ((log10(matric_potential)-log10(762500.d0))/ &
     &                     ( log10(1500000.d0 ) -  log10(762500.d0 ) ) )
      endif
      if (matric_potential .gt. 1500000.d0) then
         f_moisture_humus = 0.0d0
      endif
! --- microbial respiration rate at soil surface [kg/m3/d]
      r_microbial_z0 = specific_resp_humus * carbon_humuspools *        &
     &             ( q10_microbial**( (t_humus-298.d0) / 10.d0 ) )      &
     &             * f_moisture_humus
      return
      end

      subroutine MICRO (c_mroot,w_root,f_senes,q10_root,                &
     &                 soil_temp,sat_water_cont,gas_filled_porosity,    &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro, resp_factor)
! --- calculate minimum oxygen concentration in the gas phase of the soil
! --- that is needed to provide all cells within a plant root with a
! --- sufficient amount of oxygen.                                                                     *)

      implicit none
      real(8) c_min_micro,c_mroot,w_root,f_senes,q10_root,soil_temp
      real(8) sat_water_cont,gas_filled_porosity
      real(8) d_o2inwater,d_root,perc_org_mat,soil_density
      real(8) specific_resp_humus,q10_microbial,depth
      real(8) shape_factor_microbialr,root_radius
      real(8) waterfilm_thickness,bunsencoeff
      real(8) resp_factor

      real(8) r_mref,r_mroot,waterfilm_porosity,d_waterfilm,lambda
      real(8) r_waterfilm_lengthroot,alpha_alpha
      real(8) f_moisture_humus,t_humus,carbon_humuspools
      real(8) c_min_micro_interphase

      real(8) r_microbial_volumetric_wf ,r_microbial_z0_wf
      real(8) pi
      parameter (pi = 3.1415926535897932d0)

! --- calc respiration per unit length of root [kg/m/d]
      r_mref = c_mroot * w_root * resp_factor

      r_mroot = f_senes * r_mref *                                      &
     &          ( q10_root ** ( (soil_temp - 298.d0 ) / 10.d0 ) )
! --- calc porosity of water film [-]
      waterfilm_porosity = ( sat_water_cont - gas_filled_porosity ) /   &
     &                ( 1 - gas_filled_porosity )
! --- calc diffusion coeff of water film [m2/d]
      d_waterfilm = d_o2inwater * ( waterfilm_porosity**(4.0d0/3.0d0) )
! --- calc ratio of d_root and d_waterfilm, input for c_min_micro_interphase
      lambda = d_root / d_waterfilm
! --- start microbial respiration rate in water film
! --- temperature humus [k]
      t_humus = soil_temp
! --- amount of carbon in humus [kg/m3]
      carbon_humuspools = 0.48d0 * ( perc_org_mat/100.d0 ) *soil_density
! --- reduction factor for moisture [-] (saturation)
      f_moisture_humus = 0.5d0
! --- microbial respiration rate at soil surface [kg/m3/d]
      r_microbial_z0_wf = specific_resp_humus * carbon_humuspools *     &
     &                 ( q10_microbial**( (t_humus-298.d0) / 10.d0 ) ) *&
     &                 f_moisture_humus
! --- volumetric microbial resp rate at depth z [kg/m3/d]
      r_microbial_volumetric_wf = r_microbial_z0_wf*                    &
     &                          exp(-depth/shape_factor_microbialr)
! --- respiration rate in water film per unit length of root [kg/m/d]
      r_waterfilm_lengthroot = pi *                                     &
     &                     ( (root_radius + waterfilm_thickness)**2 -   &
     &                     root_radius**2 )*r_microbial_volumetric_wf
! --- end microbial respiration rate in water film*)

! --- calculation procedure following de willigen & van noordwijk 1983
! --- p.220-221
! --- ratio of rhizosphere respiration to total respiration [-]
      alpha_alpha = r_waterfilm_lengthroot /                            &
     &              (r_waterfilm_lengthroot+r_mroot)
! --- volumetric respiration rate of the root + rhizosphere, but attributed
! ---   to the root (following de willigen & van noordwijk 1983) [kg/m3/d]
      c_min_micro_interphase = ((r_waterfilm_lengthroot + r_mroot)/     &
     &          (2.d0*pi*d_root)) *                                     &
     &          ( 0.5d0 + ( ( lambda - 1.d0 ) * alpha_alpha / 2.d0 ) +  &
     &          lambda * log(1.d0 + waterfilm_thickness / root_radius) -&
     &          ( lambda * alpha_alpha *                                &
     &          ( 1.d0 + waterfilm_thickness / root_radius )**2 ) *     &
     &          log(1.d0 + waterfilm_thickness / root_radius ) /        &
     &          ( ( waterfilm_thickness / root_radius ) *               &
     &          (2.d0 + waterfilm_thickness / root_radius) )  )

      c_min_micro = c_min_micro_interphase/bunsencoeff
! --- end calculation procedure following de willigen & van noordwijk 1983*)
      return
      end

      subroutine MACRO (c_macro,depth,resp_factor,                      &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)
     
! --- calculate oxygen concentration in the soil at a specific depth with respect
! --- to the soil surface. diffusion from atmosphere into soil is considered,
! --- including microbial and root respiration.

      implicit none
      real(8) c_macro, depth, resp_factor
      real(8) c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop
      real(8) shape_factor_microbialr,shape_factor_rootr,r_microbial_z0

      real(8) dum,r_mref_z0,r_mroot_z0,d_soil,l,lnew,lnew_initial,fi,   &
     &        fi_a
      integer counterMacro, counterMacroSub
      character(len=200) error_messag

! --- reference respiration at soil surface per unit volume of roots [kg/m3] *)
      r_mref_z0 = c_mroot * w_root_z0 * resp_factor

! --- respiration at soil surface per unit volume of roots [kg/m3]
      r_mroot_z0 = f_senes * r_mref_z0 *                                &
     &             ( q10_root ** ( (soil_temp - 298.d0 ) / 10.d0 ) )

! --- calculate oxygen concentration at certain depth (c_macro [kg/m3]
! --- in the soil. method according to cook 1995
! --- calculate criterium (dum) to switch to specific equation (if then else)
      dum = (( shape_factor_microbialr**2) * r_microbial_z0 /d_soil)+   &
     &      ( ( shape_factor_rootr**2 ) * r_mroot_z0 / d_soil )
! --- as z goes to infinity, the oxygen concentration asymptotically approaches
! ---    a constant non-zero value
      if (dum .lt. ctop) then
         c_macro = ctop -                                               &
     &     ( (shape_factor_microbialr**2) * r_microbial_z0 / d_soil)*   &
     &     ( 1.d0- exp( - depth / shape_factor_microbialr ) ) -         &
     &     ( (shape_factor_rootr**2) * r_mroot_z0 / d_soil ) *          &
     &     ( 1.d0- exp( - depth / shape_factor_rootr ) )
! --- at z = l the oxygen concentration goes to zero. find l through
! ---    newton-raphson method
      else
         counterMacro = 0
         counterMacroSub = 0
         lnew_initial = 0.1d0 !initialize
         fi = 1.d0 !initialize
             fi_a = 1.d0 !initialize
         lnew = lnew_initial
         do while ((abs(fi).gt.1.d-8) .AND.                             &
     &                 (dabs(fi_a).gt.0.d0) .AND.                       &
     &             (lnew.gt.1.d-4))
            counterMacro = counterMacro + 1
            counterMacroSub = counterMacroSub + 1
            l = lnew
            fi = ctop -                                                 &
     &         ((shape_factor_microbialr**2)*r_microbial_z0/d_soil)*    &
     &         ( 1.d0- (l/shape_factor_microbialr)*                     &
     &         exp( - l / shape_factor_microbialr ) -                   &
     &         exp( - l / shape_factor_microbialr ) ) -                 &
     &         ( (shape_factor_rootr**2) * r_mroot_z0 / d_soil ) *      &
     &         (1.d0-(l/shape_factor_rootr)*exp(-l/shape_factor_rootr )-&
     &         exp( - l / shape_factor_rootr ) )
! --- derivative of fi to l
            fi_a = -r_microbial_z0/d_soil*l*                            &
     &             exp(-l/shape_factor_microbialr)-                     &
     &             r_mroot_z0/d_soil*l*exp(-l/shape_factor_rootr)
               
            if (dabs(fi_a) .gt. 0.d0) then
                lnew = abs(l - (fi / fi_a))
            endif
            
! --- for convergence --> depends on initial value of lnew
            if (lnew.gt.1.d3                                            &
     &           .or.                                                   &
     &           counterMacroSub.gt.100) then
               !restart do while loop with new lnew_initial value
               lnew_initial = lnew_initial + 0.1d0
               lnew = lnew_initial
               counterMacroSub = 0
            endif
! ---       fatal error if too many iterations --> 
            if (counterMacro .gt. 1.d6) then
              error_messag = '1 Too much iterations for macroscopic '   &
     &               //' oxygen diffusion.'
!D              call warn ('rootextraction',error_messag,logf,swscre)
              call fatalerr ('rootextraction',error_messag)
            endif            
         enddo
        
         if (depth .lt. l) then
            c_macro = ctop -                                            &
     &         ((shape_factor_microbialr**2)*r_microbial_z0/d_soil )*   &
     &         ( 1.d0- (depth/shape_factor_microbialr)*                 &
     &         exp( - l / shape_factor_microbialr ) -                   &
     &         exp( - depth / shape_factor_microbialr ) ) -             &
     &         ( (shape_factor_rootr**2) * r_mroot_z0 / d_soil ) *      &
     &         ( 1.d0- (depth/shape_factor_rootr)*                      &
     &         exp( - l / shape_factor_rootr ) -                        &
     &         exp( - depth / shape_factor_rootr ) )
         else
            c_macro = 0.0d0
         endif
      endif
      return
      end

      real(8) function SOLVE(xi,accuracy,                               &
     &                 c_mroot,w_root,w_root_z0,f_senes,q10_root,       &
     &                 soil_temp,                                       &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro,                                     &
     &                 c_macro,ctop,                                    &
     &                 shape_factor_rootr,r_microbial_z0,               &
     &                 max_resp_factor,                                 &
     &                 d_soil)     
     
      implicit none
      
      real(8) xiplus1,ximin1
      real(8) delta
      real(8) xi,accuracy
      real(8) c_mroot,w_root,w_root_z0,f_senes,q10_root,soil_temp
      real(8) sat_water_cont,gas_filled_porosity
      real(8) d_o2inwater,d_root,perc_org_mat,soil_density
      real(8) specific_resp_humus,q10_microbial,depth
      real(8) shape_factor_microbialr,root_radius
      real(8) waterfilm_thickness,bunsencoeff
      real(8) c_min_micro
      real(8) c_macro,ctop 
      real(8) shape_factor_rootr,r_microbial_z0
      real(8) max_resp_factor
      real(8) d_soil
      
      real(8) xx,fxi,xiplusdelta,ximindelta,fxiplusdelta,fximindelta,   &
     &        dfxi
      integer counterSolve
      character(len=200) error_messag

! --- iterative procedure to find the RESPIRATION FACTOR for which holds !RB20131106
! --- that c_macro = c_micro. reference value = c_micro
! --- Newton-Raphson method

! --- speed up simulations: cut of if ctop = 0 and if waterfilm_thickness is extremely high (due to very low gas filled porosity)     !RB20140826
      if (abs(ctop) .lt. 1.0d-6 .OR.                                    &
     & waterfilm_thickness .gt. 1.d0) then 
         xx = 0.d0
         C_macro = 0.0d0
         SOLVE = xx
         return
      endif

      delta= 1.d-8
      fxi = 100.d0
      counterSolve = 1
      do while (fxi .gt. accuracy)
         call MICRO (c_mroot,w_root,f_senes,q10_root,soil_temp,         &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro, xi)
         call MACRO (c_macro,depth,xi,                                  &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)
         fxi = abs(C_min_micro-C_macro)
         xiplusdelta = xi + delta
         ximindelta = xi - delta
         
         if (dabs(xiplusdelta-ximindelta).lt.1.d-20) then !if (dabs(xiplusdelta-ximindelta).lt.1.d-12) then !if (xiplusdelta .eq. ximindelta) then 
            xiplusdelta = xi + 1.d-6
            ximindelta = xi - 1.d-6
         endif
         
         call MICRO (c_mroot,w_root,f_senes,q10_root,soil_temp,         &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro, xiplusdelta)
         call MACRO (c_macro,depth,xiplusdelta,                         &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)
     
         fxiplusdelta = abs(C_min_micro-C_macro)
         call MICRO (c_mroot,w_root,f_senes,q10_root,soil_temp,         &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro, ximindelta)
         call MACRO (c_macro,depth,ximindelta,                          &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)    
     
         fximindelta = abs(C_min_micro-C_macro)  

         dfxi = (fxiplusdelta - fximindelta) / (xiplusdelta-ximindelta)
         
        if (dabs(dfxi) .lt. 1.d-20) then !if (dabs(dfxi) .lt. 1.d-15) then !RB20140312
          xiplus1 = xi
        else
          xiplus1 = max(1.0d-8,(xi - fxi/dfxi)) !xiplus1 = max(1.0d-6,(xi - fxi/dfxi)) !RB20131106, never lt 0
        endif

      !prevent endless sign change without conversion !RB20140827
      if (mod(counterSolve,2)>1.d-12) then !if even than store xi
        ximin1 = xi
      endif
!D      write(*,*) counterSolve,xi,ximin1,xiplus1 !DEBUG
      if (abs(xiplus1-ximin1) .lt. 1.d-6) then  !no change in xi value, than take average of both values between which is iterated
        xx = (xi+xiplus1)*0.5d0
        call MACRO (c_macro,depth,xx,                                   &
     &       c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop,         &
     &       shape_factor_microbialr,shape_factor_rootr,                &
     &       r_microbial_z0,d_soil)     
        SOLVE = xx
        return
      endif

! --- speed up simulations: cut of if RWU_factor will be >1 or <1e-6 (=0)       
         if (xi .lt. 1.0d-6 .AND. xi .gt. 0.0d0) then !RB20131106 greater than 0 is required
            xx = 0.d0
            C_macro = 0.0d0
            SOLVE = xx
            return
         endif
         if (xi .gt. max_resp_factor) then
            xx = max_resp_factor
            call MACRO (c_macro,depth,xx,                               &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)     
            SOLVE = xx
            return
         endif
      
         xi = xiplus1
         
! ---    error if max iterations reached--> 
         counterSolve = counterSolve + 1
         if (counterSolve .gt. 100) then
           error_messag = '1 Max iterations for solver oxygen'          &
     &            //' stress reached.'
            xx = xi
           SOLVE = xx
           call fatalerr ('rootextraction',error_messag)
           return          
        endif            
      enddo
      xx = xi
      SOLVE = xx
      return
      end

! ----------------------------------------------------------------------
      subroutine oxygen_dat (SwTopSub,NrStaring,OxygenSlope,            &
     &                       OxygenIntercept)
! ----------------------------------------------------------------------
!     date               : december 2009
!     purpose            : set parameter values of metafunction for oxygenstress
! ----------------------------------------------------------------------

! --- global
      integer SwTopSub,NrStaring
      real(8) OxygenSlope(6),OxygenIntercept(6)

! --- local
      integer i
      real(8) xtop(18,6),xsub(18,6),ytop(18,6),ysub(18,6)

      data (xtop(1,i),i=1,6)                                            &
     & /5.07d-03,2.40d+02,-4.39d+00,-6.31d+02,1.67d+00,9.08d+02/
      data (xtop(2,i),i=1,6)                                            &
     & /1.20d-02,3.26d+02,-8.91d+00,-9.29d+02,2.49d+00,1.64d+03/
      data (xtop(3,i),i=1,6)                                            &
     & /1.21d-02,3.64d+02,-9.09d+00,-9.90d+02,2.60d+00,1.70d+03/
      data (xtop(4,i),i=1,6)                                            &
     & /1.67d-02,4.42d+02,-1.21d+01,-1.26d+03,3.34d+00,2.20d+03/
      data (xtop(5,i),i=1,6)                                            &
     & /4.17d-03,1.93d+02,-3.72d+00,-5.28d+02,1.44d+00,7.77d+02/
      data (xtop(6,i),i=1,6)                                            &
     & /2.11d-02,5.02d+02,-1.51d+01,-1.40d+03,3.69d+00,2.71d+03/
      data (xtop(7,i),i=1,6)                                            &
     & /2.86d-02,5.57d+02,-1.97d+01,-1.67d+03,4.50d+00,3.41d+03/
      data (xtop(8,i),i=1,6)                                            &
     & /1.75d-02,5.55d+02,-1.32d+01,-1.34d+03,3.31d+00,2.48d+03/
      data (xtop(9,i),i=1,6)                                            &
     & /2.34d-02,6.00d+02,-1.70d+01,-1.40d+03,3.42d+00,3.09d+03/
      data (xtop(10,i),i=1,6)                                           &
     & /3.12d-02,6.62d+02,-2.20d+01,-1.53d+03,3.65d+00,3.91d+03/
      data (xtop(11,i),i=1,6)                                           &
     & /2.58d-02,6.42d+02,-1.83d+01,-1.61d+03,4.00d+00,3.26d+03/
      data (xtop(12,i),i=1,6)                                           &
     & /2.50d-02,6.53d+02,-1.79d+01,-1.45d+03,3.40d+00,3.21d+03/
      data (xtop(13,i),i=1,6)                                           &
     & /2.53d-02,5.97d+02,-1.79d+01,-1.72d+03,4.58d+00,3.18d+03/
      data (xtop(14,i),i=1,6)                                           &
     & /2.82d-02,7.10d+02,-2.04d+01,-1.54d+03,3.56d+00,3.71d+03/
      data (xtop(15,i),i=1,6)                                           &
     & /2.08d-02,4.72d+02,-1.46d+01,-1.37d+03,3.65d+00,2.57d+03/
      data (xtop(16,i),i=1,6)                                           &
     & /1.99d-02,4.80d+02,-1.39d+01,-1.34d+03,3.47d+00,2.45d+03/
      data (xtop(17,i),i=1,6)                                           &
     & /2.27d-02,6.31d+02,-1.62d+01,-1.58d+03,3.91d+00,2.91d+03/
      data (xtop(18,i),i=1,6)                                           &
     & /2.23d-02,6.50d+02,-1.60d+01,-1.74d+03,4.45d+00,2.89d+03/

      data (xsub(1,i),i=1,6)                                            &
     & /7.21d-04,1.76d+02,-1.59d+00,-4.33d+02,1.16d+00,4.52d+02/
      data (xsub(2,i),i=1,6)                                            &
     & /3.75d-03,2.25d+02,-3.61d+00,-5.96d+02,1.59d+00,7.91d+02/
      data (xsub(3,i),i=1,6)                                            &
     & /7.06d-03,2.86d+02,-5.91d+00,-7.52d+02,2.00d+00,1.19d+03/
      data (xsub(4,i),i=1,6)                                            &
     & /1.29d-02,3.74d+02,-9.74d+00,-1.03d+03,2.72d+00,1.82d+03/
      data (xsub(5,i),i=1,6)                                            &
     & /2.92d-03,1.86d+02,-3.06d+00,-5.07d+02,1.39d+00,6.85d+02/
      data (xsub(6,i),i=1,6)                                            &
     & /3.00d-02,5.41d+02,-2.06d+01,-1.69d+03,4.61d+00,3.55d+03/
      data (xsub(7,i),i=1,6)                                            &
     & /2.76d-02,6.63d+02,-1.95d+01,-1.67d+03,4.15d+00,3.48d+03/
      data (xsub(8,i),i=1,6)                                            &
     & /1.85d-02,4.88d+02,-1.34d+01,-1.34d+03,3.50d+00,2.43d+03/
      data (xsub(9,i),i=1,6)                                            &
     & /2.08d-02,5.46d+02,-1.50d+01,-1.50d+03,3.92d+00,2.72d+03/
      data (xsub(10,i),i=1,6)                                           &
     & /1.85d-02,5.94d+02,-1.39d+01,-1.45d+03,3.60d+00,2.60d+03/
      data (xsub(11,i),i=1,6)                                           &
     & /2.29d-02,6.33d+02,-1.67d+01,-1.65d+03,4.21d+00,3.05d+03/
      data (xsub(12,i),i=1,6)                                           &
     & /2.60d-02,5.89d+02,-1.83d+01,-1.33d+03,3.12d+00,3.26d+03/
      data (xsub(13,i),i=1,6)                                           &
     & /3.36d-02,6.27d+02,-2.29d+01,-1.40d+03,3.26d+00,3.95d+03/
      data (xsub(14,i),i=1,6)                                           &
     & /3.84d-02,6.74d+02,-2.71d+01,-1.40d+03,3.21d+00,4.83d+03/
      data (xsub(15,i),i=1,6)                                           &
     & /2.25d-02,6.16d+02,-1.66d+01,-1.37d+03,3.22d+00,3.05d+03/
      data (xsub(16,i),i=1,6)                                           &
     & /1.88d-02,4.75d+02,-1.32d+01,-1.29d+03,3.31d+00,2.33d+03/
      data (xsub(17,i),i=1,6)                                           &
     & /2.19d-02,5.40d+02,-1.53d+01,-1.50d+03,3.87d+00,2.70d+03/
      data (xsub(18,i),i=1,6)                                           &
     & /1.98d-02,5.00d+02,-1.41d+01,-1.40d+03,3.67d+00,2.52d+03/

      data (ytop(1,i),i=1,6)                                            &
     & /1.89d-04,-3.91d-01,-8.65d-02,2.11d+01,-8.61d-02,7.12d+00/
      data (ytop(2,i),i=1,6)                                            &
     & /5.02d-05,-7.56d-02,-1.01d-02,1.80d+01,-7.27d-02,-2.81d+00/
      data (ytop(3,i),i=1,6)                                            &
     & /5.36d-06,3.54d-01,1.20d-02,1.73d+01,-7.09d-02,-5.51d+00/
      data (ytop(4,i),i=1,6)                                            &
     & /-2.67d-05,4.87d-02,3.03d-02,1.76d+01,-7.02d-02,-7.90d+00/
      data (ytop(5,i),i=1,6)                                            &
     & /2.60d-04,-1.74d+00,-1.16d-01,2.28d+01,-8.99d-02,9.67d+00/
      data (ytop(6,i),i=1,6)                                            &
     & /-1.23d-04,-6.11d-02,8.42d-02,1.66d+01,-6.63d-02,-1.51d+01/
      data (ytop(7,i),i=1,6)                                            &
     & /-1.02d-04,8.54d-03,7.12d-02,1.76d+01,-7.06d-02,-1.33d+01/
      data (ytop(8,i),i=1,6)                                            &
     & /-8.09d-05,2.23d-01,5.65d-02,1.62d+01,-6.61d-02,-1.08d+01/
      data (ytop(9,i),i=1,6)                                            &
     & /-9.68d-05,1.56d-01,6.56d-02,1.67d+01,-6.76d-02,-1.21d+01/
      data (ytop(10,i),i=1,6)                                           &
     & /-1.37d-04,-4.78d-02,8.94d-02,1.81d+01,-7.09d-02,-1.55d+01/
      data (ytop(11,i),i=1,6)                                           &
     & /-1.71d-04,-2.25d-01,1.08d-01,1.68d+01,-6.46d-02,-1.78d+01/
      data (ytop(12,i),i=1,6)                                           &
     & /-1.36d-04,-5.72d-01,8.92d-02,1.02d+01,-3.91d-02,-1.54d+01/
      data (ytop(13,i),i=1,6)                                           &
     & /-1.19d-04,-8.68d-03,8.28d-02,1.77d+01,-6.97d-02,-1.54d+01/
      data (ytop(14,i),i=1,6)                                           &
     & /-9.91d-05,-9.32d-01,7.13d-02,1.35d+01,-5.10d-02,-1.35d+01/
      data (ytop(15,i),i=1,6)                                           &
     & /-7.90d-05,2.75d-02,5.85d-02,1.55d+01,-6.20d-02,-1.16d+01/
      data (ytop(16,i),i=1,6)                                           &
     & /-7.91d-05,2.32d-01,5.71d-02,1.10d+01,-4.40d-02,-1.12d+01/
      data (ytop(17,i),i=1,6)                                           &
     & /-2.21d-04,-5.23d-01,1.39d-01,1.44d+01,-5.42d-02,-2.27d+01/
      data (ytop(18,i),i=1,6)                                           &
     & /-1.16d-04,-3.60d-01,7.85d-02,1.05d+01,-3.96d-02,-1.41d+01/

      data (ysub(1,i),i=1,6)                                            &
     & /4.17d-04,-1.41d+00,-2.06d-01,2.42d+01,-9.84d-02,2.20d+01/
      data (ysub(2,i),i=1,6)                                            &
     & /2.56d-04,-4.40d-01,-1.22d-01,2.19d+01,-8.92d-02,1.16d+01/
      data (ysub(3,i),i=1,6)                                            &
     & /1.93d-04,-2.84d-01,-8.88d-02,1.96d+01,-8.08d-02,7.68d+00/
      data (ysub(4,i),i=1,6)                                            &
     & /5.69d-05,-7.44d-02,-1.48d-02,1.80d+01,-7.36d-02,-2.01d+00/
      data (ysub(5,i),i=1,6)                                            &
     & /5.21d-04,-1.81d+00,-2.61d-01,1.96d+01,-7.85d-02,3.01d+01/
      data (ysub(6,i),i=1,6)                                            &
     & /-1.10d-04,6.88d-02,7.90d-02,1.77d+01,-7.11d-02,-1.48d+01/
      data (ysub(7,i),i=1,6)                                            &
     & /-1.76d-04,-3.81d-01,1.11d-01,1.79d+01,-6.90d-02,-1.85d+01/
      data (ysub(8,i),i=1,6)                                            &
     & /-2.38d-05,5.72d-01,2.47d-02,1.64d+01,-6.82d-02,-6.50d+00/
      data (ysub(9,i),i=1,6)                                            &
     & /-4.68d-05,4.19d-01,3.84d-02,1.74d+01,-7.16d-02,-8.59d+00/
      data (ysub(10,i),i=1,6)                                           &
     & /-1.10d-04,1.03d-01,7.32d-02,1.69d+01,-6.77d-02,-1.32d+01/
      data (ysub(11,i),i=1,6)                                           &
     & /-1.48d-04,-3.36d-01,9.58d-02,1.78d+01,-6.94d-02,-1.64d+01/
      data (ysub(12,i),i=1,6)                                           &
     & /-1.58d-04,1.36d-01,9.92d-02,1.56d+01,-6.18d-02,-1.65d+01/
      data (ysub(13,i),i=1,6)                                           &
     & /-2.69d-04,-6.12d-01,1.66d-01,1.27d+01,-4.80d-02,-2.65d+01/
      data (ysub(14,i),i=1,6)                                           &
     & /-6.34d-05,-4.27d-01,5.05d-02,1.73d+01,-6.82d-02,-1.06d+01/
      data (ysub(15,i),i=1,6)                                           &
     & /3.08d-05,5.25d-02,-6.04d-03,8.44d+00,-3.54d-02,-2.00d+00/
      data (ysub(16,i),i=1,6)                                           &
     & /-6.20d-05,3.40d-01,4.76d-02,9.54d+00,-3.85d-02,-9.98d+00/
      data (ysub(17,i),i=1,6)                                           &
     & /-2.90d-05,3.99d-01,2.75d-02,8.07d+00,-3.30d-02,-6.73d+00/
      data (ysub(18,i),i=1,6)                                           &
     & /-6.76d-05,4.41d-01,5.01d-02,1.50d+01,-6.10d-02,-1.01d+01/

      if (SwTopSub .eq. 1) then
        do i = 1,6
          OxygenSlope(i) = xtop(NrStaring,i)
          OxygenIntercept(i) = ytop(NrStaring,i)
        enddo
      else
        do i = 1,6
          OxygenSlope(i) = xsub(NrStaring,i)
          OxygenIntercept(i) = ysub(NrStaring,i)
        enddo
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine OxygenReproFunction (OxygenSlope,OxygenIntercept,      &
     &   theta,thetas,tsoil,node,z,dz,rwu_factor)
! ----------------------------------------------------------------------
!     date               : January 2010
!     purpose            : Calculate oxygen stress according to reproduction function
! ----------------------------------------------------------------------
      implicit none
      include 'arrays.fi'

! --- global
      integer node,i
      real(8) OxygenSlope(6),OxygenIntercept(6),theta(macp),thetas(macp)
      real(8) tsoil(macp),z(macp),dz(macp)

! --- local
      real(8) intercept,slope,sum_porosity
      real(8) gas_filled_porosity
      real(8) soil_temp,depth_ss,mean_gas_filled_porosity
      real(8) rwu_factor

      gas_filled_porosity = thetas(node) - theta(node)
      soil_temp = tsoil(node) + 273.d0
      depth_ss = -z(node) * 0.01d0

      if (gas_filled_porosity .lt. 1.d-10) then
         rwu_factor = 0.d0 
         return
      endif

! --- mean gas filled porosity
      sum_porosity = 0.0d0
      do i = 1,node
        sum_porosity = sum_porosity +                                   &
     &                 (thetas(i) - theta(i)) * dz(node)
      enddo
      mean_gas_filled_porosity = sum_porosity /(-z(node)+0.5d0*dz(node))
             
      intercept = OxygenIntercept(1)*soil_temp**2 +                     &
     &  OxygenIntercept(2)*depth_ss**2 +                                &
     &  OxygenIntercept(3)*soil_temp +                                  &
     &  OxygenIntercept(4)*depth_ss +                                   &
     &  OxygenIntercept(5)*soil_temp*depth_ss + OxygenIntercept(6)

      slope = OxygenSlope(1)*soil_temp**2 +                             &
     &  OxygenSlope(2)*depth_ss**2 +                                    &
     &  OxygenSlope(3)*soil_temp +                                      &
     &  OxygenSlope(4)*depth_ss +                                       &
     &  OxygenSlope(5)*soil_temp*depth_ss + OxygenSlope(6)

! --- Calculate the sink term (Root Water Uptake) variable due to oxygen stress.
      rwu_factor = intercept + slope*mean_gas_filled_porosity
      if (rwu_factor .gt. 1.d0) then
         rwu_factor = 1.d0
      endif
      if (rwu_factor .lt. 0.d0) then
         rwu_factor = 0.d0
      endif
      return

      end
