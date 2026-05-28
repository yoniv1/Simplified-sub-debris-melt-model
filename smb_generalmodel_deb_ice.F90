 
!-----------------------------------------------------------------------            
! PROGRAM SMB                                        
!-----------------------------------------------------------------------        

!-----------------------------------------------------------------------
! Program: SMB model for Djankuat Glacier
!          Calculates clean-ice and debris-covered SMB using simplified
!          SEB and vertical debris heat conduction.
! Main configuration:
!   full_seb_on      = 0/1    (Include or exclude full energy balance calculation)
!   lin_temp_grad_on = 0/1    (Include or exclude internal heat storage within debris pack)
!   patchy_deb_on    = 0/1    (Include or exclude patchy debris cover for thin debris)
! Units:
!   debris thickness input: cm
!   temperatures: celsius 
!   precipitation: m 3h-1
!   insolation: W m-2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------                     
! VARIABLES                                               
!-----------------------------------------------------------------------   

!     Glacier specific parameters                                                                                                                                                                                           
!     ---------------------------                                                                                                                                                                                           

!     Number of rows/cols of DEM                                                                                                                                                                                            
      integer,parameter :: maxrow=150, maxcol=178

!     Numerical resolution                                                                                                                                                                                                  
      real,parameter    :: res=25.   

!     Timestep meteorological and insolation data                                                                                                                                                                           
      real,parameter    :: timestep_smb = 3.        ! Hours (temporal resolution forcing data)                                                                                                                              
      real,parameter    :: pi = 3.1415927           ! Redefined in subroutines                                                                                                                                                   
      real,parameter    :: sechr = 3600.            ! To convert seconds into hours                                                                                                                                              

!     SMB parameters for clean ice                                                                                                                                                                                          
      real,parameter    :: tau = 0.54               ! Transmissivity atmosphere (-)                                                                                                                                            
      real,parameter    :: lat = 43.20              ! Latitude (deg)                                                                                                                                                             
      real,parameter    :: lon = 42.76              ! Longitude (deg)                                                                                                                                                            
      real,parameter    :: altT = 2141.             ! Altitude Terskol (m ASL)                                                                                                                                                      
      real,parameter    :: precRate1 = 0.0023       ! Precipitation lapse rate (m w.e. yr-1 m-1)                                                                                                                                                                
      real,parameter    :: lapseRate_w = -0.0049    ! Temperature lapse rate winter (degC m-1)                                                                                                                                        
      real,parameter    :: lapseRate_s = -0.0067    ! Temperature lapse rate summer (degC m-1)                                                                                                                                        
      real,parameter    :: snowdepchar = 0.011      ! Characteristic snow depth (m w.e.)                                                                                                                                             
      real,parameter    :: ttresh = 2.              ! Threshold snow rain (degC)                                                                                                                                                   
      real,parameter    :: ttip = 0.                ! Threshold Ta-dependent fluxes (degC)                                                                                                                                        
      real,parameter    :: albsnow = 0.90           ! Albedo of snow (-)                                                                                                                                                       
      real,parameter    :: c0 = -45.                ! Ta-dependent fluxes intercept (W m-2)                                                                                                                                        
      real,parameter    :: c1 = 13.0                ! Ta-dependent fluxes slope (W m-2 degC-1)                                                                                                                                            
      real,parameter    :: albice = 0.35            ! Albedo of ice (-)
      real,parameter    :: albfirn = 0.50           ! Albedo of firn (-)                                                                                                                                                       
      real,parameter    :: eta = 0.15               ! Snowpack retention capacity of water (-)                                                                                                                                 
      real,parameter    :: tchar = 21.9             ! Timescale for decrease of snow albedo to firn albedo (days)                                                                                                                 

!     SMB parameters for debris-covered ice                                                                                                                                                                                 
      real,parameter     :: albdeb = 0.10            ! Albedo for debris (-)                                                                                                                                                    
      real,parameter     :: d0 = -60.                ! Ta-Ts-dependent fluxes intercept (W m-2)                                                                                                                                    
      real,parameter     :: d1 = 36.                 ! Ta-Ts-dependent fluxes slope (W m-2 degC-1)                                                                                                                                        
      integer, parameter :: N_layers = 10            ! Number of layers for vertical debris discretization (-)                                                                                                                 
      integer, parameter :: Nt = N_layers + 1        ! Vertical discretization debris pack (-)                                                                                                                                 
      real,parameter     :: vhc = 1.5E6              ! Effective volumetric heat capacity debris (J m-3 degC-1)                                                                                                                           
      real,parameter     :: ked = 1.8                ! Effective thermal conductivity debris (W m-1 degC-1)                                                                                                                               

!     SMB paramteres for full surface debris-covered ice energy balance                                                                                                                                                     
      real, parameter    :: rhoa = 1.29              ! Air density (kg m-3)                                                                                                                                                 
      real, parameter    :: cA = 1010.0              ! Specific heat capacity of air (J kg-1 degC-1)                                                                                                                           
      real, parameter    :: C_ex = 0.004             ! Bulk exchange coefficient over debris (-)                                                                                                                            
      real, parameter    :: Lv = 2.49e6              ! Latent heat of vaporization (J kg-1)                                                                                                                                 
      real, parameter    :: cW = 4184.0              ! Specific heat capacity of water (J kg-1 degC-1)                                                                                                                         
      real, parameter    :: rhow = 1000.0            ! Density of water (kg m-3)                                                                                                                                            
      real, parameter    :: em_s = 0.95              ! Debris surface emissivity (-)                                                                                                                                        
      real, parameter    :: stf_bltz = 5.6703744e-8  ! Stefan-Boltzmann constant (W m-2 K-4)                                                                                                                                

! ---------------------------------------------------------------------
! Choose the model complexity that is desired
! ---------------------------------------------------------------------

        lin_temp_grad_on = 0           ! Linear temperature gradient in the debris pack = 1, else 0
        full_seb_on = 0                ! Full energy balance calculation = 1, else 0
        patchy_deb_on = 0              ! Patchy debris cover parameterization = 1, else 0
             
!-----------------------------------------------------------------------     
!---Initialise matrices, set everything to zero                             
!-----------------------------------------------------------------------

        count_loop = 0
        shade = 0.
        insol = 0.
        snowdep = 0.
        albedo = 0.
        snowmelt = 0.
        snowret = 0.
        runoff = 0.
        ice_melt = 0.
        G = 0.
        cumbal = 0.
        yn = 0.                                                                         
        shading = 0.
        daysnow = 0.
        tsnow = 0.
        eng_flux = 0.
	snowflow = 0.
        qnet_flux_out = 0.
        qnet_deb = 0.
        ta_flx_deb = 0.
	enflux_deb = 0.
        qc_deb = 0.
        flx_deb = 0.
        dqnet_deb = 0.
        dta_flx_deb = 0.
        dqc_deb = 0.
        dflx_deb = 0.
	snowdep_deb = 0.
	albedo_deb = 0.
	snowmelt_deb = 0.
	snowret_deb = 0.
	runoff_deb = 0.
        G_deb = 0.
        ice_melt_deb = 0.
	cumbal_deb = 0.
	snowflow_deb = 0.
	daysnow_deb = 0.
        Td(:) = 273.15
        Td_past(:,:,:) = 273.15
	Td_gradient = 0.
        a_Crank(:) = 0.
        b_Crank(:) = 0.
        c_Crank(:) = 0.
	d_Crank(:) = 0.
	E_Crank(:) = 0.
	S_Crank(:) = 0.
	n_iterations = 0.
	vol_heat_cap_deb(:) = 0.
	k_eff_deb(:) = 0.
        Ts_past = 273.15
        lnet_deb = 0.
        shf_deb = 0.
        lhf_deb = 0.
        qrain_deb = 0.
        dlnet_deb = 0.
        dshf_deb = 0.
        dlhf_deb = 0.
        dqrain_deb = 0.
        qZ = 0.
        qS = 0.

!-----------------------------------------------------------------------            
! Run the SMB model loop over the years
!-----------------------------------------------------------------------   

       write(*,*), 'SMB calculations starting in year:, ', year_data(1,1)

       it = 1

       do while (it.lt.(dims+1))

           !-----------------------------------------------------------------------               
	   ! Calculate the solar geometry (incl. shading)
	   !-----------------------------------------------------------------------    

           h = hour_data(it,1)        ! Hour of the day
           day = day_data(it,1)       ! Number day of the year
           y = year_data(it,1)        ! Year
           
           call SOL(y,day,h,e_val,a_val,height_angle,decl,phi,eot,tbc,lambda,r_value)  ! Call subroutine to calculate sun-related variables                                                                         
           if (e_val.lt.0.) then
              e_val=0.
           end if
           
           call shaderoutine(solaralt,azimuth,HT,shade) ! Call subroutine to calculate if there is shade in each grid cell from topographic shading       
                     
           !-----------------------------------------------------------------------             
	   ! Start the loop over the area
	   !-----------------------------------------------------------------------   

           do j=1,maxcol
              do i=1,maxrow
                 
                 !-----------------------------------------------------------------------      
                 ! Start spatial SMB loop only over ice-covered pixels                    
                 !-----------------------------------------------------------------------    

                 if (mask_new(i,j).eq.1) then
                 
                 ! Initialize for new loop

                 insol=0.
                 melt=0.
		 melt_deb=0.
                 runoff=0.
		 runoff_deb=0.
                 prec=0.
                 Psolid=0.

                ! Get full surface energy balance variables                                                                                                                                                                                
                if(full_seb_on.eq.1)then
                    qsin = qsin_data(it,1)   ! In W m-2
                    pres = pres_data(it,1)   ! In hPa
                    uwind = uwind_data(it,1) ! In m/s
                    skyem = skyem_data(it,1) 
                    rha = rha_data(it,1)     ! In % (0-100) 
                 endif
                    
                 !--------------------------------------------------------------------              
                 ! Calculate insolation part of the energy balance
		 !--------------------------------------------------------------------

                    ! Calculate cos(incidence angle)
                    cos_theta = sin(declination) * sin(lat_rad) * cos(slope_rad(i,j)) - &
                         sin(declination) * cos(lat_rad) * sin(slope_rad(i,j)) * cos(aspect_inc(i,j)) + &
                         cos(declination) * cos(lat_rad) * cos(slope_rad(i,j)) * cos(hourangle) + &
                         cos(declination) * sin(lat_rad) * sin(slope_rad(i,j)) * cos(aspect_inc(i,j)) * cos(hourangle) + &
                         cos(declination) * sin(slope_rad(i,j)) * sin(aspect_inc(i,j)) * sin(hourangle)

                    ! Calculate the incidence angle in radians
                    inc_angle = acos(cos_theta)

                    ! Calculation solar radiation on sloping surface
		    ! if shade in grid cell and solar elevation angle > 0 : Qdif 
                    if ((e_val.gt.0.).and.(shade(i,j).eq.1).or.(inc_angle.ge.(pi_value/2.0))) then
                       insol= Qdif
                    ! Else if no shade and solar elevation angle > 0 : Qdir with cos_theta
                    elseif ((e_val.gt.0.).and.(shade(i,j).eq.0).and.(inc_angle.lt.(pi_value/2.0)))  then
                       insol= Qdir
                    endif
                    ! Else, solar elevation (e_val) <0 : Q = 0                                                                                                                                                         
                    if (e_val.le.0.or.insol.le.0.)then
                       insol=0.                                                                                                                                                 
                    endif

                    !--------------------------------------------------------------------       
		    ! Calculate air temperature
		    !--------------------------------------------------------------------  

                    ! Temperature 3-hourly values, depending on altitude of cell and season
                    if ((count_lpyear.eq.3).and.(day.ge.92.).and.(day.le.274.))then
                        Tair = temp_data(it,1)+lapseRate_s*(ELEV_OBS(i,j,idx_noyear)-altT)
                    elseif (count_lpyear.eq.3)then
                        Tair = temp_data(it,1)+lapseRate_w*(ELEV_OBS(i,j,idx_noyear)-altT)
                    elseif ((count_lpyear.lt.3).and.(day.ge.91.).and.(day.le.273.))then
                        Tair = temp_data(it,1)+lapseRate_s*(ELEV_OBS(i,j,idx_noyear)-altT)
                    elseif (count_lpyear.lt.3)then
                        Tair = temp_data(it,1)+lapseRate_w*(ELEV_OBS(i,j,idx_noyear)-altT)
                    endif

                    !--------------------------------------------------------------------  
		    ! Calculate precipitation
		    !--------------------------------------------------------------------  

                    ! Precipitation 3-hourly values, depending on altitude of cell         
                    if (ELEV_OBS(i,j,idx_noyear).gt.0.) then
			if (prec_yearly(it,1) > 0.0) then
   				Pratio = (prec_yearly(it,1) + precRate*(ELEV_OBS(i,j,idx_noyear)-altT)) / prec_yearly(it,1)
			else
   				Pratio = 1.0
			endif                       
			prec = prec_data(it,1)*Pratio
                    endif   

		    ! For safety, shouldn't happen

		    if (prec.lt.0)then
                        prec = 0
                    endif	
                    
                    !-------------------------------------------------------------------- 
		    ! Calculate solid precipitation (snow)
		    !--------------------------------------------------------------------  

                    if (Tair.lt.ttresh)then
                       Psolid = prec
                    else
                       Psolid = 0
                    endif

                    ! Set tsnow for snow and firn albedo

                    if (Psolid.gt.0)then
                       tsnow(i,j) = 0
                    endif
                    
		    ! For safety, shouldn't happen

		    if (Psolid.lt.0)then
                        Psolid = 0
                    endif	

                    !-------------------------------------------------------------------- 
		    ! Calculate albedo for snow, clean ice and debris-covered ice
		    !-------------------------------------------------------------------- 

		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    !!!!!! For snow and clean ice !!!!!!
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    albsnow_par = albfirn + (albsnow - albfirn)*exp(-tsnow(i,j)/tchar)
                    albedo = albsnow_par+(albice-albsnow_par)*exp(-snowdep(i,j)/snowdepchar)

		   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		   !!!!!! For sub-debris ice !!!!!!
		   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		    if (mask_deb(i,j).eq.1) then 
		    	albedo_deb = albsnow_par+(albdeb-albsnow_par)*exp(-snowdep_deb(i,j)/snowdepchar)
	            endif 

                    !-----------------------------------------------------------------------------                  
		    ! Get debris thickness-related variables and allocate matrices if appropriate
		    !-----------------------------------------------------------------------------
	
      	            if (mask_deb(i,j).eq.1)then
			if (a_deb(i,j) > 1.0e-6) then
   				debris_thickness = (th_deb(i,j)/100.0) / a_deb(i,j)
			else
   				debris_thickness = 0.0
			endif		      
		      h_deb = debris_thickness / N_layers
                      if (it.eq.starth) then
                        ! Initialize surface/internal debris properties
                        Td(:)=273.15
                        Td_past(:,:,:)=273.15
                      endif
		      ! Debris vertical properties
	              do jt = 1, Nt-1
                        vol_heat_cap_deb(jt) = vhc
                        k_eff_deb(jt) = ked
                      end do
                      ! Internal debris temperature of previous time step
                      if(it.eq.starth)then
                         do jt = 1, Nt
                            Td_past(i,j,jt)=273.15
                         end do
                      endif
		    end if

                    !--------------------------------------------------------------------------       
		    ! Calculate surface energy flux for snow, clean ice and debris-covered ice
		    !-------------------------------------------------------------------------- 

		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    !!!!!! For snow and clean ice !!!!!!
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                                           
                    if (Tair.ge.ttip) then
                       enflux = (1.-albedo)*tau*insol+c1*Tair+c0
                       ta_flx = c1*Tair+c0
                    else
                       enflux = (1.-albedo)*tau*insol+c0
                       ta_flx = c0
                    endif
                       qnet = (1.-albedo)*tau*insol

		   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		   !!!!!! For sub-debris ice !!!!!!
		   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		   if (mask_deb(i,j) .eq. 1) then  ! Check if in the debris region

    			if (snowdep_deb(i,j) .gt. 0) then  ! Snow present on debris
        
        			! Energy flux calculation based on air temperature
        			if (Tair .ge. ttip) then
            				enflux_deb = (1.0 - albedo_deb) * tau * insol + c1 * Tair + c0
            				ta_flx_deb = c1 * Tair + c0
        			else
            				enflux_deb = (1.0 - albedo_deb) * tau * insol + c0
            				ta_flx_deb = c0
        			end if
        
        			! Net shortwave radiation energy flux to debris
                                qnet_deb = (1.0 - albedo_deb) * tau * insol
           
                                ! For snow-covered conditions over debris, debris temperature is 0 degrees through the entire profile
                                Td(1:Nt)=273.15

                        elseif (snowdep_deb(i,j) .le. 0) then  ! No snow on debris (debris-covered ice surface)

        			! Initialize surface temperature iteration variables
        			n_iterations = 0
                                Ts_past = 273.15
                                Td(1:Nt-1) = 273.15 
        			Td(Nt) = 273.15   ! Ice temperature is 0 degrees Celsius

        			! Set initial guess for surface temperature
                                if (it .eq. starth)then
                                   Td(1) = 273.15
                                endif
                                
                                if(it.gt.starth.and.daysnow_deb(i,j) .eq. 1) then
            				Td(1) = Tair+273.15
        			elseif(it.gt.starth.and.daysnow_deb(i,j) .ne. 1)then
            				Td(1) = Td_past(i,j,1)
        			end if

        			! Heat diffusion by conduction: initial guess
        			if (it .eq. starth .or. daysnow_deb(i,j) .eq. 1) then
            				! Initial condition: linear temperature profile in debris
            				Td_gradient = (Td(1) - Td(Nt)) / debris_thickness
            				do jt = 2, Nt-1
                				Td(jt) = Td(1) - (jt * h_deb) * Td_gradient
            				end do

        			else
            			! Crank-Nicholson scheme for heat diffusion for first guess vertical temperature profile
            				do jt = 2, Nt-1

                                                ! Crank-Nicholson coefficients
                                                C_deb(jt) = k_eff_deb(jt) * (timestep_smb * sechr) / (2.0 * vol_heat_cap_deb(jt) * h_deb**2)
                				a_Crank(jt) = C_deb(jt)
                				b_Crank(jt) = 2.0 * C_deb(jt) + 1.0
                				c_Crank(jt) = C_deb(jt)

                				! Equation handling based on the current index
                				if (jt.eq.2) then
                    					d_Crank(jt) = C_deb(jt) * (Td(1) + Td_past(i,j,1) + Td_past(i,j,jt+1)) + &
                                    			(1.0 - 2.0 * C_deb(jt)) * Td_past(i,j,jt)
                				elseif (jt.lt.Nt-1) then
                    					d_Crank(jt) = C_deb(jt) * (Td_past(i,j,jt-1) + Td_past(i,j,jt+1)) + &
                                                        (1.0 - 2.0 * C_deb(jt)) * Td_past(i,j,jt)
                                                endif
                				if (jt.eq.Nt-1) then
                    					d_Crank(jt) = 2.0 * C_deb(jt) * Td(Nt) + &
                                     			C_deb(jt) * Td_past(i,j,Nt-2) + (1.0 - 2.0 * C_deb(jt)) * Td_past(i,j,Nt-1)
                				end if

                				! Forward substitution in Crank-Nicholson scheme
                				if (jt.eq.2) then
                    					E_Crank(jt) = b_Crank(jt)
                    					S_Crank(jt) = d_Crank(jt)
                				else
                    					E_Crank(jt) = b_Crank(jt) - a_Crank(jt) / E_Crank(jt-1) * c_Crank(jt-1)
                    					S_Crank(jt) = d_Crank(jt) + a_Crank(jt) / E_Crank(jt-1) * S_Crank(jt-1)
                				end if
            				end do

            				! Backward substitution
            				do jt = Nt-1, 2, -1
                				if (jt.eq.Nt-1) then
                    					Td(jt) = S_Crank(jt) / E_Crank(jt)
                				else
                    					Td(jt) = (S_Crank(jt) + c_Crank(jt) * Td(jt+1)) / E_Crank(jt)
                				end if
            				end do
        			end if

        			! Surface energy flux calculations for debris (first guess)
                                if(full_seb_on.eq.1)then
                                   qnet_deb = (1.0 - albedo_deb) * qsin ! (tau already included in input)
                                   lnet_deb = skyem*(stf_bltz*(Tair+273.15)**4) - em_s*(stf_bltz*(Td(1))**4)
                                   shf_deb = rhoa*cA*(C_ex)*uwind*(Tair+273.15-(Td(1)))
                                   qZ = (0.622 * (rha/100) * 6.112 * exp((17.67 * Tair) / (Tair + 243.5))) / &
                                        (pres - 0.378 * (rha/100) * 6.112 * exp((17.67 * Tair) / (Tair + 243.5)))
                                   if(prec.gt.0.and.Tair.gt.ttresh)then
                                      qS = (0.622 * 6.112 * exp((17.67 * (Td(1)-273.15)) / ((Td(1)-273.15) + 243.5))) / &
                                           (pres - 0.378 * 6.112 * exp((17.67 * (Td(1)-273.15)) / ((Td(1)-273.15) + 243.5)))
                                   else
                                      qS = qZ * (Td(1)/(Tair+273.15)) 
                                   endif
                                   lhf_deb = rhoa*Lv*(C_ex)*uwind*(qZ-qS)
                                   if(prec.gt.0.and.Tair.gt.ttresh)then
                                      qrain_deb = cW*rhow*((prec)/(timestep_smb*sechr))*((Tair+273.15) - Td(1))
                                   else
                                      qrain_deb = 0
                                   endif
                                else
                                  qnet_deb = (1.0 - albedo_deb) * tau * insol
                                  ta_flx_deb = d0 + d1 * ((Tair+273.15) - Td(1))
                                endif
                                if (lin_temp_grad_on.eq.1) then
                                   qc_deb = k_eff_deb(1) * (273.15 - Td(1)) / debris_thickness
                                else
                                   qc_deb = k_eff_deb(1) * (Td(2) - Td(1)) / h_deb
                                endif
                                if (full_seb_on.eq.1) then
                                   flx_deb = qnet_deb + lnet_deb + shf_deb + lhf_deb + qrain_deb + qc_deb
                                else
                                   flx_deb = qnet_deb + ta_flx_deb + qc_deb
                                endif

        			! Derivative of energy fluxes with respect to surface temperature (first guess)
                                if(full_seb_on.eq.1)then
                                   dqnet_deb = 0.
                                   dlnet_deb = -4*em_s*stf_bltz*Td(1)**3
                                   dshf_deb = -rhoa*cA*C_ex*uwind
                                   if (prec.gt.0.and.Tair.gt.ttresh) then
                                      dlhf_deb = -rhoa*Lv*C_ex*uwind*(0.622*pres*(6.112 * exp((17.67*(Td(1)-273.15))/((Td(1)-273.15)+243.5)) * &
                                           (17.67*243.5)/((Td(1)-273.15)+243.5)**2)) / (pres - 0.378 * 6.112 * exp((17.67*(Td(1)-273.15))/((Td(1)-273.15)+243.5)))**2
                                   else
                                      dlhf_deb = -rhoa*Lv*C_ex*uwind*(qZ/(Tair+273.15))
                                   endif
                                   if (prec.gt.0.and.Tair.gt.ttresh) then
                                      dqrain_deb = -cW*rhow*((prec)/(timestep_smb*sechr))
                                   else
                                      dqrain_deb = 0.
                                   endif
                                else
                                   dqnet_deb = 0.0
                                   dta_flx_deb = -d1
                                endif
                                if (lin_temp_grad_on.eq.1) then
                                   dqc_deb = -k_eff_deb(1) / debris_thickness
                                else
                                   dqc_deb = -k_eff_deb(1) / h_deb
                                endif
                                if (full_seb_on.eq.1) then
                                   dflx_deb = dqnet_deb + dlnet_deb + dshf_deb + dlhf_deb + dqrain_deb + dqc_deb
                                else
                                   dflx_deb = dqnet_deb + dta_flx_deb + dqc_deb
                                endif

                                ! Check for false skipping of iterations
                                if((it.gt.starth).and.(abs(Td(1) - Ts_past).le.0.01))then
                                     Ts_past = Ts_past + 0.03
                                endif

        			! Newton-Raphson method for solving surface temperature
        			do while (abs(Td(1) - Ts_past) > 0.01 .and. n_iterations < 100)
            				n_iterations = n_iterations + 1
            				Ts_past = Td(1)
            				Td(1) = Ts_past - (flx_deb / dflx_deb)

            				! Apply step-size limits to avoid large changes in surface temperature
            				if (Td(1) - Ts_past > 1.0) Td(1) = Ts_past + 1.0
                                        if (Td(1) - Ts_past < -1.0) Td(1) = Ts_past - 1.0

            				! Recompute temperature profile if needed
              				if (it .eq. starth)then
                                                Td(1) = 273.15
                                                Td_gradient = (Td(1) - Td(Nt)) / debris_thickness
                                                do jt = 2, Nt-1
                                                        Td(jt) = Td(1) - (jt * h_deb) * Td_gradient
                                                end do
                                        endif
                                             
                                        if (it.gt.starth.and.daysnow_deb(i,j) .eq. 1) then
                                                Td(1) = Tair+273.15
                				Td_gradient = (Td(1) - Td(Nt)) / debris_thickness
                				do jt = 2, Nt-1
                    					Td(jt) = Td(1) - (jt * h_deb) * Td_gradient
                				end do
            				else
                				! Crank-Nicholson recalculations
                				do jt = 2, Nt-1

                                                        ! Crank-Nicholson coefficients
                                                        C_deb(jt) = k_eff_deb(jt) * (timestep_smb * sechr) / (2.0 * vol_heat_cap_deb(jt) * h_deb**2)
                    					a_Crank(jt) = C_deb(jt)
                    					b_Crank(jt) = 2.0 * C_deb(jt) + 1.0
                                                        c_Crank(jt) = C_deb(jt)

                    					if (jt.eq.2) then
                        					d_Crank(jt) = C_deb(jt) * (Td(1) + Td_past(i,j,1) + Td_past(i,j,jt+1)) + &
                                         			(1.0 - 2.0 * C_deb(jt)) * Td_past(i,j,jt)
                    					elseif (jt.lt.Nt-1) then
                        					d_Crank(jt) = C_deb(jt) * (Td_past(i,j,jt-1) + Td_past(i,j,jt+1)) + &
                                                                (1.0 - 2.0 * C_deb(jt)) * Td_past(i,j,jt)
                                                        endif
                                                        if (jt.eq.Nt-1) then
                                                                d_Crank(jt) = 2.0 * C_deb(jt) * Td(Nt) + &
                                                                C_deb(jt) * Td_past(i,j,Nt-2) + (1.0 - 2.0 * C_deb(jt)) * Td_past(i,j,Nt-1)
                                                        end if

                    					! Forward substitution
                    					if (jt.eq.2) then
                        					E_Crank(jt) = b_Crank(jt)
                        					S_Crank(jt) = d_Crank(jt)
                    					else
                        					E_Crank(jt) = b_Crank(jt) - a_Crank(jt) / E_Crank(jt-1) * c_Crank(jt-1)
                        					S_Crank(jt) = d_Crank(jt) + a_Crank(jt) / E_Crank(jt-1) * S_Crank(jt-1)
                    					end if
                				end do

                				! Backward substitution
                				do jt = Nt-1, 2, -1
                    					if (jt.eq.Nt-1) then
                        					Td(jt) = S_Crank(jt) / E_Crank(jt)
                    					else
                        					Td(jt) = (S_Crank(jt) + c_Crank(jt) * Td(jt+1)) / E_Crank(jt)
                    					end if
                				end do
            				end if

					! Surface energy flux calculations for debris
                                        if(full_seb_on.eq.1)then
                                           qnet_deb = (1.0 - albedo_deb) * qsin ! (tau already included in input)
                                           lnet_deb = skyem*(stf_bltz*(Tair+273.15)**4) - em_s*(stf_bltz*(Td(1))**4)
                                           shf_deb = rhoa*cA*(C_ex)*uwind*(Tair+273.15-(Td(1)))
                                           qZ = (0.622*(rha/100)*6.112*exp((17.67 * Tair) / (Tair + 243.5))) / &
                                                (pres - 0.378 * (rha/100) * 6.112 * exp((17.67 * Tair) / (Tair + 243.5)))
                                           if(prec.gt.0.and.Tair.gt.ttresh)then
                                              qS = (0.622 * 6.112 * exp((17.67 * (Td(1)-273.15)) / ((Td(1)-273.15) + 243.5))) / &
                                                   (pres - 0.378 * 6.112 * exp((17.67 * (Td(1)-273.15)) / ((Td(1)-273.15) + 243.5)))
                                           else
                                              qS = qZ*(Td(1)/(Tair+273.15)) 
                                           endif
                                           lhf_deb = rhoa*Lv*(C_ex)*uwind*(qZ-qS)
                                           if(prec.gt.0.and.Tair.gt.ttresh)then
                                              qrain_deb = cW*rhow*((prec)/(timestep_smb*sechr))*((Tair+273.15) - Td(1))
                                           else
                                              qrain_deb = 0
                                           endif
                                        else
                                           qnet_deb = (1.0 - albedo_deb) * tau * insol
                                           ta_flx_deb = d0 + d1 * ((Tair+273.15) - Td(1))
                                        endif
                                        if (lin_temp_grad_on.eq.1) then
                                           qc_deb = k_eff_deb(1) * (273.15 - Td(1)) / debris_thickness
                                        else
                                           qc_deb = k_eff_deb(1) * (Td(2) - Td(1)) / h_deb
                                        endif
                                        if (full_seb_on.eq.1) then
                                           flx_deb = qnet_deb + lnet_deb + shf_deb + lhf_deb + qrain_deb + qc_deb
                                        else
                                           flx_deb = qnet_deb + ta_flx_deb + qc_deb
                                        endif

        				! Derivative of energy fluxes with respect to surface temperature
                                        if(full_seb_on.eq.1)then
                                           dqnet_deb = 0.
                                           dlnet_deb = -4*em_s*stf_bltz*Td(1)**3
                                           dshf_deb = -rhoa*cA*C_ex*uwind
                                           if (prec.gt.0.and.Tair.gt.ttresh) then
                                              dlhf_deb = -rhoa*Lv*C_ex*uwind*(0.622*pres*(6.112 * exp((17.67*(Td(1)-273.15))/((Td(1)-273.15)+243.5)) * &
                                                   (17.67*243.5)/((Td(1)-273.15)+243.5)**2)) / (pres - 0.378 * 6.112 * exp((17.67*(Td(1)-273.15))/ & 
						   ((Td(1)-273.15)+243.5)))**2
                                           else
                                              dlhf_deb = -rhoa*Lv*C_ex*uwind*(qZ/(Tair + 273.15))
                                           endif
                                           if (prec.gt.0.and.Tair.gt.ttresh) then
                                              dqrain_deb = -cW*rhow*((prec)/(timestep_smb * sechr))
                                           else
                                              dqrain_deb = 0.
                                           endif
                                        else
                                           dqnet_deb = 0.0
                                           dta_flx_deb = -d1
                                        endif
                                        if (lin_temp_grad_on.eq.1) then
                                           dqc_deb = -k_eff_deb(1) / debris_thickness
                                        else
                                           dqc_deb = -k_eff_deb(1) / h_deb
                                        endif
                                        if (full_seb_on.eq.1) then
                                           dflx_deb = dqnet_deb + dlnet_deb + dshf_deb + dlhf_deb + dqrain_deb + dqc_deb
                                        else
                                           dflx_deb = dqnet_deb + dta_flx_deb + dqc_deb
                                        endif

        				! Set maximum iterations to 100 for Newton-Raphson scheme
        				if (n_iterations.eq.100) then
            					Td(1) = (Td(1) + Ts_past) / 2.0
        				end if

                                 end do ! End of Newton-Raphson Crank-Nicholson loop

                                 enflux_deb = flx_deb
                                 Td_gradient = (Td(1) - Td(Nt)) / debris_thickness
				 ! Check false values
                                 do jt = 1, Nt
                                    if(Td(jt).eq.0)then
				       write(*,*) 'WARNING: Td = 0 at i,j,jt,it = ', i, j, jt, it
                                       Td(jt) = 273.15
                                    endif
                                end do

    			end if  ! End of snow condition over debris check

                   else ! If pixel is not debris-covered

                        enflux_deb = 0.
                        qnet_deb = 0.
                        ta_flx_deb = 0.
                        qc_deb = 0.
                        flx_deb = 0.
                        Td_gradient = 0.
                        Ts_past = 0.
                        Td(:)=273.15
                        lnet_deb = 0.
                        shf_deb = 0.
                        lhf_deb = 0.
                        qrain_deb = 0.
                        qS = 0.
                        qZ = 0.
       
		    end if  ! End of mask_deb check

                    !--------------------------------------------------------------------  
		    ! Melt for snow and (sub-debris) ice
		    !-------------------------------------------------------------------- 

		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    !!!!!! For snow and clean ice !!!!!!
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    
		    melt = max(0.0, (timestep_smb*sechr*enflux)/(rhowater*xLm))

		    ! For safety, shouldn't happen

                    if (melt.le.0)then
		       melt = 0
		       snowmelt = 0
                       runoff = 0
                    endif

		    ! Calculate meltwater production and potential retention in snow pack

		    if(snowdep(i,j).gt.0)then

                        ! Snow melt
                        snowmelt = min(snowdep(i,j), melt)

                        ! Check if snowmelt exceeds the remaining retention capacity
                        if(snowret(i,j).le.(eta*snowdep(i,j)))then
                           snowflow = 0
                           snowret(i,j) = snowret(i,j) + snowmelt - snowflow
                           runoff = snowflow
                        elseif(snowret(i,j).gt.eta*snowdep(i,j)) then
                           snowflow = snowmelt - (eta*snowdep(i,j)-snowret(i,j))
                           snowret(i,j) = snowret(i,j) + snowmelt - snowflow
                           runoff = snowflow
                        endif
                           
			! For safety, shouldn't happen

                         if (snowret(i,j).lt.0)then
                            snowret(i,j) = 0
                         endif
                         if (snowflow.lt.0)then
                            snowflow = 0
                         endif
                         if(eta.eq.0)then
                            runoff = snowmelt
                         endif
                         if (runoff.lt.0)then
                           runoff = 0.
                         endif

                      elseif(snowdep(i,j).le.0)then

                         ! Outflow of remaining retained water from recently disappeared snowpack     
                         if(daysnow(i,j).eq.1.and.eta.ne.0)then
                            snowflow = snowret(i,j)
                            runoff = runoff + snowflow
                            snowret(i,j) = 0.
                         else
                            snowflow = 0
                            runoff = 0
                         endif

                        ! All snow-related variables are 0
			snowmelt = 0
			snowret(i,j) = 0

                        ! Runoff is just ice melt
                        runoff = runoff + max(0.0, melt)
                        ice_melt(i,j) = ice_melt(i,j) + max(0.0, melt)

			! For safety, shouldn't happen

			 if (runoff.lt.0)then
                           runoff = 0.
                         endif

                      endif

                   ! Some minor adjustments to ensure no false values occur

		   if(melt.le.0.)then
                     ! No melt
		       melt = 0
		       snowmelt = 0
                       runoff = 0
                    endif
                    
		  ! Update snow depth over ice

		   snowdep(i,j) = snowdep(i,j) + Psolid - snowmelt

		  ! For safety, shouldn't happen 
                   if (snowdep(i,j).lt.0)then
                      snowdep(i,j)=0. 
                   endif
                    
		  ! Calculate mass balance for clean ice
			
		  cumbal(i,j)=cumbal(i,j)+Psolid-runoff     ! <---- SMB FOR CLEAN ICE

		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  !!!!!! For sub-debris ice !!!!!!
		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  if (mask_deb(i,j) .eq. 1) then  ! Check if in the debris region

		    ! Calculate meltwater production

		    if(snowdep_deb(i,j).gt.0)then

		         melt_deb = max(0.0, (timestep_smb*sechr*enflux_deb)/(rhowater*xLm))

		         ! For safety, shouldn't happen

                         if (melt_deb.le.0)then
		            melt_deb = 0
		            snowmelt_deb = 0
                            runoff_deb = 0
                         endif

                        ! Snow melt
                        snowmelt_deb = min(snowdep_deb(i,j), melt_deb)

                        ! Check if snowmelt exceeds the remaining retention capacity
                        if(snowret_deb(i,j).le.(eta*snowdep_deb(i,j)))then
                           snowflow_deb = 0
                           snowret_deb(i,j) = snowret_deb(i,j) + snowmelt_deb - snowflow_deb
                           runoff_deb = snowflow_deb
                        elseif(snowret_deb(i,j).gt.eta*snowdep_deb(i,j)) then
                           snowflow_deb = snowmelt_deb - (eta*snowdep_deb(i,j)-snowret_deb(i,j))
                           snowret_deb(i,j) = snowret_deb(i,j) + snowmelt_deb - snowflow_deb
                           runoff_deb = snowflow_deb
                        endif
                           
			! For safety, shouldn't happen

                         if (snowret_deb(i,j).lt.0)then
                            snowret_deb(i,j) = 0
                         endif
                         if (snowflow_deb.lt.0)then
                            snowflow_deb = 0
                         endif
                         if(eta.eq.0)then
                            runoff_deb = snowmelt_deb
                         endif
                         if (runoff_deb.lt.0)then
                           runoff_deb = 0.
                         endif

                    elseif(snowdep_deb(i,j).le.0)then

                       enflux_deb = k_eff_deb(Nt-1) * (Td(Nt-1) - Td(Nt)) / h_deb
                       if(lin_temp_grad_on.eq.1)then
                          enflux_deb = k_eff_deb(Nt-1) * ((Td(1) - 273.15) / debris_thickness)     ! Linear temperature gradient
                       endif

                         ! Include fractional debris-covered area for sub-debris ice melt
                         melt_deb = (1-a_deb(i,j))*max(0.0, (timestep_smb*sechr*enflux)/(rhowater*xLm)) & 
				+ a_deb(i,j)*max(0.0, (timestep_smb*sechr*enflux_deb)/(rhowater*xLm))
                       
		         ! For safety, shouldn't happen

                         if (melt_deb.le.0)then
		            melt_deb = 0
		            snowmelt_deb = 0
                            runoff_deb = 0
                         endif

                         ! Outflow of remaining retained water from recently disappeared snowpack     
                         if(daysnow_deb(i,j).eq.1.and.eta.ne.0)then
                            snowflow_deb = snowret_deb(i,j)
                            runoff_deb = runoff_deb + snowflow_deb
                            snowret_deb(i,j) = 0.
                         else
                            snowflow_deb = 0
                            runoff_deb = 0
                         endif

                        ! All snow-related variables are 0
			snowmelt_deb = 0
			snowret_deb(i,j) = 0

                        ! Runoff is just sub-debris ice melt due to conduction
                        runoff_deb = runoff_deb + max(0.0, melt_deb)
                        ice_melt_deb(i,j) = ice_melt_deb(i,j) + max(0.0, melt_deb)

			! For safety, shouldn't happen

			 if (runoff_deb.lt.0)then
                           runoff_deb = 0.
                         endif

                      endif

                     ! Some minor adjustments to ensure no false values occur

		     if(melt_deb.le.0.)then
                       ! No melt
		         melt_deb = 0
		         snowmelt_deb = 0
                         runoff_deb = 0
                      endif

		    ! Update snow depth over debris

		     snowdep_deb(i,j) = snowdep_deb(i,j) + Psolid - snowmelt_deb

                    ! For safety, shouldn't happen
                     if (snowdep_deb(i,j).lt.0)then
                        snowdep_deb(i,j)=0. 
                    endif

		    ! Calculate mass balance for debris-covered ice
			
		    cumbal_deb(i,j)=cumbal_deb(i,j)+Psolid-runoff_deb     ! <---- SMB FOR SUB-DEBRIS ICE

		  else

                    snowdep_deb(i,j)=0.
                    melt_deb = 0.
                    runoff_deb = 0.
		    cumbal_deb(i,j)=0

		  endif 

                  !--------------------------------------------------------------------  
		  ! Save cumulative mass balance
	          !-------------------------------------------------------------------- 

                    G(i,j) = cumbal(i,j)
		    G_deb(i,j) = cumbal_deb(i,j)

                  !--------------------------------------------------------------------  
		  ! Calculate time since snow event or snow cover presence
	          !-------------------------------------------------------------------- 

                    ! Set daysnow + timestep since last present snow cover over clean ice (for snow meltwater retention)
                    if (snowdep(i,j).gt.0)then
                       daysnow(i,j) = 0
                    elseif (snowdep(i,j).le.0)then
                       daysnow(i,j) = daysnow(i,j) + 1
                    endif

		    ! Set days now for debris + timestep since last present snow cover over debris (for snow meltwater retention and temperature profile inside debris)
		    if (mask_deb(i,j).eq.1) then
                        if (snowdep_deb(i,j).gt.0)then
                           daysnow_deb(i,j) = 0
                        elseif (snowdep_deb(i,j).le.0)then
                           daysnow_deb(i,j) = daysnow_deb(i,j) + 1
                        endif
		    endif
                    
                    ! Set tsnow + time since last snow event (for snow and firn albedo)
                    tsnow(i,j) = tsnow(i,j) + (timestep_smb/24.0)
      
                    ! Re-insert internal debris temperature from previous timestep
                    if (mask_deb(i,j).eq.1) then
                         if(it.gt.starth)then
                              do jt = 1, Nt
                                 Td_past(i,j,jt)=Td(jt)
                              end do
                         endif
                    endif

                 elseif (mask_new(i,j).eq.0)then

                   G(i,j) = 0
                   G_deb(i,j) = 0
                   TMA(i,j) = 0
                   ALB(i,j) = 0
                   cumsnow(i,j) = 0
                   shading(i,j) = 0
                   SOLAR(i,j) = 0
                   eng_flux(i,j) = 0
                   qnet_flux_out(i,j) = 0
                    
                 endif                 
                 
              enddo                                                                                                                    
           enddo

           write(*,*),'Loop over the area done for iteration ',it, '3-hours'

           !----------------------------------------------------------------------- 
	   ! If at the end of the year
	   !-----------------------------------------------------------------------  

           if ((it-starth).eq.(stoph-starth)) then

              yn = yn + 1
                            
              ! Set values back to 0 for calculations of next balance year; for all (i,j)
              
              cumbal=0.
              snowdep=0.
	      snowflow = 0.
              snowret=0.
              shade=0.
              TMA = 0.
              sir = 0.
              ALB = 0.
              cumsnow = 0.
              shading = 0.
              daysnow = 0.
              SOLAR = 0.
              eng_flux = 0.
              qnet_flux_out = 0.
              ice_melt = 0.
              cumbal_deb=0.
              snowdep_deb=0.
	      snowflow_deb = 0.
              snowret_deb=0.
	      daysnow_deb=0.
              Td_gradient = 0.
              Td(:) = 273.15
              Td_past(:,:,:) = 273.15
              a_Crank(:) = 0.
              b_Crank(:) = 0.
              c_Crank(:) = 0.
              d_Crank(:) = 0.
              E_Crank(:) = 0.
              S_Crank(:) = 0.
              n_iterations = 0.
              phi_deb(:) = 0.
              vol_heat_cap_deb(:) = 0.
              k_eff_deb(:) = 0.
              Ts_past = 273.15
              ice_melt_deb = 0.
              lnet_deb = 0.
              shf_deb = 0.
              lhf_deb = 0.
              qrain_deb = 0.
              dlnet_deb = 0.
              dshf_deb = 0.
              dlhf_deb = 0.
              dqrain_deb = 0.
              qZ = 0.
              qS = 0.

              write(*,*), 'Values cleared'

	   endif
                           
!-----------------------------------------------------------------------      
!  If iterations over the years are done, end program
!-----------------------------------------------------------------------
              
	if (it.eq.dims) then
   	   write(*,*) 'We reached the end of the mass balance time series'
   	   stop
	else
   	   it = it + 1
   	   write(*,*) 'Restart the iteration'
	endif
           
      end do

      end program
         
