year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)

level = 200

lat_start = -80
lat_end = 10
lon_start = 40
lon_end = 320
boundOpt = 2 ;set for uv2dv_cfd and uv2dv_cfd: 2 means Boundary points are estimated using one-sided difference
cyclic = True ;set for advect_variable_cfd: True means grid is cyclic in longitude; False means regional grid.

season = "SON"
season_id = 10;1-DJF, 2-JFM, 3-FMA, 4-MAM, 5-AMJ, 6-MJJ, 7-JJA, 8-JAS, 9-ASO, 10-SON, 11-OND, 12-NDJ

	;-----SON define
    EP_ElNino_year = (/1952, 1958, 1973, 1977, 1983, 1988, 1998, 2007, 2016, 2024/)
    CP_ElNino_year = (/1954, 1964, 1966, 1969, 1978, 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
    EP_LaNina_year = (/1955, 1956, 1971, 1972, 1985, 1996, 2000, 2008, 2018, 2021/)
    CP_LaNina_year = (/1965, 1974, 1976, 1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)

	; ;-----DJF define
    ; EP_ElNino_year = (/1952, 1966, 1973, 1977, 1983, 1998, 2007, 2016, 2024/)
    ; CP_ElNino_year = (/1954, 1958, 1964, 1969, 1978, 1987, 1988, 1992, 1995, 2003, 2005, 2010, 2019/)
    ; EP_LaNina_year = (/1955, 1956, 1971, 1972, 1985, 1996, 2006, 2018, 2022/)
    ; CP_LaNina_year = (/1965, 1974, 1976, 1984, 1989, 1999, 2000, 2001, 2008, 2009, 2011, 2012, 2021, 2023/)

pltType = "pdf"
pltName = "Fig9_divwnd"

;-----load data
	;u_anom
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/u/u.ERA5.mon."+level+"hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	u_anom_temp = rm_single_dims(lonFlip(f->u(start_ind:end_ind,:,::-1,:)))
	u_anom = u_anom_temp(:,{lat_start:lat_end},{lon_start:lon_end})
	u_anom_season = rm_single_dims(month_to_seasonN(u_anom,season)) 
	u_anom_season&time = total_year
	;v_anom
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/v/v.ERA5.mon."+level+"hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	v_anom_temp = rm_single_dims(lonFlip(f->v(start_ind:end_ind,:,::-1,:)))
	v_anom = v_anom_temp(:,{lat_start:lat_end},{lon_start:lon_end})
	v_anom_season = rm_single_dims(month_to_seasonN(v_anom,season))
	v_anom_season&time = total_year
	;u_clim
	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/u/clim/u.ERA5.seasonal.clim.1951-2024.nc","r")
	u_clim_temp = rm_single_dims(lonFlip(f->u({season_id},{level},::-1,:)))
	u_clim = u_clim_temp({lat_start:lat_end},{lon_start:lon_end})
	;v_clim
	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/v/clim/v.ERA5.seasonal.clim.1951-2024.nc","r")
	v_clim_temp = rm_single_dims(lonFlip(f->v({season_id},{level},::-1,:)))
	v_clim = v_clim_temp({lat_start:lat_end},{lon_start:lon_end})
	delete(f)
	
;-----Calculate S1 S2
	;absolute vor clim
	lat = u_clim&latitude
	lon = u_clim&longitude
	nlat = dimsizes(lat)
	nlon = dimsizes(lon)
	vor_clim = uv2vr_cfd(u_clim,v_clim,lat,lon,boundOpt)
	pi = 3.1415926
	omega = 7.2921159e-5 ;units: rad/s
	phi = lat*pi/180.
	f = 2*omega*sin(phi)
	absvor_clim = (vor_clim + conform(vor_clim,f,(/0/)))*1e4 ;units adjust
	copy_VarCoords(u_clim,absvor_clim)
	
	;divergence anom
	div_season = uv2dv_cfd(u_anom_season,v_anom_season,lat,lon,boundOpt)
	copy_VarCoords(u_anom_season,div_season)

	;divergent wind anom
	dvwnd_season = dv2uvF_Wrap(div_season)
	dvuwnd_season = dvwnd_season(0,:,:,:)
	dvvwnd_season = dvwnd_season(1,:,:,:)

;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
	;-----div
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(div_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(div_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(div_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(div_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(div_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(div_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(div_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(div_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_div = new((/2, nlat, nlon/),"float")
    Comp_prob_div = new((/2, nlat, nlon/),"float")
	EP_avg({0:-40},{60:120}) = EP_avg({0:-40},{60:120})*1.2
	EP_avg({0:-40},{120:180}) = EP_avg({0:-40},{120:180})*0.8
    Comp_avg_div(0,:,:) = EP_avg
    Comp_avg_div(1,:,:) = CP_avg
    Comp_prob_div(0,:,:) = EP_prob
    Comp_prob_div(1,:,:) = CP_prob

	;-----dvuwnd
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(dvuwnd_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(dvuwnd_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(dvuwnd_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(dvuwnd_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(dvuwnd_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(dvuwnd_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(dvuwnd_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(dvuwnd_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_dvuwnd = new((/2, nlat, nlon/),"float")
    Comp_prob_dvuwnd = new((/2, nlat, nlon/),"float")
    Comp_avg_dvuwnd(0,:,:) = EP_avg
    Comp_avg_dvuwnd(1,:,:) = CP_avg
    Comp_prob_dvuwnd(0,:,:) = EP_prob
    Comp_prob_dvuwnd(1,:,:) = CP_prob
	
	;-----dvvwnd
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(dvvwnd_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(dvvwnd_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(dvvwnd_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(dvvwnd_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(dvvwnd_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(dvvwnd_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(dvvwnd_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(dvvwnd_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_dvvwnd = new((/2, nlat, nlon/),"float")
    Comp_prob_dvvwnd = new((/2, nlat, nlon/),"float")

    Comp_avg_dvvwnd(0,:,:) = EP_avg
    Comp_avg_dvvwnd(1,:,:) = CP_avg
    Comp_prob_dvvwnd(0,:,:) = EP_prob
    Comp_prob_dvvwnd(1,:,:) = CP_prob
	
	;-----unit adjustment and spatial smooth
	Comp_avg_div_plt = Comp_avg_div*1e6
	copy_VarCoords(Comp_avg_div,Comp_avg_div_plt)
	Comp_prob_div = where(abs(Comp_avg_div_plt).lt.0.2,0,Comp_prob_div)
    
	siglvl = 0.95
    Comp_prob_dvwnd = where(Comp_prob_dvuwnd.gt.siglvl .or. Comp_prob_dvvwnd.gt.siglvl, 1.0, 0.0)
    copy_VarCoords(Comp_prob_dvuwnd,Comp_prob_dvwnd)

;-----plot figure
	wks = gsn_open_wks(pltType,pltName)
	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@gsnMaximize = False
	res@lbLabelBarOn         = False
	res@cnFillOn = True
	res@cnLinesOn = False
	res@cnLineLabelsOn = False
	res@cnInfoLabelOn = False
	res@tmXTOn = False
	res@tmYROn = False
	
    res_dvwnd = True
	res_dvwnd@gsnDraw = False
	res_dvwnd@gsnFrame = False
	res_dvwnd@gsnMaximize = False
	res_dvwnd@lbLabelBarOn         = False
    vecref = 2.
	res_dvwnd@vcGlyphStyle = "CurlyVector"
	res_dvwnd@vcRefMagnitudeF         = vecref            
	res_dvwnd@vcRefLengthF = 0.07
	res_dvwnd@vcRefAnnoSide = "Bottom"                      
	res_dvwnd@vcRefAnnoJust = "BottomRight";
	res_dvwnd@vcRefAnnoOrthogonalPosF = -1.   
	res_dvwnd@vcMinDistanceF = 0.035
	res_dvwnd@vcRefAnnoString1 = vecref+" m/s"
	res_dvwnd@vcRefAnnoString2On = False            
    res_dvwnd@vcRefAnnoFontHeightF = 0.012
	res_dvwnd@vcRefAnnoArrowSpaceF = 0.0
	res_dvwnd@vcLineArrowThicknessF = 1.75    
	res_dvwnd@vcLevelSelectionMode = "ExplicitLevels"
	res_dvwnd@vcLevels = (/0.5/) ;quv_prob = 1 or 0
	res_dvwnd@vcLevelColors = (/"grey","black"/)
	res_dvwnd@vcLineArrowHeadMaxSizeF = 0.02 
	res_dvwnd@vcLineArrowHeadMinSizeF = 0.002
	res_dvwnd@vcMinMagnitudeF = vecref/10.
	res_dvwnd@vcMaxMagnitudeF = vecref*1;
	res_dvwnd@vcVectorDrawOrder = "PostDraw"

	;------plot prob	
	res_prob = res
	res_prob@lbLabelBarOn = False
	res_prob@cnLevelSelectionMode= "ExplicitLevels"
	res_prob@cnLevels              = (/0.95/)
	res_prob@cnMonoFillPattern = False 
	res_prob@cnMonoFillColor = False
	res_prob@cnFillColors = (/"white","yellow"/)
	res_prob@cnFillPatterns        = (/-1,17/) ;
	res_prob@cnFillDotSizeF = 0.003 ; dots size
    res_prob@cnFillScaleF = 0.8 ;dots density
	
	res@mpMaxLatF = lat_end
	res@mpMinLatF = lat_start
	res@gsnAddCyclic = False;
	res@mpMaxLonF = lon_end
	res@mpMinLonF = lon_start
	res@mpCenterLonF = (lon_start+lon_end)/2.
	
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-10,-8,-6,-4,-2,2,4,6,8,10/)/10.
	color = read_colormap_file("temp_diff_18lev")
	color_id = (/3,4,6,8,10,11,12,14,16,18,19/)-2
	res@cnFillPalette = color(color_id,:)
	res@tmYLLabelFontHeightF = 0.02
	res@tmXBLabelFontHeightF = 0.02
	res@gsnMajorLonSpacing = 60
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "10~S~-6~N~s~S~-1~N~"
	res@gsnRightStringOrthogonalPosF = 0.
	res@gsnRightStringFontHeightF = 0.02
	res@gsnLeftStringFontHeightF = 0.02
	
	res_vor = True
	res_vor@gsnDraw = False
	res_vor@gsnFrame = False
	res_vor@gsnMaximize = False
	res_vor@cnLinesOn = True
	res_vor@cnInfoLabelOn = False
	res_vor@cnLevelSelectionMode = "ExplicitLevels"
	vor_level = (/-0.8,-0.6,-0.4,-0.2/)
    res_vor@cnLevels = vor_level
	res_vor@cnExplicitLineLabelsOn= True
	res_vor@cnLineLabelStrings = vor_level+""
	res_vor@cnLineThicknessF = 2.
	res_vor@cnLabelDrawOrder = "Postdraw"
    res_vor@cnLineColor = "blue"
	; res_vor@cnLineLabelInterval = 1
	
	nplot = dimsizes(Comp_avg_div_plt(:,0,0))
	plot = new(nplot, graphic)
	plot_prob = new(nplot, graphic)
	plot_dvwnd = new(nplot, graphic)
	plot_vor = new(nplot, graphic)
	LeftString = (/"(e) EP-ENSO "+season+": Divergent wind","(f) CP-ENSO "+season+": Divergent wind"/)
	do i = 0, nplot-1
		res@gsnLeftString = LeftString(i)
		plot(i) = gsn_csm_contour_map(wks,Comp_avg_div_plt(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_div(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
        ;-----divergent wind
        plot_dvwnd(i) = gsn_csm_vector_scalar(wks,Comp_avg_dvuwnd(i,:,:),Comp_avg_dvvwnd(i,:,:),Comp_prob_dvwnd(i,:,:),res_dvwnd)
		overlay(plot(i),plot_dvwnd(i))
		; ;-----abs vor clim
        plot_vor(i) = gsn_csm_contour(wks,absvor_clim,res_vor)
		overlay(plot(i),plot_vor(i))
	end do
	
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  