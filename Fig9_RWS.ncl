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
pltName = "Fig9_RWS"

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
	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/u/clim/u.ERA5.seasonal.clim.1979-2024.nc","r")
	u_clim_temp = rm_single_dims(lonFlip(f->u({season_id},{level},::-1,:)))
	u_clim = u_clim_temp({lat_start:lat_end},{lon_start:lon_end})
	;v_clim
	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/v/clim/v.ERA5.seasonal.clim.1979-2024.nc","r")
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
	absvor_clim = vor_clim + conform(vor_clim,f,(/0/))
	copy_VarCoords(u_clim,absvor_clim)
	
	;divergence anom
	div_anom = uv2dv_cfd(u_anom_season,v_anom_season,lat,lon,boundOpt)
	copy_VarCoords(u_anom_season,div_anom)

	;divergent wind anom
	dvwnd = dv2uvF_Wrap(div_anom)
	dvuwnd = dvwnd(0,:,:,:)
	dvvwnd = dvwnd(1,:,:,:)

	longName = "S1 term"
	units     = "s^-2"
	opt   = 0      ; Return only the advected variable
	S1_season = -1*advect_variable_cfd(dvuwnd,dvvwnd,conform(dvuwnd,absvor_clim,(/1,2/)),lat,lon,cyclic,longName,units,opt)
	copy_VarCoords(u_anom_season,S1_season)
	
	S2_season = -1*conform(div_anom,absvor_clim,(/1,2/))*div_anom
	copy_VarCoords(u_anom_season,S2_season)

;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
    ;-----S1
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(S1_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(S1_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(S1_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(S1_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)	
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(S1_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(S1_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(S1_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(S1_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_S1 = new((/2, nlat, nlon/),"float")
    Comp_prob_S1 = new((/2, nlat, nlon/),"float")
    Comp_avg_S1(0,:,:) = EP_avg
    Comp_avg_S1(1,:,:) = CP_avg
    Comp_prob_S1(0,:,:) = EP_prob
    Comp_prob_S1(1,:,:) = CP_prob
    
    ;-----S2
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(S2_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(S2_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(S2_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(S2_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)	
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(S2_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(S2_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(S2_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(S2_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_S2 = new((/2, nlat, nlon/),"float")
    Comp_prob_S2 = new((/2, nlat, nlon/),"float")
    Comp_avg_S2(0,:,:) = EP_avg
    Comp_avg_S2(1,:,:) = CP_avg
    Comp_prob_S2(0,:,:) = EP_prob
    Comp_prob_S2(1,:,:) = CP_prob
	
	;-----unit adjustment and spatial smooth
	Comp_avg_S1_plt = Comp_avg_S1*1e10
	copy_VarCoords(Comp_avg_S1,Comp_avg_S1_plt)
	Comp_avg_S2_plt = Comp_avg_S2*1e10
	copy_VarCoords(Comp_avg_S2,Comp_avg_S2_plt)
	
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
	res@cnSmoothingOn = True
    
	;------S1
	res_S1 = res
	res_S1@cnFillOn = False
	res_S1@cnLinesOn = True
	res_S1@cnLevelSelectionMode = "ExplicitLevels"
    res_S1@cnLevels = (/-0.6,-0.4,-0.2,0.2,0.4,0.6/)
	res_S1@cnLineThicknessF = 2.
    res_S1@cnLineLabelsOn = True
	res_S1@cnLineLabelPlacementMode = "Computed" 
	res_S1@cnLineLabelDensityF = 1.5
	res_S1@cnLineLabelInterval = 1 
	res_S1@cnExplicitLineLabelsOn= True
	res_S1@cnLineLabelStrings = res_S1@cnLevels+""

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
	res@gsnAddCyclic = False;True
	res@mpMaxLonF = lon_end
	res@mpMinLonF = lon_start
	res@mpCenterLonF = (lon_start+lon_end)/2.
	
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1/)
	color = read_colormap_file("temp_diff_18lev")
	color_id = (/3,4,6,8,10,11,12,14,16,18,19/)-2
	res@cnFillPalette = color(color_id,:)
	res@tmYLLabelFontHeightF = 0.02
	res@tmXBLabelFontHeightF = 0.02
	res@gsnMajorLonSpacing = 60
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "10~S~-10~N~s~S~-2~N~"
	res@gsnRightStringOrthogonalPosF = 0.
	res@gsnRightStringFontHeightF = 0.02
	res@gsnLeftStringFontHeightF = 0.02
    
	nplot = 2
	plot = new(nplot, graphic)
	plot_prob = new(nplot, graphic)
	plot_S1 = new(nplot, graphic)
	
	LeftString = (/"(c) EP-ENSO "+season+": RWS","(d) CP-ENSO "+season+": RWS"/)
	do i = 0, nplot-1
		res@gsnLeftString = LeftString(i)
		res@gsnRightString = "10~S~-10~N~s~S~-2~N~"
		plot(i) = gsn_csm_contour_map(wks,Comp_avg_S2_plt(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_S2(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
		plot_S1(i) = gsn_csm_contour(wks,Comp_avg_S1_plt(i,:,:),res_S1)
        plot_S1(i) = ColorNegDashZeroPosContour(plot_S1(i),"blue","transparent","red")
		overlay(plot(i),plot_S1(i))
	end do
	
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  