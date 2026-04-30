load "/home/561/ls3248/NCL/WAF_TN2001_lysun.ncl"
lat_start = -80
lat_end = 10
lon_start = 40
lon_end = 320
year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)
level = 200
rCyclic_lon = True ;True for 0-360 or -180-180; otherwise False for regional data
season = "SON"
season_id = 10 ;1-DJF, 2-JFM, 3-FMA, 4-MAM, 5-AMJ, 6-MJJ, 7-JJA, 8-JAS, 9-ASO, 10-SON, 11-OND, 12-NDJ

	;-----SON define
    EP_ElNino_year = (/1952, 1958, 1973, 1977, 1983, 1988, 1998, 2007, 2016, 2024/)
    CP_ElNino_year = (/1954, 1964, 1966, 1969, 1978, 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
    EP_LaNina_year = (/1955, 1956, 1971, 1972, 1985, 1996, 2000, 2008, 2018, 2021/)
    CP_LaNina_year = (/1965, 1974, 1976, 1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)

pltType = "pdf"
pltName = "Fig9_Z200_WAF"

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/hgt/hgt.ERA5.mon."+level+"hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	hgt_temp = lonFlip(rm_single_dims(f->hgt(start_ind:end_ind,:,:,:)))
	hgt = hgt_temp(:,{lat_start:lat_end},{lon_start:lon_end})
	hgt_season = rm_single_dims(month_to_seasonN(hgt,season))
	hgt_season&time = total_year
    nlat = dimsizes(hgt&latitude)
    nlon = dimsizes(hgt&longitude)
    
;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(hgt_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(hgt_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(hgt_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(hgt_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(hgt_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(hgt_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(hgt_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(hgt_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)

    ;-----combine
    Comp_avg_hgt = new((/2, nlat, nlon/),"float")
    Comp_prob_hgt = new((/2, nlat, nlon/),"float")
    Comp_avg_hgt(0,:,:) = EP_avg
    Comp_avg_hgt(1,:,:) = CP_avg
    Comp_prob_hgt(0,:,:) = EP_prob
    Comp_prob_hgt(1,:,:) = CP_prob

;-----WAF calculation
	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/u/clim/u.ERA5.seasonal.clim.1951-2024.nc","r")
	u_temp = lonFlip(f->u({season_id},{level},{lat_start:lat_end},:))
	u = conform(Comp_avg_hgt,u_temp(:,{lon_start:lon_end}),(/1,2/))
	copy_VarCoords(Comp_avg_hgt,u)
	
	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/v/clim/v.ERA5.seasonal.clim.1951-2024.nc","r")
	v_temp = lonFlip(f->v({season_id},{level},{lat_start:lat_end},:))
	v = conform(Comp_avg_hgt,v_temp(:,{lon_start:lon_end}),(/1,2/))
	copy_VarCoords(Comp_avg_hgt,v)
    
    WAF = WAF_TN2001_lysun(Comp_avg_hgt, u, v, level, rCyclic_lon)
	WAF_x = WAF[0]
    WAF_y = WAF[1]
	
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
	res@gsnAddCyclic = False
	res@mpMaxLonF = lon_end
	res@mpMinLonF = lon_start
	res@mpCenterLonF = (lon_start+lon_end)/2.
	
	res@cnLevelSelectionMode = "ManualLevels"
	res@cnLevelSpacingF = 5
	res@cnMinLevelValF = -40
	res@cnMaxLevelValF = 40
	res@cnFillPalette = "BlueDarkRed18"
	res@tmYLLabelFontHeightF = 0.02
	res@tmXBLabelFontHeightF = 0.02
	res@gsnMajorLonSpacing = 60
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "gpm"
	res@gsnRightStringOrthogonalPosF = 0.
	res@gsnRightStringFontHeightF = 0.02
	res@gsnLeftStringFontHeightF = 0.02

	;-----plot WAF
	res_WAF = True
	res_WAF@gsnDraw = False
	res_WAF@gsnFrame = False
	res_WAF@vcGlyphStyle = "FillArrow"     ; turn on curly vectors
	res_WAF@vcFillArrowEdgeThicknessF = 3.
	res_WAF@vcFillArrowWidthF = 0.07
	res_WAF@vcFillArrowHeadYF = 0.1
	res_WAF@vcFillArrowHeadXF = 0.32
	res_WAF@vcFillArrowHeadInteriorXF = 0.25
	res_WAF@vcLineArrowThicknessF = 3.0               ; change vector thickness
    res_WAF@vcMinDistanceF = 0.025
	RefMag = 1.
	res_WAF@vcRefMagnitudeF  =  RefMag            ; define vector ref mag
	res_WAF@vcRefAnnoString1 = RefMag+" m~S~2~N~/s~S~2~N~"
	res_WAF@vcRefAnnoArrowSpaceF = 0.0 ;minimize the height of RefVec box
    res_WAF@vcRefLengthF = 0.09           ; define length of vec ref 
	res_WAF@vcRefAnnoArrowSpaceF = 0.
	res_WAF@vcRefAnnoFontHeightF = 0.012
	res_WAF@vcRefAnnoString2On = False
	res_WAF@vcRefAnnoOrthogonalPosF = -1.
	res_WAF@vcMinMagnitudeF = RefMag/10.
	res_WAF@vcMaxMagnitudeF = RefMag*1;
	
	nplot = dimsizes(Comp_avg_hgt(:,0,0))
	plot = new(nplot, graphic)
	plot_prob = new(nplot, graphic)
	plot_WAF = new(nplot, graphic)
	LeftString = (/"(a) EP-ENSO "+season+": Z200 WAF","(b) CP-ENSO "+season+": Z200 WAF"/)
	do i = 0, nplot-1
		res@gsnLeftString = LeftString(i)
		;---hgt
		plot(i) = gsn_csm_contour_map(wks,Comp_avg_hgt(i,:,:),res)
		;---prob
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_hgt(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
		;---WAF
		plot_WAF(i) = gsn_csm_vector(wks,WAF_x(i,:,:),WAF_y(i,:,:),res_WAF)
		overlay(plot(i),plot_WAF(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  