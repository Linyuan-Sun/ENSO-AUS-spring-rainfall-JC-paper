year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)

level = 200

lat_start = -60
lat_end = 0
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
pltName = "Fig10_RWS"

;-----load data
	;RWS
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/RWS/RWS.mon."+level+"hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	RWS_temp = rm_single_dims(f->RWS(start_ind:end_ind,::-1,:))
	RWS = RWS_temp(:,{lat_start:lat_end},{lon_start:lon_end})
	RWS_season = rm_single_dims(month_to_seasonN(RWS,season)) 
	RWS_season&time = total_year

;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
    ;-----RWS
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(RWS_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(RWS_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(RWS_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(RWS_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)	
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(RWS_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(RWS_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(RWS_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(RWS_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
	nlat = dimsizes(RWS_season&lat)
	nlon = dimsizes(RWS_season&lon)
    Comp_avg_RWS = new((/2, nlat, nlon/),"double")
    Comp_prob_RWS = new((/2, nlat, nlon/),"double")

    Comp_avg_RWS(0,:,:) = EP_avg
    Comp_avg_RWS(1,:,:) = CP_avg
    Comp_prob_RWS(0,:,:) = EP_prob
    Comp_prob_RWS(1,:,:) = CP_prob
   
	;-----unit adjustment and spatial smooth
	Comp_avg_RWS_plt = Comp_avg_RWS*1e10
	copy_VarCoords(Comp_avg_RWS,Comp_avg_RWS_plt)
	Comp_prob_RWS = where(abs(Comp_avg_RWS_plt).lt.0.2,0,Comp_prob_RWS)

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

	;------plot prob	
	res_prob = res
	res_prob@lbLabelBarOn = False
	res_prob@cnLevelSelectionMode= "ExplicitLevels"
	res_prob@cnLevels              = (/0.95/)
	res_prob@cnMonoFillPattern = False 
	res_prob@cnMonoFillColor = False
	res_prob@cnFillColors = (/"white","yellow"/)
	res_prob@cnFillPatterns        = (/-1,17/) ;
	res_prob@cnFillDotSizeF = 0.0035 ; dots size
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
	res@tmYLLabelFontHeightF = 0.015
	res@tmXBLabelFontHeightF = 0.015
	res@gsnMajorLonSpacing = 60
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "10~S~-10~N~s~S~-2~N~"
	res@gsnRightStringOrthogonalPosF = 0.
	res@gsnRightStringFontHeightF = 0.015
	res@gsnLeftStringFontHeightF = 0.015
    
	nplot = 2
	plot = new(nplot, graphic)
	plot_prob = new(nplot, graphic)
	plot_RWS = new(nplot, graphic)
	line = new(nplot, graphic)
	
	LeftString = (/"(b) EP-ENSO "+season+": Total RWS anom.","(c) CP-ENSO "+season+": Total RWS anom."/)
	do i = 0, nplot-1
		res@gsnLeftString = LeftString(i)
		res@gsnRightString = "10~S~-10~N~s~S~-2~N~"
		plot(i) = gsn_csm_contour_map(wks,Comp_avg_RWS_plt(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_RWS(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))

		plres = True
		plres@gsLineColor = "black"
		plres@gsLineThicknessF = 2.
		plres@gsLineDashPattern = 1     
		line(i) = gsn_add_polyline(wks, plot(i), (/0,120,240,360/),(/-30,-30,-30,-30/), plres)
	end do
	
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.65
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/2,1/),resP)  