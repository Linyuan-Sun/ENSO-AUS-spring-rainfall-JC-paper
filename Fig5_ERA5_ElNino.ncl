year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)

season = "SON"

latS = -44
latN = -10
lonW = 112
lonE = 156

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
pltName = "Fig5_ERA5_ElNino"

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/moisture_budget/moisture_budget.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	mtpr_temp = lonFlip(f->mtpr(start_ind:end_ind,:,:))
	mtpr = mtpr_temp*30*24*3600. ;/1000*1000 (convert from kg/m2/s to mm/month)
	copy_VarCoords(mtpr_temp,mtpr)
	precip_season = rm_single_dims(month_to_seasonN(mtpr(:,{latS:latN},{lonW:lonE}),season))
	precip_season&time = total_year
    lat = precip_season&latitude
    lon = precip_season&longitude
    nlat = dimsizes(lat)
    nlon = dimsizes(lon)

;-----Composite and T-test
	total_avg = dim_avg_n_Wrap(precip_season,0)
	total_var = dim_variance_n_Wrap(precip_season,0)
	total_num = dimsizes(total_year)
	iflag    = True      
	tval_opt = False 
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(precip_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(precip_season({EP_LaNina_year},:,:),0)
	EP_avg = EP_ElNino_avg;(EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(precip_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(precip_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(precip_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(precip_season({CP_LaNina_year},:,:),0)
	CP_avg = CP_ElNino_avg;(CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(precip_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(precip_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)

	
	; f = addfile("ERA5.nc", "c")
	; f->EP_avg = EP_avg
	; f->CP_avg = CP_avg
	; exit
	; delete(precip_season)

    ;-----combine
    Comp_avg_precip = new((/2, nlat, nlon/),"float")
    Comp_prob_precip = new((/2, nlat, nlon/),"float")
    Comp_avg_precip(0,:,:) = EP_avg
    Comp_avg_precip(1,:,:) = CP_avg
    Comp_prob_precip(0,:,:) = EP_prob
    Comp_prob_precip(1,:,:) = CP_prob
	
	Comp_prob_precip = where(abs(Comp_avg_precip).lt.5,Comp_prob_precip@_FillValue,Comp_prob_precip)
	delete([/EP_avg, CP_avg, EP_prob, CP_prob/])
    
;-----plot figure
	wks = gsn_open_wks(pltType,pltName)
	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@gsnAddCyclic = False
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
	res_prob@cnLevels              = (/0.90,0.95/)
	res_prob@cnMonoFillPattern = False 
	res_prob@cnMonoFillColor = False
	res_prob@cnFillColors = (/"white","grey","yellow"/)
	res_prob@cnFillPatterns        = (/-1,17,17/) ;
	res_prob@cnFillDotSizeF = 0.0035 ; dots size
    res_prob@cnFillScaleF = 1. ;dots density

	; res@mpShapeMode = "FreeAspect"
	; res@vpWidthF = 1.0
	; res@vpHeightF = 0.25
	; res@mpCenterLonF = 180
    res@gsnMajorLonSpacing = 15
    res@mpMinLatF = lat(0)
    res@mpMaxLatF = lat(dimsizes(lat)-1)
    res@mpMinLonF = lon(0)
    res@mpMaxLonF = lon(dimsizes(lon)-1)
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = (/-20,-15,-10,-5,5,10,15,20/)
	res@cnFillPalette = "MPL_BrBG"
	res@tmYLLabelFontHeightF = 0.02
	res@tmXBLabelFontHeightF = 0.02
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = ""
	res@gsnRightStringOrthogonalPosF = 0.
	res@gsnRightStringFontHeightF = 0.02
	res@gsnLeftStringFontHeightF = 0.02
    
    res@mpFillOn = False
    res@mpDataSetName = "Earth..4"
    res@mpOutlineBoundarySets = "AllBoundaries"
	res@mpProvincialLineColor  = "grey"
    res@cnFillDrawOrder = "PreDraw"
	res_prob@cnFillDrawOrder = "PreDraw"
    res@mpFillDrawOrder = "Draw"
    res@mpDataBaseVersion  = "MediumRes"
    res@mpOutlineDrawOrder = "Draw"
	res@mpGeophysicalLineColor = "grey30"
	res@mpGeophysicalLineThicknessF = 2.
	
	plot = new(2, graphic)
	plot_prob = new(2, graphic)
	LeftString = (/"(a) EP-ElNino SON: ERA5 rainfall","(b) CP-ElNino SON: ERA5 rainfall"/)
	do i = 0, 1
		res@gsnLeftString = LeftString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_avg_precip(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_precip(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.6
	resP@pmLabelBarHeightF = 0.05
	resP@lbLabelFontHeightF = 0.015
	gsn_panel(wks,plot,(/1,2/),resP)  