year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)

; lat_start = -25
; lat_end = -65
lat_start = 0
lat_end = -90
; lon_start = 80
; lon_end = 240
lon_start = 40
lon_end = 320


level = 200

season = "SON"

	;-----SON define
    EP_ElNino_year = (/1952, 1958, 1973, 1977, 1983, 1988, 1998, 2007, 2016, 2024/)
    CP_ElNino_year = (/1954, 1964, 1966, 1969, 1978, 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
    EP_LaNina_year = (/1955, 1956, 1971, 1972, 1985, 1996, 2000, 2008, 2018, 2021/)
    CP_LaNina_year = (/1965, 1974, 1976, 1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)
	
pltType = "pdf"
pltName = "Fig12"

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/Feddy/"+season+"/Feddy.ERA5."+level+"hPa."+season+".anom_slidingClim.nc","r")
	Feddy = f->Feddy({year_start-1:year_end-1},:,:)
	Feddy&time = total_year

	lmbda = Feddy(:,0,0)
	lmbda = 0.0
	pi = 3.1415926
	omega = 2*pi/24/60/60.
	phi = conform(Feddy,(Feddy&lat)*pi/180.,(/1/))
	g = 9.8
	FeddyGT_temp = tofloat(ilapsF(Feddy,lmbda)*2*omega*sin(phi)/g*24*60*60.) ;units: m/day
	copy_VarCoords(Feddy,FeddyGT_temp)
	FeddyGT = FeddyGT_temp(:,{lat_start:lat_end},{lon_start:lon_end})

    nlat = dimsizes(FeddyGT&lat)
    nlon = dimsizes(FeddyGT&lon)

;-----adjust EP
    EP_ElNino_avg_temp = dim_avg_n_Wrap(Feddy({EP_ElNino_year},:,:),0)
	EP_LaNina_avg_temp = dim_avg_n_Wrap(Feddy({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg_temp - EP_LaNina_avg_temp)/2
	copy_VarCoords(EP_ElNino_avg_temp,EP_avg)

	delete([/lmbda,phi/])
	lmbda = 0.0
	phi = conform(EP_avg,(EP_avg&lat)*pi/180.,(/0/))
	FeddyGT_EP = tofloat(ilapsF(EP_avg,lmbda)*2*omega*sin(phi)/g*24*60*60.) ;units: m/day
	copy_VarCoords(EP_avg,FeddyGT_EP)
	delete(EP_avg)
;-----Composite and T-test
	iflag    = True      
	tval_opt = False
	
;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(FeddyGT({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(FeddyGT({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(FeddyGT({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(FeddyGT({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(FeddyGT({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(FeddyGT({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(FeddyGT({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(FeddyGT({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)

    ;-----combine
    Comp_avg_FeddyGT = new((/2, nlat, nlon/),"float")
    Comp_prob_FeddyGT = new((/2, nlat, nlon/),"float")
    Comp_avg_FeddyGT(0,:,:) = FeddyGT_EP({lat_start:lat_end},{lon_start:lon_end});EP_avg

    Comp_avg_FeddyGT(1,:,:) = CP_avg
    Comp_prob_FeddyGT(0,:,:) = EP_prob
    Comp_prob_FeddyGT(1,:,:) = CP_prob

	Comp_prob_FeddyGT = where(abs(Comp_avg_FeddyGT).lt.1,0,Comp_prob_FeddyGT)
	delete([/EP_avg,CP_avg,EP_ElNino_avg,EP_LaNina_avg,EP_ElNino_var,EP_LaNina_var,CP_ElNino_avg,CP_LaNina_avg,CP_ElNino_var,CP_LaNina_var/])
;-----Z200
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
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(hgt_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(hgt_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    ;-----combine
    Comp_avg_hgt = new((/2, nlat, nlon/),"float")
    Comp_avg_hgt(0,:,:) = EP_avg
    Comp_avg_hgt(1,:,:) = CP_avg
	Comp_avg_hgt(:,{0:-25},:) = Comp_avg_hgt@_FillValue ;mask low-latitude

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

	res_contour = res
	res_contour@cnLevelSelectionMode = "ExplicitLevels"	
	; res_contour@cnLevels = (/-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10/)
	res_contour@cnLevels = (/-50,-40,-30,-20,15,25,35,45,55/)
	res_contour@cnLinesOn = True
	res_contour@cnLineLabelsOn = True
	res_contour@cnLineLabelInterval = 1
	res_contour@cnLineThicknessF = 3.
	res_contour@gsnContourNegLineDashPattern = 1
	res_contour@cnFillOn = False
	
	res@mpMaxLatF = lat_start
	res@mpMinLatF = lat_end
	res@gsnAddCyclic = False
	res@mpMaxLonF = lon_end
	res@mpMinLonF = lon_start
	res@mpCenterLonF = (lon_start+lon_end)/2.
	; res@cnLevelSelectionMode = "ManualLevels"
	; res@cnLevelSpacingF = 1.5
	; res@cnMinLevelValF = -8
	; res@cnMaxLevelValF = 8
    ; res@cnLevelSelectionMode = "ExplicitLevels"
    ; res@cnLevels = (/-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10/)
	; res@cnFillPalette = "BlueDarkRed18"
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = (/-8,-6,-4,-3,-2,-1,1,2,3,4,6,8/)
	color = read_colormap_file("temp_diff_18lev")
	color_id = (/3,4,6,8,10,11,12,14,16,18,19/)-2
	res@cnFillPalette = color(color_id,:)
	cmap = read_colormap_file("MPL_BrBG")
	res@mpLandFillColor = cmap(65-2,:)
	res@mpGeophysicalLineColor = "black"
	res@tmYLLabelFontHeightF = 0.02
	res@tmXBLabelFontHeightF = 0.02
	res@gsnMajorLatSpacing = 30
	res@gsnMajorLonSpacing = 60
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "gpm/day"
	res@gsnRightStringOrthogonalPosF = 0.
	res@gsnRightStringFontHeightF = 0.02
	res@gsnLeftStringFontHeightF = 0.02
		
	nplot = dimsizes(Comp_avg_FeddyGT(:,0,0))
	plot = new(nplot, graphic)
	plot_prob = new(nplot, graphic)
	plot_FeddyGT = new(nplot, graphic)
	plot_contour = new(nplot, graphic)
	LeftString = (/"(a) EP-ENSO "+season+": Eddy-induced Geo. H. tendency","(b) CP-ENSO "+season+": Eddy-induced Geo. H. Tendency"/)
	do i = 0, nplot-1
		res@gsnLeftString = LeftString(i)
		plot(i) = gsn_csm_contour_map(wks,Comp_avg_FeddyGT(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_FeddyGT(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
		plot_contour(i) = gsn_csm_contour(wks,Comp_avg_hgt(i,:,:),res_contour)
		overlay(plot(i),plot_contour(i))
		
	end do
	
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  