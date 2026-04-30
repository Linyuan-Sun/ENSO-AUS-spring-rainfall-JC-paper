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

level = 500
level_u = 200

season = "SON"

	;-----SON define
    EP_ElNino_year = (/1952, 1958, 1973, 1977, 1983, 1988, 1998, 2007, 2016, 2024/)
    CP_ElNino_year = (/1954, 1964, 1966, 1969, 1978, 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
    EP_LaNina_year = (/1955, 1956, 1971, 1972, 1985, 1996, 2000, 2008, 2018, 2021/)
    CP_LaNina_year = (/1965, 1974, 1976, 1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)
	
pltType = "pdf"
pltName = "Fig11_EGR"

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/EGR/backup/EGR.ERA5.mon."+level+"hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	EGR_temp = f->EGR(start_ind:end_ind,{lat_start:lat_end},{lon_start:lon_end})
	
	EGR = EGR_temp/24/60/60.*1e6 ;day**-1 -> 1e-6 s**-1
	copy_VarCoords(EGR_temp,EGR)
	
	EGR_season = rm_single_dims(month_to_seasonN(EGR,season))
	EGR_season&time = total_year
    nlat = dimsizes(EGR_season&lat)
    nlon = dimsizes(EGR_season&lon)
    delete(time)
    
;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/u/u.ERA5.mon."+level_u+"hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	u_temp = lonFlip(f->u(start_ind:end_ind,{level_u},{lat_start:lat_end},:))
	u = u_temp(:,:,{lon_start:lon_end})
	u_season = rm_single_dims(month_to_seasonN(u,season))
	u_season&time = total_year
    nlat_u = dimsizes(u_season(0,:,0))
    nlon_u = dimsizes(u_season(0,0,:))
	
;-----Composite and T-test
	iflag    = True      
	tval_opt = False

;-----EGR500
;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(EGR_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(EGR_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(EGR_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(EGR_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(EGR_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(EGR_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(EGR_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(EGR_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_EGR = new((/2, nlat, nlon/),"float")
    Comp_prob_EGR = new((/2, nlat, nlon/),"float")
    Comp_avg_EGR(0,:,:) = EP_avg
    Comp_avg_EGR(1,:,:) = CP_avg
    Comp_prob_EGR(0,:,:) = EP_prob
    Comp_prob_EGR(1,:,:) = CP_prob
    delete([/EP_ElNino_avg,EP_LaNina_avg,EP_avg,CP_ElNino_avg,CP_LaNina_avg,CP_avg/])
	Comp_prob_EGR = where(abs(Comp_avg_EGR).lt.0.2,Comp_prob_EGR@_FillValue,Comp_prob_EGR)

;-----u200    
;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(u_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(u_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(u_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(u_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    ;-----combine
    Comp_avg_u = new((/2, nlat_u, nlon_u/),"float")
    Comp_avg_u(0,:,:) = EP_avg
    Comp_avg_u(1,:,:) = CP_avg

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
	
	res@mpMaxLatF = lat_start
	res@mpMinLatF = lat_end
	res@gsnAddCyclic = False
	res@mpMaxLonF = lon_end
	res@mpMinLonF = lon_start
	res@mpCenterLonF = (lon_start+lon_end)/2.
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = (/-1.,-0.8,-0.6,-0.4,-0.3,-0.2,0.2,0.3,0.4,0.6,0.8,1./)
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
	res@gsnRightString = "10~S~-6~N~s~S~-1~N~"
	res@gsnRightStringOrthogonalPosF = 0.
	res@gsnRightStringFontHeightF = 0.02
	res@gsnLeftStringFontHeightF = 0.02
	
	res_u = True
	res_u@gsnDraw = False
	res_u@gsnFrame = False
	res_u@gsnMaximize = False
	res_u@lbLabelBarOn         = False
	res_u@cnFillOn = False 
	res_u@cnLinesOn = True
	res_u@cnLineLabelsOn = True
	res_u@cnInfoLabelOn = False
	res_u@cnLevelSelectionMode = "ExplicitLevels"
	res_u@cnLevels = (/-6,-4,-3,-2,-1,1,2,4,6,8/)
	res_u@cnLineLabelPlacementMode = "Computed" 
	res_u@cnLineLabelDensityF = 0.5
	res_u@cnLineLabelInterval = 1 
    res_u@cnLineThicknessF = 2.
    res_u@gsnContourNegLineDashPattern = 1
		
	nplot = dimsizes(Comp_avg_EGR(:,0,0))
	plot = new(nplot, graphic)
	plot_prob = new(nplot, graphic)
	plot_EGR = new(nplot, graphic)
	plot_u = new(nplot, graphic)
	LeftString = (/"(a) EP-ENSO "+season+": EGR500 U200","(b) CP-ENSO "+season+": EGR500 U200"/)
	do i = 0, nplot-1
		res@gsnLeftString = LeftString(i)
		plot(i) = gsn_csm_contour_map(wks,Comp_avg_EGR(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_EGR(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
		plot_u(i) = gsn_csm_contour(wks,Comp_avg_u(i,:,:),res_u)
		overlay(plot(i),plot_u(i))
	end do
	
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  