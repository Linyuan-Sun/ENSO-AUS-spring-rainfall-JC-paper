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
pltName = "Fig11_vu"

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/vu/vu.ERA5."+level+"hPa."+season+".anom_slidingClim.nc","r")
	vu = f->vu({year_start-1:year_end-1},{lat_start:lat_end},{lon_start:lon_end})
	vu&time = total_year

    nlat = dimsizes(vu&lat)
    nlon = dimsizes(vu&lon)
	
	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/vu/vu.ERA5."+season+"."+level+"hPa.nc","r")
	vu_clim= dim_avg_n_Wrap(f->vu,0)
	
;-----Composite and T-test
	iflag    = True      
	tval_opt = False
	
;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(vu({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(vu({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(vu({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(vu({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(vu({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(vu({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(vu({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(vu({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)

    ;-----combine
    Comp_avg_vu = new((/2, nlat, nlon/),"float")
    Comp_prob_vu = new((/2, nlat, nlon/),"float")
    Comp_avg_vu(0,:,:) = EP_avg
    Comp_avg_vu(1,:,:) = CP_avg
    Comp_prob_vu(0,:,:) = EP_prob
    Comp_prob_vu(1,:,:) = CP_prob
    Comp_prob_vu = where(abs(Comp_avg_vu).lt.2,0,Comp_prob_vu)

	Comp_avg_vu!0 = "type"
	Comp_avg_vu!1 = "lat"
	Comp_avg_vu!2 = "lon"	
	Comp_avg_vu&lat = EP_avg&lat
	Comp_avg_vu&lon = EP_avg&lon
	copy_VarCoords(Comp_avg_vu,Comp_prob_vu)

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
	res@cnLevels = (/-6,-5,-4,-3,-2,-1,1,2,3,4,5,6/)
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
	res@gsnRightString = "m~S~2~N~s~S~-2~N~"
	res@gsnRightStringOrthogonalPosF = 0.
	res@gsnRightStringFontHeightF = 0.02
	res@gsnLeftStringFontHeightF = 0.02
	
	res_clim = True
	res_clim@gsnDraw = False
	res_clim@gsnFrame = False
	res_clim@gsnMaximize = False
	res_clim@lbLabelBarOn         = False
	res_clim@cnFillOn = False 
	res_clim@cnLinesOn = True
	res_clim@cnLineLabelsOn = True
	res_clim@cnInfoLabelOn = False
	res_clim@cnLevelSelectionMode = "ExplicitLevels"
	res_clim@cnLevels = (/-30,-22.5,-15,-10,-5,2,4,6,8/)
	res_clim@cnLineThicknessF = 2.
	res_clim@gsnContourNegLineDashPattern = 1
	res_clim@cnLineLabelPlacementMode = "Computed" 
	res_clim@cnLineLabelDensityF = 0.5 
	res_clim@cnLineLabelInterval = 1 
    res_clim@gsnLeftString = ""
		
	nplot = dimsizes(Comp_avg_vu(:,0,0))
	plot = new(nplot, graphic)
	plot_prob = new(nplot, graphic)
	plot_vu = new(nplot, graphic)
	plot_vu_clim = new(nplot, graphic)
	LeftString = (/"(e) EP-ENSO "+season+": vu","(f) CP-ENSO "+season+": vu"/)
	do i = 0, nplot-1
		res@gsnLeftString = LeftString(i)
		plot(i) = gsn_csm_contour_map(wks,Comp_avg_vu(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_vu(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
		plot_vu_clim(i) = gsn_csm_contour(wks,vu_clim,res_clim)
		overlay(plot(i),plot_vu_clim(i))
	end do
	
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  