year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)

year_start_olr = 1976
year_end_olr = 2024
total_year_olr = ispan(year_start_olr,year_end_olr,1)

latS = -35
latN = 35

season = "SON"

	;-----SON define
    EP_ElNino_year = (/1952, 1958, 1973, 1977, 1983, 1988, 1998, 2007, 2016, 2024/)
    CP_ElNino_year = (/1954, 1964, 1966, 1969, 1978, 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
    EP_LaNina_year = (/1955, 1956, 1971, 1972, 1985, 1996, 2000, 2008, 2018, 2021/)
    CP_LaNina_year = (/1965, 1974, 1976, 1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)

    EP_ElNino_year_olr = (/1977, 1983, 1988, 1998, 2007, 2016, 2024/)
    CP_ElNino_year_olr = (/1978, 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
    EP_LaNina_year_olr = (/1985, 1996, 2000, 2008, 2018, 2021/)
    CP_LaNina_year_olr = (/1976, 1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)

	; ;-----DJF define
    ; EP_ElNino_year = (/1952, 1966, 1973, 1977, 1983, 1998, 2007, 2016, 2024/)
    ; CP_ElNino_year = (/1954, 1958, 1964, 1969, 1978, 1987, 1988, 1992, 1995, 2003, 2005, 2010, 2019/)
    ; EP_LaNina_year = (/1955, 1956, 1971, 1972, 1985, 1996, 2006, 2018, 2022/)
    ; CP_LaNina_year = (/1965, 1974, 1976, 1984, 1989, 1999, 2000, 2001, 2008, 2009, 2011, 2012, 2021, 2023/)
	
pltType = "pdf"
pltName = "Fig2_sym"
FontHeightF = 0.015

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/HadISST/HadISST.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	sst = f->sst(start_ind:end_ind,{latS:latN},:)
	sst_season = rm_single_dims(month_to_seasonN(sst,season))
	sst_season&time = total_year
    nlat_sst = dimsizes(sst&latitude)
    nlon_sst = dimsizes(sst&longitude)
	delete(time)
	
	f = addfile("/g/data/if69/ls3248/data_HF/OLR/olr.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start_olr-1)*100+1))
	end_ind = ind(time.eq.((year_end_olr-1)*100+12))
	olr = lonFlip(f->olr(start_ind:end_ind,{latS:latN},:)) ;make lon consistent with SST
	olr_season = rm_single_dims(month_to_seasonN(olr,season))
	olr_season&time = total_year_olr
    nlat_olr = dimsizes(olr&lat)
    nlon_olr = dimsizes(olr&lon)
	delete(time)
	
	f = addfile("/g/data/if69/ls3248/data_HF/AGCD/precip/precip.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	precip = f->precip(start_ind:end_ind,:,:)
	precip_season = rm_single_dims(month_to_seasonN(precip,season))
	precip_season&time = total_year
    lat_precip = precip&lat
    lon_precip = precip&lon
    nlat_precip = dimsizes(lat_precip)
    nlon_precip = dimsizes(lon_precip)
	delete(precip)
    
;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
	
	;-----sst
	;-----EP ENSO
    EP_ElNino_avg_sst = dim_avg_n_Wrap(sst_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg_sst = dim_avg_n_Wrap(sst_season({EP_LaNina_year},:,:),0)
	EP_avg_sst = (EP_ElNino_avg_sst - EP_LaNina_avg_sst)/2
	copy_VarCoords(EP_ElNino_avg_sst,EP_avg_sst)
    EP_ElNino_var_sst = dim_variance_n_Wrap(sst_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num_sst = dimsizes(EP_ElNino_year)
    EP_LaNina_var_sst = dim_variance_n_Wrap(sst_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num_sst = dimsizes(EP_LaNina_year)
    EP_prob_sst = 1. - ttest(EP_ElNino_avg_sst, EP_ElNino_var_sst, EP_ElNino_num_sst, EP_LaNina_avg_sst, EP_LaNina_var_sst, EP_LaNina_num_sst, iflag, tval_opt)
	copy_VarCoords(EP_avg_sst,EP_prob_sst)
	;-----CP ENSO
    CP_ElNino_avg_sst = dim_avg_n_Wrap(sst_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg_sst = dim_avg_n_Wrap(sst_season({CP_LaNina_year},:,:),0)
	CP_avg_sst = (CP_ElNino_avg_sst - CP_LaNina_avg_sst)/2
	copy_VarCoords(CP_ElNino_avg_sst,CP_avg_sst)
    CP_ElNino_var_sst = dim_variance_n_Wrap(sst_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num_sst = dimsizes(CP_ElNino_year)
    CP_LaNina_var_sst = dim_variance_n_Wrap(sst_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num_sst = dimsizes(CP_LaNina_year)
    CP_prob_sst = 1. - ttest(CP_ElNino_avg_sst, CP_ElNino_var_sst, CP_ElNino_num_sst, CP_LaNina_avg_sst, CP_LaNina_var_sst, CP_LaNina_num_sst, iflag, tval_opt)
	copy_VarCoords(CP_avg_sst,CP_prob_sst)
    ;-----combine
    Comp_avg_sst = new((/2, nlat_sst, nlon_sst/),"float")
    Comp_prob_sst = new((/2, nlat_sst, nlon_sst/),"float")
    Comp_avg_sst(0,:,:) = EP_avg_sst
    Comp_avg_sst(1,:,:) = CP_avg_sst
    Comp_prob_sst(0,:,:) = EP_prob_sst
    Comp_prob_sst(1,:,:) = CP_prob_sst

	;-----olr
	EP_ElNino_avg_olr = dim_avg_n_Wrap(olr_season({EP_ElNino_year_olr},:,:),0)
	CP_ElNino_avg_olr = dim_avg_n_Wrap(olr_season({CP_ElNino_year_olr},:,:),0)
	EP_LaNina_avg_olr = dim_avg_n_Wrap(olr_season({EP_LaNina_year_olr},:,:),0)
	CP_LaNina_avg_olr = dim_avg_n_Wrap(olr_season({CP_LaNina_year_olr},:,:),0)
	EP_avg_olr = (EP_ElNino_avg_olr - EP_LaNina_avg_olr)/2
	copy_VarCoords(EP_ElNino_avg_olr,EP_avg_olr)
	CP_avg_olr = (CP_ElNino_avg_olr - CP_LaNina_avg_olr)/2
	copy_VarCoords(CP_ElNino_avg_olr,CP_avg_olr)
	;-----combine
    Comp_avg_olr = new((/2, nlat_olr, nlon_olr/),"float")
    Comp_avg_olr(0,:,:) = EP_avg_olr
    Comp_avg_olr(1,:,:) = CP_avg_olr
	
	;-----precip
	;-----EP ENSO
    EP_ElNino_avg_precip = dim_avg_n_Wrap(precip_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg_precip = dim_avg_n_Wrap(precip_season({EP_LaNina_year},:,:),0)
	EP_avg_precip = (EP_ElNino_avg_precip - EP_LaNina_avg_precip)/2
	copy_VarCoords(EP_ElNino_avg_precip,EP_avg_precip)
    EP_ElNino_var_precip = dim_variance_n_Wrap(precip_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num_precip = dimsizes(EP_ElNino_year)
    EP_LaNina_var_precip = dim_variance_n_Wrap(precip_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num_precip = dimsizes(EP_LaNina_year)
    EP_prob_precip = 1. - ttest(EP_ElNino_avg_precip, EP_ElNino_var_precip, EP_ElNino_num_precip, EP_LaNina_avg_precip, EP_LaNina_var_precip, EP_LaNina_num_precip, iflag, tval_opt)
	copy_VarCoords(EP_avg_precip,EP_prob_precip)
	;-----CP ENSO
    CP_ElNino_avg_precip = dim_avg_n_Wrap(precip_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg_precip = dim_avg_n_Wrap(precip_season({CP_LaNina_year},:,:),0)
	CP_avg_precip = (CP_ElNino_avg_precip - CP_LaNina_avg_precip)/2
	copy_VarCoords(CP_ElNino_avg_precip,CP_avg_precip)
    CP_ElNino_var_precip = dim_variance_n_Wrap(precip_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num_precip = dimsizes(CP_ElNino_year)
    CP_LaNina_var_precip = dim_variance_n_Wrap(precip_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num_precip = dimsizes(CP_LaNina_year)
    CP_prob_precip = 1. - ttest(CP_ElNino_avg_precip, CP_ElNino_var_precip, CP_ElNino_num_precip, CP_LaNina_avg_precip, CP_LaNina_var_precip, CP_LaNina_num_precip, iflag, tval_opt)
	copy_VarCoords(CP_avg_precip,CP_prob_precip)
    ;-----combine
    Comp_avg_precip = new((/2, nlat_precip, nlon_precip/),"float")
    Comp_prob_precip = new((/2, nlat_precip, nlon_precip/),"float")
    Comp_avg_precip(0,:,:) = EP_avg_precip
    Comp_avg_precip(1,:,:) = CP_avg_precip
    Comp_prob_precip(0,:,:) = EP_prob_precip
    Comp_prob_precip(1,:,:) = CP_prob_precip
	
	Comp_prob_precip = where(abs(Comp_avg_precip).lt.5,Comp_prob_precip@_FillValue,Comp_prob_precip)
	Comp_prob_sst = where(abs(Comp_avg_sst).lt.0.1,Comp_prob_sst@_FillValue,Comp_prob_sst)

;-----plot figure
	wks = gsn_open_wks(pltType,pltName)
	res_sst = True
	res_sst@gsnDraw = False
	res_sst@gsnFrame = False
	res_sst@gsnAddCyclic = True
	res_sst@gsnMaximize = False
	res_sst@lbLabelBarOn         = False
	res_sst@cnFillOn = True
	res_sst@cnLinesOn = False
	res_sst@cnLineLabelsOn = False
	res_sst@cnInfoLabelOn = False
	res_sst@tmXTOn = False
	res_sst@tmYROn = False

	;------plot olr
	res_olr = res_sst
	res_olr@cnLinesOn = True
	res_olr@cnFillOn = False
	res_olr@cnLevelSelectionMode = "ExplicitLevels"
	olr_level = (/-20,-16,-12,-8,-4,4,8,12,16,20/)
	res_olr@cnLevels = olr_level
	;res_olr@gsnContourNegLineDashPattern = 14       ; sets negative contours to dash pattern 1
	res_olr@cnLineThicknessF = 3.
	res_olr@cnSmoothingOn = True
	res_olr@cnLineLabelsOn = True
	res_olr@cnExplicitLineLabelsOn= True
	res_olr@cnNoDataLabelOn  = False
	res_olr@cnLineLabelPlacementMode = "Computed"
	res_olr@cnLineLabelDensityF  = 2.             ; increase the number of line labels/line
	res_olr@cnLineLabelInterval = 1
	res_olr@cnLineLabelFontHeightF = FontHeightF
	res_olr@cnLineLabelStrings = olr_level+""
	res_olr@cnLabelDrawOrder = "PostDraw"
	
	;------plot prob
	res_sst_prob = res_sst
	res_sst_prob@lbLabelBarOn = False
	res_sst_prob@cnLevelSelectionMode= "ExplicitLevels"
	res_sst_prob@cnLevels              = (/0.95/)
	res_sst_prob@cnMonoFillPattern = False 
	res_sst_prob@cnMonoFillColor = False
	res_sst_prob@cnFillColors = (/"white","yellow"/)
	res_sst_prob@cnFillPatterns        = (/-1,17/) ;
	res_sst_prob@cnFillDotSizeF = 0.003 ; dots size
    res_sst_prob@cnFillScaleF = 1. ;dots density

	res_sst@mpShapeMode = "FreeAspect"
	res_sst@mpCenterLonF = 180
    res_sst@mpMinLatF = latS
    res_sst@mpMaxLatF = latN
	cmap = read_colormap_file("MPL_BrBG")
	res_sst@mpLandFillColor = cmap(65-2,:)
	res_sst@mpGeophysicalLineColor = "grey30"
	res_sst@mpOutlineDrawOrder = "Draw"
	res_sst@cnFillDrawOrder = "PreDraw"
	res_sst@cnLevelSelectionMode = "ExplicitLevels"
	res_sst@cnLevels = (/-1.5,-1.,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,1.,1.5/)
	res_sst@cnFillPalette = "NCV_blu_red"
	res_sst@tmYLLabelFontHeightF = FontHeightF
	res_sst@tmXBLabelFontHeightF = FontHeightF
	res_sst@gsnMajorLonSpacing = 60
	res_sst@gsnLeftStringOrthogonalPosF = 0.
	res_sst@gsnRightString = ""
	res_sst@gsnRightStringOrthogonalPosF = 0.
	res_sst@gsnRightStringFontHeightF = FontHeightF
	res_sst@gsnLeftStringFontHeightF = FontHeightF
	
	res_sst@gsnLeftString = "(a) EP-ENSO "+season+": SST OLR"
	res_sst@vpXF                 = 0.06                  ; position and size
	res_sst@vpYF                 = 0.85                  ; for contour plots
	res_sst@vpWidthF             = 0.58
	res_sst@vpHeightF            = 0.25
	plot0 = gsn_csm_contour_map_ce(wks,Comp_avg_sst(0,:,:),res_sst)
	plot_prob0 = gsn_csm_contour(wks,Comp_prob_sst(0,:,:),res_sst_prob)
	overlay(plot0,plot_prob0)
	plot_olr0 = gsn_csm_contour(wks,Comp_avg_olr(0,:,:),res_olr) ;OLR
	plot_olr0 = ColorNegDashZeroPosContour(plot_olr0,"blue","transparent","magenta") 
	overlay(plot0,plot_olr0) 
	
	res_sst@gsnLeftString = "(c) CP-ENSO "+season+": SST OLR"
	res_sst@vpYF                 = 0.52
	res_sst@lbLabelBarOn         = True
	res_sst@pmLabelBarWidthF = 0.55
	res_sst@pmLabelBarHeightF = 0.05
	res_sst@pmLabelBarOrthogonalPosF = 0.15
	res_sst@lbLabelFontHeightF = 0.012
	plot1 = gsn_csm_contour_map_ce(wks,Comp_avg_sst(1,:,:),res_sst)
	plot_prob1 = gsn_csm_contour(wks,Comp_prob_sst(1,:,:),res_sst_prob)
	overlay(plot1,plot_prob1)
	plot_olr1 = gsn_csm_contour(wks,Comp_avg_olr(1,:,:),res_olr) ;OLR
	plot_olr1 = ColorNegDashZeroPosContour(plot_olr1,"blue","transparent","magenta")
	overlay(plot1,plot_olr1) 
	
	;-----precip
	res_precip = True
	res_precip@gsnDraw = False
	res_precip@gsnFrame = False
	res_precip@gsnAddCyclic = False
	res_precip@gsnMaximize = False
	res_precip@lbLabelBarOn         = False
	res_precip@cnFillOn = True
	res_precip@cnLinesOn = False
	res_precip@cnLineLabelsOn = False
	res_precip@cnInfoLabelOn = False
	res_precip@tmXTOn = False
	res_precip@tmYROn = False
	;------plot prob
	res_precip_prob = res_precip
	res_precip_prob@lbLabelBarOn = False
	res_precip_prob@cnLevelSelectionMode= "ExplicitLevels"
	res_precip_prob@cnLevels              = (/0.95,0.99/)
	res_precip_prob@cnMonoFillPattern = False 
	res_precip_prob@cnMonoFillColor = False
	res_precip_prob@cnFillColors = (/"white","grey","yellow"/)
	res_precip_prob@cnFillPatterns        = (/-1,17,17/) ;
	res_precip_prob@cnFillDotSizeF = 0.0035 ; dots size
    res_precip_prob@cnFillScaleF = 1. ;dots density
	
    res_precip@gsnMajorLonSpacing = 15
    res_precip@mpMinLatF = lat_precip(0)
    res_precip@mpMaxLatF = lat_precip(dimsizes(lat_precip)-1)
    res_precip@mpMinLonF = lon_precip(0)
    res_precip@mpMaxLonF = lon_precip(dimsizes(lon_precip)-1)
	res_precip@cnLevelSelectionMode = "ExplicitLevels"
	res_precip@cnLevels = (/-20,-15,-10,-5,5,10,15,20/)
	res_precip@cnFillPalette = "MPL_BrBG"
	res_precip@tmYLLabelFontHeightF = FontHeightF
	res_precip@tmXBLabelFontHeightF = FontHeightF
	res_precip@gsnLeftStringOrthogonalPosF = 0.
	; res_precip@gsnRightString = "mm/mon"
	; res_precip@gsnRightStringOrthogonalPosF = 0.
	; res_precip@gsnRightStringFontHeightF = FontHeightF
	res_precip@gsnLeftStringFontHeightF = FontHeightF
    res_precip@mpFillOn = True
    res_precip@mpDataSetName = "Earth..4"
    res_precip@mpOutlineBoundarySets = "AllBoundaries"
	res_precip@mpProvincialLineColor  = "grey"
    res_precip@cnFillDrawOrder = "PreDraw"
	res_precip_prob@cnFillDrawOrder = "PreDraw"
    res_precip@mpOceanFillColor = "white" 
    res_precip@mpFillDrawOrder = "Draw"
    res_precip@mpLandFillColor = "transparent"
    res_precip@mpDataBaseVersion  = "MediumRes"
    res_precip@mpOutlineDrawOrder = "Draw"
	res_precip@mpGeophysicalLineColor = "grey30"
	res_precip@mpShapeMode = "FreeAspect"
	
	res_precip@gsnLeftString = "(b) EP-ENSO "+season+": Rainfall"
	res_precip@vpXF                 = 0.71                  ; position and size
	res_precip@vpYF                 = 0.85                  ; for contour plots
	res_precip@vpWidthF             = 0.28
	res_precip@vpHeightF            = 0.25
	plot_precip0 = gsn_csm_contour_map_ce(wks,Comp_avg_precip(0,:,:),res_precip)
	plot_precip_prob0 = gsn_csm_contour(wks,Comp_prob_precip(0,:,:),res_precip_prob)
	overlay(plot_precip0,plot_precip_prob0)	

	res_precip@gsnLeftString = "(d) CP-ENSO "+season+": Rainfall"
	res_precip@lbLabelBarOn         = True
	res_precip@pmLabelBarWidthF = 0.28
	res_precip@pmLabelBarHeightF = 0.05
	res_precip@pmLabelBarOrthogonalPosF = 0.15
	res_precip@lbLabelFontHeightF = 0.012
	res_precip@vpYF                 = 0.52                  ; for contour plots
	plot_precip1 = gsn_csm_contour_map_ce(wks,Comp_avg_precip(1,:,:),res_precip)
	plot_precip_prob1 = gsn_csm_contour(wks,Comp_prob_precip(1,:,:),res_precip_prob)
	overlay(plot_precip1,plot_precip_prob1)	
	    
	; drawNDCGrid(wks) 
	
	maximize_output(wks,True)  