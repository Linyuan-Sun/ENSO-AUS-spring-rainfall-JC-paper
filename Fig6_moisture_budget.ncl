year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)

latS = -55
latN = -5
lonW = 90
lonE = 190

; latS = -44.5
; latN = -10
; lonW = 112
; lonE = 156

season = "SON"

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
pltName = "Fig6_moisture_budget"

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/moisture_budget/moisture_budget.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	
	mtpr_temp = lonFlip(f->mtpr(start_ind:end_ind,:,:))
	mtpr = mtpr_temp*1000000.
	copy_VarCoords(mtpr_temp,mtpr)
	
	mvimd_temp = lonFlip(f->mvimd(start_ind:end_ind,:,:))
	mvimd = mvimd_temp*-1000000.
	copy_VarCoords(mvimd_temp,mvimd)
	
	mer_temp = lonFlip(f->mer(start_ind:end_ind,:,:))
	mer = mer_temp*-1000000.
	copy_VarCoords(mer_temp,mer)

	mtpr_season = rm_single_dims(month_to_seasonN(mtpr(:,{latS:latN},{lonW:lonE}),season))
	mtpr_season&time = total_year
	
	mvimd_season = rm_single_dims(month_to_seasonN(mvimd(:,{latS:latN},{lonW:lonE}),season))
	mvimd_season&time = total_year
	
	mer_season = rm_single_dims(month_to_seasonN(mer(:,{latS:latN},{lonW:lonE}),season))
	mer_season&time = total_year
	
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/quv/quv.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	qu = lonFlip(f->qu(start_ind:end_ind,{latS:latN},:))
	qu_season = rm_single_dims(month_to_seasonN(qu(:,:,{lonW:lonE}),season))
	qu_season&time = total_year
    
	qv = lonFlip(f->qv(start_ind:end_ind,{latS:latN},:))
	qv_season = rm_single_dims(month_to_seasonN(qv(:,:,{lonW:lonE}),season))
	qv_season&time = total_year
	
    lat = mtpr_season&latitude
    lon = mtpr_season&longitude
    nlat = dimsizes(lat)
    nlon = dimsizes(lon)

;-----Composite and T-test
	iflag    = True      
	tval_opt = False
	
	;-----mtpr
	;-----EP ENSO
    EP_ElNino_avg_mtpr = dim_avg_n_Wrap(mtpr_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg_mtpr = dim_avg_n_Wrap(mtpr_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg_mtpr - EP_LaNina_avg_mtpr)/2
	copy_VarCoords(EP_ElNino_avg_mtpr,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(mtpr_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(mtpr_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg_mtpr, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg_mtpr, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg_mtpr = dim_avg_n_Wrap(mtpr_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg_mtpr = dim_avg_n_Wrap(mtpr_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg_mtpr - CP_LaNina_avg_mtpr)/2
	copy_VarCoords(CP_ElNino_avg_mtpr,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(mtpr_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(mtpr_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg_mtpr, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg_mtpr, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_mtpr = new((/2, nlat, nlon/),"float")
    Comp_prob_mtpr = new((/2, nlat, nlon/),"float")
    Comp_avg_mtpr(0,:,:) = EP_avg
    Comp_avg_mtpr(1,:,:) = CP_avg
    Comp_prob_mtpr(0,:,:) = EP_prob
    Comp_prob_mtpr(1,:,:) = CP_prob
	Comp_prob_mtpr = where(abs(Comp_avg_mtpr).lt.0.5,Comp_prob_mtpr@_FillValue,Comp_prob_mtpr)

	;-----mvimd
	;-----EP ENSO
    EP_ElNino_avg_mvimd = dim_avg_n_Wrap(mvimd_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg_mvimd = dim_avg_n_Wrap(mvimd_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg_mvimd - EP_LaNina_avg_mvimd)/2
	copy_VarCoords(EP_ElNino_avg_mvimd,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(mvimd_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(mvimd_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg_mvimd, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg_mvimd, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg_mvimd = dim_avg_n_Wrap(mvimd_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg_mvimd = dim_avg_n_Wrap(mvimd_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg_mvimd - CP_LaNina_avg_mvimd)/2
	copy_VarCoords(CP_ElNino_avg_mvimd,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(mvimd_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(mvimd_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg_mvimd, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg_mvimd, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_mvimd = new((/2, nlat, nlon/),"float")
    Comp_prob_mvimd = new((/2, nlat, nlon/),"float")
    Comp_avg_mvimd(0,:,:) = EP_avg
    Comp_avg_mvimd(1,:,:) = CP_avg
    Comp_prob_mvimd(0,:,:) = EP_prob
    Comp_prob_mvimd(1,:,:) = CP_prob
	Comp_prob_mvimd = where(abs(Comp_avg_mvimd).lt.0.5,Comp_prob_mvimd@_FillValue,Comp_prob_mvimd)
	
	;-----mer
	;-----EP ENSO
    EP_ElNino_avg_mer = dim_avg_n_Wrap(mer_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg_mer = dim_avg_n_Wrap(mer_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg_mer - EP_LaNina_avg_mer)/2
	copy_VarCoords(EP_ElNino_avg_mer,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(mer_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(mer_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg_mer, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg_mer, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg_mer = dim_avg_n_Wrap(mer_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg_mer = dim_avg_n_Wrap(mer_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg_mer - CP_LaNina_avg_mer)/2
	copy_VarCoords(CP_ElNino_avg_mer,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(mer_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(mer_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg_mer, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg_mer, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_mer = new((/2, nlat, nlon/),"float")
    Comp_prob_mer = new((/2, nlat, nlon/),"float")
    Comp_avg_mer(0,:,:) = EP_avg
    Comp_avg_mer(1,:,:) = CP_avg
    Comp_prob_mer(0,:,:) = EP_prob
    Comp_prob_mer(1,:,:) = CP_prob
	Comp_prob_mer = where(abs(Comp_avg_mer).lt.0.5,Comp_prob_mer@_FillValue,Comp_prob_mer)
	
	;-----qu
	;-----EP ENSO
    EP_ElNino_avg_qu = dim_avg_n_Wrap(qu_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg_qu = dim_avg_n_Wrap(qu_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg_qu - EP_LaNina_avg_qu)/2
	copy_VarCoords(EP_ElNino_avg_qu,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(qu_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(qu_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg_qu, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg_qu, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg_qu = dim_avg_n_Wrap(qu_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg_qu = dim_avg_n_Wrap(qu_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg_qu - CP_LaNina_avg_qu)/2
	copy_VarCoords(CP_ElNino_avg_qu,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(qu_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(qu_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg_qu, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg_qu, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_qu = new((/2, nlat, nlon/),"float")
    Comp_prob_qu = new((/2, nlat, nlon/),"float")
    Comp_avg_qu(0,:,:) = EP_avg
    Comp_avg_qu(1,:,:) = CP_avg
    Comp_prob_qu(0,:,:) = EP_prob
    Comp_prob_qu(1,:,:) = CP_prob
	Comp_prob_qu = where(abs(Comp_avg_qu).lt.0.5,Comp_prob_qu@_FillValue,Comp_prob_qu)
	
	;-----qv
	;-----EP ENSO
    EP_ElNino_avg_qv = dim_avg_n_Wrap(qv_season({EP_ElNino_year},:,:),0)
	EP_LaNina_avg_qv = dim_avg_n_Wrap(qv_season({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg_qv - EP_LaNina_avg_qv)/2
	copy_VarCoords(EP_ElNino_avg_qv,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(qv_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(qv_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg_qv, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg_qv, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP ENSO
    CP_ElNino_avg_qv = dim_avg_n_Wrap(qv_season({CP_ElNino_year},:,:),0)
	CP_LaNina_avg_qv = dim_avg_n_Wrap(qv_season({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg_qv - CP_LaNina_avg_qv)/2
	copy_VarCoords(CP_ElNino_avg_qv,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(qv_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(qv_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg_qv, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg_qv, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_qv = new((/2, nlat, nlon/),"float")
    Comp_prob_qv = new((/2, nlat, nlon/),"float")
    Comp_avg_qv(0,:,:) = EP_avg
    Comp_avg_qv(1,:,:) = CP_avg
    Comp_prob_qv(0,:,:) = EP_prob
    Comp_prob_qv(1,:,:) = CP_prob
	Comp_prob_qv = where(abs(Comp_avg_qv).lt.0.5,Comp_prob_qv@_FillValue,Comp_prob_qv)
	
	siglvl = 0.95
    Comp_prob_quv = where(Comp_prob_qu.gt.siglvl .or. Comp_prob_qv.gt.siglvl, 1.0, 0.0)
    copy_VarCoords(Comp_prob_qu,Comp_prob_quv)
    
;-----plot figure
	wks = gsn_open_wks(pltType,pltName)
	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@gsnAddCyclic = False
	res@gsnMaximize = False
	res@lbLabelBarOn         = False
	res@tmXTOn = False
	res@tmYROn = False
    
    res_quv = res
    vecref = 20.
	res_quv@vcGlyphStyle = "CurlyVector"
	res_quv@vcRefMagnitudeF         = vecref            
	res_quv@vcRefLengthF = 0.08
	res_quv@vcRefAnnoSide = "Bottom"                      
	res_quv@vcRefAnnoJust = "BottomRight";
	res_quv@vcMinDistanceF = 0.05
	res_quv@vcRefAnnoString1 = vecref+" kg/m/s"
	res_quv@vcRefAnnoString2On = False            
    res_quv@vcRefAnnoFontHeightF = 0.015
	res_quv@vcRefAnnoArrowSpaceF = 0.0
	res_quv@vcLineArrowThicknessF = 2.5    
	res_quv@vcLevelSelectionMode = "ExplicitLevels"
	res_quv@vcLevels = (/0.5/) ;quv_prob = 1 or 0
	res_quv@vcLevelColors = (/"black","blue"/)
	res_quv@vcLineArrowHeadMaxSizeF = 0.05 
	res_quv@vcLineArrowHeadMinSizeF = 0.005
	res_quv@vcRefAnnoOrthogonalPosF = -1.    
	
    res@cnFillOn = True
	res@cnLinesOn = False
	res@cnLineLabelsOn = False
	res@cnInfoLabelOn = False
    
	;------plot_quv_div prob
	res_prob = res
	res_prob@lbLabelBarOn = False
	res_prob@cnLevelSelectionMode= "ExplicitLevels"
	res_prob@cnLevels              = (/0.95/)
	res_prob@cnMonoFillPattern = False 
	res_prob@cnMonoFillColor = False
	res_prob@cnFillColors = (/"white","yellow"/)
	res_prob@cnFillPatterns        = (/-1,17/) ;
	res_prob@cnFillDotSizeF = 0.003 ; dots size
    res_prob@cnFillScaleF = 1. ;dots density

	; res@mpShapeMode = "FreeAspect"
	; res@vpWidthF = 1.0
	; res@vpHeightF = 0.35
    res@mpMinLatF = latS
    res@mpMaxLatF = latN
    res@mpMinLonF = lonW
    res@mpMaxLonF = lonE
    res@mpCenterLonF = (lonE+lonW)/2
    res@mpDataBaseVersion  = "HighRes"
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-20,-15,-10,-5,-3,-2,-1,1,2,3,5,10,15,20/)
    res@cnFillPalette = "MPL_BrBG"
	res@tmYLLabelFontHeightF = 0.025
	res@tmXBLabelFontHeightF = 0.025
	res@gsnMajorLonSpacing = 30
	res@gsnMajorLatSpacing = 15
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnLeftStringFontHeightF = 0.025
	; res@gsnRightString = "10~S~-6~N~kg/(m~S~2~N~s)"
	; res@gsnRightStringOrthogonalPosF = 0.
    
    res@cnFillDrawOrder = "PreDraw"  
    res_prob@cnFillDrawOrder = "PreDraw"
    res@mpFillDrawOrder        = "Draw"
    res@mpLandFillColor = "transparent"
    ; res@mpOceanFillColor        = "white"      
	res@mpGeophysicalLineColor = "black"
	res@mpGeophysicalLineThicknessF = 2
	
	res@gsnLeftString = "(a) EP-ENSO "+season+": Precipitation rate"
	plot_mtpr0 = gsn_csm_contour_map_ce(wks,Comp_avg_mtpr(0,:,:),res)
	plot_mtpr_prob0 = gsn_csm_contour(wks,Comp_prob_mtpr(0,:,:),res_prob)
	overlay(plot_mtpr0,plot_mtpr_prob0)		

	res@gsnLeftString = "(b) CP-ENSO "+season+": Precipitation rate"
	plot_mtpr1 = gsn_csm_contour_map_ce(wks,Comp_avg_mtpr(1,:,:),res)
	plot_mtpr_prob1 = gsn_csm_contour(wks,Comp_prob_mtpr(1,:,:),res_prob)
	overlay(plot_mtpr1,plot_mtpr_prob1)		
	
	res@gsnLeftString = "(c) EP-ENSO "+season+": Moisture flux convergence"
	plot_mvimd0 = gsn_csm_contour_map_ce(wks,Comp_avg_mvimd(0,:,:),res)
	plot_mvimd_prob0 = gsn_csm_contour(wks,Comp_prob_mvimd(0,:,:),res_prob)
	overlay(plot_mvimd0,plot_mvimd_prob0)	
	plot_quv0 = gsn_csm_vector_scalar(wks,Comp_avg_qu(0,:,:),Comp_avg_qv(0,:,:),Comp_prob_quv(0,:,:),res_quv)
	overlay(plot_mvimd0,plot_quv0)
	
	res@gsnLeftString = "(d) CP-ENSO "+season+": Moisture flux convergence"
	plot_mvimd1 = gsn_csm_contour_map_ce(wks,Comp_avg_mvimd(1,:,:),res)
	plot_mvimd_prob1 = gsn_csm_contour(wks,Comp_prob_mvimd(1,:,:),res_prob)
	overlay(plot_mvimd1,plot_mvimd_prob1)
	plot_quv1 = gsn_csm_vector_scalar(wks,Comp_avg_qu(1,:,:),Comp_avg_qv(1,:,:),Comp_prob_quv(1,:,:),res_quv)
	overlay(plot_mvimd1,plot_quv1)
	
	res@gsnLeftString = "(e) EP-ENSO "+season+": Evaporation rate"
	plot_mer0 = gsn_csm_contour_map_ce(wks,Comp_avg_mer(0,:,:),res)
	plot_mer_prob0 = gsn_csm_contour(wks,Comp_prob_mer(0,:,:),res_prob)
	overlay(plot_mer0,plot_mer_prob0)	
	
	res@gsnLeftString = "(f) CP-ENSO "+season+": Evaporation rate"
	plot_mer1 = gsn_csm_contour_map_ce(wks,Comp_avg_mer(1,:,:),res)
	plot_mer_prob1 = gsn_csm_contour(wks,Comp_prob_mer(1,:,:),res_prob)
	overlay(plot_mer1,plot_mer_prob1)
	
	plot = new(6,graphic)
	plot(0) = plot_mtpr0
	plot(1) = plot_mtpr1
	plot(2) = plot_mvimd0
	plot(3) = plot_mvimd1
	plot(4) = plot_mer0
	plot(5) = plot_mer1
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.75
	resP@pmLabelBarHeightF = 0.05
	resP@lbLabelFontHeightF = 0.015
	gsn_panel(wks,plot,(/3,2/),resP)  