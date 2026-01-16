year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)

; latS = -60
; latN = 0
; lonW = 70
; lonE = 210

latS = -60
latN = 0
lonW = 90
lonE = 210

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
pltName = "Fig8_SLP_wnd850_LaNina"

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/u/u.ERA5.mon.850hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	u = lonFlip(f->u(start_ind:end_ind,0,{latS:latN},:))
	u_season = rm_single_dims(month_to_seasonN(u(:,:,{lonW:lonE}),season))
	u_season&time = total_year
    
    f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/v/v.ERA5.mon.850hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	v = lonFlip(f->v(start_ind:end_ind,0,{latS:latN},:))
	v_season = rm_single_dims(month_to_seasonN(v(:,:,{lonW:lonE}),season))
	v_season&time = total_year
	delete(time)
    f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/SLP/SLP.mon.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	SLP_lonFlip = lonFlip(f->SLP(start_ind:end_ind,{latS:latN},:))
    SLP_lonFlip_temp = SLP_lonFlip(:,:,{lonW:lonE})
    SLP = SLP_lonFlip_temp/100. ;units: Pa -> hPa
    copy_VarCoords(SLP_lonFlip_temp,SLP)
	SLP_season = rm_single_dims(month_to_seasonN(SLP,season))
	SLP_season&time = total_year
    
    nlat = dimsizes(SLP_season&latitude)
    nlon = dimsizes(SLP_season&longitude)
    
;-----Composite and T-test
	iflag    = True      
	tval_opt = False 

	;-----SLP
	;-----EP
	EP_LaNina_avg = dim_avg_n_Wrap(SLP_season({EP_LaNina_year},:,:),0)
	EP_avg = EP_LaNina_avg
    EP_LaNina_var = dim_variance_n_Wrap(SLP_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
	total_avg = dim_avg_n_Wrap(SLP_season,0)
    total_var = dim_variance_n_Wrap(SLP_season,0)
    total_num = dimsizes(total_year)
    EP_prob = 1. - ttest(EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP
	CP_LaNina_avg = dim_avg_n_Wrap(SLP_season({CP_LaNina_year},:,:),0)
	CP_avg = CP_LaNina_avg
    CP_LaNina_var = dim_variance_n_Wrap(SLP_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
	total_avg = dim_avg_n_Wrap(SLP_season,0)
    total_var = dim_variance_n_Wrap(SLP_season,0)
    total_num = dimsizes(total_year)
    CP_prob = 1. - ttest(CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_SLP = new((/2, nlat, nlon/),"float")
    Comp_prob_SLP = new((/2, nlat, nlon/),"float")
    Comp_avg_SLP(0,:,:) = EP_avg
    Comp_avg_SLP(1,:,:) = CP_avg
    Comp_prob_SLP(0,:,:) = EP_prob
    Comp_prob_SLP(1,:,:) = CP_prob

	;-----u
	;-----EP
	EP_LaNina_avg = dim_avg_n_Wrap(u_season({EP_LaNina_year},:,:),0)
	EP_avg = EP_LaNina_avg
    EP_LaNina_var = dim_variance_n_Wrap(u_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
	total_avg = dim_avg_n_Wrap(u_season,0)
    total_var = dim_variance_n_Wrap(u_season,0)
    total_num = dimsizes(total_year)
    EP_prob = 1. - ttest(EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP
	CP_LaNina_avg = dim_avg_n_Wrap(u_season({CP_LaNina_year},:,:),0)
	CP_avg = CP_LaNina_avg
    CP_LaNina_var = dim_variance_n_Wrap(u_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
	total_avg = dim_avg_n_Wrap(u_season,0)
    total_var = dim_variance_n_Wrap(u_season,0)
    total_num = dimsizes(total_year)
    CP_prob = 1. - ttest(CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_u = new((/2, nlat, nlon/),"float")
    Comp_prob_u = new((/2, nlat, nlon/),"float")
    Comp_avg_u(0,:,:) = EP_avg
    Comp_avg_u(1,:,:) = CP_avg
    Comp_prob_u(0,:,:) = EP_prob
    Comp_prob_u(1,:,:) = CP_prob

	;-----v
	;-----EP
	EP_LaNina_avg = dim_avg_n_Wrap(v_season({EP_LaNina_year},:,:),0)
	EP_avg = EP_LaNina_avg
    EP_LaNina_var = dim_variance_n_Wrap(v_season({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
	total_avg = dim_avg_n_Wrap(v_season,0)
    total_var = dim_variance_n_Wrap(v_season,0)
    total_num = dimsizes(total_year)
    EP_prob = 1. - ttest(EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	;-----CP
	CP_LaNina_avg = dim_avg_n_Wrap(v_season({CP_LaNina_year},:,:),0)
	CP_avg = CP_LaNina_avg
    CP_LaNina_var = dim_variance_n_Wrap(v_season({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
	total_avg = dim_avg_n_Wrap(v_season,0)
    total_var = dim_variance_n_Wrap(v_season,0)
    total_num = dimsizes(total_year)
    CP_prob = 1. - ttest(CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)
    ;-----combine
    Comp_avg_v = new((/2, nlat, nlon/),"float")
    Comp_prob_v = new((/2, nlat, nlon/),"float")
    Comp_avg_v(0,:,:) = EP_avg
    Comp_avg_v(1,:,:) = CP_avg
    Comp_prob_v(0,:,:) = EP_prob
    Comp_prob_v(1,:,:) = CP_prob

	siglvl = 0.95
    Comp_prob_uv = where(Comp_prob_u.gt.siglvl .or. Comp_prob_v.gt.siglvl, 1.0, 0.0)
    copy_VarCoords(Comp_prob_u,Comp_prob_uv)
    
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
    
    res_uv = res
    vecref = 1.
	res_uv@vcGlyphStyle = "CurlyVector"
	res_uv@vcRefMagnitudeF         = vecref            
	res_uv@vcRefLengthF = 0.06
	res_uv@vcRefAnnoSide = "Bottom"                      
	res_uv@vcRefAnnoJust = "BottomRight";
	res_uv@vcMinDistanceF = 0.05
	res_uv@vcRefAnnoString1 = vecref+"m/s"
	res_uv@vcRefAnnoString2On = False            
    res_uv@vcRefAnnoFontHeightF = 0.015
	res_uv@vcRefAnnoArrowSpaceF = 0.0
	res_uv@vcLineArrowThicknessF = 2.
	res_uv@vcLevelSelectionMode = "ExplicitLevels"
	res_uv@vcLevels = (/0.5/) ;uv_prob = 1 or 0
	res_uv@vcLevelColors = (/"black","blue"/)
	res_uv@vcLineArrowHeadMaxSizeF = 0.03 
	res_uv@vcLineArrowHeadMinSizeF = 0.005
	res_uv@vcRefAnnoOrthogonalPosF = -1.    
	
    res@cnFillOn = True
	res@cnLinesOn = False
	res@cnLineLabelsOn = False
	res@cnInfoLabelOn = False
    
	;------plot_SLP prob
	res_prob = res
	res_prob@lbLabelBarOn = False
	res_prob@cnLevelSelectionMode= "ExplicitLevels"
	res_prob@cnLevels              = (/0.95/)
	res_prob@cnMonoFillPattern = False 
	res_prob@cnMonoFillColor = False
	res_prob@cnFillColors = (/"white","yellow"/)
	res_prob@cnFillPatterns        = (/-1,17/) ;
	res_prob@cnFillDotSizeF = 0.0035 ; dots size
    res_prob@cnFillScaleF = 1. ;dots density

	res@mpShapeMode = "FreeAspect"
	res@vpHeightF         = 0.5                    ; change aspect ratio of plot
	res@vpWidthF          = 1.
    res@mpMinLatF = latS
    res@mpMaxLatF = latN
    res@mpMinLonF = lonW
    res@mpMaxLonF = lonE
    res@mpCenterLonF = (lonE+lonW)/2
    ; res@mpDataBaseVersion  = "HighRes"
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = fspan(-2.4,2.4,17)
    res@cnFillPalette = "BlueDarkRed18"
	res@tmYLLabelFontHeightF = 0.025
	res@tmXBLabelFontHeightF = 0.025
	res@gsnMajorLonSpacing = 30
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "hPa"
	res@gsnRightStringOrthogonalPosF = 0.
    
    ; res@cnFillDrawOrder = "PreDraw"  
    ; res_prob@cnFillDrawOrder = "PreDraw"
    ; res@mpFillDrawOrder        = "Draw"
    ; res@mpLandFillColor = "transparent"
    ; ; res@mpOceanFillColor        = "white"      
	; res@mpGeophysicalLineColor = "blue"
	res@mpGeophysicalLineThicknessF = 1.5
	
	nplot = dimsizes(Comp_avg_SLP(:,0,0))
	plot_SLP = new(nplot, graphic)
	plot_SLP_prob = new(nplot, graphic)
	plot_uv = new(nplot, graphic)
	LeftString = (/"(g) EP-La Nina "+season+": SLP wind850","(h) CP-La Nina "+season+": SLP wind850"/)
	do i = 0, nplot-1
		res@gsnLeftString = LeftString(i)
		
		;-----SLP & prob
		plot_SLP(i) = gsn_csm_contour_map_ce(wks,Comp_avg_SLP(i,:,:),res)
		plot_SLP_prob(i) = gsn_csm_contour(wks,Comp_prob_SLP(i,:,:),res_prob)
		overlay(plot_SLP(i),plot_SLP_prob(i))
        ;-----uv
        plot_uv(i) = gsn_csm_vector_scalar(wks,Comp_avg_u(i,:,:),Comp_avg_v(i,:,:),Comp_prob_uv(i,:,:),res_uv)
		overlay(plot_SLP(i),plot_uv(i))
		
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot_SLP,(/1,2/),resP)  