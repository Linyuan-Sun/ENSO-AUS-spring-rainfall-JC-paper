year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)
level = 200

season = "SON"

	;-----SON define
    EP_ElNino_year = (/1952, 1958, 1973, 1977, 1983, 1988, 1998, 2007, 2016, 2024/)
    CP_ElNino_year = (/1954, 1964, 1966, 1969, 1978, 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
    
pltType = "pdf"
pltName = "Fig8_Z200_ElNino"

;-----load data
	f = addfile("/g/data/if69/ls3248/data_HF/ERA5/monthly/hgt/hgt.ERA5.mon."+level+"hPa.anom_slidingClim.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	hgt = rm_single_dims(f->hgt(start_ind:end_ind,:,:,:))
	hgt_season = rm_single_dims(month_to_seasonN(hgt,season))
    nlat = dimsizes(hgt&latitude)
    nlon = dimsizes(hgt&longitude)
    hgt_season&time = total_year
    
;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
	;-----EP
	EP_ElNino_avg = dim_avg_n_Wrap(hgt_season({EP_ElNino_year},:,:),0)
	EP_avg = EP_ElNino_avg
    EP_ElNino_var = dim_variance_n_Wrap(hgt_season({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
	total_avg = dim_avg_n_Wrap(hgt_season,0)
    total_var = dim_variance_n_Wrap(hgt_season,0)
    total_num = dimsizes(total_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	
	;-----CP
	CP_ElNino_avg = dim_avg_n_Wrap(hgt_season({CP_ElNino_year},:,:),0)
	CP_avg = CP_ElNino_avg
    CP_ElNino_var = dim_variance_n_Wrap(hgt_season({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
	total_avg = dim_avg_n_Wrap(hgt_season,0)
    total_var = dim_variance_n_Wrap(hgt_season,0)
    total_num = dimsizes(total_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, total_avg, total_var, total_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)

    ;-----combine
    Comp_avg_hgt = new((/2, nlat, nlon/),"float")
    Comp_prob_hgt = new((/2, nlat, nlon/),"float")
    Comp_avg_hgt(0,:,:) = EP_avg
    Comp_avg_hgt(1,:,:) = CP_avg
    Comp_prob_hgt(0,:,:) = EP_prob
    Comp_prob_hgt(1,:,:) = CP_prob
    
;-----plot figure
	wks = gsn_open_wks(pltType,pltName)
	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@gsnAddCyclic = True
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
	res_prob@cnFillDotSizeF = 0.0025 ; dots size
    res_prob@cnFillScaleF = 0.6 ;dots density
	
	res@mpCenterLonF = 180
	res@mpShapeMode = "FreeAspect"
	res@vpHeightF         = 0.5                    ; change aspect ratio of plot
	res@vpWidthF          = 1.
	res@cnLevelSelectionMode = "ManualLevels"
	res@cnLevelSpacingF = 5
	res@cnMinLevelValF = -40
	res@cnMaxLevelValF = 40
	res@cnFillPalette = "BlueDarkRed18"
	res@tmYLLabelFontHeightF = 0.025
	res@tmXBLabelFontHeightF = 0.025
	res@gsnMajorLonSpacing = 60
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "gpm"
	res@gsnRightStringOrthogonalPosF = 0.
	
	plot = new(2, graphic)
	plot_prob = new(2, graphic)
	LeftString = (/"(a) EP-El Nino "+season+": Z200","(b) CP-El Nino "+season+": Z200"/)
	do i = 0, 1
		res@gsnLeftString = LeftString(i)
		plot(i) = gsn_csm_contour_map_ce(wks,Comp_avg_hgt(i,:,:),res)
		plot_prob(i) = gsn_csm_contour(wks,Comp_prob_hgt(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  