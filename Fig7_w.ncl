year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)

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
pltName = "Fig7_w"

;-----load data
	f = addfile("../Fig4_w_out.nc","r")
	time = cd_calendar(f->time,1)
	start_ind = ind(time.eq.((year_start-1)*100+1))
	end_ind = ind(time.eq.((year_end-1)*100+12))
	w_temp = f->w(start_ind:end_ind,:,:,:)
	w = w_temp*100. ;Pa/s -> 10^-2Pa/s
	copy_VarCoords(w_temp,w)
	w_season = rm_single_dims(month_to_seasonN(w,season))
    nlevel = dimsizes(w&level)
    nlon = dimsizes(w&longitude)
    w_season&time = total_year

	w_season_mm = dim_avg_n_Wrap(w_season,2)
    
;-----Composite and T-test
	iflag    = True      
	tval_opt = False 
	;-----EP ENSO
    EP_ElNino_avg = dim_avg_n_Wrap(w_season_mm({EP_ElNino_year},:,:),0)
	EP_LaNina_avg = dim_avg_n_Wrap(w_season_mm({EP_LaNina_year},:,:),0)
	EP_avg = (EP_ElNino_avg - EP_LaNina_avg)/2
	copy_VarCoords(EP_ElNino_avg,EP_avg)
    EP_ElNino_var = dim_variance_n_Wrap(w_season_mm({EP_ElNino_year},:,:),0)
    EP_ElNino_num = dimsizes(EP_ElNino_year)
    EP_LaNina_var = dim_variance_n_Wrap(w_season_mm({EP_LaNina_year},:,:),0)
    EP_LaNina_num = dimsizes(EP_LaNina_year)
    EP_prob = 1. - ttest(EP_ElNino_avg, EP_ElNino_var, EP_ElNino_num, EP_LaNina_avg, EP_LaNina_var, EP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(EP_avg,EP_prob)
	
	;-----CP ENSO
    CP_ElNino_avg = dim_avg_n_Wrap(w_season_mm({CP_ElNino_year},:,:),0)
	CP_LaNina_avg = dim_avg_n_Wrap(w_season_mm({CP_LaNina_year},:,:),0)
	CP_avg = (CP_ElNino_avg - CP_LaNina_avg)/2
	copy_VarCoords(CP_ElNino_avg,CP_avg)
    CP_ElNino_var = dim_variance_n_Wrap(w_season_mm({CP_ElNino_year},:,:),0)
    CP_ElNino_num = dimsizes(CP_ElNino_year)
    CP_LaNina_var = dim_variance_n_Wrap(w_season_mm({CP_LaNina_year},:,:),0)
    CP_LaNina_num = dimsizes(CP_LaNina_year)
    CP_prob = 1. - ttest(CP_ElNino_avg, CP_ElNino_var, CP_ElNino_num, CP_LaNina_avg, CP_LaNina_var, CP_LaNina_num, iflag, tval_opt)
	copy_VarCoords(CP_avg,CP_prob)

    ;-----combine
    Comp_avg_w = new((/2, nlevel, nlon/),"float")
    Comp_prob_w = new((/2, nlevel, nlon/),"float")
    Comp_avg_w(0,:,:) = EP_avg
    Comp_avg_w(1,:,:) = CP_avg
    Comp_prob_w(0,:,:) = EP_prob
    Comp_prob_w(1,:,:) = CP_prob
    
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
	res@tmYRMode             = "Automatic"
	res@tiYAxisOn = False
	
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
    res_prob@cnFillScaleF = 1.2 ;dots density
	
	res@vpHeightF         = 0.35                    ; change aspect ratio of plot
	res@vpWidthF          = 1.
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = (/-1.5,-1.2,-0.9,-0.6,-0.3,-0.15,0,0.15,0.3,0.6,0.9,1.2,1.5/)
	res@cnFillPalette = "BlueDarkRed18"
	res@tmYLLabelFontHeightF = 0.03
	res@tmXBLabelFontHeightF = 0.03
	res@gsnMajorLonSpacing = 30
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = "10~S~-2~N~Pa/s"
	res@gsnRightStringOrthogonalPosF = 0.
	res@tmYLMode = "Explicit"
	res@tmYLValues = (/100,200,300,400,500,700,850,1000/)
	res@tmYLLabels = res@tmYLValues+""
	
	plot = new(2, graphic)
	plot_prob = new(2, graphic)
	LeftString = (/"(c) EP-ENSO "+season+": omega","(d) CP-ENSO "+season+": omega"/)
	do i = 0, 1
		res@gsnLeftString = LeftString(i)
		plot(i) = gsn_csm_pres_hgt(wks,Comp_avg_w(i,:,:),res)
		plot_prob(i) = gsn_csm_pres_hgt(wks,Comp_prob_w(i,:,:),res_prob)
		overlay(plot(i),plot_prob(i))
	end do
	resP                  = True   
	resP@gsnMaximize = True
	resP@gsnPanelLabelBar = True   
	resP@pmLabelBarWidthF = 0.5
	resP@pmLabelBarHeightF = 0.045
	resP@lbLabelFontHeightF = 0.012
	gsn_panel(wks,plot,(/1,2/),resP)  