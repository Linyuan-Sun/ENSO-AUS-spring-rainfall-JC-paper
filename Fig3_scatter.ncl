        year_start = 1951
        year_end = 2024
        total_year = ispan(year_start,year_end,1)
        total_year@_FillValue = 9999
        season = "SON"
        seasonIndex = 10

        sigLev = 0.95
        corLev = 0.42

        EP_ElNino_year = (/1952, 1958, 1973, 1977, 1983, 1988, 1998, 2007, 2016, 2024/)
        CP_ElNino_year = (/1954, 1964, 1966, 1969, 1978, 1987, 1992, 1995, 2003, 2005, 2010, 2019/)
        EP_LaNina_year = (/1955, 1956, 1971, 1972, 1985, 1996, 2000, 2008, 2018, 2021/)
        CP_LaNina_year = (/1965, 1974, 1976, 1984, 1989, 1999, 2001, 2011, 2012, 2022, 2023/)

        EP_ENSO_year = array_append_record(EP_ElNino_year,EP_LaNina_year,0)
        CP_ENSO_year = array_append_record(CP_ElNino_year,CP_LaNina_year,0)
        nonEP_ENSO_year = total_year(get1Dindex_Exclude(total_year,EP_ENSO_year))
        nonCP_ENSO_year = total_year(get1Dindex_Exclude(total_year,CP_ENSO_year))

        city_names  = (/ "Sydney", "Melbourne", "Brisbane", "Perth", "Adelaide", "Canberra", "Hobart", "Darwin"/)
        pltType = "pdf"
        pltType@wkPaperHeightF  = 11
        pltType@wkPaperWidthF = 11
        pltName = "Fig3_scatter"
        FontHeightF = 0.04
        plot_EP_nonEPENSO = new(dimsizes(city_names),graphic)
        wks = gsn_open_wks(pltType,pltName)
        
        f = addfile("/g/data/if69/ls3248/data/ONI/ONI.nc","r")
        ONI = f->ONI({year_start-1:year_end-1},{seasonIndex})
        ONI&year = total_year

        f = addfile("/g/data/if69/ls3248/data_HF/HadISST/HadISST.mon.anom_slidingClim.nc","r")
        time = cd_calendar(f->time,1)
        start_ind = ind(time.eq.((year_start-1)*100+1))
        end_ind = ind(time.eq.((year_end-1)*100+12))
        sst = lonFlip(f->sst(start_ind:end_ind,:,:))
        Nino3 = wgt_areaave_Wrap(sst(:,{-5:5},{210:270}),1.,1.,0)
        Nino4 = wgt_areaave_Wrap(sst(:,{-5:5},{160:210}),1.,1.,0)
        Nino3_season = rm_single_dims(month_to_seasonN(Nino3,season))
        Nino4_season = rm_single_dims(month_to_seasonN(Nino4,season))
        delete(time)
        
        ; ;-----EP CP index (Ren & Jin 2011)
        ; EP_index_season = Nino3_season
        ; CP_index_season = Nino3_season
        ; do i = 0, dimsizes(Nino3_season)-1
        ; if Nino3_season(i)*Nino4_season(i) .gt. 0 then
        ; EP_index_season(i) = Nino3_season(i) - 2/5.*Nino4_season(i)
        ; CP_index_season(i) = Nino4_season(i) - 2/5.*Nino3_season(i)
        ; else
        ; EP_index_season(i) = Nino3_season(i)
        ; CP_index_season(i) = Nino4_season(i)
        ; end if
        ; end do
        ; EP_index_season&time = total_year
        ; CP_index_season&time = total_year
		
        ;-----EP CP index (Sullivan 2016)
        EP_index_season = dim_standardize(Nino3_season,0) - 0.5*dim_standardize(Nino4_season,0)
		CP_index_season = dim_standardize(Nino4_season,0) - 0.5*dim_standardize(Nino3_season,0)
		EP_index_season!0 = "time"
        EP_index_season&time = total_year
		CP_index_season!0 = "time"
        CP_index_season&time = total_year

do k = 0,dimsizes(city_names)-1

        cityIndex = k+1
        cityName = city_names(k)
        ;-----load data
        f = addfile("/g/data/if69/ls3248/data_HF/AGCD/precip/precip.mon.anom_slidingClim.city.new.nc","r")
        time = cd_calendar(f->time,1)
        start_ind = ind(time.eq.((year_start-1)*100+1))
        end_ind = ind(time.eq.((year_end-1)*100+12))
        precip = f->precip(start_ind:end_ind,{cityIndex})
        precip_season = rm_single_dims(month_to_seasonN(precip,season))
        precip_season&time = total_year
        delete(time)

        ;-----regline
        rc    = regline(EP_index_season({EP_ENSO_year}),precip_season({EP_ENSO_year}))
        EP_index_season_new = fspan(min(EP_index_season),max(EP_index_season),20)
        precip_season_EP_regline = rc*(EP_index_season_new-rc@xave) + rc@yave
        df   = rc@nptxy-2
        prob_EP = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
        cor_EP = escorc(EP_index_season({EP_ENSO_year}),precip_season({EP_ENSO_year}))

        rc    = regline(CP_index_season({CP_ENSO_year}),precip_season({CP_ENSO_year}))
        CP_index_season_new = fspan(min(CP_index_season),max(CP_index_season),20)
        precip_season_CP_regline = rc*(CP_index_season_new-rc@xave) + rc@yave
        df   = rc@nptxy-2
        prob_CP = (1 - betainc(df/(df+rc@tval^2), df/2.0, 0.5) )
        cor_CP = escorc(CP_index_season({CP_ENSO_year}),precip_season({CP_ENSO_year}))

        ;-----plot figure
        res                     = True                   ; plot mods desired
        res@gsnDraw            = False       ; Don't draw plot
        res@gsnFrame           = False       ; Don't advance frame
        res@gsnMaximize         = False                   ; maximize plot in frame
        res@xyMarkLineMode    = "Markers"
        res@xyMarkers           = 16                     ; choose type of marker 
        res@xyMarkerSizeF       = 0.015                 ; Marker size (default 0.01)
        res@trYMinF            =  -80               ; min value on y-axis
        res@trYMaxF            =  80               ; max value on y-axis
        res@trXMinF            =  -2.8                  ; min value on x-axis
        res@trXMaxF            =  2.8                 ; max value on x-axis
        res@xyLineThicknessF = 4.
        res@tmBorderThicknessF = 1.
        res@tmXTOn = False
        res@tmXBMajorOutwardLengthF = 0.0
        res@tmXBMinorOutwardLengthF = 0.0
        res@tmXBLabelFontHeightF = FontHeightF
        res@tmYLLabelFontHeightF = FontHeightF
        res@tiYAxisString = "";"SON rainfall anomaly (mm/mon)"
        res@gsnLeftString = "";cityName
		res@gsnXRefLine = 0.0
		res@gsnYRefLine = 0.0
		res@gsnXRefLineColor = "grey"
		res@gsnYRefLineColor = "grey"
		res@gsnXRefLineThicknessF = 1.5
		res@gsnYRefLineThicknessF = 1.5
		res@xyCurveDrawOrder = "PostDraw"
        ; create plot
		res@tmYLOn = True
        ; res@tmXBLabelFontHeightF = 0.1
        ; res@tmYLLabelFontHeightF = 0.1

        ;adjust the tmXB
        if k.eq.dimsizes(city_names)-1 .or. k.eq.dimsizes(city_names)-2 then
            res@tmXBLabelsOn = True
            res@tiXAxisString = "EP Nino index"
            res@tiXAxisFontHeightF = 0.1
            res@tiXAxisOffsetYF = 0.01
        else
            res@tmXBLabelsOn = False
            res@tiXAxisString = ""
        end if

        res@tiYAxisString = "Rain. anom. (mm/mon)"
        res@tiYAxisFontHeightF = 0.1
        res@tiYAxisOffsetXF = 0.01

        res@xyMarkerColor       = "grey" 
        plot_EP_nonEPENSO(k)  = gsn_csm_xy (wks,EP_index_season({nonEP_ENSO_year}),precip_season({nonEP_ENSO_year}),res) 
        res@xyMarkerColor       = "red"                  ; Marker color
        plot_EP_ElNino  = gsn_csm_xy (wks,EP_index_season({EP_ElNino_year}),precip_season({EP_ElNino_year}),res)      
        overlay(plot_EP_nonEPENSO(k),plot_EP_ElNino)
        res@xyMarkerColor       = "blue" 
        plot_EP_LaNina  = gsn_csm_xy (wks,EP_index_season({EP_LaNina_year}),precip_season({EP_LaNina_year}),res) 
        overlay(plot_EP_nonEPENSO(k),plot_EP_LaNina)
        res@xyMarkLineMode   = "Lines"
        if prob_EP.gt.sigLev then
            res@xyDashPattern = 0
        else
            res@xyDashPattern = 16
        end if
        plot_EP_regline  = gsn_csm_xy (wks,EP_index_season_new,precip_season_EP_regline,res) 
        overlay(plot_EP_nonEPENSO(k),plot_EP_regline)

        txtres = True
        txtres@txFontHeightF = 0.035;FontHeightF
        if abs(cor_EP) .gt.corLev then
            text = "Cor: "+decimalPlaces(cor_EP ,2,False)+"~S~**~N~"
        else
            text = "Cor: "+decimalPlaces(cor_EP ,2,False)
        end if
        textEP = gsn_add_text(wks,plot_EP_nonEPENSO(k),text,1.37,68,txtres)

        ;adjust the tmXB
        if k.eq.dimsizes(city_names)-1 .or. k.eq.dimsizes(city_names)-2 then
            res@tmXBLabelsOn = True
            res@tiXAxisString = "CP Nino index"
            res@tiXAxisFontHeightF = 0.1
            res@tiXAxisOffsetYF = 0.01
        else
            res@tmXBLabelsOn = False
            res@tiXAxisString = ""
        end if
        res@tiYAxisString = ""
        res@xyMarkLineMode   = "Markers"
        res@xyMarkerColor       = "grey" 
        plot_CP_nonCPENSO  = gsn_csm_xy (wks,CP_index_season({nonCP_ENSO_year}),precip_season({nonCP_ENSO_year}),res) 
        res@xyMarkerColor       = "red"                  ; Marker color
        plot_CP_ElNino  = gsn_csm_xy (wks,CP_index_season({CP_ElNino_year}),precip_season({CP_ElNino_year}),res)      
        overlay(plot_CP_nonCPENSO,plot_CP_ElNino)
        res@xyMarkerColor       = "blue" 
        plot_CP_LaNina  = gsn_csm_xy (wks,CP_index_season({CP_LaNina_year}),precip_season({CP_LaNina_year}),res) 
        overlay(plot_CP_nonCPENSO,plot_CP_LaNina)
        res@xyMarkLineMode   = "Lines"
        if prob_CP.gt.sigLev then
            res@xyDashPattern = 0
        else
            res@xyDashPattern = 16
        end if
        plot_CP_regline  = gsn_csm_xy (wks,CP_index_season_new,precip_season_CP_regline,res) 
        overlay(plot_CP_nonCPENSO,plot_CP_regline)

        if abs(cor_CP) .gt.corLev then
            text = "Cor: "+decimalPlaces(cor_CP ,2,False)+"~S~**~N~"
        else
            text = "Cor: "+decimalPlaces(cor_CP ,2,False)
        end if
        textCP = gsn_add_text(wks,plot_CP_nonCPENSO,text,1.37,68,txtres)
        
        amid = gsn_attach_plots(plot_EP_nonEPENSO(k),(/plot_CP_nonCPENSO/),True,True)
        print(cityName+" is OK.")
end do

    citytxres               = True
    citytxres@txFontHeightF = 0.015
    gsn_text_ndc(wks,"Sydney",0.28,0.99,citytxres) 
    gsn_text_ndc(wks,"Melbourne",0.78,0.99,citytxres) 

	gsn_text_ndc(wks,"Brisbane",0.28,0.75,citytxres)
	gsn_text_ndc(wks,"Perth",0.78,0.75,citytxres) 

    gsn_text_ndc(wks,"Adelaide",0.28,0.51,citytxres) 
    gsn_text_ndc(wks,"Canberra",0.78,0.51,citytxres) 
    
    gsn_text_ndc(wks,"Hobart",0.28,0.27,citytxres) 
    gsn_text_ndc(wks,"Darwin",0.78,0.27,citytxres) 
    
    resP                  = True   
    resP@gsnMaximize = True
    ; drawNDCGrid(wks)
    resP@gsnPanelYF = (/0.98,0.98,0.74,0.74,0.5,0.5,0.26,0.26/)      ; Adjust bottommost plots.
    resP@gsnPanelXF = (/0.07,0.57,0.07,0.57,0.07,0.57,0.07,0.57/)
	resP@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h"/)
	resP@amJust = "TopLeft"
	resP@gsnPanelFigureStringsFontHeightF = 0.015
	; resP@gsnPanelFigureStringsPerimOn = False
	gsn_panel(wks,plot_EP_nonEPENSO,(/4,2/),resP)  



