year_start = 1951
year_end = 2024
total_year = ispan(year_start,year_end,1)
level = 200
lat_start = -60
lat_end = 0
lon_start = 40
lon_end = 320
season = "SON"
season_ind = 10

pltType = "pdf"
pltName = "Fig10_U_Ks"

;-----load data
	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/Ks/Ks.ERA5."+level+"hPa.nc","r")
	Ks_temp = lonFlip(f->Ks({season_ind},:,:))
	Ks = Ks_temp({lat_start:lat_end},{lon_start:lon_end})

	f = addfile("/g/data/if69/ls3248/data/ERA5/monthly/u/clim/u.ERA5.seasonal.clim.1951-2024.nc","r")
	u_temp = rm_single_dims(lonFlip(f->u({season_ind},{level},::-1,:)))
    u = u_temp({lat_start:lat_end},{lon_start:lon_end})

;-----plot figure
	wks = gsn_open_wks(pltType,pltName)
	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@gsnMaximize = False
	res@cnInfoLabelOn = False
	res@tmXTOn = False
	res@tmYROn = False
	
	res@mpMaxLatF = lat_end
	res@mpMinLatF = lat_start
	res@gsnAddCyclic = False
	res@mpMaxLonF = lon_end
	res@mpMinLonF = lon_start
	res@mpCenterLonF = (lon_start+lon_end)/2.
    
	res@tmYLLabelFontHeightF = 0.015
	res@tmXBLabelFontHeightF = 0.015
	res@gsnMajorLonSpacing = 60
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnLeftString = "(a) "+season+" climatology: U200 Ks"
	res@gsnLeftStringFontHeightF = 0.015
	res@gsnLeftStringOrthogonalPosF = 0.
	res@gsnRightString = ""

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/0,1,2,3,4,5,6,7,8/)
	res@cnFillPalette = "MPL_rainbow";"wh-bl-gr-ye-re"
	res@cnFillOn = True
	res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@pmLabelBarOrthogonalPosF = 0.2
    plot0 = gsn_csm_contour_map(wks,Ks,res)
    
    delete([/res@mpMaxLatF, res@mpMinLatF, res@mpMaxLonF, res@mpMinLonF, res@mpCenterLonF, res@gsnMajorLonSpacing, res@cnLevels/])
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/30,32.5,35,37.5,40/)
	res@cnFillOn = False
	res@cnLinesOn = True
    res@cnLineLabelsOn = True
    res@cnLineThicknessF = 2.
	res@cnLineLabelInterval = 1
    res@cnLineLabelFontHeightF = 0.01
    res@cnMonoLineDashPattern = False
    res@cnLineDashPatterns = (/0,0,0,0,0,0,0/)
	res@cnLineColor = "blue"
	res@gsnLeftString = ""
    plot1 =  gsn_csm_contour(wks,u,res)
    
    overlay(plot0,plot1)
	plres = True
	plres@gsLineColor = "black"
	plres@gsLineThicknessF = 2.
	plres@gsLineDashPattern = 1     
	line = gsn_add_polyline(wks, plot0, (/0,120,240,360/),(/-30,-30,-30,-30/), plres)

    draw(plot0)
    frame(wks)
