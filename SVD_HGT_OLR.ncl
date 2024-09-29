load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"

;olr 区域
     latS_olr   =  -10
     latN_olr   =  30
     lonL_olr   =  0
     lonR_olr   =  360

;hgt
     latS_hgt   =  30
     latN_hgt   = 70
     lonL_hgt   =  0
     lonR_hgt   = 180

    yrStrt = 1979
    yrLast = 2021
    season = "JJA"
    years  = ispan(yrStrt,yrLast,1)
    nyears = dimsizes(years)



  olr_in        = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/olr.mon.mean.nc", "r")
  time_olr          = olr_in->time
  XXXX1          = cd_calendar(time_olr,-1 )
  XXXX_start1        = ind(XXXX1.eq.197901)
  XXXX_end1        = ind(XXXX1.eq.202112)
  olr_trpcl           = olr_in->olr(XXXX_start1:XXXX_end1,{latS_olr:latN_olr},{lonL_olr:lonR_olr});截取时间
  olr_trpcl           = lonFlip(olr_trpcl)
  printVarSummary(olr_trpcl)
  olr_JJA      = month_to_season(olr_trpcl, "JJA")
  printVarSummary(olr_JJA)
  olr_JJA_1   = dim_standardize_n_Wrap(olr_JJA, 1, 0)
  printVarSummary(olr_JJA_1)



     hgt_in       = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_hgt_1959_2022_JJA_360180.nc", "r")  
     print(hgt_in)
     time_hgt     = hgt_in->time
     XXXX2        = cd_calendar(time_hgt,-1 )
     XXXX_start2  = ind(XXXX2.eq.197906)
     XXXX_end2    = ind(XXXX2.eq.202108)     
     hgt          = hgt_in->z 
     hgt         := hgt(XXXX_start2:XXXX_end2,{200},{latS_hgt:latN_hgt},{lonL_hgt:lonR_hgt}) 
     printVarSummary(hgt)  
     hgt_flt      = short2flt(hgt)
     printVarSummary(hgt_flt)
     hgt_JJA      = new((/nyears,40,181/),"float")    ;40181  17 73
     hgt_JJA!0    ="time"
     hgt_JJA&time = ispan(1979,2021,1)  
     do nn        =0,42,1
     hgt_JJA(nn,:,:)= (hgt_flt(nn*3,:,:)+hgt_flt((nn*3)+1,:,:)+hgt_flt((nn*3)+2,:,:))/(3.0*9.8)
     end do
     copy_VarCoords(hgt_flt(0,:,:),hgt_JJA(0,:,:))
     printVarSummary(hgt_JJA)
     hgt_JJA_1 = dim_standardize_n_Wrap(hgt_JJA, 1, 0)
     printVarSummary(hgt_JJA_1)




    time_olr  := years
    lat_olr   := olr_JJA&lat
    lon_olr   := olr_JJA&lon
    ntime_olr := dimsizes(time_olr)
    nlat_olr  := dimsizes(lat_olr)
    nlon_olr  := dimsizes(lon_olr)

    olr_JJA := reshape(olr_JJA_1,(/nyears,nlat_olr*nlon_olr/))
  

    time_hgt := years
    lat_hgt  := hgt_JJA&lat
    lon_hgt  := hgt_JJA&lon
    ntime_hgt = dimsizes(time_hgt)
    nlat_hgt  = dimsizes(lat_hgt)
    nlon_hgt  = dimsizes(lon_hgt)

    hgt_JJA := reshape(hgt_JJA_1,(/ntime_hgt,nlat_hgt*nlon_hgt/))


  olr_JJA!0 = "time"
  olr_JJA!1 = "col"
  hgt_JJA!0     = "time"
  hgt_JJA!1     = "col"

  ntime   = nyears                         ; # time steps
  ncols   = dimsizes(olr_JJA(0,:))          ; # columns (stations or grid pts) for S
  ncolz   = dimsizes(hgt_JJA(0,:))          ; # columns (stations or grid pts) for Z
  nsvd    = 2                          ; # svd patterns to calculate ;  [nsvd <= min(ncols, ncolz) ]
  xmsg    = -999.9                     ; missing value

  print(ncols)
  print(ncolz)

  homlft = new((/nsvd,ncols/),float)
  hetlft = new((/nsvd,ncols/),float)
  homrgt = new((/nsvd,ncolz/),float)
  hetrgt = new((/nsvd,ncolz/),float)
   printVarSummary(homlft)
   printVarSummary(hetrgt)
   printVarSummary(olr_JJA)
   printVarSummary(hgt_JJA)
   ;exit

  x = svdstd(olr_JJA(col|:,time|:),hgt_JJA(col|:,time|:),nsvd,homlft,hetlft,homrgt,hetrgt)
  print("svdcov: percent variance= " + x)

  printVarSummary(x)

  ak   = onedtond(x@ak,(/nsvd,ntime/))
  bk   = onedtond(x@bk,(/nsvd,ntime/))
  ak!0 = "sv"
  ak!1 = "time"
  bk!0 = "sv"
  bk!1 = "time"

  ccr1 = escorc(ak(0,:), bk(0,:))
  print(ccr1)
  ccr2 = escorc(ak(1,:), bk(1,:))
  print(ccr2)
  ; exit
  ; print(ak)
  ; print(bk)

  ak_std = dim_standardize_Wrap(ak,1)
  bk_std = dim_standardize_Wrap(bk,1)
  ak_std = ak_std 
  bk_std = bk_std 
  printVarSummary(ak_std)
  printVarSummary(bk_std)
  printVarSummary(homlft)
  printVarSummary(hetrgt)


  homlft     := reshape(homlft, (/nsvd,nlat_olr, nlon_olr/))
  hetrgt     := reshape(hetrgt, (/nsvd, nlat_hgt, nlon_hgt/))
  homlft     := homlft 
  hetrgt     := hetrgt  
  copy_VarCoords(olr_trpcl(0,:,:), homlft(0,:,:))
  copy_VarCoords(hgt_flt(0,:,:), hetrgt(0,:,:))
  printVarSummary(homlft)
  printVarSummary(hetrgt)



     u_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_u_1959_2022_JJA.nc", "r")
     ;print(u_in)
     u       = u_in->u
     u      := u(XXXX_start2:XXXX_end2,{200},{latS_hgt:latN_hgt},{lonL_hgt:lonR_hgt}) 
     u_flt  = short2flt(u)
     printVarSummary(u_flt)
     u_JJA = new((/nyears,161,720/),"float")
     u_JJA!0="time"
     u_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     u_JJA(nn,:,:)= (u_flt(nn*3,:,:)+u_flt((nn*3)+1,:,:)+u_flt((nn*3)+2,:,:))/3.0
     end do
     copy_VarMeta(u_flt(0,:,:),u_JJA(0,:,:))
     printVarSummary(u_JJA)

    ;---------------------------------------------------------
     v_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_v_1959_2022_JJA.nc", "r")
     ;print(v_in)
     v       = v_in->v
     v      := v(XXXX_start2:XXXX_end2,{200},{latS_hgt:latN_hgt},{lonL_hgt:lonR_hgt}) 
     v_flt  = short2flt(v)
     printVarSummary(v_flt)
     v_JJA = new((/nyears,161,720/),"float")
     v_JJA!0="time"
     v_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     v_JJA(nn,:,:)= (v_flt(nn*3,:,:)+v_flt((nn*3)+1,:,:)+v_flt((nn*3)+2,:,:))/3.0
     end do
     copy_VarMeta(v_flt(0,:,:),v_JJA(0,:,:))
     printVarSummary(v_JJA)


    u_reg_HGT1 = regCoef_n(bk_std(0,:),u_JJA(latitude|:,longitude|:,time|:),0,2)     
    u_reg_HGT1!0 = "lat"
    u_reg_HGT1!1 = "lon"
    u_reg_HGT1&lat = u_JJA&latitude
    u_reg_HGT1&lon = u_JJA&longitude
    ;printVarSummary(u_reg_HGT1)

    u_reg_HGT2 = regCoef_n(bk_std(1,:),u_JJA(latitude|:,longitude|:,time|:),0,2)     
    u_reg_HGT2!0 = "lat"
    u_reg_HGT2!1 = "lon"
    u_reg_HGT2&lat = u_JJA&latitude
    u_reg_HGT2&lon = u_JJA&longitude
    ;printVarSummary(u_reg_HGT2)

    v_reg_HGT1 = regCoef_n(bk_std(0,:),v_JJA(latitude|:,longitude|:,time|:),0,2)     
    v_reg_HGT1!0 = "lat"
    v_reg_HGT1!1 = "lon"
    v_reg_HGT1&lat = v_JJA&latitude
    v_reg_HGT1&lon = v_JJA&longitude
    ;printVarSummary(v_reg_HGT1)

    v_reg_HGT2 = regCoef_n(bk_std(1,:),v_JJA(latitude|:,longitude|:,time|:),0,2)     
    v_reg_HGT2!0 = "lat"
    v_reg_HGT2!1 = "lon"
    v_reg_HGT2&lat = v_JJA&latitude
    v_reg_HGT2&lon = v_JJA&longitude
    ;printVarSummary(v_reg_HGT2)


u_reg_HGT1=where(sqrt((/u_reg_HGT1/)^2+(/v_reg_HGT1/)^2).lt. 0.25, u_reg_HGT1@_FillValue, u_reg_HGT1)
v_reg_HGT1=where(sqrt((/u_reg_HGT1/)^2+ (/v_reg_HGT1/)^2).lt.0.25, v_reg_HGT1@_FillValue, v_reg_HGT1)

u_reg_HGT2=where(sqrt((/u_reg_HGT2/)^2+(/v_reg_HGT2/)^2).lt. 0.35, u_reg_HGT2@_FillValue, u_reg_HGT2)
v_reg_HGT2=where(sqrt((/u_reg_HGT2/)^2+ (/v_reg_HGT2/)^2).lt.0.35, v_reg_HGT2@_FillValue, v_reg_HGT2)

; ==============================================================
; Set the figure parameters
; ==============================================================

    wks  = gsn_open_wks("eps","SVD")
    plot = new(6,graphic)     ;shading

    gsn_define_colormap(wks,"cb_9step")

    res                  = True
    res@gsnDraw          = False
    res@gsnFrame         = False      ;dont advance frame
    res@gsnAddCyclic     = False     ;EOF data is not cyclic
    res@cnFillOn         = True         ; turn on color fill
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False      ;turn off cn line labels
    res@cnInfoLabelOn    = False      ;turn off contour information label

    res@tmXTOn = False
    res@tmYROn = False

    res@mpOutlineBoundarySets       ="Geophysical"
    res@mpGeophysicalLineColor      = "gray50"
    res@mpGeophysicalLineDashPattern= 0
    res@mpGeophysicalLineDashSegLenF= 0.2  
    res@mpGeophysicalLineThicknessF = 1.0  ;地图边界的粗细
    res@mpLandFillColor             = "white"
    ;res@mpOutlineDrawOrder          = "PreDraw"
    
    res@tmXBMinorOn               = True
    res@tmXBMinorThicknessF       = 0.
    res@tmXBMinorLengthF          = 0.
    res@tmXBMinorOutwardLengthF   = 0.
    res@tmXBMajorThicknessF       = 3.0
    res@tmXBMajorLengthF          = 0.003
    res@tmXBMajorOutwardLengthF   = 0.003
    res@tmXBLabelDeltaF           = -0.6
    res@tmXBLabelFontHeightF        = 0.0075
    res@tmXBLabelFontThicknessF     = 0.03    
    res@tmBorderThicknessF        = 3.0     ;边框粗细
    res@tmXBMode          = "Explicit"
    res@tmXBValues        = (/ 0, 30, 60, 90, 120 , 150, 180/)
    res@tmXBLabels        = (/"0", "30~S~o~N~E","60~S~o~N~E", "90~S~o~N~E","120~S~o~N~E","150~S~o~N~E", "180~S~o~N~"/)

    res@tmYLMinorOn               = True
    res@tmYLMinorThicknessF       = 0.
    res@tmYLMinorLengthF          = 0.
    res@tmYLMinorOutwardLengthF   = 0.
    res@tmYLMajorThicknessF       = 3.0
    res@tmYLMajorLengthF          = 0.003
    res@tmYLMajorOutwardLengthF   = 0.003
    res@tmYLLabelDeltaF           = -0.6
    res@tmYLLabelFontHeightF        = 0.0075
    res@tmYLLabelFontThicknessF     = 0.03    
    res@tmBorderThicknessF        = 3.0     ;边框粗细
    res@tmYLMode          = "Explicit"
    res@tmYLValues        = (/30,40,50,60,70/)
    res@tmYLLabels        = (/"30~S~o~N~N","40~S~o~N~N","50~S~o~N~N","60~S~o~N~N","70~S~o~N~N"/)


    res@mpShapeMode = "FreeAspect"
    res@vpHeightF  = 0.12
    res@vpWidthF   = 0.33
    res@txFont   = "helvetica"

    res@cnFillPalette             = "MPL_PuOr"  ;"NCV_blue_red"  
    res@cnLevelSelectionMode   = "ExplicitLevels"  
    res@cnFillColors           =  (/118,106,95,85,76,68, -1,-1, 60,50,40,35,30,25/) 
    res@cnSmoothingOn          = True
    ;res@pmTickMarkDisplayMode  = "Always"

    res@mpFillOn         = True        ; turn off map fill
    res@mpMinLatF                = latS_hgt
    res@mpMaxLatF                = latN_hgt
    res@mpMinLonF                =  lonL_hgt
    res@mpMaxLonF                =  lonR_hgt 
    res@mpCenterLonF             = ( lonL_hgt + lonR_hgt )/2

    res@lbLabelBarOn                 = False
    res@lbBoxLinesOn                 = False
    res@lbOrientation                = "Vertical"
    res@lbBoxEndCapStyle             = "TriangleBothEnds"     
    res@pmLabelBarWidthF             =  0.02       ;色标宽度
    res@pmLabelBarHeightF            =  0.12             ;   色标高度
    res@lbLabelStride                =  2
    res@lbLabelOffsetF               =  0.1              ;色标与标签之间的距离
    res@lbLabelFontHeightF           =  0.005           ;色标标签字体大小    
    res@pmLabelBarOrthogonalPosF     = 0.01
    ;res2@pmLabelBarParallelPosF       = 0.5 

    res@gsnLeftStringFontHeightF  = 0.01
    res@gsnLeftStringOrthogonalPosF = 0.023  
    res@cnLevels              = (/-0.4,-0.3,-0.2,-0.15,-0.1,-0.05,0,0.1,0.2,0.3,0.4,0.5,0.6/)    
    res@gsnLeftString = "~F22~MCA1_200Z"           
    plot(0) = gsn_csm_contour_map(wks,hetrgt(0,:,:),res)


    res@lbLabelBarOn                 = True
    res@gsnLeftString = "~F22~MCA2_200Z"  
    plot(1) = gsn_csm_contour_map(wks,hetrgt(1,:,:),res)


  ;--------------------200hpa回归风场叠加--------------
  resv                      = True
  resv@vcVectorDrawOrder    = "PostDraw"
  resv@gsnDraw              = False;
  resv@gsnFrame             = False;
  resv@gsnAddCyclic         = False
  resv@gsnLeftString        = ""
  resv@gsnRightString       = ""
  resv@vcGlyphStyle              = "CurlyVector"  ;其余三种选项为“LineArrow”、“WindBarb” 、“CurlyVector”
  resv@vcLineArrowThicknessF     = 1.2
  resv@txFont   = "helvetica"

  ;****设定参考箭头****
  resv@vcRefAnnoOn               = True  
  resv@vcRefLengthF              = 0.035   ;标准长度箭头在单位平方中的大小
  resv@vcRefAnnoBackgroundColor  = 0 ;背景颜色     
  resv@vcRefAnnoPerimOn          = True     ;关闭边框                                       
  resv@vcRefAnnoFontHeightF      = 0.01  ;参考箭头标签字体大小      

  resv@vcRefAnnoString1On     = False   ;设定参考箭头上、下的字符        
  resv@vcRefAnnoString2On     = True    ; 这里仅设定其下方的字符
  resv@vcRefAnnoSide            = "Bottom" ; 参考箭头放至图形上方
  resv@vcRefAnnoOrthogonalPosF  = -0.176
  resv@vcRefMagnitudeF        = 1.5    ;标准长度箭头所表示的大小
  resv@vcRefAnnoString2       = 1.5 
  resv@vcMinDistanceF            = 0.03
  plot_200uv_HGT1=gsn_csm_vector(wks,u_reg_HGT1,v_reg_HGT1,resv)
  overlay(plot(0),plot_200uv_HGT1)

  plot_200uv_HGT2=gsn_csm_vector(wks,u_reg_HGT2,v_reg_HGT2,resv)
  overlay(plot(1),plot_200uv_HGT2)

  ;---------------------AC气旋反气旋标志--------------------
  res_text_AC                                =  True
  res_text_AC@gsnLeftString                  = ""
  res_text_AC@gsnRightString                 = "" 
  res_text_AC@vpHeightF  = 0.12
  res_text_AC@vpWidthF   = 0.33       
  res_text_AC@gsnDraw                        = False;
  res_text_AC@gsnFrame                       = False;
  res_text_AC@cnFillOn                       = False
  res_text_AC@gsnAddCyclic                   = True 
  res_text_AC@cnInfoLabelOn                  = False 
  res_text_AC@txAngleF                       = 0
  res_text_AC@txFontHeightF                  = 0.02
  res_text_AC@txFontThicknessF               = 3.0
  res_text_AC@txFontColor                    = "brown1"
  res_text_AC@txBackgroundFillColor          = -1
  res_text_AC@txFont                         = "helvetica"

  text200Z1_1 = gsn_add_text(wks,plot(0),"~F22~C",65,40,res_text_AC)

  text200Z2_1 = gsn_add_text(wks,plot(1),"~F22~C",19,56,res_text_AC)
  text200Z2_2 = gsn_add_text(wks,plot(1),"~F22~C",104,55,res_text_AC)

  res_text_AC@txFontColor                    = "blue1"
  text200Z1_2 = gsn_add_text(wks,plot(0),"~F22~A",35,50,res_text_AC)
  text200Z1_3 = gsn_add_text(wks,plot(0),"~F22~A",97,42,res_text_AC)

  text200Z2_3 = gsn_add_text(wks,plot(1),"~F22~A",63,56,res_text_AC)
  text200Z2_4 = gsn_add_text(wks,plot(1),"~F22~A",126,41,res_text_AC)



 ;------------------------左场OLR--------------------------------------------
    res2                          = True
    res2@gsnDraw          = False
    res2@gsnFrame         = False      ;dont advance frame
    res2@gsnAddCyclic     = True      ;EOF data is not cyclic
    res2@cnFillOn         = True         ; turn on color fill
    res2@cnLinesOn        = False
    res2@cnLineLabelsOn   = False      ;turn off cn line labels
    res2@cnInfoLabelOn    = False      ;turn off contour information label
    

    res2@tmXTOn = False
    res2@tmYROn = False

    res2@mpOutlineBoundarySets       ="Geophysical"
    res2@mpGeophysicalLineColor      = "gray50"
    res2@mpGeophysicalLineDashPattern= 0
    res2@mpGeophysicalLineDashSegLenF= 0.2  
    res2@mpGeophysicalLineThicknessF = 1.0  ;地图边界的粗细
    res2@mpLandFillColor             = "white"

    res2@tmXBMinorOn               = True
    res2@tmXBMinorThicknessF       = 0.
    res2@tmXBMinorLengthF          = 0.
    res2@tmXBMinorOutwardLengthF   = 0.
    res2@tmXBMajorThicknessF       = 3.0
    res2@tmXBMajorLengthF          = 0.003
    res2@tmXBMajorOutwardLengthF   = 0.003
    res2@tmXBLabelDeltaF           = -0.6
    res2@tmXBLabelFontHeightF        = 0.0075
    res2@tmXBLabelFontThicknessF     = 0.03    
    res2@tmBorderThicknessF        = 3.0     ;边框粗细
    res2@tmXBMode          = "Explicit"
    res2@tmXBValues        = (/-240, -180,-120, -60, 0, 60/)
    res2@tmXBLabels        = (/"120~S~o~N~E","180~S~o~N~", "120~S~o~N~W","60~S~o~N~W", "0", "60~S~o~N~E"/)

    res2@tmYLMinorOn               = True
    res2@tmYLMinorThicknessF       = 0.
    res2@tmYLMinorLengthF          = 0.
    res2@tmYLMinorOutwardLengthF   = 0.
    res2@tmYLMajorThicknessF       = 3.0
    res2@tmYLMajorLengthF          = 0.003
    res2@tmYLMajorOutwardLengthF   = 0.003
    res2@tmYLLabelDeltaF           = -0.6
    res2@tmYLLabelFontHeightF        = 0.0075
    res2@tmYLLabelFontThicknessF     = 0.03    
    res2@tmBorderThicknessF        = 3.0     ;边框粗细
    res2@tmYLMode          = "Explicit"
    res2@tmYLValues        = (/-10, 0, 10, 20 ,30/)
    res2@tmYLLabels        = (/"10~S~o~N~S","0~S~o~N~","10~S~o~N~N","20~S~o~N~N" ,"30~S~o~N~N"/)

    res2@mpShapeMode = "FreeAspect"
    res2@vpHeightF  = 0.08
    res2@vpWidthF   = 0.33
    res2@txFont     = "helvetica"

    res2@gsnAddCyclic             = True
    res2@mpMinLatF                = latS_olr-1
    res2@mpMaxLatF                = latN_olr
    res2@mpMinLonF                   = -270
    res2@mpMaxLonF                   = 90
    res2@mpCenterLonF                = (-270+90)/2.
     
    res2@cnSmoothingOn          = True
    res2@cnFillOn              = True  

    res2@lbLabelBarOn                 = False
    res2@lbBoxLinesOn                 = False
    res2@cnLevelSelectionMode         = "ExplicitLevels"     
    res2@cnLevels                     = (/-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,0,0.2,0.3,0.4,0.5,0.6,0.7/)   
    res2@cnFillPalette         = "MPL_BrBG"
    res2@cnFillColors          = (/109,100,92,85,79,74,-1,-1,55,50,44,37,29,20/)
    ;res2@pmTickMarkDisplayMode  = "Always"

    res2@lbOrientation                = "Vertical"
    res2@lbBoxEndCapStyle             = "TriangleBothEnds"     
    res2@pmLabelBarWidthF             =  0.02       ;色标宽度
    res2@pmLabelBarHeightF            =  0.08             ;   色标高度
    res2@lbLabelStride                =  2
    res2@lbLabelOffsetF               =  0.1              ;色标与标签之间的距离
    res2@lbLabelFontHeightF           =  0.005           ;色标标签字体大小    
    res2@pmLabelBarOrthogonalPosF     = 0.01
    ;res2@pmLabelBarParallelPosF      = 0.5 

    res2@gsnLeftStringFontHeightF  = 0.01
    res2@gsnLeftStringOrthogonalPosF = 0.023
    res2@gsnLeftString = "~F22~MCA1_OLR"       
    plot(2) = gsn_csm_contour_map(wks,homlft(0,:,:),res2) 

    res2@lbLabelBarOn          =  True
    res2@gsnLeftString = "~F22~MCA2_OLR"       
    plot(3) = gsn_csm_contour_map(wks,homlft(1,:,:),res2) 


    ;-------------------------时间序列画图设置--------------------------------------------------
    res3          = True
    res3@gsnDraw  = False       ; dont draw yet
    res3@gsnFrame = False       ; dont advance frame yet

    res3@tmXTOn = False
    res3@tmYROn = False

    res3@vpHeightF  = 0.08
    res3@vpWidthF   = 0.33

    res3@tmXBMinorOn               = True
    res3@tmXBMinorThicknessF       = 0.
    res3@tmXBMinorLengthF          = 0.
    res3@tmXBMinorOutwardLengthF   = 0. 
    res3@tmXBMajorThicknessF       = 3.0
    res3@tmXBMajorLengthF          = 0.003
    res3@tmXBMajorOutwardLengthF   = 0.003
    res3@tmXBLabelDeltaF           = -0.6
    res3@tmXBLabelFontHeightF        = 0.0075
    res3@tmXBLabelFontThicknessF     = 0.03    
    res3@tmBorderThicknessF        = 3.0     ;边框粗细

    res3@tmYLMinorOn               = True
    res3@tmYLMinorThicknessF       = 0.
    res3@tmYLMinorLengthF          = 0.
    res3@tmYLMinorOutwardLengthF   = 0.
    res3@tmYLMajorThicknessF       = 3.0
    res3@tmYLMajorLengthF          = 0.003
    res3@tmYLMajorOutwardLengthF   = 0.003
    res3@tmYLLabelDeltaF           = -0.6
    res3@tmYLLabelFontHeightF        = 0.0075
    res3@tmYLLabelFontThicknessF     = 0.03    
    res3@tmBorderThicknessF        = 3.0     ;边框粗细
    res3@tmYLMode          = "Explicit"
    res3@tmYLValues        = (/-3.0,-1.5,0,1.5,3.0/)
    res3@tmYLLabels        = (/"-3.0","-1.5","0.0","1.5","3.0"/)

    res3@gsnYRefLine       = 0.              ; reference line
    res3@gsnYRefLineDashPattern = 5
    res3@gsnYRefLineThicknessF = 2.0
    res3@gsnYRefLineColor  = "gray50"
    res3@xyLineColors      = (/"dodgerblue2","brown1"/) ; colors chosen
    res3@xyLineThicknesses = (/3.5,3.5/)      ; line thicknesses
    res3@xyDashPatterns    = (/0.,0./)      ; make all lines solid

    res3@trXMinF          = yrStrt  ; leave a margin for legend
    res3@trXMaxF          = yrLast
    res3@txFont     = "helvetica"

    res3@tmXBMode         = "Manual"
    res3@tmXBTickStartF   = yrStrt
    res3@tmXBTickEndF     = yrLast
    res3@tmXBTickSpacingF = 5

    res3@trYMinF      = -3.4            ; min value on x-axis
    res3@trYMaxF      = 3.4            ; max value on x-axis

    ;res3@tiYAxisString             = "Standardized"         ; yaxis label
    res3@gsnLeftStringFontHeightF  = 0.01
    res3@gsnLeftStringOrthogonalPosF = 0.023

    res3@gsnLeftString = "~F22~TC1 (scf="+decimalPlaces(x(0), 0, True)+ "%)"
    plot(4) = gsn_csm_xy(wks,years,data1,res3)

    res3@gsnLeftString = "~F22~TC2 (scf="+decimalPlaces(x(1), 0, True)+ "%)"
    plot(5) = gsn_csm_xy(wks,years,data2,res3) 

    ;--------------------------折线图叠加图例---------------------------------
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsLineThicknessF  = 3.5                   ; line thicker

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.008                ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification
    res_text@txFont     = "helvetica"

    res_lines@gsLineColor = "dodgerblue2"
    yy                    = (/2.7,2.7/)

    xx1                    = (/2002,2004/)
    text_pos1              = 2005
    xx2                    = (/1990,1992/)
    text_pos2              = 1993


    dum11 = gsn_add_polyline(wks,plot(4),xx1,yy,res_lines)              ; add polyline
    dum12 = gsn_add_text(wks,plot(4),"~F21~TC1_OLR",text_pos1,yy(0),res_text); add text

    dum13 = gsn_add_polyline(wks,plot(5),xx1,yy,res_lines)              ; add polyline
    dum14 = gsn_add_text(wks,plot(5),"~F21~TC2_OLR",text_pos1,yy(0),res_text); add text

    res_lines@gsLineColor = "brown1"
    dum21 = gsn_add_polyline(wks,plot(4),xx2,yy,res_lines)              ; add polyline
    dum22 = gsn_add_text(wks,plot(4),"~F21~TC1_200Z",text_pos2,yy(0),res_text); add text

    dum22 = gsn_add_polyline(wks,plot(5),xx2,yy,res_lines)              ; add polyline
    dum24 = gsn_add_text(wks,plot(5),"~F21~TC2_200Z",text_pos2,yy(0),res_text); add text

    dum31 = gsn_add_text(wks,plot(4),"~F21~r="+decimalPlaces(ccr1,2,True),2013,-2.8,res_text); add text   +decimalPlaces(ccr1,2,True)
    dum31 = gsn_add_text(wks,plot(5),"~F21~r=0.80" ,2013,-2.8,res_text); add text +decimalPlaces(ccr2,2,True)     
    ;dum31 = gsn_add_text(wks,plot(5),"~F21~r="+decimalPlaces(ccr2,2,True),2013,-2.8,res_text); add text +decimalPlaces(ccr2,2,True)     


    res_lines@gsLineColor       = "firebrick2"
    res_lines@gsLineDashPattern = 11
    dum_box1= gsn_add_polyline(wks,plot(2), (/220,360/) ,(/25,25/),res_lines)  
    dum_box2= gsn_add_polyline(wks,plot(2), (/0,40/) ,(/25,25/),res_lines)  
    dum_box3= gsn_add_polyline(wks,plot(2), (/220,360/) ,(/-10,-10/),res_lines)
    dum_box4= gsn_add_polyline(wks,plot(2), (/0,40/) ,(/-10,-10/),res_lines) 
    dum_box5= gsn_add_polyline(wks,plot(2), (/220,220/) ,(/-10,25/),res_lines)  
    dum_box6= gsn_add_polyline(wks,plot(2), (/40,40/) ,(/-10,25/),res_lines)  

    dum_box11= gsn_add_polyline(wks,plot(3), (/220,360/) ,(/25,25/),res_lines)  
    dum_box22= gsn_add_polyline(wks,plot(3), (/0,40/) ,(/25,25/),res_lines)  
    dum_box33= gsn_add_polyline(wks,plot(3), (/220,360/) ,(/-10,-10/),res_lines)
    dum_box44= gsn_add_polyline(wks,plot(3), (/0,40/) ,(/-10,-10/),res_lines) 
    dum_box55= gsn_add_polyline(wks,plot(3), (/220,220/) ,(/-10,25/),res_lines)  
    dum_box66= gsn_add_polyline(wks,plot(3), (/40,40/) ,(/-10,25/),res_lines)  




    resP = True
    resP@gsnPanelRowSpec    = True
    resP@gsnPanelLabelBar   = False
    resP@gsnMaximize        = True
    resP@gsnPanelXWhiteSpacePercent = 10.0
    resP@gsnPanelYWhiteSpacePercent = 2.0
    resP@vpXF             = 0.01
    resP@vpYF            = 0.99
    ;drawNDCGrid(wks)
    resP@gsnPanelYF               = (/0.9,0.9,0.7,0.7,0.55,0.55/)
    resP@gsnPanelXF               = (/0.05,0.49,0.05,0.49,0.05,0.49/)
    resP@gsnPanelFigureStrings    = (/"~F22~a","~F22~d","~F22~b","~F22~e","~F22~c","~F22~f"/)
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF =  0.02  
    resP@gsnPanelFigureStringsJust    = "TopLeft"
    resP@amOrthogonalPosF         = -0.8        ;y轴方向往上方向调整
    resP@amParallelPosF           = -0.58        ;x轴方向往负方向调整

  gsn_panel(wks,plot,(/2,2,2/),resP)     ; now draw as one plot


   
