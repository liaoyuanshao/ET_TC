
    begin

     latS_olr   =  -10
     latN_olr   =  25
     lonL_olr   =  0
     lonR_olr   =  360

    yrStrt = 1979
    yrLast = 2021
    years  = ispan(yrStrt,yrLast,1)
    years!0="time"
    years@units="year"
    nyears = dimsizes(years)
    ;print(nyears)

    olr_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/olr.mon.mean.nc", "r")
    time_olr = olr_in->time
    XXXX_OLR     = cd_calendar(time_olr,-1 )
    t_start1 = ind(XXXX_OLR.eq.197901)
    t_end1   = ind(XXXX_OLR.eq.202112)
    olr      = olr_in->olr(t_start1:t_end1,{latS_olr:latN_olr},{lonL_olr:lonR_olr});截取时间
    olr      = lonFlip(olr)
    printVarSummary(olr)
    olr_JJA      = month_to_season(olr, "JJA")
    printVarSummary(olr_JJA)
    olr_JJA_dtrd   = dtrend_n(olr_JJA,True, 0)
    copy_VarCoords(olr_JJA,olr_JJA_dtrd)
    printVarSummary(olr_JJA_dtrd)
     

    lonR_olr1 = 40
    lonL_olr1 = -140
    olr_JJA_1  = olr_JJA(:,:,{lonL_olr1:lonR_olr1})
    printVarSummary(olr_JJA_1)

    dims = dimsizes(olr_JJA_1)
    print(dims)

    ;>--------------------   计算权重 -------------------------<
    w  = sqrt(cos(0.01745329 * olr_JJA_1&lat))
    wp = olr_JJA_1 * conform(olr_JJA_1, w, 1)
    copy_VarCoords(olr_JJA_1, wp)
    ;print(w)
    printVarSummary(wp)

    ;>--------------------- EOF分解----------------------<
    neof   = 10
    eof    = eofunc_n_Wrap(wp, neof, False, 0) 
    eof_ts = eofunc_ts_n_Wrap (wp, eof, False, 0)

    tzh = eof@pcvar

    lam = tzh * sqrt(2. / dims(0))                ;计算解释方差误差范围

    eof_ts_stddev = dim_stddev(eof_ts)            ;计算时间系数标准差
    eof_ts = dim_standardize_n_Wrap(eof_ts,0,1)   ;时间序列标准化

    print(tzh)
    print(lam)
    printVarSummary(eof)
    printVarSummary(eof_ts)
    printVarSummary(eof_ts_stddev)
    ntim = dims(0)
    sig_pcv = eofunc_north(eof@pcvar,ntim,True)
    print(sig_pcv)

    ;----------------------------------------
    soizs_1 = eof_ts(0,:)
    soizs_2 = eof_ts(1,:)
    print(soizs_1)
    print(soizs_2)
 

    ;---------------------原始序列数据: data_pc -----
    data_pc = new((/2,nyears/),"float")
    data_pc(0,:) = soizs_1  
    data_pc(1,:) = soizs_2
    printVarSummary(data_pc) 

    ;--------------pc1&pc2 滑动平均 滤除年代际信号-----------
    data_run_ave_9yr = runave_n_Wrap(data_pc(:,:),9,0,1)
    printVarSummary(data_run_ave_9yr) 


    ;>---------------------------------------------------------------------------------------<
    ;; Normalize time series: Sum spatial weights over the area of used
    ;>---------------------------------------------------------------------------------------<
    dimx   = dimsizes( wp )
    mln    = dimx(2)
    sumWgt = mln*sum( w )
    eof_ts = eof_ts/sumWgt
    print("eof_ts: min="+min(eof_ts)+"   max="+max(eof_ts) )

    do i = 0,neof-1
        eof(i,:,:) = (/eof(i,:,:)/)* (/eof_ts_stddev(i) /)  
    end do
    printVarSummary(eof)





     u_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_u_1959_2022_JJA.nc", "r")
     ;print(u_in)
     time = u_in->time
     XXXX     = cd_calendar(time,-1 )
     ;print(XXXX)
     t_start = ind(XXXX.eq.197906)
     t_end   = ind(XXXX.eq.202108)     
     u       = u_in->u
     printVarSummary(u)  
     u      := u(t_start:t_end,:,{latS_olr:latN_olr},:) 
     u_flt  = short2flt(u)
     printVarSummary(u_flt)
     u_JJA = new((/nyears,3,141,1440/),"float")
     u_JJA!0="time"
     u_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     u_JJA(nn,:,:,:)= (u_flt(nn*3,:,:,:)+u_flt((nn*3)+1,:,:,:)+u_flt((nn*3)+2,:,:,:))/3.0
     end do
     copy_VarMeta(u_flt(0,:,:,:),u_JJA(0,:,:,:))
     printVarSummary(u_JJA)

      u850_JJA    = u_JJA(:,{850},:,:)
      copy_VarMeta(u_JJA(:,0,:,:),u850_JJA(:,:,:))
      printVarSummary(u850_JJA)

    ;---------------------------------------------------------
     v_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_v_1959_2022_JJA.nc", "r")
     ;print(v_in)
     v       = v_in->v
     printVarSummary(v)  
     v      := v(t_start:t_end,:,{latS_olr:latN_olr},:) 
     v_flt  = short2flt(v)
     printVarSummary(v_flt)
     v_JJA = new((/nyears,3,141,1440/),"float")
     v_JJA!0="time"
     v_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     v_JJA(nn,:,:,:)= (v_flt(nn*3,:,:,:)+v_flt((nn*3)+1,:,:,:)+v_flt((nn*3)+2,:,:,:))/3.0
     end do
     copy_VarMeta(v_flt(0,:,:,:),v_JJA(0,:,:,:))
     printVarSummary(v_JJA)

      v850_JJA    = v_JJA(:,{850},:,:)
      copy_VarMeta(v_JJA(:,0,:,:),v850_JJA(:,:,:))
      printVarSummary(v850_JJA)

    ;------------------------------ 计算850UUUUVVVVVV回归系数并做检验=---------------------------
    u_reg_pc = regCoef_n(data_pc,u850_JJA(latitude|:,longitude|:,time|:),1,2)   
    u_reg_pc!0 = "sv"
    u_reg_pc!1 = "lat"
    u_reg_pc!2 = "lon"
    u_reg_pc&lat = u850_JJA&latitude
    u_reg_pc&lon = u850_JJA&longitude
    ;printVarSummary(u_reg_pc)

    u_prob_pc_1 = student_t(u_reg_pc@tval,u_reg_pc@nptxy-2)
    u_prob_pc   = onedtond(u_prob_pc_1,dimsizes(u_reg_pc))
    copy_VarCoords(u_reg_pc, u_prob_pc)
    ;printVarSummary(u_prob_pc) 

    v_reg_pc = regCoef_n(data_pc,v850_JJA(latitude|:,longitude|:,time|:),1,2)     
    v_reg_pc!0 = "sv"
    v_reg_pc!1 = "lat"
    v_reg_pc!2 = "lon"
    v_reg_pc&lat = u850_JJA&latitude
    v_reg_pc&lon = u850_JJA&longitude
    ;printVarSummary(v_reg_pc)

    v_prob_pc_1 = student_t(v_reg_pc@tval,v_reg_pc@nptxy-2)
    v_prob_pc   = onedtond(v_prob_pc_1,dimsizes(v_reg_pc))
    copy_VarCoords(v_reg_pc, v_prob_pc)
    ;printVarSummary(v_prob_pc) 

    u_reg_pc1 = u_reg_pc(0,:,:)
    v_reg_pc1 = v_reg_pc(0,:,:)   
    u_reg_pc2 = u_reg_pc(1,:,:)
    v_reg_pc2 = v_reg_pc(1,:,:)

    u_reg_pc1=where(sqrt((/u_reg_pc1/)^2+(/v_reg_pc1/)^2).lt. 0.15, u_reg_pc1@_FillValue, u_reg_pc1)   ;mask风场
    v_reg_pc1=where(sqrt((/u_reg_pc1/)^2+ (/v_reg_pc1/)^2).lt.0.15, v_reg_pc1@_FillValue, v_reg_pc1)

    u_reg_pc2=where(sqrt((/u_reg_pc2/)^2+(/v_reg_pc2/)^2).lt. 0.1, u_reg_pc2@_FillValue, u_reg_pc2)
    v_reg_pc2=where(sqrt((/u_reg_pc2/)^2+ (/v_reg_pc2/)^2).lt.0.1, v_reg_pc2@_FillValue, v_reg_pc2)


    ;-----------------SVD序列OLR1 OLR2--------------------------------------
    SVD_file_in=addfile("/nuist/scratch/zhuzhiwei/shaoly/code3/SVD_data/HGT_OLR_79_21_JJA_SVD/undtrend/SVD_pc.nc","r")
    SVD_in = SVD_file_in->ak_std
    OLR1     = SVD_in(0,:)
    OLR2     = SVD_in(1,:)

    PC1_OLR1_COR = esccr(soizs_1,OLR1,0)
    print(PC1_OLR1_COR)
    PC2_OLR2_COR = esccr(soizs_2,OLR2,0)
    print(PC2_OLR2_COR)  

    ;----------------------散点图计算回归线--------------------
    rc_pc1_OLR1    = regline(soizs_1,OLR1)                     
    ;print(rc_pc1_OLR1)
    pltarry_pc1_OLR1   = new ( (/2,nyears/), "float")
    pltarry_pc1_OLR1(0,:) = OLR1                              ; use markers
    pltarry_pc1_OLR1(1,:) = rc_pc1_OLR1*soizs_1 + rc_pc1_OLR1@yintercept              ; use solid line
    ;printVarSummary(pltarry_pc1_OLR1)

    rc_pc2_OLR2    = regline(soizs_2,OLR2)                     
    ;print(rc_pc2_OLR2)
    pltarry_pc2_OLR2   = new ( (/2,nyears/), "float")
    pltarry_pc2_OLR2(0,:) = OLR2                              ; use markers
    pltarry_pc2_OLR2(1,:) = rc_pc2_OLR2*soizs_2 + rc_pc2_OLR2@yintercept              ; use solid line
    ;printVarSummary(pltarry_pc2_OLR2)







;>------------------------------------------wks1 EOF出图 ------------------------------------------------------------------------------------------------------<
;>-------------------------------------------wks1 EOF出图 -----------------------------------------------------------------------------------------------------<
    wks1 = gsn_open_wks("eps","olr-eof-value")
    plot_eof = new(4, "graphic")
    gsn_define_colormap(wks1,"cb_9step")

    res                      = True         
    res@gsnDraw              = False        ; don't draw yet
    res@gsnFrame             = False        ; don't advance frame yet
    res@gsnAddCyclic         = False        ; plotted dataa are not cyclic
    res@mpFillOn             = False        ; turn off map fill
    res@mpOutlineOn          = True         ; turn the map outline on
    ; res@pmTickMarkDisplayMode= "Always"

    res@mpMaxLonF            = lonR_olr1
    res@mpMinLonF            = lonL_olr1
    res@mpMaxLatF            = latN_olr
    res@mpMinLatF            = latS_olr
    res@mpCenterLonF         = (lonR_olr1+lonL_olr1)/2.
    res@mpShapeMode = "FreeAspect"
    res@vpHeightF  = 0.13
    res@vpWidthF   = 0.3

    res@tmXBMinorOn = True
    res@tmYLMinorOn = True
    res@mpOutlineBoundarySets ="Geophysical"
    res@mpGeophysicalLineColor      = "gray50"
    res@mpGeophysicalLineDashPattern= 0
    res@mpGeophysicalLineDashSegLenF= 0.2  
    res@mpGeophysicalLineThicknessF = 1.2  ;地图边界的粗细

    ;>------设置经纬度数值大小------
    res@tmXBMinorOn               = True
    res@tmXBMinorThicknessF       = 2.5
    res@tmXBMinorLengthF          = 0.000
    res@tmXBMinorOutwardLengthF   = 0.000 
    res@tmXBMajorThicknessF       = 4.0
    res@tmXBMajorLengthF          = 0.003
    res@tmXBMajorOutwardLengthF   = 0.003
    res@tmXBLabelDeltaF           = -0.6
    res@tmXBLabelFontHeightF        = 0.008
    res@tmXBLabelFontThicknessF     = 0.03    
    res@tmBorderThicknessF        = 4.0     ;边框粗细
    res@tmXBMode          = "Explicit"
    res@tmXBValues        = (/-120,-90,-60,-30, 0, 30/)
    res@tmXBLabels        = (/"120~S~o~N~W","90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0", "30~S~o~N~E" /)

    res@tmYLMinorOn               = True
    res@tmYLMinorThicknessF       = 2.5
    res@tmYLMinorLengthF          = 0.000
    res@tmYLMinorOutwardLengthF   = 0.000
    res@tmYLMajorThicknessF       = 4.0
    res@tmYLMajorLengthF          = 0.003
    res@tmYLMajorOutwardLengthF   = 0.003
    res@tmYLLabelDeltaF           = -0.6
    res@tmYLLabelFontHeightF        = 0.008
    res@tmYLLabelFontThicknessF     = 0.03    
    res@tmBorderThicknessF        = 4.0     ;边框粗细
    res@tmYLMode          = "Explicit"
    res@tmYLValues        = (/-10,0,10,20/)
    res@tmYLLabels        = (/"10~S~o~N~S","0","10~S~o~N~N","20~S~o~N~N"/)

    res@tmYLValues           = (/-10,0,10,20/)
    res@tmXBValues           = (/-120,-90,-60,-30,0,30/)
   
    res@cnFillOn             = True         ; turn on color fill
    res@cnLinesOn            = False       ; True is default
    res@cnLineLabelsOn       = False        ; True is default

    res@cnSmoothingOn        = True
    res@tmXTOn               = False
    res@tmYROn               = False

    res@gsnLeftStringFontHeightF      = 0.01
    res@gsnLeftStringOrthogonalPosF   = 0.01

    res@cnLevelSelectionMode        = "ExplicitLevels"
    res@cnLevels                    = (/-8,-7,-5,-3,-2,-1, 1,2,3,5,7,8/)   
    res@cnFillPalette               = "MPL_BrBG"    ; 绿棕配色
    res@cnFillColors                = (/114,110,104,92,80,72,-1,60,50,40,30,25,19/)
    res@lbLabelBarOn                =  False
    res@gsnLeftString               = "~F22~EOF1                                                               "+ decimalPlaces(eof@pcvar(0), 0, True)+"%"
    ;res@gsnLeftString               = "~F22~(a) EOF1                                                  "+ decimalPlaces(eof@pcvar(0), 0, True)+"%"
    plot_eof(0) = gsn_csm_contour_map(wks1, eof(0,:,:), res)


    res@lbLabelBarOn                =  True
    res@lbBoxLinesOn                = False
    res@lbBoxEndCapStyle            = "TriangleBothEnds"  
    res@lbOrientation               = "Horizontal"  ;Vertical
    res@lbLabelStride               = 1
    res@lbLabelOffsetF              = 0.1
    res@lbLabelFontHeightF          = 0.006        ;字体大小
    res@pmLabelBarWidthF            = 0.24
    res@pmLabelBarHeightF           = 0.025       ;调整高度
    res@pmLabelBarOrthogonalPosF    = 0.1  ;调整色标垂直方向位置(越大越上)       ;越小越左
    res@pmLabelBarParallelPosF      = 0.5
    res@gsnLeftString               = "~F22~EOF2                                                               "+ decimalPlaces(eof@pcvar(1), 0, True)+"%"
    ;res@gsnLeftString               = "~F22~(d) EOF2                                                  "+ decimalPlaces(eof@pcvar(1), 0, True)+"%"
    plot_eof(3) = gsn_csm_contour_map(wks1, eof(1,:,:), res)

    ;--------------------850hpa回归风场叠加--------------
    res_uv                      = True
    res_uv@vcVectorDrawOrder    = "PostDraw"
    res_uv@gsnDraw              = False;
    res_uv@gsnFrame             = False;
    res_uv@gsnAddCyclic         = True
    res_uv@vcGlyphStyle              = "CurlyVector"  ;其余三种选项为“LineArrow”、“WindBarb” 、“CurlyVector”
    res_uv@vcLineArrowThicknessF     = 3.0
    res_uv@vpHeightF  = 0.13
    res_uv@vpWidthF   = 0.3
    res_uv@vcLineArrowHeadMaxSizeF    = 0.0075
    res_uv@vcLineArrowHeadMinSizeF    = 0.003
    res_uv@vcLineArrowColor            = "gray20"

    ;****设定参考箭头****
    res_uv@vcRefAnnoOn               = True  
    res_uv@vcRefLengthF              = 0.04   ;标准长度箭头在单位平方中的大小
    res_uv@vcRefAnnoBackgroundColor  = 0 ;背景颜色     
    res_uv@vcRefAnnoPerimOn          = True     ;关闭边框                                       
    res_uv@vcRefAnnoFontHeightF      = 0.005 ;参考箭头标签字体大小      

    res_uv@vcRefAnnoString1On     = False   ;设定参考箭头上、下的字符        
    res_uv@vcRefAnnoString2On     = True    ; 这里仅设定其下方的字符
    res_uv@vcRefAnnoSide            = "Bottom" ; 参考箭头放至图形上方
    res_uv@vcRefMagnitudeF        = 0.5    ;标准长度箭头所表示的大小
    res_uv@vcRefAnnoString2       = 0.5
    res_uv@vcMinDistanceF            = 0.01

    res_uv@vcRefAnnoOrthogonalPosF  = -0.145    
    plot_850uv_pc1=gsn_csm_vector(wks1,u_reg_pc1,v_reg_pc1,res_uv)
    overlay(plot_eof(0),plot_850uv_pc1)

    res_uv@vcRefAnnoOrthogonalPosF  = -0.435
    plot_850uv_pc2=gsn_csm_vector(wks1,u_reg_pc2,v_reg_pc2,res_uv)
    overlay(plot_eof(3),plot_850uv_pc2)


    ;----------------------------------------------------------散点图--------------------------------------------------------------------------------------- 
    rts          = True
    rts@gsnDraw  = False       ; dont draw yet
    rts@gsnFrame = False       ; dont advance frame yet

    rts@tmXTOn = False
    rts@tmYROn = False
    rts@vpHeightF  = 0.13
    rts@vpWidthF   = 0.13

    rts@tmXBMinorOn               = True
    rts@tmXBMinorThicknessF       = 2.5
    rts@tmXBMinorLengthF          = 0.000
    rts@tmXBMinorOutwardLengthF   = 0.000 
    rts@tmXBMajorThicknessF       = 4.0
    rts@tmXBMajorLengthF          = 0.003
    rts@tmXBMajorOutwardLengthF   = 0.003
    rts@tmXBLabelDeltaF           = -0.6
    rts@tmXBLabelFontHeightF        = 0.008
    rts@tmXBLabelFontThicknessF     = 0.03    
    rts@tmBorderThicknessF        = 4.0     ;边框粗细

    rts@tmYLMinorOn               = True
    rts@tmYLMinorThicknessF       = 2.5
    rts@tmYLMinorLengthF          = 0.000
    rts@tmYLMinorOutwardLengthF   = 0.000
    rts@tmYLMajorThicknessF       = 4.0
    rts@tmYLMajorLengthF          = 0.003
    rts@tmYLMajorOutwardLengthF   = 0.003
    rts@tmYLLabelDeltaF           = -0.6
    rts@tmYLLabelFontHeightF        = 0.008
    rts@tmYLLabelFontThicknessF     = 0.03    
    rts@tmBorderThicknessF        = 4.0     ;边框粗细

    rts@trXMinF     = -3
    rts@trXMaxF     = 3.1
    rts@trYMinF     = -3
    rts@trYMaxF     = 3

    rts@tmXBTickStartF        = -3
    rts@tmXBTickEndF          = 3
    rts@tmYLTickStartF        = -3
    rts@tmYLTickEndF          = 3

    rts@gsnXRefLine            = 0          ;X轴参考线
    rts@gsnXRefLineColor       = "gray60"
    rts@gsnXRefLineDashPattern = 11
    rts@gsnXRefLineThicknessF  = 4.0
    rts@gsnYRefLine            = 0          ;Y轴参考线
    rts@gsnYRefLineColor       = "gray60"
    rts@gsnYRefLineDashPattern = 11
    rts@gsnYRefLineThicknessF  = 4.0

    rts@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
    rts@xyMarkers           = 16                     ; choose type of marker 
    rts@xyMarkerColor       = "red1"                  ; Marker color
    rts@xyMarkerSizeF       = 0.0035                 ; Marker size (default 0.01)
    rts@xyDashPatterns      = 1                      ; solid line 
    rts@xyLineThicknesses   = (/1,5/)                ; set second line to 2
    rts@xyLineColors        = (/"red1","red3"/)

    rts@tiXAxisOn             = True         ;X轴标题启用
    rts@tiXAxisFontHeightF    = 0.008        ;X轴标题字体大小
    rts@tiXAxisOffsetXF       = 0.00         ;X轴标题X方向移动
    rts@tiXAxisOffsetYF       = 0.00         ;X轴标题Y方向移动
    rts@tiYAxisOn             = True          ;y轴标题启用
    rts@tiYAxisFontHeightF    = 0.008         ;y轴标题字体大小
    rts@tiYAxisOffsetXF       = 0.005         ;y轴标题X方向移动
    rts@tiYAxisOffsetYF       = 0.00          ;y轴标题y方向移动
    rts@xyCurveDrawOrder      = "PostDraw"

    rts@gsnLeftStringFontHeightF      = 0.01
    rts@gsnLeftStringOrthogonalPosF   = -0.02

    rts@tiXAxisString       = "TC1_OLR"
    rts@tiYAxisString       = "PC1"    
    rts@gsnLeftString       = "~F22~" ;"~F22~(b) "
    plot_eof(1)   = gsn_csm_xy(wks1,pltarry_pc1_OLR1,soizs_1, rts)   ;(wks,X变量,Y变量,rts)  

    rts@tiXAxisString       = "TC2_OLR"  
    rts@tiYAxisString       = "PC2"
    rts@gsnLeftString       = "~F22~" ;"~F22~(e) "
    plot_eof(4)   = gsn_csm_xy(wks1,pltarry_pc2_OLR2,soizs_2, rts)


    rts_text                    = True                  ; Cor图例文字
    rts_text@txFontHeightF      = 0.009                ; change text size
    rts_text@txJust             = "CenterLeft"          ; text justification
    rts_text@txFontColor        = "red1"  
    dum1 = gsn_add_text(wks1,plot_eof(1),"~F21~r="+decimalPlaces(rc_pc1_OLR1, 2, True),-2.8,2.6,rts_text) 
    dum2 = gsn_add_text(wks1,plot_eof(4),"~F21~r="+decimalPlaces(rc_pc2_OLR2, 2, True),-2.8,2.6,rts_text) 


   

    ;------------------------------------------wks1组图---------------------------------------------------------------------------------------------------
    resP = True
    resP@gsnPanelRowSpec    = True
    resP@gsnPanelLabelBar   = False
    resP@gsnMaximize        = True
    resP@gsnPanelXWhiteSpacePercent = 3.8
    resP@gsnPanelYWhiteSpacePercent = 7.0
    resP@vpXF             = 0.01
    resP@vpYF            = 0.99
    resP@gsnPanelYF               = (/0.9,0.90.725,0.725/)
    resP@gsnPanelXF               = (/0.05,0.38,0.525,0.05/)
    resP@gsnPanelFigureStrings    = (/"~F22~a","~F22~b","~F22~c","~F22~d"/)
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF =  0.018  
    resP@gsnPanelFigureStringsJust    = "TopLeft"
    resP@amOrthogonalPosF         = -0.72        ;y轴方向往上方向调整
    resP@amParallelPosF           = -0.58        ;x轴方向往负方向调整
    gsn_panel(wks1, plot_eof, (/2,2/), resP)

    end

