    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
    ;-----------------------------定义TN-WAF函数-----------------------
    undef("TN_WAF")
    function TN_WAF(hgt_anom,uwnd_clim,vwnd_clim,lev)
    local zavar,cuvar,cvvar,lat,lon,nlat,nlon,ga,re,pi,f,ilat,coslat,coslattmp,\
    ftmp,cumag,psidev,dpsidevdlon,ddpsidevdlonlon,dpsidevdlat,ddpsidevdlonlat,\
    ddpsidevdlatlat,xuterm,xvterm,yvterm,Fx,Fy
    begin

    zavar     := hgt_anom       ;- (nlat,nlon), geopotential height anomaly
    cuvar     := uwnd_clim      ;- (nlat,nlon), climatology zonal wind
    cvvar     := vwnd_clim      ;- (nlat,nlon), climatology meriditional wind 

    lat       := zavar&lat
    lon       := zavar&lon
    nlat      := dimsizes(lat)
    nlon      := dimsizes(lon)

    ;  Gravitational acceleration
    ga        := 9.80665
    ;  Radius of the earth
    re        := 6378388
    ; pi
    pi        := atan(1.0)*4.

    zavar1  = doubletofloat(zavar)
    ; Coriolis parameter
    f            :=  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
    f!0           = "lat"
    f&lat         = lat
    f@_FillValue  = zavar1@_FillValue

    ; missing for 10S - 10N
    do ilat = 0, nlat-1
    if (abs(lat(ilat) ).lt. 20..or.lat(ilat).le.0. ) then
    f(ilat)= f@_FillValue
    end if
    end do

    ; cosine
    coslat           := cos(lat(:)*pi/180.)

    ; 1-D -> 2-D
    ;leveltmp := conform_dims(dimsizes(zavar),lev,1)
    coslattmp        := conform_dims(dimsizes(zavar),coslat,0)
    ftmp             := conform_dims(dimsizes(zavar),f,0)

    ; magnitude of climatological wind
    cumag            := sqrt(cuvar^2 + cvvar^2)
    cumag@_FillValue := cuvar@_FillValue
    cumag            := where(cumag .gt. 0, cumag, cumag@_FillValue)

    ; QG steam function for anomaly
    psidev           := zavar*ga/ftmp

    ;dpsidev/dlon
    dpsidevdlon      :=  center_finite_diff_n(psidev,lon*pi/180.,True,0,1)
    ;ddpsidev/dlonlon
    ddpsidevdlonlon  :=  center_finite_diff_n(dpsidevdlon,lon*pi/180.,True,0,1)
    ;dpsidev/dlat
    dpsidevdlat      := center_finite_diff_n(psidev, lat*pi/180., False,0,0)
    ;ddpsidev/dlonlat
    ddpsidevdlonlat  :=  center_finite_diff_n(dpsidevdlon,lat*pi/180.,False,0,0)
    ;ddpsidev/dlatdlat
    ddpsidevdlatlat  := center_finite_diff_n(dpsidevdlat, lat*pi/180.,False,0,0)

    xuterm           := (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
    xvterm           := (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)
    ;yuterm := xvterm
    yvterm           := (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)

    ; Mask out where westerlies is small or negative (less than 5 ms).
    ;  by using mask

    ;x-component of (38)
    Fx := mask(lev/1000./(2.*cumag*re*re)*( cuvar/coslattmp * xuterm + cvvar * xvterm), cuvar.lt.3,False)
    ;ycomponent 
    Fy := mask(lev/1000./(2.*cumag*re*re)*( cuvar*xvterm + coslattmp*cvvar*yvterm),  cuvar.lt.3,False)
    copy_VarMeta(hgt_anom, Fx)
    copy_VarMeta(hgt_anom, Fy)
    Fx@units = "m^2/s^2"
    Fx@units = "m^2/s^2"
    return([/Fx,Fy/])
    end



    begin
    latS      = -30
    latN      =  90
    lonL      = 0
    lonR      = 360

    yrStrt = 1979
    yrLast = 2021
    years  = ispan(yrStrt,yrLast,1)
    nyears = dimsizes(years)
    ;print(nyears)
    ;**************************************读取数据 ****************************************************
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;读取时间序列;;;;;;;;;;;;;;;;;;;;;;;
     soizs_1      =asciiread("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/pc1.txt", -1, "float")   
     soizs_2      =asciiread("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/pc2.txt", -1, "float")     
     soizs        = new((/2,nyears/),"float")
     soizs(0,:)   = soizs_1
     soizs(1,:)   = soizs_2
     soizs!0      = "sv"
     soizs!1      = "time"
     printVarSummary(soizs)  

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;读高度场数据;;;;;;;;;;;;;;;;
     hgt_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_hgt_1959_2022_JJA_14473.nc", "r")
     print(hgt_in)
     time = hgt_in->time
     XXXX     = cd_calendar(time,-1 )
     ;print(XXXX)
     t_start = ind(XXXX.eq.197906)
     t_end   = ind(XXXX.eq.202108)     
     hgt       = hgt_in->z
     printVarSummary(hgt)  
     hgt      := hgt(t_start:t_end,:,{latS:latN},:) 
     hgt_flt  = short2flt(hgt)
     printVarSummary(hgt_flt)
     hgt_JJA = new((/nyears,3,49,144/),"float")
     hgt_JJA!0="time"
     hgt_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     hgt_JJA(nn,:,:,:)= (hgt_flt(nn*3,:,:,:)+hgt_flt((nn*3)+1,:,:,:)+hgt_flt((nn*3)+2,:,:,:))/(3.0*9.8)
     end do
     copy_VarMeta(hgt_flt(0,:,:,:),hgt_JJA(0,:,:,:))
     printVarSummary(hgt_JJA)
    ;z200_JJA_clm = dim_avg_n_Wrap(hgt_flt(:,{200},{35:85},:),0)  ;{35:85}
    ;printVarSummary(z200_JJA_clm)  

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;读uuuuu场数据;;;;;;;;;;;;;;;;
     u_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_u_1959_2022_JJA_14473.nc", "r")
     ;print(u_in)
     time = u_in->time
     XXXX     = cd_calendar(time,-1 )
     ;print(XXXX)
     u       = u_in->u
     printVarSummary(u)  
     u      := u(t_start:t_end,:,{latS:latN},:) 
     u_flt  = short2flt(u)
     printVarSummary(u_flt)
     u_JJA = new((/nyears,3,49,144/),"float")
     u_JJA!0="time"
     u_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     u_JJA(nn,:,:,:)= (u_flt(nn*3,:,:,:)+u_flt((nn*3)+1,:,:,:)+u_flt((nn*3)+2,:,:,:))/3.0
     end do
     copy_VarMeta(u_flt(0,:,:,:),u_JJA(0,:,:,:))
     printVarSummary(u_JJA)

    u200_JJA    = u_JJA(:,{200},:,:)
    copy_VarMeta(u_JJA(:,0,:,:),u200_JJA(:,:,:))
    printVarSummary(u200_JJA)
    u200_JJA_clm = dim_avg_n_Wrap(u_flt(:,{200},{30:85},:),0) 
    printVarSummary(u200_JJA_clm)   


    ;---------------------------------------------------------
     v_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_v_1959_2022_JJA_14473.nc", "r")
     ;print(v_in)
     time = v_in->time
     XXXX     = cd_calendar(time,-1 )
     ;print(XXXX)
     v       = v_in->v
     printVarSummary(v)  
     v      := v(t_start:t_end,:,{latS:latN},:) 
     v_flt  = short2flt(v)
     printVarSummary(v_flt)
     v_JJA = new((/nyears,3,49,144/),"float")
     v_JJA!0="time"
     v_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     v_JJA(nn,:,:,:)= (v_flt(nn*3,:,:,:)+v_flt((nn*3)+1,:,:,:)+v_flt((nn*3)+2,:,:,:))/3.0
     end do
     copy_VarMeta(v_flt(0,:,:,:),v_JJA(0,:,:,:))
     printVarSummary(v_JJA)

      v200_JJA    = v_JJA(:,{200},:,:)
      copy_VarMeta(v_JJA(:,0,:,:),v200_JJA(:,:,:))
      printVarSummary(v200_JJA)
      v200_JJA_clm = dim_avg_n_Wrap(v_flt(:,{200},{30:85},:),0)
      printVarSummary(v200_JJA_clm)   


    ;-----------------------precip------------------------------------------------
    precip_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/precip.mon.mean.nc", "r")
    ;print(precip_in)                                ;GPCP数据197901-202212月
    time_precip = precip_in->time
    XXXX2     = cd_calendar(time_precip,-1 )
    ;print(XXXX2)
    precip      = precip_in->precip(:,{latS:latN},:)
    printVarSummary( precip)
    precip_JJA = new((/nyears,48,144/),"float")
    precip_JJA!0="time"
    precip_JJA&time = ispan(1979,2021,1)  
    do nn=0,42,1                     
    precip_JJA(nn,:,:)= (precip((nn*12)+5,:,:)+precip((nn*12)+6,:,:)+precip((nn*12)+7,:,:))/3.0
    end do
    copy_VarCoords(precip(0,:,:),precip_JJA(0,:,:))
    printVarSummary(precip_JJA)

    ;----------------------sst------------------------------------------------
    sst_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/HadISST_sst.nc", "r")
    ;print(sst_in)
    time_sst = sst_in->time
    XXXX3     = cd_calendar(time_sst,-1 )
    t_start1 = ind(XXXX3.eq.197801)
    t_end1   = ind(XXXX3.eq.202212)
    sst      = sst_in->sst(t_start1:t_end1,{latS:latN},:)
    printVarSummary( sst)
    time_sst1 = sst&time
    XXXX33     = cd_calendar(time_sst1,-1 )     ;197801-202212月
    ;print(XXXX33)                   
    sst_JJA = new((/nyears,120,360/),"float")
    sst_JJA!0="time"
    sst_JJA&time = ispan(1979,2021,1)  
    do nn=0,42,1                      
    sst_JJA(nn,:,:)= (sst((nn*12)+17,:,:)+sst((nn*12)+18,:,:)+sst((nn*12)+19,:,:))/3.0
    end do
    copy_VarCoords(sst(0,:,:),sst_JJA(0,:,:))
    printVarSummary(sst_JJA)

    ;---------------------tem------------------------------------------------
    tem_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/cru_ts4.06.1901.2021.tmp.dat.nc", "r")
    ;print(tem_in)
    time_tem = tem_in->time
    XXXX4     = cd_calendar(time_tem,-1 )
    t_start1 = ind(XXXX4.eq.197901)
    t_end1   = ind(XXXX4.eq.202112)
    tem      = tem_in->tmp(t_start1:t_end1,{latS:latN},:)
    ;printVarSummary(tem)
    tem_JJA    = month_to_season(tem, "JJA")
    printVarSummary(tem_JJA )



    rc_Z = regCoef_n(soizs(sv|:,time|:),hgt_JJA(level|:,lat|:,lon|:,time|:),1,3)    
    rc_Z!0 = "sv"
    rc_Z!1 = "level"
    rc_Z!2 = "lat"
    rc_Z!3 = "lon"
    rc_Z&level = hgt_JJA&level
    rc_Z&lat = hgt_JJA&lat
    rc_Z&lon = hgt_JJA&lon
    printVarSummary(rc_Z)

    prob_Z_1 = student_t(rc_Z@tval,rc_Z@nptxy-2)
    prob_Z   = onedtond(prob_Z_1,dimsizes(rc_Z))
    copy_VarCoords(rc_Z,prob_Z)
    printVarSummary(prob_Z) 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 计算200UUUU回归系数并做检验
    rc_U = regCoef_n(soizs(sv|:,time|:),u_JJA(level|:,lat|:,lon|:,time|:),1,3)    
    rc_U!0 = "sv"
    rc_U!1 = "level"
    rc_U!2 = "lat"
    rc_U!3 = "lon"
    rc_U&level = u_JJA&level
    rc_U&lat = u_JJA&lat
    rc_U&lon = u_JJA&lon
    printVarSummary(rc_U)
    rc_U := smth9_Wrap(rc_U, 0.50, 0.25, True)

    prob_U_1 = student_t(rc_U@tval,rc_U@nptxy-2)
    prob_U   = onedtond(prob_U_1,dimsizes(rc_U))
    copy_VarCoords(rc_U,prob_U)
    printVarSummary(prob_U) 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 计算200VVVV回归系数并做检验
    rc_V = regCoef_n(soizs(sv|:,time|:),v_JJA(level|:,lat|:,lon|:,time|:),1,3)    
    rc_V!0 = "sv"
    rc_V!1 = "level"
    rc_V!2 = "lat"
    rc_V!3 = "lon"
    rc_V&level = u_JJA&level
    rc_V&lat = u_JJA&lat
    rc_V&lon = u_JJA&lon
    printVarSummary(rc_V)
    rc_V := smth9_Wrap(rc_V, 0.50, 0.25, True)

    prob_V_1 = student_t(rc_V@tval,rc_V@nptxy-2)
    prob_V   = onedtond(prob_V_1,dimsizes(rc_V))
    copy_VarCoords(rc_V,prob_V)
    printVarSummary(prob_V) 


    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;计算PRECIP归系数并做检验
    rc_pre = regCoef_n(soizs(sv|:,time|:),precip_JJA(lat|:,lon|:,time|:),1,2)     ;返回值为tval(t值统计量tval=rc@tval)  nptxy(自由度df=rc@nptxy-2）
    rc_pre!0 = "sv"
    rc_pre!1 = "lat"
    rc_pre!2 = "lon"
    rc_pre&lat = precip_JJA&lat
    rc_pre&lon = precip_JJA&lon
    printVarSummary(rc_pre)

    prob_pre_1 = student_t(rc_pre@tval,rc_pre@nptxy-2)
    prob_pre = onedtond(prob_pre_1,dimsizes(rc_pre))
    copy_VarCoords(rc_pre, prob_pre)
    ;printVarSummary(prob_pre) 
    prob_pre = where(abs(rc_pre).ge.0.2,prob_pre,prob_pre@_Fillvalue)
    copy_VarCoords(rc_pre, prob_pre)
    printVarSummary(prob_pre)   


    ;;;;;;;;;;;;;;;;;;;;;; ; 计算sst-IO回归系数并做检验
    rc_sst = regCoef_n(soizs(sv|:,time|:),sst_JJA(latitude|:,longitude|:,time|:),1,2)     
    rc_sst!0 = "sv"
    rc_sst!1 = "lat"
    rc_sst!2 = "lon"
    rc_sst&lat = sst_JJA&latitude
    rc_sst&lon = sst_JJA&longitude
    printVarSummary(rc_sst)

    prob_sst_1 = student_t(rc_sst@tval,rc_sst@nptxy-2)
    prob_sst   = onedtond(prob_sst_1,dimsizes(rc_sst))
    copy_VarCoords(rc_sst, prob_sst)
    ;printVarSummary(prob_sst) 
    rc_sst_test = where(abs(rc_sst).lt.5,rc_sst,rc_sst@_Fillvalue)
    copy_VarCoords(rc_sst, rc_sst_test)
    printVarSummary(rc_sst_test)

    ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;计算tem-IO回归系数并做检验
    rc_tem = regCoef_n(soizs(sv|:,time|:),tem_JJA(lat|:,lon|:,time|:),1,2)     
    rc_tem!0 = "sv"
    rc_tem!1 = "lat"
    rc_tem!2 = "lon"
    rc_tem&lat = tem_JJA&lat
    rc_tem&lon =tem_JJA&lon
    printVarSummary(rc_tem)

    prob_tem_1 = student_t(rc_tem@tval,rc_tem@nptxy-2)
    prob_tem   = onedtond(prob_tem_1,dimsizes(rc_tem))
    copy_VarCoords(rc_tem, prob_tem)
    printVarSummary(prob_tem) 

    ;----------------------计算TN-WAF部分-----------------------------------
    lev = 200
    u_clim = u200_JJA_clm
    v_clim = v200_JJA_clm
    ;--------计算pc1 TN-WAF部分--------
    regh_PC1   = rc_Z(0,{200},{30:85},:)
    WAF_PC1 = TN_WAF(regh_PC1,u_clim,v_clim,lev)   ;-- hgt_anom is anomaly based on climatology
    Fx_PC1  = WAF_PC1[0]
    Fy_PC1  = WAF_PC1[1]
    copy_VarCoords(u200_JJA_clm,Fx_PC1)
    copy_VarCoords(u200_JJA_clm,Fy_PC1)  
    Fx_PC1 = where(sqrt(Fx_PC1^2+Fy_PC1^2).le.0.003, Fx_PC1@_FillValue, Fx_PC1) 
    Fy_PC1 = where(sqrt(Fx_PC1^2+Fy_PC1^2).le.0.003, Fy_PC1@_FillValue, Fy_PC1)      
    Fx_PC1 := smth9_Wrap(Fx_PC1, 0.50, 0.25, True)
    Fy_PC1 := smth9_Wrap(Fy_PC1, 0.50, 0.25, True)
    ;printVarSummary(Fx_PC1)
    ;printVarSummary(Fy_PC1) 


    ;--------计算pc2 TN-WAF部分-----
    regh_PC2   = rc_Z(1,{200},{30:85},:)
    WAF_PC2 = TN_WAF(regh_PC2,u_clim,v_clim,lev)   ;-- hgt_anom is anomaly based on climatology
    Fx_PC2  = WAF_PC2[0]
    Fy_PC2  = WAF_PC2[1]
    copy_VarCoords(u200_JJA_clm,Fx_PC2)
    copy_VarCoords(u200_JJA_clm,Fy_PC2)
    Fx_PC2 := smth9_Wrap(Fx_PC2, 0.50, 0.25, True)
    Fy_PC2 := smth9_Wrap(Fy_PC2, 0.50, 0.25, True)
    Fx_PC2 = where(sqrt(Fx_PC2^2+Fy_PC2^2).le.0.007, Fx_PC2@_FillValue, Fx_PC2) 
    Fy_PC2 = where(sqrt(Fx_PC2^2+Fy_PC2^2).le.0.007, Fy_PC2@_FillValue, Fy_PC2) 
    printVarSummary(Fx_PC2)
    printVarSummary(Fy_PC2)



    ;------------------------         200pha 风场          --------------------
    rc_U_200vc = rc_U(:,{200},:,:)
    printVarSummary(rc_U_200vc)
    rc_V_200vc = rc_V(:,{200},:,:)
    printVarSummary(rc_V_200vc)
    rc_U_200vc=where(sqrt((/rc_U_200vc/)^2+(/rc_V_200vc/)^2).lt.0.17, rc_U_200vc@_FillValue, rc_U_200vc)
    rc_V_200vc=where(sqrt((/rc_U_200vc/)^2+ (/rc_V_200vc/)^2).lt.0.17, rc_V_200vc@_FillValue, rc_V_200vc)

    ;------------------------        500pha 风场          --------------------
    rc_U_500vc = rc_U(:,{500},:,:)
    printVarSummary(rc_U_500vc)
    rc_V_500vc = rc_V(:,{500},:,:)
    printVarSummary(rc_V_500vc)
    rc_U_500vc=where(sqrt((/rc_U_500vc/)^2+(/rc_V_500vc/)^2).lt.0.1, rc_U_500vc@_FillValue, rc_U_500vc)
    rc_V_500vc=where(sqrt((/rc_U_500vc/)^2+ (/rc_V_500vc/)^2).lt.0.1, rc_V_500vc@_FillValue, rc_V_500vc)

    ;------------------------       850pha 风场          --------------------
    rc_U_850vc = rc_U(:,{850},:,:)
    printVarSummary(rc_U_850vc)
    rc_V_850vc = rc_V(:,{850},:,:)
    printVarSummary(rc_V_850vc)
    rc_U_850vc=where(sqrt((/rc_U_850vc/)^2+(/rc_V_850vc/)^2).lt.0.022, rc_U_850vc@_FillValue, rc_U_850vc)
    rc_V_850vc=where(sqrt((/rc_U_850vc/)^2+ (/rc_V_850vc/)^2).lt.0.022, rc_V_850vc@_FillValue, rc_V_850vc)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;画图参数设置;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    wks = gsn_open_wks("eps","Reg_PC_JJA")
    base = new(6,"graphic") 
    gsn_define_colormap(wks,"cb_rainbow")
    res                       = True
    res@gsnDraw               = False
    res@gsnFrame              = False   
    res@gsnAddCyclic          = True

    res@cnFillOn              = True             ; turn on color
    res@cnLinesOn             = False             ; turn off contour lines
    res@cnLineLabelsOn        = False            ; turn off contour linelabels
    res@cnLevelSelectionMode  = "ManualLevels"   ; set manual contour levels
    ;res@pmTickMarkDisplayMode = "Always"


    res@cnSmoothingOn        = True
    res@lbLabelFont          = 0.25

    res@tmXTOn = False
    res@tmYROn = False

    res@mpOutlineOn                 = True       ; turn off outline
    res@mpOutlineBoundarySets       ="Geophysical"
    res@mpGeophysicalLineColor      = 206
    res@mpGeophysicalLineDashPattern= 0
    res@mpGeophysicalLineDashSegLenF= 0.2  
    res@mpGeophysicalLineThicknessF = 1.5  ;地图边界的粗细
    res@mpFillOn                    = False    ; turn off map fill    
    ;res@mpOutlineDrawOrder          = "PreDraw"

    res@tmXBMinorOn               = True
    res@tmXBMinorThicknessF       = 2.5
    res@tmXBMinorLengthF          = 0.000
    res@tmXBMinorOutwardLengthF   = 0.000 
    res@tmXBMajorThicknessF       = 3.0
    res@tmXBMajorLengthF          = 0.003
    res@tmXBMajorOutwardLengthF   = 0.003
    res@tmXBLabelDeltaF           = -0.6
    res@tmXBLabelFontHeightF        = 0.008
    res@tmXBLabelFontThicknessF     = 0.03    
    res@tmBorderThicknessF        = 3.0     ;边框粗细
    res@tmXBMode          = "Explicit"
    res@tmXBValues        = (/-135,-90,-45, 0, 45, 90, 135 , 180/)
    res@tmXBLabels        = (/"135~S~o~N~W","90~S~o~N~W","45~S~o~N~W", "0", "45~S~o~N~E", "90~S~o~N~E", "135~S~o~N~E", "180~S~o~N~"/)

    res@tmYLMinorOn               = True
    res@tmYLMinorThicknessF       = 2.5
    res@tmYLMinorLengthF          = 0.000
    res@tmYLMinorOutwardLengthF   = 0.000
    res@tmYLMajorThicknessF       = 3.0
    res@tmYLMajorLengthF          = 0.003
    res@tmYLMajorOutwardLengthF   = 0.003
    res@tmYLLabelDeltaF           = -0.6
    res@tmYLLabelFontHeightF        = 0.008
    res@tmYLLabelFontThicknessF     = 0.03    
    res@tmBorderThicknessF        = 3.0     ;边框粗细
    res@tmYLMode          = "Explicit"
    res@tmYLValues        = (/0,25,50,75/)
    res@tmYLLabels        = (/"0","25~S~o~N~N","50~S~o~N~N","75~S~o~N~N"/)

    res@gsnMajorLatSpacing          = 25
    res@gsnMinorLatSpacing          = 25
    res@tmYLLabelStride             = 1

    res@gsnMajorLonSpacing         = 45
    res@gsnMinorLonSpacing         = 45
    res@tmXBLabelStride            = 1

    res@mpShapeMode = "FreeAspect"
    res@vpHeightF  = 0.14
    res@vpWidthF   = 0.33

    res@gsnRightString              =   " "
    res@gsnRightString              =   " "
    res@gsnLeftStringFontHeightF      = 0.01
    res@gsnLeftStringOrthogonalPosF   = 0.01 
    res@lbLabelBarOn                =  False
    res@cnFillOn                    =  False
    res@cnInfoLabelOn               = False

    res@mpMinLatF                   = -20.
    res@mpMaxLatF                   = 90.
    res@mpMinLonF                   = -150
    res@mpMaxLonF                   = 210.
    res@mpCenterLonF                = (-150+210)/2.

    res@gsnLeftString          =   "~F22~Reg 200UV&WAF onto PC1"
    base(0) = gsn_csm_contour_map(wks,rc_pre(0,:,:),res) 

    res@gsnLeftString          =   "~F22~Reg 200UV&WAF onto PC2"
    base(1) = gsn_csm_contour_map(wks,rc_pre(0,:,:),res) 

    res@gsnLeftString          =   "~F22~Reg 500ZUV&T2m onto PC1"
    base(2) = gsn_csm_contour_map(wks,rc_sst(0,:,:),res) 

    res@gsnLeftString          =   "~F22~Reg 500ZUV&T2m onto PC2"
    base(3) = gsn_csm_contour_map(wks,rc_pre(0,:,:),res) 

    res@gsnLeftString          =   "~F22~Reg 850ZUV&PRE onto PC1"
    base(4) = gsn_csm_contour_map(wks,rc_pre(0,:,:),res) 

    res@gsnLeftString          =   "~F22~Reg 850ZUV&PRE onto PC2"
    base(5) = gsn_csm_contour_map(wks,rc_sst(0,:,:),res) 

    ;----------显著性检验打点-----------------------------------------------------------
    res2                      = True
    res2@gsnDraw              = False;
    res2@gsnFrame             = False;
    res2@gsnAddCyclic         = True
    res2@cnInfoLabelOn        = False;
    res2@cnFillOn             = True
    res2@cnLinesOn            = False;
    res2@cnLineLabelsOn       = False
    res2@lbLabelBarOn         = False; 
    res2@cnMonoFillPattern    = False
    res2@gsnLeftString        = ""
    res2@gsnRightString       = ""
    res2@cnLevelSelectionMode ="ExplicitLevels"
    res2@cnFillPatterns       = (/10,-1/)                     ; ;填充形状类型
    res2@cnFillColors         = (/"white","Transparent"/)     ;填充颜色-1为透明
    res2@cnFillScaleF         =   0.45                           ;越小越密集 
    res2@cnFillDotSizeF       =  0.0009                            
    res2@cnLevels             = (/0.1/)  
    plot_pre_prob_PC1    = gsn_csm_contour(wks,prob_pre(0,:,:),res2)     ;t检验pre显著性
    plot_sst_prob_PC1    = gsn_csm_contour(wks,prob_sst(0,:,:),res2)     ;t检验sst显著性
    plot_tem_prob_PC1    = gsn_csm_contour(wks,prob_tem(0,:,:),res2)     ;t检验tem显著性

    plot_pre_prob_PC2    = gsn_csm_contour(wks,prob_pre(1,:,:),res2)     ;t检验pre显著性
    plot_sst_prob_PC2    = gsn_csm_contour(wks,prob_sst(1,:,:),res2)     ;t检验sst显著性
    plot_tem_prob_PC2    = gsn_csm_contour(wks,prob_tem(1,:,:),res2)     ;t检验tem显著性

    ;--------------------------------------风场----------------------------------------
    res3                      = True
    res3@vcVectorDrawOrder    = "PostDraw"
    res3@gsnDraw             	 = False;
    res3@gsnFrame            	 = False;
    res3@gsnAddCyclic       	  = True
    res3@gsnLeftString      	  = ""
    res3@gsnRightString     	   = ""
    res3@vcGlyphOpacityF         = 1.0
    res3@vcMinFracLengthF       = 1.0
    res3@vcGlyphStyle                = "CurlyVector" 
    res3@vcLineArrowThicknessF     = 1.5

    res3@vcMinFracLengthF          = 1.0
    res3@vcLineArrowHeadMaxSizeF   = 0.012
    res3@vcLineArrowHeadMinSizeF   = 0.008

    ;****设定参考箭头****
    res3@vcRefAnnoOn               = True  
    res3@vcRefAnnoBackgroundColor  = 0 ;背景颜色     
    res3@vcRefAnnoPerimOn          =  True     ;关闭边框                                       
    res3@vcRefAnnoString1On        = False   ;设定参考箭头上、下的字符        
    res3@vcRefAnnoString2On        = True    ; 这里仅设定其下方的字符
    res3@vcRefAnnoSide             = "Bottom" ; 参考箭头放至图形上方
    res3@vcRefAnnoFontHeightF      = 0.008   ;参考箭头标签字体大小    
    res3@vcRefAnnoParallelPosF     = 0.998  
    res3@vcRefAnnoOrthogonalPosF   = -0.118

    ;--------------------200hpa-----------------
    res3@vcLineArrowColor            = "gray20"
    res3@vcRefMagnitudeF        =  2.0     ;标准长度箭头所表示的大小
    res3@vcRefAnnoString2       =  "2.0"

    res3@vcRefLengthF              = 0.03   ;标准长度箭头在单位平方中的大小
    res3@vcMinDistanceF            = 0.026
    plot_200uv_PC1 =gsn_csm_vector(wks,rc_U_200vc(0,:,:),rc_V_200vc(0,:,:),res3)
    plot_200uv_PC2 =gsn_csm_vector(wks,rc_U_200vc(1,:,:),rc_V_200vc(1,:,:),res3)

    ;--------------------500hpa--------------------
    res3@vcRefMagnitudeF        = 0.6  ;标准长度箭头所表示的大小
    res3@vcRefAnnoString2       = 0.6
    plot_500uv_PC1 =gsn_csm_vector(wks,rc_U_500vc(0,:,:),rc_V_500vc(0,:,:),res3)
    plot_500uv_PC2 =gsn_csm_vector(wks,rc_U_500vc(1,:,:),rc_V_500vc(1,:,:),res3)

    ;--------------------850hpa------------------
    res3@vcRefMagnitudeF        = 0.3    ;标准长度箭头所表示的大小
    res3@vcRefAnnoString2       = 0.3

    plot_850uv_PC1 =gsn_csm_vector(wks,rc_U_850vc(0,:,:),rc_V_850vc(0,:,:),res3)
    plot_850uv_PC2 =gsn_csm_vector(wks,rc_U_850vc(1,:,:),rc_V_850vc(1,:,:),res3)


    ;--------------------------------------高度场------------------------------------------
    res4                                =  True
    res4@gsnLeftString                  = ""
    res4@gsnRightString                 = ""   
    res4@cnLineDrawOrder                = "PostDraw"     
    res4@gsnDraw                        = False;
    res4@gsnFrame                       = False;
    res4@cnFillOn                       = False
    res4@gsnAddCyclic                   = True 
    res4@cnInfoLabelOn                  = False 
    res4@cnLinesOn                      = True      
    res4@cnLevelSelectionMode           = "ExplicitLevels"    ;ManualLevels   
    res4@cnLineLabelsOn                 = False
    res4@cnInfoLabelOn                  = False;
    res4@cnLineThicknessF               = 1.5  

    res4@cnLineColor                    = "brown1"
    res4@cnLineDashPattern              = 0  
    res4@cnLevels                       = (/6,9,12,15/)
    plot_500Z_PC1_red = gsn_csm_contour(wks,rc_Z(0,{500},:,:),res4)

    res4@cnLevels                       = (/3,5,7,9/)  
    plot_850Z_PC1_red = gsn_csm_contour(wks,rc_Z(0,{850},:,:),res4)

    res4@cnLevels                       = (/1,5,9,13/)
    plot_500Z_PC2_red = gsn_csm_contour(wks,rc_Z(1,{500},:,:),res4)

    res4@cnLevels                       =  (/1,4,7,10/)
    plot_850Z_PC2_red = gsn_csm_contour(wks,rc_Z(1,{850},:,:),res4)

    res4@cnLineColor                    = "blue1"
    res4@cnLineDashPattern              = 1
    res4@cnLevels                       = (/-1,-2,-3,-4/)
    plot_500Z_PC1_blue = gsn_csm_contour(wks,rc_Z(0,{500},:,:),res4)

    res4@cnLevels                       = (/-2,-4,-6,-8/)
    plot_850Z_PC1_blue = gsn_csm_contour(wks,rc_Z(0,{850},:,:),res4)

    res4@cnLevels                       = (/-5,-9,-13,-17/)
    plot_500Z_PC2_blue = gsn_csm_contour(wks,rc_Z(1,{500},:,:),res4)

    res4@cnLevels                       = (/-4,-6,-8,-10/)
    plot_850Z_PC2_blue = gsn_csm_contour(wks,rc_Z(1,{850},:,:),res4)


    ;----------------------地表气温------------------------
    res6                     = True
    res6@gsnLeftString        = ""
    res6@gsnRightString       = ""
    res6@gsnDraw               = False
    res6@gsnFrame              = False   
    res6@cnFillOn              = True             
    res6@cnLinesOn             = False            
    res6@cnLineLabelsOn        = False            
    res6@gsnAddCyclic          = True
    res6@cnFillPalette         = "BlueWhiteOrangeRed"  
    res6@cnLevelSelectionMode  = "ExplicitLevels"  
    res6@cnFillColors          = (/16,34,52,70,88,-1,148,166,184,202,220/) 
    res6@cnLevels              = (/-0.6,-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.5,0.6/)
    res6@cnInfoLabelOn         = False;
    res6@lbBoxLinesOn          = False      
    res6@lbBoxEndCapStyle      = "TriangleBothEnds"      
    res6@lbOrientation         = "Vertical"   
    res6@pmLabelBarWidthF      =  0.04        ;色标宽度
    res6@pmLabelBarHeightF     =  0.6            ;   色标高度
    res6@lbLabelOffsetF        =  0.1              ;色标与标签之间的距离
    res6@lbLabelFontHeightF    =  0.01          ;色标标签字体大小 
    res6@pmLabelBarOrthogonalPosF = 0.01
    ;res6@pmLabelBarParallelPosF   = 0.0 
    res6@lbLabelBarOn          =  False
    plot_tem_PC1 = gsn_csm_contour(wks,rc_tem(0,:,:),res6) 
    res6@lbLabelBarOn          =  True
    plot_tem_PC2 = gsn_csm_contour(wks,rc_tem(1,:,:),res6) 

    ;---------------------------降水-------------------------------------
    res8                       = True
    res8@gsnLeftString        = ""
    res8@gsnRightString       = ""
    res8@gsnDraw               = False
    res8@gsnFrame              = False   
    res8@cnFillPalette         = "NEO_div_vegetation_a"
    res8@cnFillOn              = True             ; turn on color
    res8@cnLinesOn             = False            ; turn off contour lines
    res8@cnLineLabelsOn        = False            ; turn off contour linelabels
    res8@cnLevelSelectionMode  = "ExplicitLevels" 
    res8@cnFillColors          = (/52,68,84,100,116,-1,146,162,178,194,210/)
    res8@cnLevels              = (/-1.6,-1.2,-0.8,-0.4,-0.2,0.2,0.4,0.8,1.2,1.6/)  
    ;res8@lbLabelStrings        = (/"-3.00","-1.50","-0.75","-0.50","-0.15","0.15","0.50","0.75","1.50","3.00"/)
    res8@gsnAddCyclic          = True
    res8@cnInfoLabelOn         = False;
    res8@lbBoxLinesOn           = False
    res8@lbBoxEndCapStyle      ="TriangleBothEnds"     
    res8@lbOrientation         = "Vertical"   
    res8@pmLabelBarWidthF      =  0.04        ;色标宽度
    res8@pmLabelBarHeightF     =  0.6            ;   色标高度
    res8@lbLabelOffsetF        =  0.1              ;色标与标签之间的距离
    res8@lbLabelFontHeightF    =  0.01          ;色标标签字体大小 
    res8@pmLabelBarOrthogonalPosF = 0.01
    ;res8@pmLabelBarParallelPosF   = 0.5 
    res8@lbLabelBarOn          =  False
    plot_pre_PC1 = gsn_csm_contour(wks,rc_pre(0,{-30:75},:),res8)
    res8@lbLabelBarOn          =  True      
    plot_pre_PC2 = gsn_csm_contour(wks,rc_pre(1,{-30:75},:),res8) 

    ;---------------------------SST海温-------------------------------------
    res9                     = True
    res9@gsnLeftString        = ""
    res9@gsnRightString       = ""  
    res9@gsnDraw               = False
    res9@gsnFrame              = False   
    res9@cnFillOn              = True             
    res9@cnLinesOn             = False            
    res9@cnLineLabelsOn        = False            
    res9@gsnAddCyclic          = True
    res9@cnFillPalette         = "MPL_seismic"   ;MPL_seismic
    res9@cnLevelSelectionMode  = "ExplicitLevels"  
    res9@cnFillColors          = (/45,50,55,60,-1,67,72,77,82/) 
    res9@cnLevels              = (/-0.4,-0.3,-0.2,-0.05,0.05,0.2,0.3,0.4/)
    res9@cnInfoLabelOn         = False;
    res9@lbBoxLinesOn           = False
    res9@lbOrientation         = "Vertical" 
    res9@lbBoxEndCapStyle      = "TriangleBothEnds"     
    res9@lbOrientation         = "Vertical"   
    res9@pmLabelBarWidthF      =  0.04        ;色标宽度
    res9@pmLabelBarHeightF     =  0.6            ;   色标高度
    res9@lbLabelOffsetF        =  0.1              ;色标与标签之间的距离
    res9@lbLabelFontHeightF    =  0.01          ;色标标签字体大小 
    res9@pmLabelBarOrthogonalPosF = 0.01
    ;res9@pmLabelBarParallelPosF   = 0.5 
    res9@lbLabelBarOn          =  False    
    plot_sst_PC1 = gsn_csm_contour(wks,rc_sst_test(0,:,:),res9)
    res9@lbLabelBarOn          =  True
    plot_sst_PC2 = gsn_csm_contour(wks,rc_sst_test(1,:,:),res9)

    ;---------------------------WAF波活动通量-------------------------------------
    res11    = True
    res11@vcVectorDrawOrder    = "PostDraw"
    res11@gsnDraw              = False;
    res11@gsnFrame             = False;
    res11@gsnAddCyclic         = True
    res11@gsnLeftString             = ""
    res11@gsnRightString            = ""     
    res11@vcGlyphStyle              = "CurlyVector"
    res11@vcLineArrowThicknessF     = 2.5   
    res11@vcLineArrowHeadMaxSizeF   = 0.024      ;箭头最大尺寸 默认0.05
    res11@vcLineArrowHeadMinSizeF   = 0.008    ;箭头最小尺寸 默认0.005
    res11@vcMinFracLengthF          = 0.9
    res11@vcLineArrowColor          = "mediumvioletred";"brown3" 

    ;****设定参考箭头****
    res11@vcRefAnnoOn               = True  
    res11@vcRefAnnoBackgroundColor  = 0   
    res11@vcRefAnnoPerimOn          =  True                                        
    res11@vcRefAnnoFontHeightF      = 0.008   
    res11@vcRefAnnoString1On        = False        
    res11@vcRefAnnoString2On        = True    
    res11@vcRefAnnoSide             = "Bottom"   
    res11@vcRefAnnoOrthogonalPosF   = -0.118
    res11@vcRefAnnoParallelPosF     = 0.952
    res11@vcRefLengthF              = 0.04    
    res11@vcMinDistanceF            = 0.028
    res11@vcRefMagnitudeF           = 0.05   ;同样长度表示的实际大小
    res11@vcRefAnnoString2          = 0.05
    plot_waf_PC1   = gsn_csm_vector(wks,Fx_PC1,Fy_PC1,res11)
    res11@vcRefMagnitudeF           = 0.1   ;同样长度表示的实际大小
    res11@vcRefAnnoString2          = 0.1    
    plot_waf_PC2   = gsn_csm_vector(wks,Fx_PC2,Fy_PC2,res11)

    ;----------------------------200hpa位势高度-------------------------------
    res12                     = True
    res12@gsnLeftString        = ""
    res12@gsnRightString       = ""
    res12@gsnDraw               = False
    res12@gsnFrame              = False   
    res12@cnFillOn              = True             
    res12@cnLinesOn             = False            
    res12@cnLineLabelsOn        = False            
    res12@gsnAddCyclic          = True
    res12@cnFillDrawOrder       = "PreDraw"
    res12@cnLevelSelectionMode  = "ExplicitLevels"  
    ;res12@cnFillPalette         = "NCV_blue_red"   ;MPL_seismic
    ;res12@cnFillColors          = (/87,100,110,120,-1,135,145,155,168/) 
    res12@cnFillPalette             = "MPL_PuOr"  
    res12@cnFillColors          =  (/110,105,95,85,72, -1, 58,45,35,25,20/)
    res12@cnLevels              =  (/-25,-20,-15,-10,-5,5,10,15,20,25/)
    ;res12@cnLevels              =  (/-3.,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.,3./)
    res12@cnInfoLabelOn         = False;
    res12@lbLabelBarOn          =  True
    res12@lbBoxLinesOn           = False
    res12@lbOrientation         = "Vertical"   
    res12@lbBoxEndCapStyle      = "TriangleBothEnds"  
    res12@pmLabelBarWidthF      =  0.04        ;色标宽度
    res12@pmLabelBarHeightF     =  0.6            ;   色标高度
    res12@lbLabelOffsetF        =  0.1              ;色标与标签之间的距离
    res12@lbLabelFontHeightF    =  0.01          ;色标标签字体大小 
    res12@pmLabelBarOrthogonalPosF = 0.01
    ;res12@pmLabelBarParallelPosF   = 0.5 

    res12@lbLabelBarOn          =  False
    plot_200Z_PC1 = gsn_csm_contour(wks,rc_Z(0,{200},:,:),res12)
    ;plot_200Z_PC1 = gsn_csm_contour(wks,sf_pc(0,:,:,{200}),res12)

    res12@lbLabelBarOn          =  True
    plot_200Z_PC2 = gsn_csm_contour(wks,rc_Z(1,{200},:,:),res12)
    ;plot_200Z_PC2 = gsn_csm_contour(wks,sf_pc(1,:,:,{200}),res12)
   


    ;overlay(base(0),plot_200Z_PC1)      
    overlay(base(0),plot_200uv_PC1) 
    overlay(base(0),plot_waf_PC1)

    overlay(base(2),plot_tem_PC1) 
    overlay(base(2),plot_tem_prob_PC1)
    overlay(base(2),plot_500Z_PC1_blue)
    overlay(base(2),plot_500Z_PC1_red)
    overlay(base(2),plot_500uv_PC1) 

    overlay(base(4),plot_pre_PC1)
    overlay(base(4),plot_pre_prob_PC1) 
    overlay(base(4),plot_850Z_PC1_red)
    overlay(base(4),plot_850Z_PC1_blue)
    overlay(base(4),plot_850uv_PC1)

    ;overlay(base(1),plot_200Z_PC2)      
    overlay(base(1),plot_200uv_PC2) 
    overlay(base(1),plot_waf_PC2)

    overlay(base(3),plot_tem_PC2) 
    overlay(base(3),plot_tem_prob_PC2)
    overlay(base(3),plot_500Z_PC2_blue)
    overlay(base(3),plot_500Z_PC2_red)
    overlay(base(3),plot_500uv_PC2) 

    overlay(base(5),plot_pre_PC2)
    overlay(base(5),plot_pre_prob_PC2) 
    overlay(base(5),plot_850Z_PC2_red)
    overlay(base(5),plot_850Z_PC2_blue)
    overlay(base(5),plot_850uv_PC2) 

    ;---------------------AC气旋反气旋标志--------------------
    res_AC                                =  True
    res_AC@gsnLeftString                  = ""
    res_AC@gsnRightString                 = ""      
    res_AC@gsnDraw                        = False;
    res_AC@gsnFrame                       = False;
    res_AC@cnFillOn                       = False
    res_AC@gsnAddCyclic                   = True 
    res_AC@cnInfoLabelOn                  = False 
    res_AC@txAngleF                       = 0
    res_AC@txFontHeightF                  = 0.015
    res_AC@txFontThicknessF               = 3.0
    res_AC@txFontColor                    = "firebrick1"
    res_AC@txBackgroundFillColor          = -1
    ;-------------------回归场PC1  CCCCCCCC-----------------------;
    text850_PC1_c = gsn_add_text(wks,base(4),"~F22~C",-30,19,res_AC)

    text200_PC1_1 = gsn_add_text(wks,base(0),"~F22~C",-18,57,res_AC)
    text200_PC1_2 = gsn_add_text(wks,base(0),"~F22~C",65,44,res_AC)

    text500_PC1_1 = gsn_add_text(wks,base(2),"~F22~C",-15,55,res_AC)
    text500_PC1_2 = gsn_add_text(wks,base(2),"~F22~C",62,42,res_AC)

    text850_PC1_1 = gsn_add_text(wks,base(4),"~F22~C",-15,55,res_AC)
    text850_PC1_2 = gsn_add_text(wks,base(4),"~F22~C",68,42,res_AC)
    ;-------------------回归场PC2  CCCCCCCCC----------------------;
    text850_PC2_c = gsn_add_text(wks,base(5),"~F22~C",-94,22,res_AC)

    text200_PC2_1 = gsn_add_text(wks,base(1),"~F22~C",-63,72,res_AC)
    text200_PC2_2 = gsn_add_text(wks,base(1),"~F22~C",15,49,res_AC)
    text200_PC2_3 = gsn_add_text(wks,base(1),"~F22~C",100,53,res_AC)

    text500_PC2_1 = gsn_add_text(wks,base(3),"~F22~C",-63,67,res_AC)
    text500_PC2_2 = gsn_add_text(wks,base(3),"~F22~C",15,49,res_AC)
    text500_PC2_3 = gsn_add_text(wks,base(3),"~F22~C",96,58,res_AC)

    text850_PC2_1 = gsn_add_text(wks,base(5),"~F22~C",-59,67,res_AC)
    text850_PC2_2 = gsn_add_text(wks,base(5),"~F22~C",23,55,res_AC)
    text850_PC2_3 = gsn_add_text(wks,base(5),"~F22~C",104,57,res_AC)



    res_AC@txFontColor                    = "blue1"
    ;-------------------回归场PC1 AAAAAAAAAA-----------------------;
    text200_PC1_4 = gsn_add_text(wks,base(0),"~F22~A",-67,44,res_AC)
    text200_PC1_5 = gsn_add_text(wks,base(0),"~F22~A",32,50,res_AC)
    text200_PC1_6 = gsn_add_text(wks,base(0),"~F22~A",105,44,res_AC)

    text500_PC1_4 = gsn_add_text(wks,base(2),"~F22~A",-55,43,res_AC)
    text500_PC1_5 = gsn_add_text(wks,base(2),"~F22~A",43,55,res_AC)
    text500_PC1_6 = gsn_add_text(wks,base(2),"~F22~A",110,45,res_AC)

    text850_PC1_4 = gsn_add_text(wks,base(4),"~F22~A",-57,38,res_AC)
    text850_PC1_5 = gsn_add_text(wks,base(4),"~F22~A",51,57,res_AC)
    text850_PC1_6 = gsn_add_text(wks,base(4),"~F22~A",120,45,res_AC)

    ;-------------------回归场PC2 AAAAAAAAAA-----------------------;
    text200_PC2_4 = gsn_add_text(wks,base(1),"~F22~A",-89,50,res_AC)
    text200_PC2_5 = gsn_add_text(wks,base(1),"~F22~A",-22,55,res_AC)
    text200_PC2_6 = gsn_add_text(wks,base(1),"~F22~A",61,53,res_AC)
    text200_PC2_7 = gsn_add_text(wks,base(1),"~F22~A",130,38,res_AC)

    text500_PC2_4 = gsn_add_text(wks,base(3),"~F22~A",-90,48,res_AC)
    text500_PC2_5 = gsn_add_text(wks,base(3),"~F22~A",-21,53,res_AC)
    text500_PC2_6 = gsn_add_text(wks,base(3),"~F22~A",57,59,res_AC)
    text500_PC2_7 = gsn_add_text(wks,base(3),"~F22~A",134,37,res_AC)

    text850_PC2_4 = gsn_add_text(wks,base(5),"~F22~A",-88,43,res_AC)
    text850_PC2_5 = gsn_add_text(wks,base(5),"~F22~A",-21,49,res_AC)
    text850_PC2_6 = gsn_add_text(wks,base(5),"~F22~A",58,61,res_AC)
    text850_PC2_7 = gsn_add_text(wks,base(5),"~F22~A",139,36,res_AC)

    ;----------------------------850hpa位青藏高原叠加-------------------------------
    shptp="/nuist/scratch/zhuzhiwei/shaoly/Map/TP/TP_polygon.shp"    ;修改为你的shp文件存放路径
    tpres        = True
    tpres@gsFillColor = "grey"
    shp_plot1     = gsn_add_shapefile_polygons(wks,base(4),shptp,tpres)        ;base为底图
    shp_plot2     = gsn_add_shapefile_polygons(wks,base(5),shptp,tpres)        ;base为底图

    resP = True                        ; 绘制panel图      
    resP@gsnPanelRowSpec = True
    resP@gsnMaximize = True
    resP@vpXF             = 0.01
    resP@vpYF            = 0.99
    resP@gsnPanelXWhiteSpacePercent = 5.0
    resP@gsnPanelYWhiteSpacePercent = 1.0    
    resP@gsnPanelYF               =  (/0.94,0.94,0.7,0.7,0.46,0.46/)
    resP@gsnPanelXF               = (/0.05,0.53,0.05,0.53,0.05,0.53/)
    resP@gsnPanelFigureStrings    = (/"~F22~a","~F22~d","~F22~b","~F22~e","~F22~c","~F22~f"/)
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF =  0.02    
    resP@gsnPanelFigureStringsJust    = "TopLeft"
    resP@amOrthogonalPosF         = -0.72        ;y轴方向往上方向调整
    resP@amParallelPosF           = -0.56        ;x轴方向往负方向调整
    gsn_panel(wks,base,(/2,2,2/),resP) 
    frame(wks)

    end