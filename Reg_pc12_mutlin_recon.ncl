load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
;---------------------------------------脚本编写--------------------------------------------

begin
    latS      = -90
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
    pc1 =asciiread("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/pc1.txt", -1, "float")
    AMO =asciiread("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/AMO.txt", -1, "float")
    CEP_PC1 =asciiread("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ZDY_CEP_PC1.txt", -1, "float")

    pc2 =asciiread("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/pc2.txt", -1, "float")
    IOD =asciiread("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ZDY_IOD.txt", -1, "float")
    CEP_PC2 =asciiread("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ZDY_CEP_PC2.txt", -1, "float")

     np   = 2

    xp1   = new( (/43,np/), "float" )  
    xp1(:,0) = CEP_PC1
    xp1(:,1) = AMO
         
    opt = True
    opt@print_anova = True     
    opt@print_data  = True
    beta_pc1   =  reg_multlin_stats(pc1,xp1,opt) 
    ;print(beta_pc1)

    pc1_reconstru = new((/nyears/),"float")
    do i = 0,nyears-1,1
    pc1_reconstru(i)=beta_pc1(1)*CEP_PC1(i)+beta_pc1(2)*AMO(i)
    end do
    ;print(pc1_reconstru)


    xp2   = new( (/43,np/), "float" )  
    xp2(:,0) = CEP_PC2
    xp2(:,1) = IOD

    beta_pc2   =  reg_multlin_stats(pc2,xp2,opt) 
    ;print(beta_pc2)

    pc2_reconstru = new((/nyears/),"float")
    do i = 0,nyears-1,1
    pc2_reconstru(i)=beta_pc2(1)*CEP_PC2(i)+beta_pc2(2)*IOD(i)
    end do
    ;print(pc2_reconstru)



     soizs        = new((/4,nyears/),"float")
     soizs(0,:)   = pc1
     soizs(1,:)   = pc2
     soizs(2,:)   = pc1_reconstru
     soizs(3,:)   = pc2_reconstru
     soizs!0      = "sv"
     soizs!1      = "time"
     ;printVarSummary(soizs) 

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;读uuuuu场数据;;;;;;;;;;;;;;;;
     u_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_u_1959_2022_JJA_14473.nc", "r")
     ;print(u_in)
     time = u_in->time
     XXXX     = cd_calendar(time,-1 )
     ;print(XXXX)
     t_start = ind(XXXX.eq.197906)
     t_end   = ind(XXXX.eq.202108)     
     u       = u_in->u
     ;printVarSummary(u)  
     u      := u(t_start:t_end,:,{latS:latN},:) 
     u_flt  = short2flt(u)
     ;printVarSummary(u_flt)
     u_JJA = new((/nyears,3,73,144/),"float")
     u_JJA!0="time"
     u_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     u_JJA(nn,:,:,:)= (u_flt(nn*3,:,:,:)+u_flt((nn*3)+1,:,:,:)+u_flt((nn*3)+2,:,:,:))/3.0
     end do
     copy_VarMeta(u_flt(0,:,:,:),u_JJA(0,:,:,:))
     printVarSummary(u_JJA)

    ;---------------------------------------------------------
     v_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/ERA5_monthly_v_1959_2022_JJA_14473.nc", "r")
     ;print(v_in)
     time = v_in->time
     XXXX     = cd_calendar(time,-1 )
     ;print(XXXX) 
     v       = v_in->v
     ;printVarSummary(v)  
     v      := v(t_start:t_end,:,{latS:latN},:) 
     v_flt  = short2flt(v)
     ;printVarSummary(v_flt)
     v_JJA = new((/nyears,3,73,144/),"float")
     v_JJA!0="time"
     v_JJA&time = ispan(1979,2021,1)  
     do nn=0,42,1
     v_JJA(nn,:,:,:)= (v_flt(nn*3,:,:,:)+v_flt((nn*3)+1,:,:,:)+v_flt((nn*3)+2,:,:,:))/3.0
     end do
     copy_VarMeta(v_flt(0,:,:,:),v_JJA(0,:,:,:))
     printVarSummary(v_JJA)
 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    precip_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/precip.mon.mean.nc", "r")
    ;print(precip_in)                                ;GPCP数据197901-202212月
    time_precip = precip_in->time
    XXXX2     = cd_calendar(time_precip,-1 )
    ;print(XXXX2)
    precip      = precip_in->precip(:,{-30:90},:)
    ;printVarSummary( precip)
    precip_JJA = new((/nyears,48,144/),"float")
    precip_JJA!0="time"
    precip_JJA&time = ispan(1979,2021,1)  
    do nn=0,42,1                      
    precip_JJA(nn,:,:)= (precip((nn*12)+5,:,:)+precip((nn*12)+6,:,:)+precip((nn*12)+7,:,:))/3.0
    end do
    copy_VarCoords(precip(0,:,:),precip_JJA(0,:,:))
    ;printVarSummary(precip_JJA)

    sst_in   = addfile("/nuist/scratch/zhuzhiwei/shaoly/data/shuju/HadISST_sst.nc", "r")
    ;print(sst_in)
    time_sst = sst_in->time
    XXXX3     = cd_calendar(time_sst,-1 )
    t_start1 = ind(XXXX3.eq.197801)
    t_end1   = ind(XXXX3.eq.202212)
    sst      = sst_in->sst(t_start1:t_end1,{-30:90},:)
    ;printVarSummary(sst)
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
    ;printVarSummary(sst_JJA)

    lat = u_JJA&lat
    lon = u_JJA&lon
    vrdv = uv2dv_cfd(u_JJA, v_JJA, lat, lon, 2)
    printVarSummary(vrdv)
    dv = vrdv(:,:,:,:)      ; 散度场
    dv := smth9_Wrap(dv, 0.50, 0.25, True)
    copy_VarCoords(v_JJA,dv)
    printVarSummary(dv)
     printMinMax(dv, 0)

    uvd = dv2uvF_Wrap(dv)      ;辐散风
    printVarSummary(uvd)
    u_uvd = uvd(0,:,:,:,:)
    printVarSummary(u_uvd)
    v_uvd = uvd(1,:,:,:,:)
    printVarSummary(v_uvd)


    u_uvd_pc = regCoef_n(soizs,u_uvd(lat|:,lon|:,level|:,time|:),1,3) 
    u_uvd_pc!0 = "sv"
    u_uvd_pc!1 = "lat"
    u_uvd_pc!2 = "lon"
    u_uvd_pc!3 = "level"    
    u_uvd_pc&lat = u_uvd&lat
    u_uvd_pc&lon = u_uvd&lon
    u_uvd_pc&level = u_uvd&level
    printVarSummary(u_uvd_pc)   ;[sv | 2] x [lat | 49] x [lon | 144] x [level | 3]
    printMinMax(u_uvd_pc, 0)


    v_uvd_pc = regCoef_n(soizs,v_uvd(lat|:,lon|:,level|:,time|:),1,3) 
    v_uvd_pc!0 = "sv"
    v_uvd_pc!1 = "lat"
    v_uvd_pc!2 = "lon"
    v_uvd_pc!3 = "level"    
    v_uvd_pc&lat = v_uvd&lat
    v_uvd_pc&lon = v_uvd&lon
    v_uvd_pc&level = v_uvd&level
    printVarSummary(v_uvd_pc)   ;[sv | 2] x [lat | 49] x [lon | 144] x [level | 3]
    printMinMax(v_uvd_pc, 0)



;***********************************************************计算回归 ****************************************************
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 计算200hpa回归A0的回归系数并做检验;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    U_pc = regCoef_n(soizs,u_JJA(lat|:,lon|:,level|:,time|:),1,3) 
    U_pc!0 = "sv"
    U_pc!1 = "lat"
    U_pc!2 = "lon"
    U_pc!3 = "level"    
    U_pc&lat = u_JJA&lat
    U_pc&lon = u_JJA&lon
    U_pc&level = u_JJA&level
    ;printVarSummary(U_pc)   ;[sv | 2] x [lat | 49] x [lon | 144] x [level | 3]

    U_prob_pc_1 = student_t(U_pc@tval,U_pc@nptxy-2)
    U_prob_pc   = onedtond(U_prob_pc_1,dimsizes(U_pc))
    copy_VarCoords(U_pc,U_prob_pc)
    ;printVarSummary(U_prob_pc)  ;[sv | 2] x [lat | 49] x [lon | 144] x [level | 3]

    V_pc = regCoef_n(soizs,v_JJA(lat|:,lon|:,level|:,time|:),1,3)    
    V_pc!0 = "sv"
    V_pc!1 = "lat"
    V_pc!2 = "lon"
    V_pc!3 = "level"   
    V_pc&lat = v_JJA&lat
    V_pc&lon = v_JJA&lon
    V_pc&level = v_JJA&level
    ;printVarSummary(V_pc)  ;[sv | 2] x [lat | 49] x [lon | 144] x [level | 3]

    V_prob_pc_1 = student_t(V_pc@tval,V_pc@nptxy-2)
    V_prob_pc   = onedtond(V_prob_pc_1,dimsizes(V_pc))
    copy_VarCoords(V_pc, V_prob_pc)
    ;printVarSummary(V_prob_pc)     ;[sv | 2] x [lat | 49] x [lon | 144] x [level | 3]

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;计算PRECIP-IO回归系数并做检验
    pre_pc = regCoef_n(soizs,precip_JJA(lat|:,lon|:,time|:),1,2)     ;返回值为tval(t值统计量tval=rc@tval)  nptxy(自由度df=rc@nptxy-2）
    pre_pc!0 = "sv"
    pre_pc!1 = "lat"
    pre_pc!2 = "lon"
    pre_pc&lat = precip_JJA&lat
    pre_pc&lon = precip_JJA&lon
    ;printVarSummary(pre_pc)  ;[sv | 2] x [lat | 48] x [lon | 144]

    pre_prob_pc_1 = student_t(pre_pc@tval,pre_pc@nptxy-2)
    pre_prob_pc = onedtond(pre_prob_pc_1,dimsizes(pre_pc))
    copy_VarCoords(pre_pc, pre_prob_pc)
    ;printVarSummary(pre_prob_pc)   ;[sv | 2] x [lat | 48] x [lon | 144]


   ;-------------------------------------------------
   sst_pc = regCoef_n(soizs,sst_JJA(latitude|:,longitude|:,time|:),1,2)   
   sst_pc!0 = "sv"
   sst_pc!1 = "lat"
   sst_pc!2 = "lon"
   sst_pc&lat = sst_JJA&latitude
   sst_pc&lon = sst_JJA&longitude
   ;printVarSummary(sst_pc)   ;[sv | 2] x [lat | 120] x [lon | 360]

   sst_prob_pc_1 = student_t(sst_pc@tval,sst_pc@nptxy-2)
   sst_prob_pc   = onedtond(sst_prob_pc_1,dimsizes(sst_pc))
   copy_VarCoords(sst_pc, sst_prob_pc)
   ;printVarSummary(sst_prob_pc)   ; [sv | 2] x [lat | 120] x [lon | 360]
   sst_pc_test = where(abs(sst_pc).lt.5,sst_pc,sst_pc@_Fillvalue)
   copy_VarCoords(sst_pc, sst_pc_test)
   ;printVarSummary(sst_pc_test)    ;[sv | 2] x [lat | 120] x [lon | 360]

  u_200_pc = U_pc(:,:,:,{200})
  v_200_pc = V_pc(:,:,:,{200})


  u_850_pc = U_pc(:,:,:,{850})
  v_850_pc = V_pc(:,:,:,{850})
  u_850_pc := smth9_Wrap(u_850_pc, 0.50, 0.25, True)
  v_850_pc := smth9_Wrap(v_850_pc, 0.50, 0.25, True)

  u_200_pc=where(sqrt((/u_200_pc/)^2+(/v_200_pc/)^2).lt.0.4, u_200_pc@_FillValue, u_200_pc)
  v_200_pc=where(sqrt((/u_200_pc/)^2+ (/v_200_pc/)^2).lt.0.4, v_200_pc@_FillValue, v_200_pc)

  u_850_pc=where(sqrt((/u_850_pc/)^2+(/v_850_pc/)^2).lt.0.02, u_850_pc@_FillValue, u_850_pc)
  v_850_pc=where(sqrt((/u_850_pc/)^2+ (/v_850_pc/)^2).lt.0.02, v_850_pc@_FillValue, v_850_pc)

  u_850_pc(1,{-10:16},{-10:50}) = U_pc(1,{-10:16},{-10:50},{850})       ;第二模态北非东风很重要   不要mask
  v_850_pc(1,{-10:16},{-10:50}) = V_pc(1,{-10:16},{-10:50},{850})




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;画图参数设置;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
wks = gsn_open_wks("eps","Reg_pc12_mutlin_recon")
base = new(6,"graphic") 
gsn_define_colormap(wks,"cb_9step")
 res                       = True
   res@gsnDraw               = False
   res@gsnFrame              = False   
   res@gsnAddCyclic          = True

   res@cnFillOn              = True             ; turn on color
   res@cnLinesOn             = False             ; turn off contour lines
   res@cnLineLabelsOn        = False            ; turn off contour linelabels
   res@cnLevelSelectionMode  = "ManualLevels"   ; set manual contour levels
   ;res@pmTickMarkDisplayMode = "Always"

    res@mpOutlineOn                 = True       ; turn off outline
    res@mpFillOn                    = False    ; turn off map fill   
    res@mpMinLatF                   = -30.
    res@mpMaxLatF                   = 75.
    res@mpMinLonF                   = -200
    res@mpMaxLonF                   = 160
    res@mpCenterLonF                = (-200+160)/2.
    res@cnSmoothingOn        = True
    res@lbLabelFont          = 0.25

    res@tmXTOn = False
    res@tmYROn = False

    res@mpOutlineOn                  = True       ; turn off outline
    res@mpFillOn                     = False    ; turn off map fill 
    res@mpOutlineBoundarySets        ="Geophysical"
    res@mpGeophysicalLineColor       = "gray50"
    res@mpGeophysicalLineDashPattern = 0
    res@mpGeophysicalLineDashSegLenF = 0.2  
    res@mpGeophysicalLineThicknessF  = 1.5  ;地图边界的粗细

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
    res@tmXBValues        = (/-180, -120, -60, 0, 60, 120 /)
    res@tmXBLabels        = (/"~F21~180~S~o~N~", "~F21~120~S~o~N~W", "~F21~60~S~o~N~W", "~F21~0", "~F21~60~S~o~N~E", "~F21~120~S~o~N~E"/)


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
    res@tmYLValues        = (/-20,0,20,40,60/)
    res@tmYLLabels        = (/"~F21~20~S~o~N~S","~F21~0","~F21~20~S~o~N~N","~F21~40~S~o~N~N","~F21~60~S~o~N~N"/)

    res@mpShapeMode = "FreeAspect"
    res@vpHeightF  = 0.14
    res@vpWidthF   = 0.33
    res@txFont   = "helvetica"

    res@gsnRightString              =   " "
    res@gsnLeftStringFontHeightF      = 0.01
    res@gsnLeftStringOrthogonalPosF   = 0.015
    res@lbLabelBarOn                =  False
    res@cnFillOn                    =  False
    res@cnInfoLabelOn               = False

   res@gsnLeftString          =   "~F22~Reg 850Div Winds&SST onto PC1"
   base(0) = gsn_csm_contour_map(wks,pre_pc(0,:,:),res) 

   res@gsnLeftString          =   "~F22~Reg 850Div Winds&SST onto PC2"
   base(1) = gsn_csm_contour_map(wks,pre_pc(0,:,:),res)    

   res@gsnLeftString          =   "~F22~Reg 850Div Winds&SST onto PC1_Rec"   
   base(4) = gsn_csm_contour_map(wks,pre_pc(0,:,:),res)   

   res@gsnLeftString          =   "~F22~Reg 850Div Winds&SST onto PC2_Rec"   
   base(5) = gsn_csm_contour_map(wks,pre_pc(0,:,:),res) 



;---------------------显著性打点---------------
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
   res2@cnFillColors         = (/0,-1/)                        ;填充颜色-1为透明;
   res2@cnFillScaleF         =   0.6
   res2@cnFillDotSizeF       =  0.002
   res2@cnLevels             = (/0.1/)  

   plot_pre_prob_pc1 = gsn_csm_contour(wks,pre_prob_pc(0,:,:),res2)       
   plot_pre_prob_pc2 = gsn_csm_contour(wks,pre_prob_pc(1,:,:),res2)  

   plot_sst_prob_pc1 = gsn_csm_contour(wks,sst_prob_pc(0,:,:),res2)    
   plot_sst_prob_pc2 = gsn_csm_contour(wks,sst_prob_pc(1,:,:),res2) 
   plot_sst_prob_pc1_rec = gsn_csm_contour(wks,sst_prob_pc(2,:,:),res2)    
   plot_sst_prob_pc2_rec = gsn_csm_contour(wks,sst_prob_pc(3,:,:),res2) 

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
    res3@vcLineArrowColor            = "black"

    res3@vcMinFracLengthF          = 1.0
    res3@vcLineArrowHeadMaxSizeF   = 0.015
    res3@vcLineArrowHeadMinSizeF   = 0.01

   ;****设定参考箭头****
    res3@vcRefAnnoOn               = True  
    res3@vcRefLengthF              = 0.03  ;标准长度箭头在单位平方中的大小
    res3@vcRefAnnoBackgroundColor  = 0 ;背景颜色     
    res3@vcRefAnnoPerimOn          = True     ;关闭边框                                       
    res3@vcRefAnnoFontHeightF      = 0.008   ;参考箭头标签字体大小      
    
    res3@vcRefAnnoString1On     = False   ;设定参考箭头上、下的字符        
    res3@vcRefAnnoString2On     = True    ; 这里仅设定其下方的字符
    res3@vcRefAnnoSide            = "Top" ; 参考箭头放至图形上方
    res3@vcRefAnnoOrthogonalPosF  = -0.225
    res3@vcRefMagnitudeF        =  1.2     ;标准长度箭头所表示的大小
    res3@vcRefAnnoString2       = "~F21~1.2"
    res3@vcMinDistanceF            = 0.023
    plot_200uv_pc1 = gsn_csm_vector(wks,u_200_pc(0,{-30:30},:),v_200_pc(0,{-30:30},:),res3)
    plot_200uv_pc2 = gsn_csm_vector(wks,u_200_pc(1,{-30:30},:),v_200_pc(1,{-30:30},:),res3)

    res3@vcRefMagnitudeF        = 0.5    ;标准长度箭头所表示的大小
    res3@vcRefAnnoString2       = "~F21~0.5"
    res3@vcMinDistanceF            = 0.028
    ;plot_850uv_pc1 = gsn_csm_vector(wks,u_850_pc(0,{-30:30},:),v_850_pc(0,{-30:30},:),res3)
    ;plot_850uv_pc2 = gsn_csm_vector(wks,u_850_pc(1,{-30:30},:),v_850_pc(1,{-30:30},:),res3)
    plot_850uv_pc1 = gsn_csm_vector(wks,u_uvd_pc(0,{-30:35},:,{850}),v_uvd_pc(0,{-30:35},:,{850}),res3)
    plot_850uv_pc2 = gsn_csm_vector(wks,u_uvd_pc(1,{-30:35},:,{850}),v_uvd_pc(1,{-30:35},:,{850}),res3)

    ;plot_850uv_pc1_rec = gsn_csm_vector(wks,u_850_pc(2,{-30:30},:),v_850_pc(2,{-30:30},:),res3)
    ;plot_850uv_pc2_rec = gsn_csm_vector(wks,u_850_pc(3,{-30:30},:),v_850_pc(3,{-30:30},:),res3)
    plot_850uv_pc1_rec = gsn_csm_vector(wks,u_uvd_pc(2,{-30:35},:,{850}),v_uvd_pc(2,{-30:35},:,{850}),res3)
    plot_850uv_pc2_rec = gsn_csm_vector(wks,u_uvd_pc(3,{-30:35},:,{850}),v_uvd_pc(3,{-30:35},:,{850}),res3)
  ;-------------------降水------------------------
    res8                       = True
    res8@gsnLeftString        = ""
    res8@gsnRightString       = ""
    res8@vpHeightF  = 0.14
    res8@vpWidthF   = 0.33
    res8@gsnDraw               = False
    res8@gsnFrame              = False      
    res8@gsnAddCyclic          = True
    res8@cnFillPalette         = "NEO_div_vegetation_a"
    res8@cnFillOn              = True             ; turn on color
    res8@cnLinesOn             = False            ; turn off contour lines
    res8@cnLineLabelsOn        = False            ; turn off contour linelabels
    res8@cnLevelSelectionMode  = "ExplicitLevels" 
    res8@cnFillColors          = (/52,68,84,100,116,-1,146,162,178,194,210/)
    res8@cnLevels              = (/-1.6,-1.2,-0.8,-0.4,-0.2,0.2,0.4,0.8,1.2,1.6/)  
    res8@cnInfoLabelOn         = False;
    res8@lbLabelBarOn          = False 
    res8@lbBoxLinesOn          = False
    res8@lbOrientation         = "Vertical"
    res8@lbBoxEndCapStyle      ="TriangleBothEnds"     
    res8@pmLabelBarWidthF      =   0.02        ;色标宽度
    res8@pmLabelBarHeightF     =  0.14             ;   色标高度
    res8@lbLabelOffsetF        =  0.1              ;色标与标签之间的距离
    res8@lbLabelFontHeightF    =  0.005           ;色标标签字体大小    
    res8@pmLabelBarOrthogonalPosF = 0.01
    ;res8@pmLabelBarParallelPosF   = 0.5 
    plot_pre_pc1 = gsn_csm_contour(wks,pre_pc(0,:,:),res8)
    res8@lbLabelBarOn          = True      
    plot_pre_pc2 = gsn_csm_contour(wks,pre_pc(1,:,:),res8) 

  ;-------------------SST海温-----------------------
   res9                     = True
   res9@gsnLeftString        = ""
   res9@gsnRightString       = ""  
   res9@vpHeightF  = 0.14
   res9@vpWidthF   = 0.33
   res9@gsnDraw               = False
   res9@gsnFrame              = False   
   res9@cnFillOn              = True             
   res9@cnLinesOn             = False            
   res9@cnLineLabelsOn        = False            
   res9@gsnAddCyclic          = True
   res9@cnLevelSelectionMode  = "ExplicitLevels"  
   res9@cnFillPalette         = "MPL_RdBu"  
   res9@cnFillColors          = (/114,103,92,81,70,-1,55,44,33,22,11/) 
   res9@cnLevels              = (/-0.5,-0.4,-0.3,-0.15,-0.05,0.05,0.15,0.3,0.4,0.5/)
   res9@cnLevels              = (/-0.5,-0.4,-0.3,-0.15,-0.05,0.05,0.15,0.3,0.4,0.5/)
   res9@cnInfoLabelOn         = False;
   res9@lbLabelBarOn          =  False 
   res9@lbBoxLinesOn          = False
   res9@lbBoxEndCapStyle      = "TriangleBothEnds"      
   res9@lbOrientation         = "Vertical"   
   res9@pmLabelBarWidthF      =   0.02         ;色标宽度
   res9@pmLabelBarHeightF     =  0.14            ;   色标高度
   res9@lbLabelOffsetF        =  0.1              ;色标与标签之间的距离
   res9@lbLabelFontHeightF    =  0.005           ;色标标签字体大小 
   res9@pmLabelBarOrthogonalPosF = 0.01
   ;res9@pmLabelBarParallelPosF   = 0.01 
   plot_sst_pc1 = gsn_csm_contour(wks,sst_pc_test(0,:,:),res9)
   plot_sst_pc1_rec = gsn_csm_contour(wks,sst_pc_test(2,:,:),res9) 

   res9@lbLabelBarOn          =  True 
   plot_sst_pc2 = gsn_csm_contour(wks,sst_pc_test(1,:,:),res9)
   plot_sst_pc2_rec = gsn_csm_contour(wks,sst_pc_test(3,:,:),res9)



 

    ;----------------------------850hpa位青藏高原叠加-------------------------------
    shptp="/nuist/scratch/zhuzhiwei/shaoly/Map/TP/TP_polygon.shp"    ;修改为你的shp文件存放路径
    tpres        = True
    tpres@gsFillColor = "grey"
    shp_plot1     = gsn_add_shapefile_polygons(wks,base(0),shptp,tpres)        ;base为底图
    shp_plot2     = gsn_add_shapefile_polygons(wks,base(1),shptp,tpres)        ;base为底图
    shp_plot3     = gsn_add_shapefile_polygons(wks,base(4),shptp,tpres)        ;base为底图
    shp_plot4     = gsn_add_shapefile_polygons(wks,base(5),shptp,tpres)        ;base为底图

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
      res_text_AC@txFontHeightF                  = 0.015
      res_text_AC@txFontThicknessF               = 2.0
      res_text_AC@txFontColor                    = "brown1"
      res_text_AC@txBackgroundFillColor          = -1
      res_text_AC@txFont                         = "helvetica"
      ;text850_PC1_c = gsn_add_text(wks,base(0),"~F22~C",-30,19,res_text_AC)
      ;text850_PC2_c = gsn_add_text(wks,base(1),"~F22~C",-94,22,res_text_AC)
      ;text850_PC1_rec_c = gsn_add_text(wks,base(4),"~F22~C",-30,19,res_text_AC)
      ;text850_PC2_rec_c = gsn_add_text(wks,base(5),"~F22~C",-94,22,res_text_AC)

    overlay(base(0),plot_sst_pc1) 
    overlay(base(0),plot_sst_prob_pc1) 
    overlay(base(0),plot_850uv_pc1) 

    overlay(base(1),plot_sst_pc2) 
    overlay(base(1),plot_sst_prob_pc2) 
    overlay(base(1),plot_850uv_pc2) 


    overlay(base(4),plot_sst_pc1_rec) 
    overlay(base(4),plot_sst_prob_pc1_rec) 
    overlay(base(4),plot_850uv_pc1_rec) 


    overlay(base(5),plot_sst_pc2_rec) 
    overlay(base(5),plot_sst_prob_pc2_rec) 
    overlay(base(5),plot_850uv_pc2_rec)



    pc1_TC = new((/3,nyears/), "float")
    pc1_TC(0,:) = pc1_reconstru
    pc1_TC(1,:) = AMO
    pc1_TC(2,:) = CEP_PC1
    ;printVarSummary(pc1_TC)
    r_PC1_AMO = esccr(pc1, AMO,0)
    r_PC1_CEP = esccr(pc1, CEP_PC1,0)
    r_pc1_recon = esccr(pc1_reconstru,pc1,0) 

    rmsd_PC1_AMO = dim_rmsd(pc1,AMO)
    rmsd_PC1_CEP = dim_rmsd(pc1,CEP_PC1)
    rmsd_PC1_recon = dim_rmsd(pc1,pc1_reconstru)

    pc2_TC = new((/3,nyears/), "float")
    pc2_TC(0,:) = pc2_reconstru
    pc2_TC(1,:) = IOD
    pc2_TC(2,:) = CEP_PC2
    ;printVarSummary(pc2_TC)
    r_PC2_IOD = esccr(pc2, IOD,0)
    r_PC2_CEP = esccr(pc2, CEP_PC2,0)
    r_pc2_recon = esccr(pc2_reconstru,pc2,0) 

    rmsd_PC2_IOD = dim_rmsd(pc2,IOD)
    rmsd_PC2_CEP = dim_rmsd(pc2,CEP_PC2)
    rmsd_PC2_recon = dim_rmsd(pc2,pc2_reconstru)



    ;------------------------时间序列 拟合序列---------------------------------------
    rts          = True
    rts@gsnDraw  = False       ; dont draw yet
    rts@gsnFrame = False       ; dont advance frame yet
    rts@tmXTOn = False
    rts@tmYROn = False
    rts@vpHeightF  = 0.14
    rts@vpWidthF   = 0.33
    rts@txFont   = "helvetica"

    rts@gsnYRefLine       = 0.              ; reference line
    rts@trXMinF          = yrStrt 
    rts@trXMaxF          = yrLast
    rts@tmXBMode         = "Manual"
    rts@tmXBTickStartF   = 1979
    rts@tmXBTickEndF     = 2021
    rts@tmXBTickSpacingF = 5
    rts@tmXBMinorPerMajor= 4  ;次刻度线数量

    rts@trYMinF          = -4.0        ; min value on x-axis
    rts@trYMaxF          = 4.0        ; max value on x-axis
    rts@tmYLTickStartF   = -4
    rts@tmYLTickEndF     = 3.0
    rts@tmYLTickSpacingF = 1.5
    rts@tmYLMinorPerMajor= 1  ;次刻度线数量

    rts@tmXBMinorOn               = True
    rts@tmXBMinorThicknessF       = 0.
    rts@tmXBMinorLengthF          = 0.
    rts@tmXBMinorOutwardLengthF   = 0.
    rts@tmXBMajorThicknessF       = 3.0
    rts@tmXBMajorLengthF          = 0.003
    rts@tmXBMajorOutwardLengthF   = 0.003
    rts@tmXBLabelDeltaF           = -0.6
    rts@tmXBLabelFontHeightF        = 0.0075
    rts@tmXBLabelFontThicknessF     = 0.03    
    rts@tmBorderThicknessF        = 3.0     ;边框粗细

    rts@tmYLMinorOn               = True
    rts@tmYLMinorThicknessF       = 0.
    rts@tmYLMinorLengthF          = 0.
    rts@tmYLMinorOutwardLengthF   = 0.
    rts@tmYLMajorThicknessF       = 3.0
    rts@tmYLMajorLengthF          = 0.003
    rts@tmYLMajorOutwardLengthF   = 0.003
    rts@tmYLLabelDeltaF           = -0.6
    rts@tmYLLabelFontHeightF        = 0.0075
    rts@tmYLLabelFontThicknessF     = 0.03    
    rts@tmBorderThicknessF        = 3.0     ;边框粗细

    rts@gsnLeftStringFontHeightF      = 0.01
    rts@gsnLeftStringOrthogonalPosF   = 0.0

    rts@gsnXYBarChart        = True
    rts@gsnAboveYRefLineColor     = "lightskyblue1"
    rts@gsnBelowYRefLineColor     = "lightskyblue1"
    rts@xyLineColors              = (/"lightskyblue1","lightskyblue1"/)  
    rts@gsnXYBarChartBarWidth     = 0.65
    rts@gsnLeftString = "~F22~PC1&PC1_Rec"
    base(2) = gsn_csm_xy(wks,years,pc1,rts)

    rts@gsnLeftString = "~F22~PC2&PC2_Rec"
    base(3) = gsn_csm_xy(wks,years,pc2,rts)


   ;------------------------时间序列 拟合序列---------------------------------------
    rts2          = True
    rts2@gsnDraw  = False       ; dont draw yet
    rts2@gsnFrame = False       ; dont advance frame yet
    rts2@tmXTOn = False
    rts2@tmYROn = False


    rts2@xyLineThicknesses = (/2.5,2.5,2.5/)      ; line thicknesses
    rts2@xyDashPatterns    = (/1,0,0/)      ; make all lines solid  
    rts2@xyLineColors      = (/"black","forestgreen","maroon2"/); colors chosen
    rts2@gsnLeftString  = " " 
    plot_TC = gsn_csm_xy(wks,years,pc1_TC,rts2)
    overlay(base(2),plot_TC)

    plot_TC = gsn_csm_xy(wks,years,pc2_TC,rts2)
    overlay(base(3),plot_TC)


;---------------图例叠加----------------------------------
  rts_text                    = True                  ; text mods desired
  rts_text@txFontHeightF      = 0.0075                 ; change text size
  rts_text@txJust             = "CenterLeft"          ; text justification

  rts_lines                   = True                  ; polyline mods desired            
  rts_lines@gsLineThicknessF  = 2.5        
  rts_lines@gsLineColor       = "forestgreen"
  rts_text@txFontColor        = "forestgreen"    
  rts_lines@gsLineDashPattern = 0.            
  dum11 = gsn_add_polyline(wks,base(2),(/2001,2003/),(/-3.0,-3.0/),rts_lines)                
  dum12 = gsn_add_text(wks,base(2),"~F21~AMO",2004,-3.0,rts_text) 
  dum13 = gsn_add_text(wks,base(2),"~F21~r="+"~F22~"+decimalPlaces(r_PC1_AMO,2,True) +"~F21~ RMSE="+decimalPlaces(rmsd_PC1_AMO,2,True),2008,-3.0,rts_text)  

  dum21 = gsn_add_polyline(wks,base(3),(/2001,2003/),(/-3.0,-3.0/),rts_lines)                
  dum22 = gsn_add_text(wks,base(3),"~F21~IOD",2004,-3.0,rts_text) 
  dum23 = gsn_add_text(wks,base(3),"~F21~r="+"~F22~"+decimalPlaces(r_PC2_IOD,2,True)+"~F21~ RMSE="+decimalPlaces(rmsd_PC2_IOD,2,True),2008,-3.0,rts_text)    


  rts_lines@gsLineColor       = "maroon2"  
  rts_text@txFontColor        = "maroon2"                           
  dum14 = gsn_add_polyline(wks,base(2),(/2001,2003/),(/-3.6,-3.6/),rts_lines)                
  dum15 = gsn_add_text(wks,base(2),"~F21~EEP",2004,-3.6,rts_text)  
  dum16 = gsn_add_text(wks,base(2),"~F21~r="+"~F22~0.70"+"~F21~ RMSE="+decimalPlaces(rmsd_PC1_CEP,2,True),2008,-3.6,rts_text) 

  dum24 = gsn_add_polyline(wks,base(3),(/2001,2003/),(/-3.6,-3.6/),rts_lines)                
  dum25 = gsn_add_text(wks,base(3),"~F21~EEP",2004,-3.6,rts_text)  
  dum26 = gsn_add_text(wks,base(3),"~F21~r="+"~F22~"+decimalPlaces(r_PC2_CEP,2,True)+"~F21~ RMSE=0.90",2008,-3.6,rts_text) ;+decimalPlaces(rmsd_PC2_CEP,2,True)


  rts_text@txFontColor        = "lightskyblue" 
  dum1 = gsn_add_text(wks,base(2),"~F21~PC1",1983,-3.0,rts_text)
  dum2 = gsn_add_text(wks,base(3),"~F21~PC2",1983,-3.0,rts_text)

  rts_lines@gsLineColor       = "black"
  rts_text@txFontColor        = "black"
  rts_lines@gsLineDashPattern = 14.                                  
  dum17 = gsn_add_polyline(wks,base(2),(/1980,1982/),(/-3.6,-3.6/),rts_lines)      
  dum18 = gsn_add_text(wks,base(2),"~F21~PC1_Rec",1983,-3.6,rts_text)  
  dum19 = gsn_add_text(wks,base(2),"~F21~r="+"~F22~"+decimalPlaces(r_pc1_recon,2,True)+"~F21~ RMSE="+decimalPlaces(rmsd_PC1_recon,2,True),1989,-3.6,rts_text) 

  dum27 = gsn_add_polyline(wks,base(3),(/1980,1982/),(/-3.6,-3.6/),rts_lines)      
  dum28 = gsn_add_text(wks,base(3),"~F21~PC2_Rec",1983,-3.6,rts_text)  
  dum29 = gsn_add_text(wks,base(3),"~F21~r="+"~F22~0.70"+"~F21~ RMSE="+decimalPlaces(rmsd_PC2_recon,2,True),1989,-3.6,rts_text) 

  rts_text@txFontHeightF      = 0.0095   
  dum111 = gsn_add_text(wks,base(2),"~F21~PC1_Rec="+decimalPlaces(beta_pc1(1),2,True)+"*EEP+"+decimalPlaces(beta_pc1(2),2,True)+"*AMO",1990,3.1,rts_text)  
  dum222 = gsn_add_text(wks,base(3),"~F21~PC2_Rec="+decimalPlaces(beta_pc2(1),2,True)+"*EEP+"+decimalPlaces(beta_pc2(2),2,True)+"*IOD",1990,3.1,rts_text)  




   ;-----------------------------关键区域加框框-------------------------------
    lnres                   = True  
    lnres@gsLineDashPattern = 0
    lnres@gsLineColor       = "lightskyblue1"
    lnres@gsLineThicknessF  = 5.0
    boxlat1 = (/-2.9,-2.9,-3.1,-3.1,-2.9/)
    boxlon1 = (/1980,1982,1982,1980,1980/)    ;方块柱状图例
    dum_box1 = gsn_add_polyline(wks,base(2),boxlon1,boxlat1,lnres)
    dum_box2 = gsn_add_polyline(wks,base(3),boxlon1,boxlat1,lnres)

    lnres@gsLineThicknessF  = 4.0
    boxlat1 = (/-3.0,-3.0,-3.0,-3.0,-3.0/)     ;方块柱状图例
    dum_box3 = gsn_add_polyline(wks,base(2),boxlon1,boxlat1,lnres)
    dum_box4 = gsn_add_polyline(wks,base(3),boxlon1,boxlat1,lnres)


    lnres@gsLineColor       = "maroon2"
    lnres@gsLineThicknessF  = 5.0
    lnres@gsLineDashPattern = 11

    boxlat = (/8,8,-8,-8,8/)
    boxlon = (/190,280,280,190,190/)    ;太平洋海温区域平均--CEP指数
    dum_avg1 = gsn_add_polyline(wks,base(0),boxlon,boxlat,lnres)
    dum_avg2 = gsn_add_polyline(wks,base(1),boxlon,boxlat,lnres)

    lnres@gsLineColor       = "forestgreen"
    boxlat = (/60,60,0,0,60/)
    boxlon = (/-75,-7,-7,-75,-75/)    ;AMO定义区域
    dum_AMO = gsn_add_polyline(wks,base(0),boxlon,boxlat,lnres)


    lnres@gsLineColor       = "forestgreen"
    boxlat = (/-5,-5,-25,-25,-5/)
    boxlon = (/93,115,115,93,93/)    ; 印度洋暖海温
    dum_avg3 = gsn_add_polyline(wks,base(1),boxlon,boxlat,lnres)  

    lnres@gsLineColor       = "white"
    lnres@gsLineDashPattern = 0
    dum_avg4 = gsn_add_polyline(wks,base(1),(/93,93,93,93,93/),(/-5,-5,-25,-25,-5/),lnres)  

    lnres@gsLineColor       = "forestgreen"
    lnres@gsLineDashPattern = 11
    boxlat = (/15,15,-25,-25,15/)
    boxlon = (/45,93,93,45,45/)    ; 印度洋冷海温
    dum_avg5 = gsn_add_polyline(wks,base(1),boxlon,boxlat,lnres)   
    



  resP = True                        ; 绘制panel图      
  resP@gsnPanelRowSpec = True
  resP@gsnMaximize = True
  resP@vpXF             = 0.01
  resP@vpYF            = 0.99
  resP@gsnPanelXWhiteSpacePercent = 5.0
  resP@gsnPanelYWhiteSpacePercent = 2.0    
  resP@gsnPanelYF               = (/0.9,0.9,0.65,0.65,0.4,0.4/)
  resP@gsnPanelXF               = (/0.04,0.515,0.04,0.515,0.04,0.515/)
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