load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
;;;;;;;;;;    ncl /cygdrive/d/xm/keti/code/Q1_draw.ncl  ;;;;;;;;;;;;;;;;;;;
begin
    latS      = -30
    latN      =  90
    lonL      = 0
    lonR      = 360

    yrStrt = 1979
    yrLast = 2021
    season = "JJA"
    years  = ispan(yrStrt,yrLast,1)
    nyears = dimsizes(years)
    ;print(nyears)
;**************************************读取数据 ****************************************************
     ;--------------------------pc1--------------------------------------
     soizs_1  = asciiread("/cygdrive/d/xm/keti/shuju/pc1.txt", -1, "float")
     ;printVarSummary(soizs_1)
     ;print(soizs_1) 

     ;--------------------------pc2--------------------------------------
     soizs_2 = asciiread("/cygdrive/d/xm/keti/shuju/pc2.txt", -1, "float")
     ;printVarSummary(soizs_2)
     ;print(soizs_2) 

     ;-----------------读入Q1数据----------------------------
     q1_in   = addfile("/cygdrive/d/xm\keti\shuju\q1_1979-2021.nc", "r")
     ;print(q1_in)
     time    = q1_in->time
     XXXX    = cd_calendar(time,-1 )
     t_start = ind(XXXX.eq.197901)
     t_end   = ind(XXXX.eq.202112)
     q1      = q1_in->q
     ;printVarSummary(q1)     
     q1_flt  = tofloat(q1)
     copy_VarCoords(q1,q1_flt)
     ;printVarSummary(q1_flt)
     q1_JJA  =  month_to_season(q1_flt, "JJA")     
     printVarSummary(q1_JJA)

     ;----------------Q1&pc1回归--------------------------------- 
     q1_reg_pc1  = regCoef_n(soizs_1,q1_JJA(level|:,lat|:,lon|:,time|:),0,3)      
     copy_VarCoords(q1_JJA(0,:,:,:),q1_reg_pc1)
     ;q1_reg_pc1  := lonFlip(q1_reg_pc1)
     printVarSummary(q1_reg_pc1)

    ;----------------Q1&pc2回归--------------------------------- 
     q1_reg_pc2  = regCoef_n(soizs_2,q1_JJA(level|:,lat|:,lon|:,time|:),0,3)      
     copy_VarCoords(q1_JJA(0,:,:,:),q1_reg_pc2) 
     ;q1_reg_pc2  := lonFlip(q1_reg_pc2)
     printVarSummary(q1_reg_pc2) 

     ;----------------pc1时北非和太平洋区域-----------------------
     q1_Pacific_pc1      = dim_avg_n_Wrap(q1_reg_pc1(:,{-5:5},{200:275}),(/1,2/))
     ;print(q1_Pacific_pc1)    
     q1_north_africa_pc1 = dim_avg_n_Wrap(q1_reg_pc1(:,{5:15},{0:40}),(/1,2/))
     ;print(q1_north_africa_pc1)

    ;----------------pc2时大西洋和太平洋区域-----------------------
     q1_Pacific_pc2      = dim_avg_n_Wrap(q1_reg_pc2(:,{-5:5},{200:275}),(/1,2/))
     ;print(q1_Pacific_pc2) 
     q1_Atlantic_pc2     = dim_avg_n_Wrap(q1_reg_pc2(:,{-5:5},{300:340}),(/1,2/))
     ;print(q1_Atlantic_pc2)

     data_q1_pc1      = new((/2,12/),float)
     data_q1_pc1(0,:) = q1_Pacific_pc1
     data_q1_pc1(1,:) = q1_north_africa_pc1
     ;print(data_q1_pc1)

     data_q1_pc2      = new((/2,12/),float)
     data_q1_pc2(0,:) = q1_Pacific_pc2
     data_q1_pc2(1,:) = q1_Atlantic_pc2
     ;print(data_q1_pc2)


    ;----------------------------------画图-Q1垂直多层画图----------------------------------------
     wks =gsn_open_wks("eps","Q1") 
    base = new(3,"graphic") 

   res                       = True
   
   res@gsnDraw               = False
   res@gsnFrame              = False   
   res@gsnAddCyclic          = True

   res@cnFillOn              = True             ; turn on color
   res@cnLinesOn             = False             ; turn off contour lines
   res@cnLineLabelsOn        = False            ; turn off contour linelabels
   res@cnLevelSelectionMode  = "ManualLevels"   ; set manual contour levels
   res@cnMaxLevelValF        = 1.5
   res@cnMinLevelValF        = -1.5
   res@cnLevelSpacingF       = 0.1
   res@pmTickMarkDisplayMode = "Always"

   res@mpOutlineBoundarySets       = "Geophysical" ;
   res@mpGeophysicalLineThicknessF = 1.5         ; double the thickness of geophysicalboundaries
   res@mpGeophysicalLineColor     = "gray"  
   res@mpNationalLineThicknessF    = 2.       ; double the thickness of national boundaries
   res@mpOutlineOn                 = True       ; turn off outline
   res@mpFillOn                    = False    ; turn off map fill   
   res@mpMinLatF                   = -30.
   res@mpMaxLatF                   = 90.
   res@mpMinLonF                   = 0
   res@mpMaxLonF                   = 360.
   res@mpCenterLonF                = 110.

   res@tmXBLabelFontHeightF        = 0.005
   res@tmXBLabelDeltaF             = -0.6 

   res@vpHeightF             = 0.15  
   res@vpWidthF              = 0.45
   res@vpXF                  = 0
   res@vpYF                  = 1   

   res@gsnLeftString          =   "200hPa |pc1"
   res@gsnRightString         =   " "
   res@gsnLeftStringOrthogonalPosF = -0.01
    base(0) = gsn_csm_contour_map(wks,q1_reg_pc1(9,:,:),res)

   res@gsnLeftString          =   "500hPa |pc1"
   res@gsnRightString         =   " "
   res@gsnLeftStringOrthogonalPosF = -0.01
    base(1) = gsn_csm_contour_map(wks,q1_reg_pc1(5,:,:),res) 
   
   res@gsnLeftString          =   "850hPa |pc1"
   res@gsnRightString         =   " "
   res@gsnLeftStringOrthogonalPosF = -0.01
   base(2) = gsn_csm_contour_map(wks,q1_reg_pc1(2,:,:),res) 

        ;-----------------------------关键区域加框框-------------------------------
    boxlat = (/5,5,-5,-5,5/)
    boxlon = (/200,275,275,200,200/) 
    lnres                   = True  
    lnres@gsLineColor       = "red"
    lnres@gsLineThicknessF  = 4.0
    lnres@gsLineDashPattern = 1
    dum1 = gsn_add_polyline(wks,base(0),boxlon,boxlat,lnres)   
    dum2 = gsn_add_polyline(wks,base(1),boxlon,boxlat,lnres)   
    dum3 = gsn_add_polyline(wks,base(2),boxlon,boxlat,lnres)   

    boxlat1 = (/15,15,5,5,15/)                  
    boxlon1 = (/0,40,40,0,0/) 
    lnres1                   = True  
    lnres1@gsLineColor       = "blue"
    lnres1@gsLineThicknessF  = 4.0
    lnres1@gsLineDashPattern = 11
    dum6 = gsn_add_polyline(wks,base(0),boxlon1,boxlat1,lnres1) 
    dum7 = gsn_add_polyline(wks,base(1),boxlon1,boxlat1,lnres1) 
    dum8 = gsn_add_polyline(wks,base(2),boxlon1,boxlat1,lnres1) 

    resP = True                        ; 绘制panel图      
   resP@gsnPanelRowSpec = True
   resP@gsnMaximize     = True
   gsn_panel(wks,base,(/1,1,1/),resP)      
   frame(wks)


  ;---------------------------------设定区域平均垂直加热廓线--------------------------------
    wks1 = gsn_open_wks("pdf","Q1_real") 
    plot = new(2,"graphic") 

    res1 = True
    res1@vpHeightF             = 0.6 
    res1@vpWidthF              = 0.45
    ;res1@vpXF                  = 0.25
    ;res1@vpYF                  = 0.9  

    res1@tmXTOn                = False
    res1@tmYROn                = False

    res1@xyLineColors           = (/"red","blue"/) ; colors chosen
    res1@xyLineThicknesses      = (/4.0,4.0/)      ; line thicknesses
    res1@xyDashPatterns         = (/0,0/)      ; make all lines solid

    res1@trXMinF                = -0.6 
    res1@trXMaxF                = 0.6
    res1@tmXBMode               = "Explicit"
    res1@tmXBValues             = (/-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6/)
    res1@tmXBLabels             = (/"-0.6","-0.4","-0.2","0.","0.2","0.4","0.6"/)
    res1@gsnXRefLine            = 0
    res1@gsnXRefLineDashPattern = 0
    res1@gsnXRefLineThicknessF  = 4.0
    res1@gsnXRefLineColors   = "gray"

    res1@trYMinF                = 4     
    res1@trYReverse             = True
    res1@tmYLMode               = "Explicit"
    res1@tmYLLabels             = (/"1000","800","600","400","200"/)
    res1@tmYLLabelDeltaF        = -0.6
    res1@tiYAxisString          = "Pres1sure/hPa"

    res1@gsnLeftString          = "(a)Q1 with pc1"
    plot(0) = gsn_csm_xy(wks1,data_q1_pc1,q1&level,res1) 

    res1@gsnLeftString          = "(b)Q1 with pc2"
    plot(1) = gsn_csm_xy(wks1,data_q1_pc2,q1&level,res1)  

  ;----------------手动添加图例-----------------------------------------------
  res_text                    = True                  ; text mods desired
  res_text@txFontHeightF      = 0.01                 ; change text size
  res_text@txJust             = "CenterLeft"          ; text justification

  res_lines                   = True                  ; polyline mods desired
  res_lines@gsLineThicknessF  = 4. 
  ;-------------q1--PC1图例-----------------------------
  res_lines@gsLineDashPattern = 0. 
  res_lines@gsLineColor       = "red"                 ; line color
  dum11 = gsn_add_polyline(wks1,plot(0),(/-0.58,-0.47/),(/815,815/),res_lines)             
  dum22 = gsn_add_text(wks1,plot(0),"Pacific",-0.46,815,res_text)
  
  res_lines@gsLineDashPattern = 1.                             
  dum33 = gsn_add_polyline(wks1,plot(0),(/-0.58,-0.47/),(/865,865/),res_lines)                
  dum44 = gsn_add_text(wks1,plot(0),"Pacific_M",-0.46,865,res_text)   

  res_lines@gsLineDashPattern = 0. 
  res_lines@gsLineColor       = "blue"                 ; line color
  dum55 = gsn_add_polyline(wks1,plot(0),(/-0.58,-0.47/),(/915,915/),res_lines)             
  dum66 = gsn_add_text(wks1,plot(0),"North africa",-0.46,915,res_text)

  res_lines@gsLineDashPattern = 1. 
  dum77 = gsn_add_polyline(wks1,plot(0),(/-0.58,-0.47/),(/965,965/),res_lines)                
  dum88 = gsn_add_text(wks1,plot(0),"North africa_M",-0.46,965,res_text)  

   ;-------------q1--PC2图例-----------------------------
  res_lines@gsLineDashPattern = 0. 
  res_lines@gsLineColor       = "red"                 ; line color
  dum111 = gsn_add_polyline(wks1,plot(1),(/-0.58,-0.47/),(/815,815/),res_lines)             
  dum222 = gsn_add_text(wks1,plot(1),"Pacific",-0.46,815,res_text)
  
  res_lines@gsLineDashPattern = 1.                             
  dum333 = gsn_add_polyline(wks1,plot(1),(/-0.58,-0.47/),(/865,865/),res_lines)                
  dum444 = gsn_add_text(wks1,plot(1),"Pacific_M",-0.46,865,res_text)   

  res_lines@gsLineDashPattern = 0. 
  res_lines@gsLineColor       = "blue"                 ; line color
  dum555 = gsn_add_polyline(wks1,plot(1),(/-0.58,-0.47/),(/915,915/),res_lines)             
  dum666 = gsn_add_text(wks1,plot(1),"Atlantic",-0.46,915,res_text)

  res_lines@gsLineDashPattern = 1. 
  dum777 = gsn_add_polyline(wks1,plot(1),(/-0.58,-0.47/),(/965,965/),res_lines)                
  dum888 = gsn_add_text(wks1,plot(1),"Atlantic_M",-0.46,965,res_text)    

  ;----------------手动设置画虚线-------------------------------------------------
  res_m_line                   = True                  
  res_m_line@gsLineThicknessF  = 4. 
  res_m_line@gsLineDashPattern = 1. 
  res_m_line@gsLineColor       = "red"            
  line_m_1 = gsn_add_polyline(wks1,plot(0),(/-0.25,-0.35/),(/850,700/),res_m_line)                      
  line_m_2 = gsn_add_polyline(wks1,plot(0),(/-0.35,-0.4/),(/700,500/),res_m_line) 
  line_m_3 = gsn_add_polyline(wks1,plot(0),(/-0.4,-0.35/),(/500,300/),res_m_line) 
  line_m_4 = gsn_add_polyline(wks1,plot(0),(/-0.35,-0.25/),(/300,200/),res_m_line)

  line_m_11 = gsn_add_polyline(wks1,plot(1),(/-0.3,-0.45/),(/850,700/),res_m_line)                      
  line_m_22 = gsn_add_polyline(wks1,plot(1),(/-0.45,-0.4/),(/700,500/),res_m_line) 
  line_m_33 = gsn_add_polyline(wks1,plot(1),(/-0.4,-0.35/),(/500,300/),res_m_line) 
  line_m_44 = gsn_add_polyline(wks1,plot(1),(/-0.35,-0.3/),(/300,200/),res_m_line)

  res_m_line@gsLineColor       = "blue"            
  line_m_5 = gsn_add_polyline(wks1,plot(0),(/0.35,0.5/),(/850,700/),res_m_line)                      
  line_m_6 = gsn_add_polyline(wks1,plot(0),(/0.5,0.55/),(/700,500/),res_m_line) 
  line_m_7 = gsn_add_polyline(wks1,plot(0),(/0.55,0.5/),(/500,300/),res_m_line) 
  line_m_8 = gsn_add_polyline(wks1,plot(0),(/0.5,0.4/),(/300,200/),res_m_line) 

  line_m_55 = gsn_add_polyline(wks1,plot(1),(/0.25,0.4/),(/850,700/),res_m_line)                      
  line_m_66 = gsn_add_polyline(wks1,plot(1),(/0.4,0.4/),(/700,500/),res_m_line) 
  line_m_77 = gsn_add_polyline(wks1,plot(1),(/0.4,0.3/),(/500,300/),res_m_line) 
  line_m_88 = gsn_add_polyline(wks1,plot(1),(/0.3,0.2/),(/300,200/),res_m_line)  
  ;---------------------------panel组图--------------------------------
   resP1 = True                        ; 绘制panel图      
   resP1@gsnPanelRowSpec = True
   resP1@gsnMaximize     = True
   gsn_panel(wks1,plot,(/1,1/),resP1)      
   frame(wks1)
   end



  