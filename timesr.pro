;+
; NAME:
;       HANTS
;
; PURPOSE:
;         This programme is based on the 2005 fortran version of HANTS programme by Wout Verhoef.
;
; INPUT:
;   ni:   Number of images = length of time series, such as 46 for MODIS 8day LST product
;   per:  Period that input time series covered, such as 365 for one year
;   ts:   Day number of each element in array y, such as [1,9,17,...]
;   nf:   Number of frequency, such as 3 (generally 3 or 4. Larger nf is, the reconstructed curve turns to include more noise)
;   pa:   Period, which is corresponding to nf. If nf is 3, pa = [361 180 120]
;   HiLo: High or low outliers. (To decide higher values or lower values will be treated as unreasonable values)
;   low:  Minimum of valid range
;   high: Maximum of valid range
;   fet:  Fit error tolerance
;   dod:  Degree of overdetermined
;   y:    Input time series curve.
;   delta = small positive number to suppress high amplitudes
;
;OUTPUT:
;   yr:   Reconstructed time series
;   amp:  Array of amplitudes, first element is the average of the curve
;   phi:  Array of phases, first element is zero
;
;AUTHOR:
;         Jing Li, Ph.D
;         E-mail: li_jing_chn@hotmail.com
; TIME:
;   September, 2008

FUNCTION HANTS,ni,per,ts,nf,pa,HiLo,low,high,fet,dod,y,looplim,delta,$
  POLYNOMIAL=polynomial
  ON_ERROR,0
  ;catch, theError
  ;if theError ne 0 then begin
  ; catch, /cancel
  ; ok = dialog_message(!error_state.msg, /error)
  ; return, {yr:yr, amp:amp, phi:phi}
  ;endif
  
  ;delta = 0.1
  
  ;output
  amp = FLTARR(nf+1)
  phi = FLTARR(nf+1)
  yr = UINTARR(ni)
  
  tagdiff = STRCMP(HiLo,'High',/FOLD_CASE)
  IF tagdiff EQ 1 THEN BEGIN
    sHiLo = -1
  ENDIF ELSE BEGIN
    sHiLo = 1
  ENDELSE
  
  nr = 2*nf+1
  
  IF N_ELEMENTS(polynomial) NE 0 THEN nr=nr+3
  
  IF nr GT ni THEN BEGIN
    nr = ni
  ENDIF
  noutmax = ni-nr-dod
  
  pi = 4*ATAN(1)
  rd = pi/180
  dg = 180/pi
  
  mat = FLTARR(nr,ni)+1
  ang = FINDGEN(per)
  ang = 2*pi*ang/per
  cs = COS(ang)
  sn = SIN(ang)
  
  ;"mat": Dimensions is ni*(2*nf+1)
  FOR i=0,nf-1 DO BEGIN
    ifr = FIX(FLOAT(per)/pa[i])
    index = ifr*(ts-1) MOD per
    ;angles(i,*) = index
    mat[2*i+1,*] = cs[index]
    mat[2*i+2,*] = sn[index]
  ENDFOR
  
  IF N_ELEMENTS(polynomial) NE 0 THEN BEGIN
    ;ifr = fix(float(per)/pa[i])
    index = cgSCALEVECTOR(ts,0,per-1)
    mat[nr-2,*]=index*0.05;/(2*pi)
    mat[nr-1,*]=LONG64(index*0.05)^2
    mat[nr-1,*]=LONG64(index*0.05)^3
  ENDIF
  ;'p': weight value!
  p = INTARR(1,ni)+1
  pos = WHERE(y GT high OR y LT low,nout)
  IF nout NE 0 THEN BEGIN
    p[pos] = 0
  ENDIF
  
  IF nout GT noutmax THEN BEGIN
    ;result = dialog_message('Samples in series are lower than needed! Reconstruction failed!',/error)
    RETURN, 0
  ENDIF
  
  nloop = 0
  ;nloopmax = 10
  ready = 0
  mat = TRANSPOSE(mat)
  
  ;cgoplot,y
  
  ;����ȥ�����������
  ;while ready eq 0 and nloop lt nloopmax and nloop lt looplim do begin
  WHILE ready EQ 0 AND nloop LT looplim DO BEGIN
    nloop = nloop + 1
    
    za = mat ## (p*y)
    pmatrix = p
    FOR ii = 1,nr-1 DO BEGIN
      pmatrix = [[pmatrix],[p]]
    ENDFOR
    A = mat*pmatrix##TRANSPOSE(mat)
    ADim = SIZE(A,/dimensions)
    Adelta = FLTARR(ADim)
    Adelta[INDGEN(ADim[0]) * (ADim[0]+1)] = 1
    Adelta = Adelta*delta
    Adelta[0,0] = 0
    A = A + Adelta
    B = INVERT(A)
    zr = B##za
    
    
    ;rewrite by wwj
    ;yr = transpose(zr)##mat
    yr = TRANSPOSE(mat)##zr
    
    diff = sHiLo*(yr-y)
    err = TRANSPOSE(p)*diff
    
    errsort_sub = SORT(err)
    maxerr=diff[errsort_sub[ni-1]]
    ready = (maxerr LE fet) OR (nout EQ noutmax)
    
    IF ready EQ 0 THEN BEGIN
      ;      i=ni-1
      ;      j=errsort_sub[i]
      ;      while p[j]*diff[j] gt maxerr*0.5 and nout lt noutmax do begin
      ;
      ;        p(j)=0
      ;        nout=nout+1
      ;        i=i-1
      ;        j=errsort_sub(i)
      ;      endwhile
      index=WHERE(err GT fet, count,COMPLEMENT=c_index)
      p[index]=p[index]-err[index]
      p_index=WHERE(p LT 0,count)
      IF count NE 0 THEN BEGIN
        p[p_index]=0
        nout=count
      ENDIF
      
      
    ENDIF
    ;index=WHERE(p LT 1, count,COMPLEMENT=c_index)
    ;y[index]=INTERPOL(y[c_index],c_index,index,/spline)
    ;  yr[index]=y[index]
    
    ;cgoplot,y
  ENDWHILE
  
  ra = zr[INDGEN(nf)*2+1]
  rb = zr[INDGEN(nf)*2+2]
  amp[1:*] = TRANSPOSE(SQRT(ra*ra+rb*rb))
  phi[1:*] = TRANSPOSE(ATAN(rb,ra)*dg)
  amp[0] = zr[0]
  phi[0] = 0
  pos = WHERE(phi LT 0, npos)
  IF npos NE 0 THEN BEGIN
    phi[pos] = phi[pos] + 360
  ENDIF
  
  RETURN, {yr:yr, amp:amp, phi:phi,valid_i:WHERE(p EQ 1)}
  
END
FUNCTION DIFF,A,Deg
  ;A=shift(D,0,1)-A
  IF N_PARAMS() EQ 1 THEN deg=1
  
  IF deg EQ 1 THEN BEGIN
    RETURN,(SHIFT(A,0,1)-A)[*,1:*]
  END
  RETURN,DIFF((SHIFT(A,0,1)-A)[*,1:*],deg-1)
END
FUNCTION WHITSMW,y, WEIGHT=w, LAMBDA=lambda, DEGREE=degree

  catch,err
  if err ne 0 then begin
    catch,/cancel
    print,y,w
    return,w*y
  endif
  IF N_ELEMENTS(degree) EQ 0 THEN degree=2
  IF N_ELEMENTS(w) EQ 0 THEN w=FLTARR(N_ELEMENTS(y))+1
  IF N_ELEMENTS(lambda) EQ 0 THEN lambda=0.5
  ;Smoothing
  m = n_elements(y);
  ;E = SPRSIN(LINDGEN(m), LINDGEN(m), REPLICATE(1.0,m), m)
  ;W = SPRSIN(LINDGEN(m), LINDGEN(m), w, m);
  wa=DIAG_MATRIX(w)
  D = IDENTITY(m)
  D=DIFF(D,degree)
  ;C = chol(W + lambda * transpose(D) * D);
  ;z = INVERT((wa + lambda * TRANSPOSE(D)## D))# (wa#y)
  A=(wa + lambda * TRANSPOSE(D)## D)
  CHOLDC, A, P
  
  
  RETURN,CHOLSOL(A, P, wa#y)
END

FUNCTION ddmat,x,d
  m=N_ELEMENTS(x)
  IF d EQ 0 THEN BEGIN
    DD=identity(m)
    ;DD=SPRSIN(indgen(m),indgen(m),REPLICATE(1.0,m),m)
  ENDIF ELSE BEGIN
    dx=x[d:m-1]-x[0:m-d-1]
    V=DIAG_MATRIX(1.0/dx)
    ;V=SPRSIN(indgen(m-d),indgen(m-d),1.0/dx,m-d)
    DD=V##diff(ddmat(x,d-1))
    ; DD=SPRSAB(V##(diff(ddmat(x,d-1))))
  ENDELSE
  RETURN,DD
END


FUNCTION whitsmddw,x, y, WEIGHT=weight, LAMBDA=lambda, DEGREE=d
  IF N_ELEMENTS(degree) EQ 0 THEN degree=2
  IF N_ELEMENTS(weight) EQ 0 THEN weight=FLTARR(N_ELEMENTS(y))+1
  IF N_ELEMENTS(lambda) EQ 0 THEN lambda=10
  m=N_ELEMENTS(y)
  e=identity(m)
  DD=ddmat(x,d)
  w=DIAG_MATRIX(weight)
  z = INVERT((W + lambda * TRANSPOSE(DD)## DD))# (w#y)
  RETURN,z
END
FUNCTION R_ha,o_s,WEIGHT=weight, PARA_INFO=para_info,$
  PER=per,TS=ts,DELTA=delta,PARR=pa,OTS=ots,OUT=out,_Extra=extrakeywords
  
  ON_ERROR,0
  ;  if n_elements(batch) ne 0 then begin
  ;     CASE (batch) OF
  ;       0: BEGIN ;all params hold the same except o_s
  ;
  ;       END
  ;       ELSE: BEGIN
  ;       END
  ;     ENDCASE
  ;
  ;  endif
  ni=N_ELEMENTS(o_s)
  
  
  IF N_ELEMENTS(delta) EQ 0 THEN delta=0.5
  IF N_ELEMENTS(weight) EQ 0 THEN weight=FLTARR(ni)+1
  
  pi = 4*ATAN(1)
  nf=N_ELEMENTS(pa)
  nr=2*nf+1
  
  mat = FLTARR(nr,ni)+1
  pha_arr=2*!PI*per/pa#ts/per
  mat[1:nr-1:2,*]=COS(pha_arr)
  mat[2:nr-1:2,*]=SIN(pha_arr)
  
  pos = WHERE(weight NE 0)
  ; p[pos] = 1
  
  w=BYTARR(ni,ni)
  w[pos,pos]=1
  
  A = TRANSPOSE(w##mat)##(mat)
  ADim = SIZE(A,/dimensions)
  Adelta = FLTARR(ADim)
  Adelta[INDGEN(ADim[0]) * (ADim[0]+1)] = 1
  Adelta = Adelta*delta
  Adelta[0,0] = 0
  A = A + Adelta
  B = INVERT(A)
  yr = mat##B##TRANSPOSE(w##mat)##o_s
  
;  M=w##mat
;  LA_SVD, M, S, U, V
;  
;  N = N_ELEMENTS(S)
;  WP = FLTARR(N, N)
;  index=where(abs(S) ge 0.1,count) 
;  ;if count ne -1 then wp[index,index]=1.0/s[index]
;  wp[indgen(N),indgen(N)]=s/(s+0.1)^2
  
;t=reform(w##o_s)
;if count ne -1 then s[index]=0
;  yr1=mat##SVSOL(U, S, V, t)   

;   A=TRANSPOSE(w##mat)##(mat)
;   CHOLDC,A,P
;   B=reform(transpose(mat)##w##o_s)
;   yr1=mat##CHOLSOL(A, P, B)
  
  
  IF N_ELEMENTS(out) ne 0  THEN BEGIN
    no=N_ELEMENTS(ots)
    mat_o = FLTARR(nr,no)+1
    pha_arr=2*!PI*per/pa#ots/per
    mat_o[1:nr-1:2,*]=COS(pha_arr)
    mat_o[2:nr-1:2,*]=SIN(pha_arr)
    
    yr=mat_o##B##TRANSPOSE(w##mat)##o_s 
    ;yr=  mat_o## V ## WP ## TRANSPOSE(U)##w##o_s
    
  ENDIF else begin
    ;yr=  mat## V ## WP ## TRANSPOSE(U)##w##o_s
    yr=mat##B##TRANSPOSE(w##mat)##o_s
  Endelse
  
  
  IF ARG_PRESENT(para_info) EQ 1 THEN BEGIN
  
    zr = B##TRANSPOSE(w##mat)##o_s
    amp = FLTARR(nf+1)
    phi = FLTARR(nf+1)
    ra = zr[INDGEN(nf)*2+1]
    rb = zr[INDGEN(nf)*2+2]
    amp[1:*] = TRANSPOSE(SQRT(ra*ra+rb*rb))
    phi[1:*] = TRANSPOSE(ATAN(rb,ra))
    amp[0] = zr[0]
    phi[0] = 0
    pos = WHERE(phi LT 0, npos)
    IF npos NE 0 THEN BEGIN
      phi[pos] = phi[pos] + 2*!PI
    ENDIF
    para_info={amp:amp,phi:phi}
  ENDIF
  
  RETURN,REFORM(yr);
END

FUNCTION GENERRIZED_RCON_MULTI,y,$
        TS=ts,OTS=ots,PER=per,DELTA=delta,PARR=pa,$
  HILO=HiLo,LOW=low,HIGH=high,FET=fet,DOD=dod,LOOPLIM=looplim,_Extra=extrakeywords
  ON_ERROR,0


  
  sHiLo = STRCMP(HiLo,'High',/FOLD_CASE) EQ 1?-1:1
  
  ;if the all the value is too small, then we don't do reconstruction.
  if max(y) le low+fet then return,1 

  y=float(y)
  
  ni=N_ELEMENTS(y)
  noutmax = (ni-dod)
  p = fltARR(1,ni)+1
  pos = WHERE(y GT high OR y LT low,nout)
  if nout ge noutmax then return, 0
  IF nout NE 0 THEN p[pos] = 0
  ;tmp_y=y
  
  ready=0
  nloop=0
  
  IF N_ELEMENTS(delta) EQ 0 THEN delta=0.5
  IF N_ELEMENTS(weight) EQ 0 THEN weight=FLTARR(ni)+1

  pi = 4*ATAN(1)
  nf=N_ELEMENTS(pa)
  nr=2*nf+1

  pos = WHERE(weight NE 0)
  ; p[pos] = 1

  w=BYTARR(ni,ni)
  w[pos,pos]=1

  mat = FLTARR(nr,ni)+1
  pha_arr=2*!PI*per/pa#ts/per
  mat[1:nr-1:2,*]=COS(pha_arr)
  mat[2:nr-1:2,*]=SIN(pha_arr)
  WHILE ready EQ 0 AND nloop LT looplim DO BEGIN
    
    
    ;**********************
    A = TRANSPOSE(w##mat)##(mat)
    ADim = SIZE(A,/dimensions)
    Adelta = FLTARR(ADim)
    Adelta[INDGEN(ADim[0]) * (ADim[0]+1)] = 1
    Adelta = Adelta*delta
    Adelta[0,0] = 0
    A = A + Adelta
    B = INVERT(A)
    yr=mat##B##TRANSPOSE(w##mat)##y
    ;************************************
  
    ;yr=WHITSMW(tmp_y,WEIGHT=p,_Extra=extrakeywords)
;    M=w##mat
;
;    LA_SVD, M, S, U, V , STATUS=status
;    if status  gt 0 then return,0
;    RES=SVSOL(U, S, V, reform(w##y))
;
;    N = N_ELEMENTS(S)
;    WP = FLTARR(N, N)
;    index=where(abs(S) ge 0.0000001,count)
;    wp[indgen(n),indgen(n)]=1.0/s
    ;if count ne -1 then wp[index,index]=1.0/s[index]
    ;wp[indgen(N),indgen(N)]=s/(s+0.1)^1
    ;yr=R_ha(y,WEIGHT=w,_Extra=extrakeywords,PARR=pa,PER=per)
   ;yr=  W##mat## V ## WP ## TRANSPOSE(U)##float(y)
;    yr=mat##res
    
    diff = sHiLo*(yr-y)
    err = w[findgen(ni),findgen(ni)]*diff
    

    e_index=WHERE(err GE fet, count)

    nout=nout+count
    ready = (count EQ 0) OR (nout EQ noutmax)
    IF ready EQ 1 THEN BREAK
    
    w[e_index,e_index]=0
    
    
  ENDWHILE
  
IF N_ELEMENTS(OTS) ne 0  THEN BEGIN
    no=N_ELEMENTS(ots)
    mat_o = FLTARR(nr,no)+1
    pha_arr=2*!PI*per/pa#ots/per
    mat_o[1:nr-1:2,*]=COS(pha_arr)
    mat_o[2:nr-1:2,*]=SIN(pha_arr)

    yr=mat_o##B##TRANSPOSE(w##mat)##y
    ;yr=  mat_o## V ## WP ## TRANSPOSE(U)##w##y
;    yr=mat##res

  ENDIF
  
  RETURN,  low>yr<high
END

;for debug use.
FUNCTION GENERRIZED_HANTS,y,$
  ITS=its,OTS=ots,$
  HILO=HiLo,LOW=low,HIGH=high,FET=fet,DOD=dod,LOOPLIM=looplim,$
  PLOT=plot,MULTI=multi,$
  SMOOTHER=smoother,_Extra=extrakeywords
  ON_ERROR,0
  ;catch, theError
  ;if theError ne 0 then begin
  ; catch, /cancel
  ; ok = dialog_message(!error_state.msg, /error)
  ; return, {yr:yr, amp:amp, phi:phi}
  ;endif
  
  ;delta = 0.1
  IF N_ELEMENTS(y) EQ 0 THEN MESSAGE,'please provide the series to be reconstructed.'
  IF N_ELEMENTS(Hilo) EQ 0 THEN Hilo='High'
  IF N_ELEMENTS(low) EQ 0 THEN low =MIN(y)
  IF N_ELEMENTS(high) EQ 0 THEN high =MAX(y)
  IF N_ELEMENTS(dod) EQ 0 THEN dod =N_ELEMENTS(y)/5
  IF N_ELEMENTS(looplim) EQ 0 THEN looplim =10
  
  sHiLo = STRCMP(HiLo,'High',/FOLD_CASE) EQ 1?-1:1
  
  
  ni=N_ELEMENTS(y)
  noutmax = N_ELEMENTS(multi) NE 0 ? 3*(ni-dod):ni-dod
  tmp_y=N_ELEMENTS(multi) NE 0 ?[y,y,y]:y
  tmp_ni=N_ELEMENTS(multi) NE 0 ? 3*ni:ni
  p = INTARR(1,tmp_ni)+1
  pos = WHERE(tmp_y GT high OR tmp_y LT low,nout)
  IF nout NE 0 THEN  p[pos] = 0
  ;tmp_y=y
  
  ready=0
  nloop=0
  WHILE ready EQ 0 AND nloop LT looplim DO BEGIN
    ;yr=CURVEFIT(ts,tmp_y,p
    ;pp=IMSL_CSSMOOTH(ts,y,weights=reform(p))
    ;yr=IMSL_SPVALUE(ts, pp)
    ;yr=SPlINEFIT(ts,Y,reform(p),ts[0:22:2],y[0:22:2],SIGYS,1.0)
    ;   X=ts
    ;          coeff = SPLINECOEFF(ts, y, lambda = 1000,sigm=reform(p))
    ;       yr = FLTARR(N_ELEMENTS(y) - 1)
    ;       x1 = ts[0:N_ELEMENTS(y)-2]+16
    ;       FOR i = 0, N_ELEMENTS(y)-2 DO yr[i] = coeff.d[I] + $
    ;                                             coeff.c[I] * (x[I+1]-x[I]) + $
    ;                                             coeff.b[I] * (x[I+1]-x[I])^2 + $
    ;                                             coeff.a[I] * (x[I+1]-x[I])^3
    
    ;    savgolFilter = SAVGOL(4, 4, 0, 4)
    ;     yr=CONVOL(y, savgolFilter, /EDGE_TRUNCATE)
    
    CASE (smoother) OF
      'HANTS': BEGIN
        yr=r_ha(tmp_y,WEIGHT=p,_Extra=extrakeywords)
      END
      'WHITSMW': BEGIN
        yr=WHITSMW(tmp_y,WEIGHT=p,_Extra=extrakeywords)
      END
      ELSE: BEGIN
      
      END
    ENDCASE
    
    
    
    diff = sHiLo*(yr-tmp_y)
    err = TRANSPOSE(p)*diff
    
    ;errsort_sub = SORT(err)
    e_index=WHERE(err GE fet, count)
    ;maxerr=max(err,max_sub)
    nout=nout+count
    ready = (count EQ 0) OR (nout EQ noutmax)
    IF ready EQ 1 THEN BREAK
    
    p[e_index]=0
    
    IF N_ELEMENTS(plot) NE 0 THEN cgplot,INDGEN(tmp_ni),yr,linestyle=1,/overplot
    IF N_ELEMENTS(plot) NE 0 THEN cgplot,(INDGEN(tmp_ni))[e_index],tmp_y[e_index],psym=21,/overplot
  ENDWHILE
  ;  IF N_ELEMENTS(combine) NE 0 THEN BEGIN
  ;    index=WHERE(p EQ 0, count)
  ;    IF count NE 0 THEN tmp_y[index]=yr[index]
  ;    yr[ni:2*ni-1]=R_ha(tmp_y[ni:2*ni-1],-3000)
  ;  ENDIF ELSE BEGIN
  ;
  ;    p[*]=1
  ;    yr=WHITSMW(yr,REFORM(p),0.5,DEGREE=3)
  ;  ENDELSE
  IF N_ELEMENTS(plot) NE 0 THEN cgplot,INDGEN(tmp_ni),yr,linestyle=0,/overplot,color='blue'
  IF  N_ELEMENTS(multi) NE 0 THEN yr=yr[ni:2*ni-1]
  
  IF N_ELEMENTS(ots) NE 0 THEN BEGIN
  
    mi=N_ELEMENTS(its)
    mo=N_ELEMENTS(ots)
    
    
    values1=REPLICATE({y:0.,w:0.,yr:0.},mo)
    
    
    values2=REPLICATE({y:0.,w:0.,yr:0.},mi)
    values2.y=yr
    values2.w=FLTARR(mi)+1
    obj_hash=HASH(ots,values1,its,values2)
    ;obj_hash[its]=values2
    
    x=(obj_hash.keys()).toarray()
    index=SORT(x)
    values_arr=((obj_hash.values()).toarray())[index]
    values_arr.yr=whitsmddw(x[index], values_arr.y, WEIGHT=values_arr.w, LAMBDA=1, DEGREE=2)
    ;values_arr.yr=r_ha(values_arr.y,WEIGHT=values_arr.w,_Extra=extrakeywords)
    obj_hash[x[index]]=values_arr
    yr_hash=obj_hash[ots]
    yr=(((yr_hash.values())[SORT((yr_hash.keys()).toarray())]).toarray()).yr
    
  ENDIF
  RETURN, yr
END
FUNCTION _conbine_ts,ots,its,I2O=i2o,O2O=o2o
  ni=N_ELEMENTS(its)
  no=N_ELEMENTS(ots)
  
  IF ARG_PRESENT (i2o) THEN BEGIN
  
    values1=INTARR(no)
    
    
    values2=INTARR(ni)+1
    
    obj_hash=HASH(ots,values1,its,values2)
    ;obj_hash[its]=values2
    
    x=(obj_hash.keys()).toarray()
    index=SORT(x)
    values_arr=((obj_hash.values()).toarray())[index]
    
    orts=x[index]
    i2o=WHERE(values_arr EQ 1)
  ENDIF
  
  IF ARG_PRESENT (o2o) THEN BEGIN
    values1=INTARR(ni)
    
    
    values2=INTARR(no)+1
    
    obj_hash=HASH(its,values1,ots,values2)
    ;obj_hash[its]=values2
    
    x=(obj_hash.keys()).toarray()
    index=SORT(x)
    values_arr=((obj_hash.values()).toarray())[index]
    
    orts=x[index]
    o2o=WHERE(values_arr EQ 1)
  ENDIF
  RETURN, orts
END
PRO TSR_SESSION,imgsetfile,ni,start_l_arr,end_l_arr,$
  tmp_packed_file,sample,line,datatype,its,ots,otype,$
  flag_mask,cur_session,n_session,pref
  
  
  on_error,0
  catch,Error_status
  if Error_status ne 0 then begin
    catch,/cancel
    
    PRINT, 'Error index: ', Error_status
    PRINT, 'Error message: ', !ERROR_STATE.MSG
    if n_elements(cur_line) gt 1 then begin
      cur_line=0
      shmunmap,'cur_line'
    endif
    MESSAGE, /REISSUE_LAST    
  endif
  
  HILO='Low';'High';
  LOW=0;11000;0;
  HIGH=10;17500;10000;1000;100;
  FET=0.5;250.0;500.0;50;5;
  DOD=10
  LOOPLIM= 10
  fillvalue=0;-3000;-3000;32767;255;
  
  if n_elements(otype) eq 0 then otype=datatype
  
  print,systime(),'Begin session'
  
  SHMMAP,'cur_line',n_session,/long
  cur_line=SHMVAR('cur_line')
  
  if flag_mask  ne 0 then begin
;    mask=read_tiff(maskfile)
    SHMMAP,'mask_arr',DIMENSION=[sample,line],/byte
    mask_arr=shmvar('mask_arr')
;    help,mask_arr
    ;mask_arr[0]=mask
 ;   flag_mask=1
  endif
  
  ;get the # of sub task
  n_sub_task=n_elements(start_l_arr)
  
  ;get the byte length of the datatype
  nblen=N_TAGS({t:MAKE_ARRAY(1,type=datatype)},/data_length)
  
  ;open the temperal output file
  openu, o_lun,tmp_packed_file,/get_lun
  
  ;constructe algrithm object
  ;intep_obj=obj_new('RECON_HA',TS=its,OTS=ots,/IOD,per=365,parr=[720,360,180,120,90],BATCH=1)
  
  for i_sub_task=0,n_sub_task-1 do begin
  
    start_line=start_l_arr[i_sub_task]
  end_line=end_l_arr[i_sub_task]
  
  ;read image data
  nsubline=end_line-start_line+1
  img_data_arr=MAKE_ARRAY(sample,nsubline,ni,type=datatype)
  img_data=MAKE_ARRAY(sample,nsubline,type=datatype)
 

  
  OPENR,i_lun,imgsetfile,/get_lun
  POINT_LUN,i_lun,start_line*nblen
  FOR i_ni=0,ni-1 DO BEGIN
    ; print,line,i_ni,start_line,sample,nblen
    ; print,(LONG64(line)*i_ni+start_line)*sample*nblen
    POINT_LUN,i_lun,(LONG64(line)*i_ni+start_line)*sample*nblen
    ;print,nblen,(LONG64(line)*i_ni+start_line)*sample*nblen
    READU,i_lun,img_data
    img_data_arr[*,*,i_ni]=img_data
  ENDFOR
  
  
  
  ;create the output array
  no=N_ELEMENTS(ots)
  o_img_data_arr=MAKE_ARRAY(sample,nsubline,no,type=otype)+fillvalue
  flag_arr=bytarr(sample,nsubline)
  
  
  
  ;reconstruction process
  ;PROFILER
  ;ROFILER, /SYSTEM
  print,systime(),'Begin reconstruction',i_sub_task
  
  
  FOR i_line=0,nsubline-1 DO BEGIN
  
    ;print,reform(mask_arr[*,start_line+i_line])
    FOR i_sample=0,sample-1 DO BEGIN
     
       if flag_mask  ne 0 then begin
         if mask_arr[i_sample,start_line+i_line] eq 0 then continue
       endif

      res=GENERRIZED_RCON_MULTI(HILO=HiLo,LOW=low,HIGH=high,FET=fet,DOD=dod,$
        float(REFORM(img_data_arr[i_sample,i_line,*])),LOOPLIM=looplim,$
        per=365,parr=[720,360,180,120],TS=its,OTS=ots)
      ; res=1
        if n_elements(res) eq 1 then begin
           o_img_data_arr[i_sample,i_line,*]= fillvalue
           continue
           ;flag_arr[i_sample,i_line]=1
        endif else begin
	   o_img_data_arr[i_sample,i_line,*]= res

	endelse 
      ;PRINT,i_sample,i_line
      flag_arr[i_sample,i_line]=1
    ENDFOR
    
    ;PRINT,i_line
    cur_line[cur_session]=start_line+i_line
  ENDFOR
  print,systime(),'end reconstruction',i_sub_task
  
  ;print a test series to check the vilidation of the reconstruction
  v_index=where(flag_arr eq 1,count) ;get all index of valid pixeles
  ;get a random index
  if count ne 0 then begin
  rt_index=v_index[0>fix(randomu(SYSTIME(1),1)*count)<(count-1)]
  ;t_sample=rt_index-rt_index/sample
  t_line=rt_index/sample
  t_sample=rt_index-t_line*sample
  ;t_sample=0>fix(randomu(SYSTIME(1),1)*sample)<(sample-1)
  ;t_line=0>fix(randomu(SYSTIME(1),1)*(end_line-start_line))<(end_line-start_line-1)
  iseries=REFORM(img_data_arr[t_sample,t_line,*])
  oseries=REFORM(o_img_data_arr[t_sample,t_line,*])
  
  tfile=string(pref,t_sample,start_line+t_line,format='(A,".S",I05,"L",I05,".txt")')
  openw,t_lun,tfile,/get_lun
  printf,t_lun,'start line: ',start_line
  printf,t_lun,'end line: ', end_line
  printf,t_lun,'series for :', t_sample,start_line+t_line
  printf,t_lun,'input time tag  ','series value'
  printf,t_lun,transpose([[its],[iseries]])
  printf,t_lun,'output time tag  ','series value'
  printf,t_lun,transpose([[ots],[oseries]])
  close,t_lun
  
  ;plot the sample series
  pixw=cgpixmap()
  ptitle=string(t_sample,start_line+t_line,format='("series for ","S",I05,"L",I05)')
  cgplot,its,iseries,color='red',psym=2,YRANGE=[low,high],title=ptitle,/window
  cgplot,ots,oseries,color='dodger blue',/addcmd,/overplot
  cd,cur=cur
  pngfile=filepath(strcompress(string(pref,t_sample,start_line+t_line,$
    randomu(systime(1))*long64(2)^20,format='(A,"-S",I05,"L",I05,"-",I08,".png")'),/remove_all),root_dir=cur)
 print,pngfile
  pixw->output,pngfile
  print,pngfile
  cgdelete,pixw
  endif
    
  img_data_arr=0
  ;cd
  ;cgplot,ts,REFORM(img_data_arr[t_sample,t_line,*]),psym=3,output=string(t_sample,start_line+t_line,format='("S",I,"L",I,".ps")')
  ;cgplot,ots,REFORM(o_img_data_arr[t_sample,t_line,*]),psym=-5
  
  
  ;write the reconstruction result to temparal file
  print,systime(),'Begin write file',i_sub_task
 ;get the byte length of the otype
  noblen=N_TAGS({t:MAKE_ARRAY(1,type=otype)},/data_length)
  for i=0,no-1 do begin
    point_lun,o_lun,sample*(long64(line)*i+start_line)*noblen
    writeu,o_lun,o_img_data_arr[*,*,i]
  endfor
  print,systime(),'end write file',i_sub_task
  o_img_data_arr=0
  
  
endfor

;obj_destroy,intep_obj
; PROFILER,/report
cur_line=0
shmunmap,'cur_line'

if flag_mask  ne 0 then begin
  ;    mask=read_tiff(maskfile)
  shmunmap,'mask_arr'
endif

print,systime(),'end session'
END



PRO imgset_recon, IMGSETFILE=imgsetfile,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=odir,$
  DATATYPE=dtype,OTYPE=otype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,MASKFILE=maskfile,PBAR=pbar,$
  PREF=pref
 ; ON_ERROR,2
  CATCH,Error_status
  IF Error_status NE 0 THEN BEGIN
    CATCH,/cancel
    PRINT, 'Error index: ', Error_status
    PRINT, 'Error message: ', !ERROR_STATE.MSG
    IF N_ELEMENTS(bridge_arr) NE 0 THEN begin
      for i=0,N_ELEMENTS(bridge_arr)-1 do begin
       if bridge_arr[i_session]->Status() EQ 1 then bridge_arr[i_session]->abort
      endfor
     OBJ_DESTROY,bridge_arr
     endif

    IF N_ELEMENTS(o_lun) NE 0 THEN FREE_LUN,o_lun
    IF N_ELEMENTS(tmp_lun) NE 0 THEN FREE_LUN,tmp_lun
    if n_elements(cur_line) gt 1 then begin
      cur_line=0
      shmunmap,'cur_line'
    endif
    MESSAGE, /REISSUE_LAST
  ENDIF
  ;ni=n_elements(its)
  IF N_ELEMENTS(ni) EQ 0 THEN MESSAGE,'please give the # of images in the packed file.'
  ;if total(file_test(imglist,/regular)) le ni then message,'some of the image files do not exist. Please check again.'
  
  
  IF N_ELEMENTS(ots) NE 0 THEN BEGIN
    IF N_ELEMENTS(its) NE ni THEN MESSAGE,'The its should has the same length as image list.'
  ENDIF
  IF N_ELEMENTS(oimgsetfiles) EQ 0 THEN MESSAGE,'please specific a full path for output files.'
  
  if N_elements(otype) eq 0 then otype=dtype
  if n_elements(pref) eq 0 then pref='Unknown'
  ;IF N_ELEMENTS(oprefix) EQ 0 THEN oprefix='ReconPrefix'
  
  ;line=208 & sample=400000 & ni =23 & no =365
  ;first determine the # of procesessors to be used based on the total task
  n_session=line/50+1 LE !CPU.HW_NCPU-1 ? line/50+1:!CPU.HW_NCPU-1
  
  ;Then determine the optimal subtasks for each processors consider the # of processor and aviliable physical memory
  ;The total pysical memory used by all processor can not exceed 70% of mochine's physical memory and each processor should hold equivalent task
  Machine_mem=long64(2)^34
  MEM_for_IDL_SESSION=Machine_mem*0.7/n_session
  no=n_elements(ots)
  ;get the byte length of the datatype
  nblen=N_TAGS({t:MAKE_ARRAY(1,type=4)},/data_length)
  max_sub_task_line=fix(MEM_for_IDL_SESSION/(sample*nblen*(ni+no)))
  sub_task_line=max_sub_task_line
  n_session_sub_task=(line / (sub_task_line*n_session))
  n_sub_task=n_session_sub_task*n_session
  if n_sub_task lt n_session then begin
    n_sub_task=n_session
    n_session_sub_task=1
    sub_task_line=line/n_session
  endif
  sub_task_line=line/n_sub_task
  
  
  start_l_arr=indgen(n_sub_task)*sub_task_line
end_l_arr=start_l_arr+sub_task_line-1

r_line=line mod (n_sub_task*sub_task_line)
if r_line ne 0 then begin
  start_l_arr[n_sub_task-r_line:n_sub_task-1]=start_l_arr[n_sub_task-r_line:n_sub_task-1]+indgen(r_line)
end_l_arr[n_sub_task-r_line:n_sub_task-1]=end_l_arr[n_sub_task-r_line:n_sub_task-1]+indgen(r_line)+1
endif

start_l_arr=reform(start_l_arr,n_session_sub_task,n_session)
end_l_arr=reform(end_l_arr,n_session_sub_task,n_session)

;print,transpose([[start_l_arr],[end_l_arr]])


;segemente the image for each session
;line=330
;n_session=1
;seg_line=line/n_session
;start_l_arr=BINDGEN(n_session)*seg_line
;end_l_arr=(BINDGEN(n_session)+1)*seg_line-1
;r_line=line MOD n_session
;start_l_arr[0:r_line]=start_l_arr[0:r_line]+1
;end_l_arr[n_session-1]=end_l_arr[n_session-1]+line MOD n_session
;print,n_session,start_l_arr,end_l_arr
cur=odir
tmpdir=FILEPATH('tmp',root_dir=cur)
IF ~file_test(tmpdir,/dir) THEN FILE_MKDIR,tmpdir



;restore save file
savfile=FILEPATH('TSR_SESSION.sav',ROOT_dir=cur,sub='tmp')
SAVE,/ROUTINES,FILENAME=savfile

;log file for each session
fname=pref+'.session'+string(INDGEN(n_session),format='(I02,".")')+string(randomu(systime(1),n_session)*long64(2)^20,format='(I08)')+'.log'
bridge_log_arr=FILEPATH(fname,ROOT_dir=cur,sub='tmp')

;temporal output file
;bridge_tmpo_arr=FILEPATH(STRING(INDGEN(n_session),FORMAT='("SESSION",I02,".dat")'),ROOT_dir=cur,sub='tmp')
tmp_packed_file=FILEPATH('recon_packed.dat',ROOT_dir=cur,sub='tmp')
;create the temperal output file
if ~file_test(tmp_packed_file) then begin
  openw,1,tmp_packed_file
  close,1
endif

;goto,tt
bridge_arr=OBJARR(n_session)



;create a shared memory to monitoring current processing state of of each session
SHMMAP,'cur_line',n_session,/long
cur_line=SHMVAR('cur_line')
cur_line[0]=reform(start_l_arr[0,*])

;if a maskfile avaliable then share the mask array between sessions
flag_mask=0
if n_elements(maskfile) ne 0 then begin
  mask=read_tiff(maskfile)
  SHMMAP,'mask_arr',DIMENSION=[sample,line],/byte
  mask_arr=shmvar('mask_arr')
  mask_arr[*,*]=mask
  flag_mask=1
endif

;;create a shared memory to retrive the sample series from each sub session
;SHMMAP,'sample_series',dimention=[ni+no,n_session], type=dtype
;sample_series=shmvar('sample_series')

 FOR i_session =0, n_session-1 DO BEGIN
  bridge_arr[i_session]=OBJ_NEW('IDL_IDLBridge',$
    OUTPUT=bridge_log_arr[i_session])

  bridge_arr[i_session]->SetVar,'imgsetfile',imgsetfile
  bridge_arr[i_session]->SetVar,'ni',ni
  bridge_arr[i_session]->SetVar,'start_l_arr',start_l_arr[*,i_session]
  bridge_arr[i_session]->SetVar,'end_l_arr',end_l_arr[*,i_session]
  ;bridge_arr[i_session]->SetVar,'tmpo',bridge_tmpo_arr[i_session]
  bridge_arr[i_session]->SetVar,'tmp_packed_file',tmp_packed_file
  bridge_arr[i_session]->SetVar,'sample',sample
  bridge_arr[i_session]->SetVar,'line',line
  bridge_arr[i_session]->SetVar,'datatype',dtype
  bridge_arr[i_session]->SetVar,'its',its
  bridge_arr[i_session]->SetVar,'ots',ots
  bridge_arr[i_session]->SetVar,'otype',otype
  bridge_arr[i_session]->SetVar,'flag_mask',flag_mask
  bridge_arr[i_session]->SetVar,'cur_session',i_session
  bridge_arr[i_session]->SetVar,'n_session',n_session
  bridge_arr[i_session]->SetVar,'pref',pref
  
  
  bridge_arr[i_session]->Execute,'COMPILE_OPT IDL2'
  bridge_arr[i_session]->Execute,'cd,"'+cur+'"'
  bridge_arr[i_session]->Execute,"restore,'"+savfile+"'"
  ;bridge_arr[i_session]->Execute,"help, /ROUTINES"
  
  cmdStr='TSR_SESSION,imgsetfile,ni,start_l_arr,end_l_arr,'+ $
  'tmp_packed_file,sample,line,datatype,its,ots,otype,'+ $
  'flag_mask,cur_session,n_session,pref'
  bridge_arr[i_session]->Execute,cmdStr,/NOWAIT
ENDFOR

;TSR_SESSION,imgsetfile,ni,start_l_arr[i_session],end_l_arr[i_session],bridge_tmpo_arr[i_session],sample,line,dtype,its,ots
flag=1
WHILE (flag) DO BEGIN

  flag=0
  FOR i_session=0,n_session-1 DO BEGIN
   ; if ~ obj_valid( bridge_arr[i_session]) then continue
   ; if bridge_arr[i_session]->Status() EQ 0 then OBJ_DESTROY,bridge_arr[i_session] & continue
    IF (bridge_arr[i_session]->Status() EQ 3) THEN MESSAGE,STRING(i_session,FORMAT='("Session",I1,"failed!")')
    flag=(bridge_arr[i_session]->Status() EQ 1) OR flag
    ; cur_line=bridge_arr[i_session]->getvar('cur_line')
    ;print,i_session,cur_line
  ENDFOR
  ;print,cur_line
  IF pbar -> CheckCancel() THEN BEGIN
    ok = Dialog_Message('The user cancelled operation.')
    FOR i_session=0,n_session-1 DO  res=bridge_arr[i_session]->abort() 

    ;destroy objects
    OBJ_DESTROY,bridge_arr

    cur_line=0
    SHMUNMAP,'cur_line'

    if flag_mask  ne 0 then begin
      ;    mask=read_tiff(maskfile)
      shmunmap,'mask_arr'
    endif
    RETURN
  ENDIF
  pbar->update,total(cur_line-start_l_arr[0,*])/float(line)*100
  print, total(cur_line-start_l_arr[0,*])/float(line)*100, format='(f4.1,"% finished!")'
  WAIT,10
ENDWHILE

;destroy objects
OBJ_DESTROY,bridge_arr

cur_line=0
SHMUNMAP,'cur_line'

if flag_mask  ne 0 then begin
  ;    mask=read_tiff(maskfile)
  shmunmap,'mask_arr'
endif


;now plot the sample series for each session
;for i_session=0,n_session-1 do begin
;  
;endfor

;  return
tt:
;after all sesssion finished, then merge result file.
;OPENW,o_lun,oimgsetfile,/GET_LUN

no=n_elements(ots)
;get the byte length of the datatype

noblen=N_TAGS({t:MAKE_ARRAY(1,type=otype)},/data_length)
o_img_arr=make_array(sample,line,TYPE=otype)
openr,o_lun,tmp_packed_file,/get_lun


;FOR i_session=0, n_session-1 DO BEGIN
;  readu,o_lun,o_img_arr
;  ;tmp_arr=MAKE_ARRAY(sample,end_l_arr[i_session]-start_l_arr[i_session]+1,ni)
;  OPENR,tmp_lun,bridge_tmpo_arr[i_session],/get_lun
;  s_o_img_arr=fltarr(sample,end_l_arr[i_session]-start_l_arr[i_session]+1,no)
;  readu,tmp_lun,s_o_img_arr
;  o_img_arr[*,start_l_arr[i_session]:end_l_arr[i_session],*]=s_o_img_arr
;  ;readu,tmp_lun,tmp_arr
;  ;POINT_LUN,o_lun,start_l_arr[i_session]*long64(nblen)
;  ;point_lun,o_lun,start_l_arr[i_session]*nblen
;  ;    FOR i_no=0,  no-1 DO BEGIN
;  ;      COPY_LUN,tmp_lun,o_lun,sample*long64(end_l_arr[i_session]-start_l_arr[i_session]+1)*nblen,TRANSFER_COUNT=n1
;  ;      point_LUN,o_lun,sample*(long64(line)*i_no+start_l_arr[i_session])*nblen    ;skip to next tempral file
;  ;      print,i_no,i_session
;  ;    ENDFOR
;  FREE_LUN,tmp_lun
;ENDFOR

;output to tiff file
FOR i_no=0,  no-1 DO BEGIN
   readu,o_lun,o_img_arr
   case otype of
    1:write_tiff,oimgsetfiles[i_no],o_img_arr,GEOTIFF=geotiff
    2:write_tiff,oimgsetfiles[i_no],o_img_arr,GEOTIFF=geotiff,/short,/signed
    3:write_tiff,oimgsetfiles[i_no],o_img_arr,GEOTIFF=geotiff,/LONG
    4:write_tiff,oimgsetfiles[i_no],o_img_arr,GEOTIFF=geotiff,/float
   endcase
  ;write_tiff,oimgsetfiles[i_no],o_img_arr,GEOTIFF=geotiff,/float;,/short,/signed
  print,i_no
ENDFOR

;close the output file
FREE_LUN,o_lun



END

PRO Pack_img, imgfilelist,outfile,sdsname=sdsname,TYPE=type
  ;+
  ; :Description:
  ;    Describe the procedure.
  ;
  ; :Params:
  ;    imgfilelist
  ;    outfile
  ;
  ;
  ;
  ; :Date: Dec 21, 2012
  ;
  ; :Author: lenovo
  ;-
  if N_elements(type) eq 0 then type='HDF'
  OPENW,lun,outfile,/GET_LUN
  case type of
    'HDF':begin
    FOR i=0, N_ELEMENTS(imgfilelist)-1 DO BEGIN
    
      hdf_id=HDF_SD_START(imgfilelist[i])
      sds_index=HDF_SD_NAMETOINDEX(hdf_id,sdsname)
      sds_id=HDF_SD_SELECT(hdf_id,sds_index)
      HDF_SD_GETDATA,sds_id,data,START=[0,0]
      
      ;write data
      WRITEU,lun,data
      
      
      HDF_SD_ENDACCESS, sds_id
      HDF_SD_END, hdf_id
    ENDFOR
  end
  'TIFF': begin
  
    FOR i=0, N_ELEMENTS(imgfilelist)-1 DO BEGIN
    
      data=read_tiff(imgfilelist[i])
      
      ;write data
      WRITEU,lun,data
      
      
      ;HDF_SD_ENDACCESS, sds_id
      ; HDF_SD_END, hdf_id
    endfor
  end
endcase

FREE_LUN,lun
END

PRO Choose_Item_Event, event
  
  if widget_info(event.id,/UNAME) eq 'OK' then begin 
  ; Get info structure pointer.
  Widget_Control, event.top, Get_UValue=ptr

  ; Get the item from the list widget and store it
  ; in the program pointer.
  listid=widget_info(event.top,FIND_BY_UNAME='list')
  Widget_Control, listid, Get_UValue=listValues
  (*ptr).answer = strjoin(listValues[widget_info(listid,/LIST_SELECT)],' ')
  (*ptr).cancel = 0
  

  ; Destroy the widget.
 Widget_Control, event.top, /Destroy
 
 endif

END ; ----------------------------------------------------


FUNCTION Choose_Item, items, Cancel=cancel, Group_Leader=group_leader

  IF N_Elements(items) EQ 0 THEN items = ['cow', 'dog', 'coyote', 'pig']

  ; Create the top-level base widget. Make it modal if you
  ; have a group_leader. Otherwise, will have to rely on this
  ; widget blocking. If so, DON'T call it from a widget program!
  IF N_Elements(group_leader) EQ 0 THEN BEGIN
    tlb = Widget_Base(Title='Animals...', Column=1 , YOFFSET=100,XOFFSET=100)
  ENDIF ELSE BEGIN
    tlb = Widget_Base(Title='Animals...', Column=1 , YOFFSET=100,XOFFSET=100, /modal,Group_Leader=group_leader)
  ENDELSE

  ; Create list widget with choices. Store the choices in the UVALUE
  ; of the list, so they are available in the event handler.
  listID = Widget_List(tlb, Value=items, UValue=items, /MULTIPLE,UNAME='list',$
   YSize=N_Elements(items) < 20, XSize=25)
  
  ok_btn=widget_button(tlb,VALUE='Processing',UNAME='OK')
  ; Create a pointer for storing the user's selection. There is a cancel
  ; flag so I can tell if user killed the widget without making a selection.
  ptr = Ptr_New({answer:"", cancel:1})

  ; Store info pointer in UVALUE of TLB.
  Widget_Control, tlb, Set_UValue=ptr

  ; Realize the widgets.
  Widget_Control, tlb, /Realize

  ; Call XManager and BLOCK the command line. This will only work
  ; if this is the FIRST blocking program. Use GROUP_LEADER and MODAL
  ; keyword to be *sure* you block.
  XMANAGER, 'choose_item', tlb

  ; Get the answer and cancel flag out of the program pointer and destroy the pointer.
  answer = (*ptr).answer
  cancel = (*ptr).cancel
  Ptr_Free, ptr

  ; Return the answer.
  ;RETURN, strsplit(answer,' ', /EXTRACT)
  return,strsplit(answer,' ')

END


PRO TIMESR

  goto,Heihe_LAI_PK;CN_LST;CN_NDVI_YD;CN_NDVI_MOD;heihe_albedo;heihe_LAI;Heihe_NDVI_YD;SN_NDVI_5days_ZB;HANTS_TIME;Test_SVD;
  ;NL_NDVI;NL_LAI;NL_NDVI;NL_NDVI;NDVISET;IMGSEET;PLOT_S;IMGSEET
  ni=46;23
  per=365
  its=INDGEN(ni)*8+1
  ots=INDGEN(46)*8+1
  nf1=4
  nf2=5
  pa1=[360,180,120]
  pa2=[720,360,180,120,90]
  HiLo='Low'
  low=11000;-2000
  high=17000;10000
  fet=250;500
  dod=10;5
  delta=0
  y=[236,97,193,217,1241,1230,1233,1313,2967,3506,6398,7721,8072,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  y=[13402,13088,13026,13253,13207,13076,13256,13252,13372,13555,13230,13704,13629,13807,14177,14489,14332,14476,14716,14963,14894,14992,15105,15224,14863,14604,14519,14660,15108,15027,15142,15008,14923,14987,14853,14587,14258,14188,14308,14212,14127,14177,13927,13457,13756,13493]

  y2=[323,264,498,58,2479,409,2651,3161,5153,5497,7238,8038,8117,8177,5341,6131,7544,5083,5183,2838,2583,2375,2331]
  y3=[0,400,900,2500,3000,2900,1600,1313,2967,3506,6398,7721,8072,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  y4=[0,400,900,2500,3000,2900,500,513,2967,3506,6398,7721,8072,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  y5=[0,400,900,2500,3200,3700,3800,3785,2967,2200,2000,3506,6398,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  y6=[0,400,900,2500,3200,3700,3800,3785,2967,900,1000,3506,6398,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  looplim= 10
  delta=0.5
  ;cgplot,y2,psym=38
  ;res1=HANTS(ni,per,ts,nf1,pa1,HiLo,low,high,fet,dod,y2,looplim,delta,/POLYNOMIAL)
  ;res2=HANTS(ni,per,ts,nf1,pa1,HiLo,low,high,fet,dod,y2,looplim,delta)
  series_arr=FLTARR(2,23)
  ;OPENR,1,'/home/rasdaman/idl/timesr/sample_series.txt'
  ;READF,1,series_arr
  ;CLOSE,1
  
  ;series_arr=[series_arr,REFORM(y2,1,23)]
  ;series_arr=[series_arr,REFORM(y3,1,23)]
  ;series_arr=[series_arr,REFORM(y4,1,23)]
  ;series_arr=[series_arr,REFORM(y5,1,23)]
  ;series_arr=[series_arr,REFORM(y6,1,23)]
  ;for i=0,8 do begin
;  intep_obj=obj_new('RECON_HA',TS=its,OTS=ots,/IOD,per=365,parr=[720,360,180,120,90],BATCH=1)
  FOR i=0,0 DO BEGIN
    ;y=REFORM(series_arr[i,*])
    
    its=INDGEN(ni)*8+1
    OTS=INDGEN(46)*8+1
    orts=_conbine_ts(ots,its,I2O=i2o,O2O=o2o)
    ;y[10:25]=0
    res4=GENERRIZED_RCON_MULTI(y,$
      TS=its,OTS=ots,$
      HILO=HiLo,LOW=low,HIGH=high,FET=fet,DOD=dod,LOOPLIM=looplim,per=365,parr=[720,360,180,80,10])
      
;      weight=fltarr(23)+1
;      weight[0:4]=0
;      res=R_ha(y,WEIGHT=weight, $
;      PER=per,TS=its,DELTA=delta,PARR=pa1,OTS=ots)
      
    ;res3=GENERRIZED_HANTS(HILO=HiLo,LOW=low,HIGH=high,FET=fet,DOD=dod,y,LOOPLIM=looplim,/MULTI,SMOOTHER='WHITSMW',$
     ; ITS=its,OTS=ots)
    ;res2=GENERRIZED_HANTS(HILO=HiLo,LOW=low,HIGH=high,FET=fet,DOD=dod,y,LOOPLIM=looplim,SMOOTHER='HANTS',$
    ; DELTA=delta,PER=per,PARR=pa1,TS=its)
    ;res1=HANTS(ni,per,ts,nf1,pa1,HiLo,low,high,fet,dod,y,looplim,delta)
    cgplot,its,y,psym=-6,COLOR='Black',linestyle=2, thick=1,XRANGE=[-10,365],yrange=[low,high]
    ;cgplot,its,res2,/overplot, Color='red',linestyle=0, thick=1
    cgplot,ots,res4,/overplot, Color='dodger blue',linestyle=2, thick=1
    ;cgplot,its,res.yr1,/overplot, Color='GREEN',linestyle=2, thick=1
    ;cgplot,ts,res1.yr,/overplot,color='green',linestyle=4,thick=2
    
    al_Legend, ['Raw', 'orignal HANTS','Modified'], PSym=[-6,-15,-16], $
      LineStyle=[1,0,2], Color=['Black','red','dodger blue'], Position=[5,-500]
    ;filename='P:\TUD\idl_projects\TimeSR\test'+STRING(i,format='(I03)')+'.jpg'
    ;  WRITE_JPEG, filename, TVRD(/TRUE), /TRUE
    ;cgdelete,/all
    ;WDELETE
  ENDFOR
  
  RETURN
  
  Test_SVD:
  ni=46;23
  per=365
  its=INDGEN(ni)*8+1
  ots=INDGEN(46)*8+1
  nf1=4
  nf2=5
  pa1=[360,180,120]
  pa2=[720,360,180,120,90]
  HiLo='Low'
  low=11000;-2000
  high=17000;10000
  fet=250;500
  dod=10;5
  delta=0.5
  y=[236,97,193,217,1241,1230,1233,1313,2967,3506,6398,7721,8072,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  y2=[323,264,498,58,2479,409,2651,3161,5153,5497,7238,8038,8117,8177,5341,6131,7544,5083,5183,2838,2583,2375,2331]
  y3=[0,400,900,2500,3000,2900,1600,1313,2967,3506,6398,7721,8072,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  y4=[0,400,900,2500,3000,2900,500,513,2967,3506,6398,7721,8072,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  y5=[0,400,900,2500,3200,3700,3800,3785,2967,2200,2000,3506,6398,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
  y6=[0,400,900,2500,3200,3700,3800,3785,2967,900,1000,3506,6398,7759,7620,7766,4953,5456,5456,2075,1866,1475,288]
 
  weight=fltarr(ni)+1
  ;weight[15:21]=0
 ; weight[0:5]=0
  
  series_arr=FLTARR(2,23)
  OPENR,1,'P:\TUD\idl_projects\TimeSR\sample_series.txt'
  READF,1,series_arr
  CLOSE,1
  
  series_arr=[series_arr,REFORM(y2,1,23)]
  series_arr=[series_arr,REFORM(y3,1,23)]
  series_arr=[series_arr,REFORM(y4,1,23)]
  series_arr=[series_arr,REFORM(y5,1,23)]
  series_arr=[series_arr,REFORM(y6,1,23)]
  
  y=series_arr[1,*]
  res=R_ha(y,WEIGHT=weight, $
    PER=per,TS=its,DELTA=delta,PARR=pa1)
  cgplot,its,y,psym=-6,COLOR='Black',linestyle=2, thick=1,XRANGE=[-10,365],yrange=[-3500,10000]
  cgplot,its,res.yr,/overplot, Color='red',linestyle=0, thick=1
  cgplot,its,res.yr1,/overplot, Color='dodger blue',linestyle=2, thick=1
  ;cgplot,its,res.yr1,/overplot, Color='GREEN',linestyle=2, thick=1
  ;cgplot,ts,res1.yr,/overplot,color='green',linestyle=4,thick=2
  
  al_Legend, ['Raw', 'orignal HANTS','Modified'], PSym=[-6,-15,-16], $
    LineStyle=[1,0,2], Color=['Black','red','dodger blue'], Position=[5,-500]
  return
  IMGSEET:
  
  ;following process intend to reconstruct MOD13A2 h25v05 NDVI dataset
  datapool_dir='H:\data\MODIS\VI\VI-16days-china\Terra\2010_china_vi'
  
  out_dir='O:\TUD\data\rcon-result\NDVI'
  h25v05_filelist=FILE_SEARCH(datapool_dir,'*h25v05*.hdf')
  
  packed_file=FILEPATH('h25v05_2010_NDVI_packed.bin',ROOT_DIR=out_dir)
  IF ~file_test(packed_file) THEN $
    Pack_img,h25v05_filelist,packed_file,sdsname='1 km 16 days NDVI'
    
    
    
  ni=23
  ;datatype=2
  sample=1200
  line=1200
  its=INDGEN(ni)*16+1
  ots=INDGEN(365)+1
  dtype=2
  
  oimgsetfiles=FILEPATH(string(ots,format='("h25v05_2010_NDVI_recon.",I03,".tif")'),ROOT_DIR=out_dir)
  imgset_recon,IMGSETFILE=packed_file,NI=23, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
    DATATYPE=dtype,SAMPLE=sample,LINE=line
    
  return
    
  ;***********************************************************************************************************************
  
  
  NDVISET:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered China
  datapool_dir='O:\Beijing\data\China_VI\Terra-MODIS\VI\mrt'
  
  out_dir='O:\TUD\data\rcon-result\NDVI\china\2010'
  china_filelist=FILE_SEARCH(datapool_dir,'Terra-MODIS_VI_2010???.1_km_16_days_NDVI.tif')
  
  packed_file=FILEPATH('china_2010_NDVI_packed.bin',ROOT_DIR=out_dir)
  IF ~file_test(packed_file) THEN $
    Pack_img,china_filelist,packed_file,TYPE='TIFF'
    
  ni=23
  ;datatype=2
  sample=4589
  line=4247
  its=INDGEN(ni)*16+1
  ots=INDGEN(365)
  dtype=2
  res=query_tiff(china_filelist[0],GEOTIFF=geotiff)
  oimgsetfiles=FILEPATH(string(ots,format='("CHINA_2010_NDVI_recon.",I03,".tif")'),ROOT_DIR=out_dir)
  imgset_recon,IMGSETFILE=packed_file,NI=23, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
    DATATYPE=dtype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff
    
  return
  
 ;************************************************************************************************************************** 
    
  PLOT_S:
  out_dir='O:\TUD\data\rcon-result\NDVI\china\2010\tmp'
  log_files=file_search(out_dir,'*.log')
  for i=0,7 do begin
    openr,1,log_files[i]
    skip_lun,1,35,/LINES
    is=fltarr(2,23)
    readf,1,is
    skip_lun,1,1,/LINES
    os=fltarr(2,365)
    readf,1,os
    pixw=cgpixmap()
    cgplot,is[0,*],is[1,*],color='red',psym=2,/window
    cgplot,os[0,*],os[1,*],color='dodger blue',/addcmd,/overplot
    pixw->output,filepath('series'+strtrim(i,2)+'.png',root_dir=out_dir)
    cgdelete,pixw
    close,1
  endfor
  
  return
  ;**********************************************************************************************************************************
  
  NL_NDVI:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered Natherlands
  datapool_dir='O:\TUD\data\netherlands-ndvi-lai\NDVI'
  
  out_dir='O:\TUD\data\netherlands-ndvi-lai\recon-ndvi'
  year=['2008','2009','2010','2011','2012']
  for i_year=0,0 do begin
  IN_filelist=FILE_SEARCH(datapool_dir,'*'+year[i_year]+'*tif')
  
  packed_file=FILEPATH('NL_'+year[i_year]+'_NDVI_packed.bin',ROOT_DIR=out_dir)
  IF ~file_test(packed_file) THEN $
    Pack_img,in_filelist,packed_file,TYPE='TIFF'
    
  ni=23
  ;datatype=2
  sample=50
  line=120
  its=INDGEN(ni)*8+1
  ots=INDGEN(365)+1
  dtype=2
  res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
  oimgsetfiles=FILEPATH('NL_'+year[i_year]+'_NDVI_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
  imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
    DATATYPE=dtype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff
    endfor
    
  return
  
  NL_LAI:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered Natherlands
  datapool_dir1='O:\Beijing\data\China_VI\aqua-MODIS\2009'
  
  out_dir='O:\TUD\data\netherlands-ndvi-lai\recon-LAI'
  year=['2008','2009','2010','2011','2012']
  for i_year=1,1 do begin
  IN_filelist=FILE_SEARCH(datapool_dir,'*'+year[i_year]+'*tif')
  
  packed_file=FILEPATH('NL_'+year[i_year]+'_LAI_packed.bin',ROOT_DIR=out_dir)
  IF ~file_test(packed_file) THEN $
    Pack_img,in_filelist,packed_file,TYPE='TIFF'
    
  ni=46
  ;datatype=2
  sample=50
  line=120
  its=INDGEN(ni)*8+1
  ots=INDGEN(365)+1
  dtype=1
  res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
  oimgsetfiles=FILEPATH('NL_'+year[i_year]+'_LAI_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
  imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
    DATATYPE=dtype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff
    endfor
  return
  
  ;****************************************************************************************************************
  CN_NDVI_YD:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered Natherlands
  ;datapool_dir_aqua='O:\Beijing\data\China_VI\aqua-MODIS'
  ;datapool_dir_terra='O:\Beijing\data\China_VI\Terra-MODIS'
  datapool_dir='P:\data\reprocessed\china\vi\mrt'
  out_dir='P:\data\reprocessed\china\vi\recon'
  
  maskfile=filepath('mask_cn_geo.tif',root_dir='P:\data\reprocessed\china')
  year=['2003','2007','2006','2005','2004','2008','2009','2010','2011','2012']
  
  for i_year=3,4 do begin
    cgProgressBar1 = Obj_New("cgProgressBar", /Cancel,$
      TEXT=year[i_year],Title='HANTS processing')
    cgProgressBar1-> Start
  IN_filelist=strarr(46)
  IN_filelist[1:45:2]=FILE_SEARCH(datapool_dir,'MYD*.A'+year[i_year]+'*tif')
  IN_filelist[0:45:2]=FILE_SEARCH(datapool_dir,'MOD*.A'+year[i_year]+'*tif')
  
  ;out_dir=filepath('',root=out_dir,sub=year[i_year])
  packed_file=FILEPATH('CN_'+year[i_year]+'_NDVI_packed.bin',ROOT_DIR=out_dir)
  IF ~file_test(packed_file) THEN $
    Pack_img,in_filelist,packed_file,TYPE='TIFF'
    
  ni=46
  ;datatype=2
  sample=6163
  line=3537
  its=INDGEN(ni)*8+1
  ots=INDGEN(365)+1
  dtype=2
  res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
  oimgsetfiles=FILEPATH('CN_'+year[i_year]+'_NDVI_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
  imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
    DATATYPE=dtype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,MASKFILE=maskfile,PBAR=cgProgressBar1
    obj_destroy,cgProgressBar1
    endfor
    
  return
  
  ;****************************************************************************************************************
  CN_NDVI_MOD:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered Natherlands
  ;datapool_dir_aqua='O:\Beijing\data\China_VI\aqua-MODIS'
  ;datapool_dir_terra='O:\Beijing\data\China_VI\Terra-MODIS'
  datapool_dir='O:\data\reprocessed\china\vi\mrt'
  out_dir='O:\data\reprocessed\china\vi\recon'
  
  maskfile=filepath('mask_cn_geo.tif',root_dir='P:\Beijing\data\China_VI\aqua-MODIS')
  year=['2001','2002','2006','2005','2004','2008','2009','2010','2011','2012']
  
  for i_year=0,1 do begin
    cgProgressBar1 = Obj_New("cgProgressBar", /Cancel,$
      TEXT=year[i_year],Title='HANTS processing')
    cgProgressBar1-> Start
    IN_filelist=strarr(23)
    ;IN_filelist[1:45:2]=FILE_SEARCH(datapool_dir,'MYD*'+year[i_year]+'*tif')
    IN_filelist=FILE_SEARCH(datapool_dir,'MOD13A2.A'+year[i_year]+'*.1_km_16_days_NDVI.tif')
    
    ;out_dir=filepath('',root=out_dir,sub=year[i_year])
    packed_file=FILEPATH('CN_'+year[i_year]+'_NDVI_packed.bin',ROOT_DIR=out_dir)
    IF ~file_test(packed_file) THEN $
      Pack_img,in_filelist,packed_file,TYPE='TIFF'
      
    ni=23
    ;datatype=2
    sample=6163
    line=3537
    its=INDGEN(ni)*16+1
    ots=INDGEN(365)+1
    dtype=2
    res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
    oimgsetfiles=FILEPATH('CN_'+year[i_year]+'_NDVI_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
    imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
      DATATYPE=dtype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,MASKFILE=maskfile,PBAR=cgProgressBar1
    obj_destroy,cgProgressBar1
  endfor
  
  return
  ;****************************************************************************************************************
  CN_LST:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered Natherlands
  ;datapool_dir_aqua='O:\Beijing\data\China_VI\aqua-MODIS'
  datapool_dir='P:\data\reprocessed\china\lst\mrt\'
  
  out_dir='P:\data\reprocessed\china\lst\recon\'
  
  maskfile=filepath('mask_cn_geo.tif',root_dir='P:\data\reprocessed\china')
  year=['2019','2012','2013','2014']
  year=string(indgen(15)+2003,format='(I4)')
  
  for i_year=7,7 do begin
    cgProgressBar1 = Obj_New("cgProgressBar", /Cancel,$
      TEXT=year[i_year],Title='HANTS processing')
    cgProgressBar1-> Start
    IN_filelist=FILE_SEARCH(filepath('',root=datapool_dir),'*'+year[i_year]+'*tif')
 
    
    ;out_dir=filepath('',root=out_dir,sub=year[i_year])
    packed_file=FILEPATH('CN_'+year[i_year]+'_LST_packed.bin',ROOT_DIR=out_dir)
    IF ~file_test(packed_file) THEN $
      Pack_img,in_filelist,packed_file,TYPE='TIFF'
      
    ni=46
    ;datatype=2
    sample=6163
    line=3537
    its=INDGEN(ni)*8+1
    ots=INDGEN(365)+1
    dtype=2
    res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
    oimgsetfiles=FILEPATH('CN_'+year[i_year]+'_LST_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
    imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
      DATATYPE=dtype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,MASKFILE=maskfile,PBAR=cgProgressBar1
    obj_destroy,cgProgressBar1
  endfor
  
  return
  
  
  Heihe_NDVI_YD:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered Natherlands
  ;datapool_dir_aqua='O:\Beijing\data\China_VI\aqua-MODIS'
  ;datapool_dir_terra='O:\Beijing\data\China_VI\Terra-MODIS'
  datapool_dir='P:\TUD\zs-comp\data\ndvi\mrt'
  out_dir='P:\TUD\zs-comp\data\ndvi\recon'
  
  ;maskfile=filepath('mask_cn_geo.tif',root_dir='O:\Beijing\data\China_VI\aqua-MODIS')
 ; year=['2007','2003','2006','2005','2004','2008','2009','2010','2011']
  
  year=['2013','2010','2009','2008']
  for i_year=0,0 do begin
    cgProgressBar1 = Obj_New("cgProgressBar", /Cancel,$
      TEXT=year[i_year],Title='HANTS processing')
    cgProgressBar1-> Start
    IN_filelist=strarr(46);FILE_SEARCH(datapool_dir,'MOD*'+year[i_year]+'*tif');
    IN_filelist[1:45:2]=FILE_SEARCH(datapool_dir,'MYD*A'+year[i_year]+'*tif')
   IN_filelist[0:45:2]=FILE_SEARCH(datapool_dir,'MOD*A'+year[i_year]+'*tif')
    
    ;out_dir=filepath('',root=out_dir,sub=year[i_year])
    packed_file=FILEPATH('heihe_'+year[i_year]+'_NDVI_packed.bin',ROOT_DIR=out_dir)
    IF ~file_test(packed_file) THEN $
      Pack_img,in_filelist,packed_file,TYPE='TIFF'
      
    ni=46
    ;datatype=2
    sample=600
    line=550
    its=INDGEN(ni)*8+1
    ots=INDGEN(365)+1
    dtype=2
    print,year[i_year]
    res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
    oimgsetfiles=FILEPATH('heihe_'+year[i_year]+'_NDVI_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
    imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
      DATATYPE=dtype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,PBAR=cgProgressBar1,PREF=year[i_year];,MASKFILE=maskfile
      obj_destroy,cgProgressBar1
  endfor
  
  return
  
  heihe_albedo:

  datapool_dir1='P:\TUD\zs-comp\data\albedo\mrt'
  
  out_dir='P:\TUD\zs-comp\data\albedo\recon'
  year=['2013']
  for i_year=0,0 do begin
    cgProgressBar1 = Obj_New("cgProgressBar", /Cancel,$
      TEXT=year[i_year],Title='HANTS processing')
    cgProgressBar1-> Start
    IN_filelist=FILE_SEARCH(datapool_dir1,'*'+year[i_year]+'*.Albedo_WSA_shortwave.tif')
    
    packed_file=FILEPATH('heihe_'+year[i_year]+'_albedo_packed.bin',ROOT_DIR=out_dir)
    IF ~file_test(packed_file) THEN $
      Pack_img,in_filelist,packed_file,TYPE='TIFF'
      
    ni=46
    ;datatype=2
    sample=600
    line=550
    its=INDGEN(ni)*8+1
    ots=INDGEN(365)+1
    dtype=2
    res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
    oimgsetfiles=FILEPATH('heihe_'+year[i_year]+'_Albedo_WSA_shortwave_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
    imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
      DATATYPE=dtype,OTYPE=otype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,PBAR=cgProgressBar1
      
    obj_destroy,cgProgressBar1
  endfor
  return
  
  heihe_LAI:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered Natherlands
  datapool_dir1='P:\TUD\zs-comp\data\LAI\mrt'
  
  out_dir='P:\TUD\zs-comp\data\LAI\recon'
  year=['2013']
  for i_year=0,0 do begin
    cgProgressBar1 = Obj_New("cgProgressBar", /Cancel,$
      TEXT=year[i_year],Title='HANTS processing')
    cgProgressBar1-> Start
    IN_filelist=FILE_SEARCH(datapool_dir1,'*'+year[i_year]+'*.LAI_1km.tif')
    
    packed_file=FILEPATH('heihe_'+year[i_year]+'_LAI_packed.bin',ROOT_DIR=out_dir)
    IF ~file_test(packed_file) THEN $
      Pack_img,in_filelist,packed_file,TYPE='TIFF'
      
    ni=46
    ;datatype=2
    sample=600
    line=550
    its=INDGEN(ni)*8+1
    ots=INDGEN(365)+1
    dtype=1
    otype=4
    res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
    oimgsetfiles=FILEPATH('heihe_'+year[i_year]+'_LAI_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
    imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
      DATATYPE=dtype,OTYPE=otype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,PBAR=cgProgressBar1
     
     obj_destroy,cgProgressBar1
  endfor
  return

  ;***********************************************************************************
  ;Following code is used to processing LAI dataset provided by Prof.Fan (Peking University) (For Dr. Chaolei Zheng 18/11/2015)
   heihe_LAI_PK:
  datapool_dir1='/mnt/shares/share1/heihe/LAI-Beida/merged'
  
  out_dir='/mnt/shares/share1/heihe/LAI-Beida/recon'
  ;year=['2013']
  year = string(indgen(13)+2000,format="(I4)")
  for i_year=2,12 do begin
    if i_year eq 1 then continue
    cgProgressBar1 = Obj_New("cgProgressBar", /Cancel,$
      TEXT=year[i_year],Title='HANTS processing')
    cgProgressBar1-> Start
    IN_filelist=FILE_SEARCH(datapool_dir1,'*'+year[i_year]+'???*.tif') 
   ; print,IN_filelist
    print,"Got the file list."
    ;break
   
    
    packed_file=FILEPATH('heihe_'+year[i_year]+'_LAI_packed.bin',ROOT_DIR=out_dir)
    IF ~file_test(packed_file) THEN $
      Pack_img,in_filelist,packed_file,TYPE='TIFF'
    
    print,"finished to pack."
   
      
    ni=46
    ;datatype=2
    sample=410
    line=561
    its=INDGEN(ni)*8+1
    ots=INDGEN(365)+1
    dtype=5
    otype=4
    res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
    oimgsetfiles=FILEPATH('heihe_'+year[i_year]+'_LAI_recon.'+string(ots,format='(I03,".tif")'),ROOT_DIR=out_dir)
    imgset_recon,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
      DATATYPE=dtype,OTYPE=otype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,PBAR=cgProgressBar1
     
     obj_destroy,cgProgressBar1
  endfor
  return  
  ;********************************************************************************************************;
  ;The code for this sucsection process the NDVI data form MODIS. These data are provided by Bo Zhong. The
  ;data spans from 2013001 to 2014176 every days. say the data for 2014 is not a complete year. Another issue
  ;is that some tiles has lost from the dataset, so during processing, we should specify concrete doy to each
  ;tile.
   
  SN_NDVI_5days_ZB:
  ;following process intend to reconstruct MOD13A2 NDVI dataset covered Natherlands
  ;datapool_dir_aqua='O:\Beijing\data\China_VI\aqua-MODIS'
  ;datapool_dir_terra='O:\Beijing\data\China_VI\Terra-MODIS'
  datapool_dir='J:\Data\MUSQ_Output\SoutheastAsia_1km\NDVI'
  out_dir='J:\Data\MUSQ_Output\SoutheastAsia_1km\NDVI_HANTS'

  tiles=["H23V04","H23V05","H24V04","H24V05","H24V06","H25V03","H25V04",$
    "H25V05","H25V06","H26V03","H26V04","H26V05","H26V06","H26V07","H27V04",$
    "H27V05","H27V06","H27V07","H27V08","H27V09","H28V04","H28V05","H28V06",$
    "H28V07","H28V08","H28V09","H29V05","H29V06","H29V07","H29V08","H29V09",$
    "H29V10","H30V07","H30V08","H30V09","H30V10","H31V08","H31V09","H31V10","H32V09"]
  

  tiles=choose_item(tiles)

 

  
  
  for i_tiles=0,n_elements(tiles)-1 do begin
    cgProgressBar1 = Obj_New("cgProgressBar", /Cancel,$
      TEXT=tiles[i_tiles],Title='HANTS processing')
    cgProgressBar1-> Start
    
    ;MuSQ.VI.1km.2013001000000.H21V03.001.h5EVI-DataSet_EVI.tif
    fpattern_str=string(tiles[i_tiles],$
      format='("MuSQ.VI.1km.201????000000.",A6,".001.h5NDVI-DataSet_NDVI.tif")')
    IN_filelist=FILE_SEARCH(datapool_dir,fpattern_str)
    
    nts=n_elements(IN_filelist)
    
    ;extract the doy infomation
    start_doy='2013001'
    end_doy='2014176'
    nodoy=365/5
    ots=[indgen(nodoy)*5+1,365+indgen(nodoy/2)*5+1]
    odoys=string((ots/365)+2013,format='(I04)')+string(ots-(ots/365)*365,format='(I03)')
    its=intarr(nts)
    basename=file_basename(IN_filelist)
    for i=0,nts-1 do begin
      year=fix(strmid(basename[i],12,4))
      doy=fix(strmid(basename[i],16,3))
      its[i]=(year-2013)*355 + doy
    endfor
    
    

    ;out_dir=filepath('',root=out_dir,sub=year[i_year])
    packed_file=FILEPATH(tiles[i_tiles]+'_NDVI_packed.bin',ROOT_DIR=out_dir)
    IF ~FILE_TEST(packed_file) THEN $
      PACK_IMG,in_filelist,packed_file,TYPE='TIFF'

    ni=nts
    ;datatype=2
    sample=1200
    line=1200
    ;its=INDGEN(ni)*8+1
    ;ots=INDGEN(365)+1
    dtype=2
    ;print,year[i_year]
    res=query_tiff(in_filelist[0],GEOTIFF=geotiff)
    oimgsetfiles=FILEPATH('MuSQ.VI.1km.NDVI.'+odoys+'.'+tiles[i_tiles]+'.001.tif',ROOT_DIR=out_dir)
    IMGSET_RECON,IMGSETFILE=packed_file,NI=ni, ITS=its, OTS=ots,OIMGSETFILES=oimgsetfiles,ODIR=out_dir,$
      DATATYPE=dtype,SAMPLE=sample,LINE=line,GEOTIFF=geotiff,PBAR=cgProgressBar1,PREF=tiles[i_tiles];,MASKFILE=maskfile
      
    ;cgProgressBar1 -> Update, i_tiles/n_elements(tiles)*100
    obj_destroy,cgProgressBar1
  endfor

  return
  ;********************************************************************************************************;
 HANTS_TIME:
  ni=26
  per=365
  its=INDGEN(26)*16+1
  ots=INDGEN(46)*8+1
  nf=3
  nf2=5
  pa1=[360,180,120]
  pa2=[720,360,180,120,90]
  HiLo='Low'
  low=0
  high=10000
  fet=500
  dod=5
  y=[236,97,193,217,1241,1230,1233,1313,2967,3506,6398,7721,8072,7759,7620,7766,4953,5456,5456,2075,1866,1475,288,323,264,498]

  looplim= 10
  delta=0.1
  
  ;for i=0,8 do begin
  ta=[1,10,100,1000,10000]
  time=fltarr(5)
  for t=0,4 do begin
    profiler,/RESET
    profiler,'HANTS'
  FOR i=0,ta[t] do BEGIN
        
    res=HANTS(ni, per, its, nf2, pa2, HiLo, low, high, fet, dod, y, looplim, delta)
;    ;
;    cgplot,its,y,psym=-6,COLOR='Black',linestyle=2, thick=1,XRANGE=[-10,365],yrange=[-3500,10000]
;    cgplot,its,res2,/overplot, Color='red',linestyle=0, thick=1
;    cgplot,ots,res4,/overplot, Color='dodger blue',linestyle=2, thick=1
;    ;cgplot,its,res.yr1,/overplot, Color='GREEN',linestyle=2, thick=1
;    ;cgplot,ts,res1.yr,/overplot,color='green',linestyle=4,thick=2
;
;    al_Legend, ['Raw', 'orignal HANTS','Modified'], PSym=[-6,-15,-16], $
;      LineStyle=[1,0,2], Color=['Black','red','dodger blue'], Position=[5,-500]
;    filename='P:\TUD\idl_projects\TimeSR\test'+STRING(i,format='(I03)')+'.jpg'
;    WRITE_JPEG, filename, TVRD(/TRUE), /TRUE
;    cgdelete,/all
;    WDELETE
  ENDFOR
  profiler,/REPORT,data=data
  time[t]=data.time
  endfor
  cgplot,alog10(ta),time,psym=-6,ytitle='Computation time for total arrays(s)',xtitle='size of input arrays(alog10)'
END
