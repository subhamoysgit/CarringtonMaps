FUNCTION disc_to_helio_map,image,year,month,day,hr,mnt

mnths=['Jan','Feb','Mar','Apr','May','Jun','July','August','Sept','Oct','Nov','Dec']

;print,strtrim(year,2)+'-'+mnths(month-1)+'-'+strtrim(day,2)
b0=get_rb0p(strtrim(year,2)+'-'+mnths(month-1)+'-'+strtrim(day,2),/b0)
;print,(180./!pi)*b0
im=fltarr(360,360)
sinb=sin(b0)
cosb=cos(b0)
rr = 256 


for i=0,359 do begin
for j=0,359 do begin
lat=-90>0.5*(j-180.)<90   
lon=-90>0.5*(i-180.)<90
x_rot=round(rr*cos(!pi*lat/180.)*sin(!pi*lon/180.)) 
y_rot=round(rr*(-cos(!pi*lat/180.)*cos(!pi*lon/180.)*sinb+sin(!pi*lat/180.)*cosb))  
pix_i=0>(x_rot+256)<511
pix_j=0>(y_rot+256)<511 
im(i,j)=image(pix_i,pix_j)
endfor
endfor

return,im
end