vth = 30.0
npart = 10000.

sgn = fltarr(npart)

;rnd = randomu(seed,npart)
;whp = where(rnd ge 0.5)
;sgn(whp) = 1.0
;whn = where(rnd lt 0.5)
;sgn(whn) = -1.0
vx = vth*sqrt(-alog(randomu(seed,npart)))*cos(!pi*randomu(seed,npart))


;rnd = randomu(seed,npart)
;whp = where(rnd ge 0.5)
;sgn(whp) = 1.0
;whn = where(rnd lt 0.5)
;sgn(whn) = -1.0
vy = vth*sqrt(-alog(randomu(seed,npart)))*cos(!pi*randomu(seed,npart))



;rnd = randomu(seed,npart)
;whp = where(rnd ge 0.5)
;sgn(whp) = 1.0
;whn = where(rnd lt 0.5)
;sgn(whn) = -1.0
vz = vth*sqrt(-alog(randomu(seed,npart)))*cos(!pi*randomu(seed,npart))


;theta = !PI*randomu(seed,npart) - !pi/2
;phi = 2*!PI*randomu(seed,npart)

;vx = v*sin(theta)*cos(phi) 
;vy = v*sin(theta)*sin(phi)
;vz = v*cos(theta)


plot,vx,vz,psym=3,xrange=[-60,60],yrange=[-60,60],/isotropic
p = plot3d(vx,vy,vz,'r.')
end

