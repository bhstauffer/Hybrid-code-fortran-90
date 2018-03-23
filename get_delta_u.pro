@get_const
uout = 200e3
rho = mp*0.18e6
a3a2 = 15/1.8
T = 6e-13
c = T*a3a2/rho
b = 2*uout

du = (-b + sqrt(b^2 - 4*c))/2

print,du

;du = T*a3a2/rho/(2*uout)
;print,du

end
