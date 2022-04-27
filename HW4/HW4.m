clear 

syms y yp ypp 

F = y * sqrt( 1 + yp^2 ); 

F_y = diff(F, y); 
F_yp = diff(F, yp); 
F_yp_y = diff(F_yp, y); 
F_yp_yp = diff(F_yp, yp); 

ddx_F_yp = F_yp_y * yp + F_yp_yp * ypp; 

EL = F_y - ddx_F_yp; 

%% 
%    x B A 
syms A B C 

y0 = C * cosh( A/C - B/C) 
y0p = diff(y0, A) 
y0pp = diff(y0p, A)

EL_subs = subs(EL, y, y0); 
EL_subs = subs(EL_subs, yp, y0p); 
EL_subs = subs(EL_subs, ypp, y0pp)

simplify(EL_subs)









