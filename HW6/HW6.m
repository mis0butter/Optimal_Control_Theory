% HW 6 

set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

%% Problem 2

% initial states 
x0 = 4; 
v0 = 0.5; 
T = 1.2; 

% ------------------------------------------------------------------------
% Problem 2a: intercept problem 

% find sigma x 
x_T = 0.2; 
sig_x = ( (v0*T + x0)/x_T - 1 ) * 6/T^3; 

% xstar, vstar, ustar trajectories 
dt = 0.01; 
t = [0 : dt : T]; 
xstar = 1/12 * sig_x * x_T * t.^3 ... 
    - 1/4 * sig_x * x_T * T * t.^2 ... 
    + v0*t + x0 ; 
vstar = 1/4 * sig_x * x_T * t.^2 ... 
    - 1/2 * sig_x * x_T * T * t ... 
    + v0; 
ustar = 1/2 * ( sig_x * x_T * t ... 
    - sig_x * x_T * T ); 

% control effort 
ustar2 = ustar.^2; 
ustar2_T = trapz(ustar2) * dt; 

ftitle = 'Intercept Problem'; 
figure('name', ftitle) 
    subplot(3,1,1) 
        plot(t, xstar); 
        title('x*(t)')
    subplot(3,1,2) 
        plot(t, vstar);
        title('v*(t)') 
    subplot(3,1,3) 
        plot(t, ustar); 
        title('u*(t)') 
        xlabel('Time') 
    sgtitle(ftitle); 

% ------------------------------------------------------------------------
% Problem 2b: rendezvous problem 

% find sigma x, v 
x_T = 0.2; 
v_T = -0.1; 
A = [   1/6*x_T*T^3     1/4*v_T*T^2 ; 
        1/4*x_T*T^2     1/2*v_T*T   ];
b = [   v0*T + x0 - x_T; 
        v0 - v_T]; 
sig_xv = A^-1*b; 
sig_x = sig_xv(1); sig_v = sig_xv(2); 

% xstar, vstar, ustar trajectories 
ustar = 1/2 * ( ( sig_x * x_T ) * t ... 
    - sig_x * x_T * T ... 
    - sig_v * v_T ); 
vstar = 1/4 * sig_x * x_T * t.^2 ... 
    - 1/2 * sig_x * x_T * T * t ... 
    - 1/2 * sig_v * v_T * t + v0; 
xstar = 1/12 * sig_x * x_T * t.^3 ... 
    - 1/4 * sig_x * x_T * T * t.^2 ... 
    - 1/4 * sig_v * v_T * t.^2 ... 
    + v0*t + x0; 

% control effort 
ustar2 = ustar.^2; 
ustar2_T = trapz(ustar2) * dt; 

ftitle = 'Rendezvous Problem'; 
figure('name', ftitle) 
    subplot(3,1,1) 
        plot(t, xstar); 
        title('x*(t)')
    subplot(3,1,2) 
        plot(t, vstar);
        title('v*(t)') 
    subplot(3,1,3) 
        plot(t, ustar); 
        title('u*(t)') 
        xlabel('Time') 
    sgtitle(ftitle); 











