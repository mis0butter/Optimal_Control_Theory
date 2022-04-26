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

x_T = 0.2; 
sig_x = ( (v0*T + x0)/x_T - 1 ) * 6/T^3; 

% xstar, vstar, ustar trajectories 
dt = 0.01; 
t = [0 : dt : T]; 
xstar = 1/12 * sig_x * x_T * t.^3 ... 
    - 1/4 * sig_x * x_T * T * t.^2 ... 
    + v0*t + x0  ; 
vstar = 1/4 * sig_x * x_T * t.^2 ... 
    - 1/2 * sig_x * x_T * T * t ... 
    + v0; 
ustar = 1/2 * ( sig_x * x_T * t ... 
    - sig_x * x_T * T ); 

% control effort 
ustar2 = ustar.^2; 
ustar2_T = trapz(ustar2) * dt; 

figure()
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

% ------------------------------------------------------------------------
% Problem 2b: rendezvous problem 


