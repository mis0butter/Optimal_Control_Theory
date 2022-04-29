% HW 6 

clear; clc 
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'defaultLineLineWidth', 2)

%% Problem 1.2

% initial states 
x0 = 4; 
v0 = 0.5; 
T = 1.2; 

% ------------------------------------------------------------------------
% Problem 1.2a: intercept problem 

% find sigma x 
x_T = 0.2; 
sig_x_i = ( (v0*T + x0)/x_T - 1 ) * 6/T^3; 
sig_v_i = 0; 

% xstar, vstar, ustar trajectories 
dt = 0.01; 
t_ir = [0 : dt : T]; 
xstar_i = 1/12 * sig_x_i * x_T * t_ir.^3 ... 
    - 1/4 * sig_x_i * x_T * T * t_ir.^2 ... 
    + v0*t_ir + x0 ; 
vstar_i = 1/4 * sig_x_i * x_T * t_ir.^2 ... 
    - 1/2 * sig_x_i * x_T * T * t_ir ... 
    + v0; 
ustar_i = 1/2 * ( sig_x_i * x_T * t_ir ... 
    - sig_x_i * x_T * T ); 

% control effort 
ustar2_i = ustar_i.^2; 
ustar2_i_T = trapz(ustar2_i) * dt; 

ftitle = 'Intercept Problem'; 
figure('name', ftitle) 
    subplot(3,1,1) 
        plot(t_ir, xstar_i); 
        title('x*(t)')
    subplot(3,1,2) 
        plot(t_ir, vstar_i);
        title('v*(t)') 
    subplot(3,1,3) 
        plot(t_ir, ustar_i); 
        title('u*(t)') 
        xlabel('Time') 
    sgtitle({ftitle; sprintf('Control effort = %.6g', ustar2_i_T)}); 

% ------------------------------------------------------------------------
% Problem 1.2b: rendezvous problem 

% find sigma x, v 
x_T = 0.2; 
v_T = -0.1; 
A = [   1/6*x_T*T^3     1/4*v_T*T^2 ; 
        1/4*x_T*T^2     1/2*v_T*T   ];
b = [   v0*T + x0 - x_T; 
        v0 - v_T]; 
sig_xv = A^-1*b; 
sig_x_r = sig_xv(1); 
sig_v_r = sig_xv(2); 

% xstar, vstar, ustar trajectories 
xstar_r = 1/12 * sig_x_r * x_T * t_ir.^3 ... 
    - 1/4 * sig_x_r * x_T * T * t_ir.^2 ... 
    - 1/4 * sig_v_r * v_T * t_ir.^2 ... 
    + v0*t_ir + x0; 
vstar_r = 1/4 * sig_x_r * x_T * t_ir.^2 ... 
    - 1/2 * sig_x_r * x_T * T * t_ir ... 
    - 1/2 * sig_v_r * v_T * t_ir + v0; 
ustar_r = 1/2 * ( ( sig_x_r * x_T ) * t_ir ... 
    - sig_x_r * x_T * T ... 
    - sig_v_r * v_T ); 

% control effort 
ustar2_r = ustar_r.^2; 
ustar2_r_T = trapz(ustar2_r) * dt; 

ftitle = 'Rendezvous Problem'; 
figure('name', ftitle) 
    subplot(3,1,1) 
        plot(t_ir, xstar_r); 
        title('x*(t)')
    subplot(3,1,2) 
        plot(t_ir, vstar_r);
        title('v*(t)') 
    subplot(3,1,3) 
        plot(t_ir, ustar_r); 
        title('u*(t)') 
        xlabel('Time') 
    sgtitle({ftitle; sprintf('Control effort = %.6g', ustar2_r_T)}); 

    
%% Problem 1.3: Finite-Horizon LQR approach: intercept 

% dynamics 
A = [0 1; 0 0 ]; 
B = [0 1]'; 

% LQR matrices 
M = [sig_x_i, 0; 0, sig_v_i]; 
Q = zeros(2); 
R = 1;
    
% intercept 
[t_i, u_i, x_i, P_i, H_i] = intRiccati(M, A, B, Q, R, T, x0, v0); 

% rendezvous 
M = [sig_x_r, 0; 0, sig_v_r]; 
[t_r, u_r, x_r, P_r, H_r] = intRiccati(M, A, B, Q, R, T, x0, v0); 



    % set ode45 params 
    rel_tol = 1e-14;         % 1e-14 accurate; 1e-6 coarse 
    abs_tol = 1e-14; 
    options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

    % solve matrix Riccati ODE (backwards) 
    Pf = M; 
    [t, Prow] = ode45(@(t, P) mRiccatiEq(t, P, A, B, Q, R), [T:-dt:0], Pf, options); 
    t = flip(t); 
    Prow = flip(Prow); 
    P0 = Prow(1,:)'; 

    % get H, x, p, u 
    xk = [x0; v0]; 
    for i = 1:length(t)

        xkm1 = xk; 
        P = Prow(i,:); 
        P = reshape(P, size(A)); 

        K(i,:) = inv(R) * B' * P; 

    end 

        tk = t; 
        X0 = [x0; v0; P0]; 
        [tk, xk] = ode45(@(t, X) mRiccatiEq_aug(t, X, A, B, Q, R), tk, X0, options);     
%  dXdt = mRiccatiEq_aug(t, X, A, B, Q, R)
 
n = 4; p = 2; 
figure()

    subplot(n,p,1) 
        plot(t_ir, ustar_i, t_i, u_i, '--')
        title({'Intercept'; 'u*(t)'}) 
        legend('classical', 'Riccati', 'location', 'best')
    subplot(n,p,3) 
        plot(t_ir, vstar_i, t_i, x_i(:,2), '--') 
        title('v*(t)') 
    subplot(n,p,5) 
        plot(t_ir, xstar_i, t_i, x_i(:,1), '--') 
        title('x*(t)') 
    subplot(n,p,7)
        plot(t_i, H_i)
        title('H*(t)')
        
    subplot(n,p,2) 
        plot(t_ir, ustar_r, t_r, u_r, '--')
        title({'Rendezvous'; 'u*(t)'}) 
        legend('classical', 'Riccati', 'location', 'best')
    subplot(n,p,4) 
        plot(t_ir, vstar_r, t_r, x_r(:,2), '--')
        title('v*(t)') 
    subplot(n,p,6) 
        plot(t_ir, xstar_r, t_r, x_r(:,1), '--', tk, xk(:,1), ':') 
        title('x*(t)') 
    subplot(n,p,8)
        plot(t_r, H_r)
        title('H*(t)')

        
        
        

%% subfunctions 

function dx = intAtilde(t, x, tk, K, A, B)

ik = find(tk == t); 
Kk = K(ik,:); 

% Atilde = [0 1; -K(1), -K(2)]; 
Atilde = A - B*Kk; 

dx = Atilde*x; 

end 

function [t, u, x, p, H] = intRiccati(M, A, B, Q, R, T, x0, v0)

    % set ode45 params 
    rel_tol = 1e-10;         % 1e-14 accurate; 1e-6 coarse 
    abs_tol = 1e-10; 
    options = odeset('reltol', rel_tol, 'abstol', abs_tol ); 

    % solve matrix Riccati ODE (backwards) 
    Pf = M; 
    [t, Prow] = ode45(@(t, P) mRiccatiEq(t, P, A, B, Q, R), [T 0], Pf, options); 
    t = flip(t); 
    Prow = flip(Prow); 

    % get H, x, p, u 
    xk = [x0; v0]; 
    for i = 2:length(t)

        xkm1 = xk; 
        P = Prow(i,:); 
        P = reshape(P, size(A)); 

        K = inv(R) * B' * P; 

        Atilde = [0 1; -K(1), -K(2)]; 
        if i == 1
            [~, xk] = ode45(@(t, x) Atilde*x, [t(i) t(i+1)], xkm1); 
        else
            [~, xk] = ode45(@(t, x) Atilde*x, [t(i-1) t(i)], xkm1); 
        end 
        xk = xk(end,:)'; 

        pk = -2*P*xk; 
        uk = - K * xk; 
        H(i,:) = pk' * [A*xk + B*uk] - xk'*Q*xk - uk'*R*uk;

        x(i,:) = xk; 
        p(i,:) = pk; 
        u(i,:) = uk; 

    end 

end 

function dPdt = mRiccatiEq(t, P, A, B, Q, R)

P    = reshape(P, size(A)); 
dPdt = - ( P*A + A'*P + Q - P*B*inv(R)*B'*P ); 
dPdt = dPdt(:); 

end 


function dXdt = mRiccatiEq_aug(t, X, A, B, Q, R)

x = X(1:2); 
P = X(3:end); 
P = reshape(P, size(A)); 

K = inv(R) * B' * P; 
Atilde = A - B*K; 
dxdt = Atilde*x; 

dPdt = - ( P*A + A'*P + Q - P*B*inv(R)*B'*P ); 
dPdt = dPdt(:); 

dXdt = [dxdt; dPdt]; 

end 






