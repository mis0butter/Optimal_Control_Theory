
%% Homework 2 Question 2

% clear all;
% close all;

% define initial condition
x0  = [0,0,1,0]';

% define constants 
M   = 2;
m   = 1;
L   = 0.5;
g   = 9.81;

% define A
A   = [0,   1,      0,              0;
       0,   0,      -(m*g)/M,       0;
       0,   0,      0,              1;
       0,   0,      g*(M+m)/(M*L),  0];

% define B
B   = [ 0;
        1/M;   
        0;  
        -1/(M*L)];

C   = [1,   0,      0,      0;
       0,   0,      1,      0];

D   = 0;

% ------------------------------------------------------------------------
% Part 2.1

% Print LQR Gain Matrix 
%  Print out the LQR gain matrix, K∞_LQR , 
%  and the eigenvalues of A ̃ = A − BK∞
%  For Q = diag(1, 1, 5, 5), R = 1

% Define Q and R
Q = diag([1, 1, 5, 5]);
R = 1;

% Find optimal gain K∞
[K_LQR, S, eVal] = lqr(A,B,Q,R);

Atilde    = (A - B*K_LQR);

disp("Part 1:");
K_LQR
disp("EigenValues of Atilde: " );
eVal

% ------------------------------------------------------------------------
% Part B

% Define Q and R
Q = diag([1, 1, 10, 10]);
R = 15;


[K_LQR, S, eVal] = lqr(A,B,Q,R);

Atilde    = (A - B*K_LQR);

disp("Part 2:");
K_LQR
disp("EigenValues of Atilde: " );
eVal


% define time interval
t    = 0:0.01:40;
XMat = [];
uMat = [];

% propagate state - since system is LTI, we can find the paticular solution
% such that x(t) = PHI(t,t0)*x0 so use the state transition matrix to
% propage the state forward and also solve for the associated control input
for ii=1:numel(t)
    xk = expm(Atilde*t(ii))*x0;
    uk = -K_LQR * xk;

    XMat = [XMat xk];
    uMat = [uMat uk];
end

ylabels = {'$x$','$\dot{x}$', '$\theta$','$\dot{\theta}$'};

figure(1);

for jj=1:size(XMat,1)
        ax = subplot(2,2,jj);
        plot(t,XMat(jj,:))
        xlabel('Time (s)')
        ylabel(ylabels{jj},'Interpreter','Latex')
        grid on;
end