
%% Problem 2

clear; 

% initial condition
x0 = [0; 0; 1; 0];

% constants 
M = 2;
m = 1;
L = 0.5;
g = 9.81;
A = [0,     1,      0,              0;
     0,     0,      -(m*g)/M,       0;
     0,     0,      0,              1;
     0,     0,      g*(M+m)/(M*L),  0];
B = [ 0; 1/M; 0; -1/(M*L)];
C = [1,     0,      0,      0;
     0,     0,      1,      0];
D = 0;

% ------------------------------------------------------------------------
% Part 2.1

% Q and R
Q = diag([1, 1, 5, 5]);
R = 1;

% K infinity 
[K_LQR_1, S_1, eVal_1] = lqr(A,B,Q,R);

Atilde = (A - B*K_LQR_1);

disp("Part 1:");
K_LQR_1
disp("EigenValues of Atilde: " );
eVal_1

% ------------------------------------------------------------------------
% Part 2.2

% Define Q and R
Q = diag([1, 1, 10, 10]);
R = 15;


[K_LQR_2, S_2, eVal_2] = lqr(A,B,Q,R);

Atilde    = (A - B*K_LQR_2);

disp("Part 2:");
K_LQR_2
disp("EigenValues of Atilde: " );
eVal_2


% define time interval
t    = 0:0.01:40;
X = [];
u = [];

% propagate state - since system is LTI, we can find the paticular solution
% such that x(t) = PHI(t,t0)*x0 so use the state transition matrix to
% propage the state forward and also solve for the associated control input
for ii=1:numel(t)
    xk = expm(Atilde*t(ii))*x0;
    uk = -K_LQR_1 * xk;

    X = [X xk];
    u = [u; uk];
end

stitle = {'x','dx', '\theta','d\theta'};

ftitle = 'Problem 3 - Infinite Horizon LQR'; 
figure('name', ftitle, 'position', [100 100 700 700]);

for jj=1:length(x0)
    ax = subplot(6,1,jj);
    plot(t,X(jj,:))
    title(stitle{jj})
end
subplot(6,1,5)
    plot(t, u); 
    xlabel('Time (s)')

subplot(6,1,6) 
    txt = { sprintf('2.1: K_{LQR} = [%.5g, %.5g, %.5g, %.5g]',  ... 
            K_LQR_1(1), K_LQR_1(2), K_LQR_1(3), K_LQR_1(4)); ... 
            sprintf('eig(Atilde) = [%.5g + j%.5g, %.5g + j%.5g, %.5g + j%.5g, %.5g + j%.5g]', ... 
            real(eVal_1(1)), imag(eVal_1(1)), real(eVal_1(2)), imag(eVal_1(2)), real(eVal_1(3)), imag(eVal_1(3)), real(eVal_1(4)), imag(eVal_1(4)) ) ; ... 
            ''; 
            sprintf('2.2: K_{LQR} = [%.5g, %.5g, %.5g, %.5g]',  ... 
            K_LQR_2(1), K_LQR_2(2), K_LQR_2(3), K_LQR_2(4)); ... 
            sprintf('eig(Atilde) = [%.5g + j%.5g, %.5g + j%.5g, %.5g + j%.5g, %.5g + j%.5g]', ... 
            real(eVal_2(1)), imag(eVal_2(1)), real(eVal_2(2)), imag(eVal_2(2)), real(eVal_2(3)), imag(eVal_2(3)), real(eVal_2(4)), imag(eVal_2(4)) ) ; ... 
            }; 
    plt_txt(txt) 

sgtitle(ftitle) 


