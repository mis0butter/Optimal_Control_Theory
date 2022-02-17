% Optimal control theory 
% HW 2
% Junette Hsin 

%% Problem 1 

clear; 

% inputs Q and b 
Q = [2 0; 0 2]; 
b = [3; 5];

% g function 
g = @(x) Q*x - b; 

% error threshold 
err   = 1e-3; 
delta = 1; 

% FIRST GUESS 
disp('FIRST GUESS')
x0 = [13; -2]; 
[x_arr, i] = min_perf(delta, err, x0, Q, b, g); 
x_arr(end,:)
i 

% SECOND GUESS 
disp('SECOND GUESS') 
x0 = [-10; 7]; 
[x_arr, i] = min_perf(delta, err, x0, Q, b, g); 
x_arr(end,:)
i 

% THIRD GUESS 
disp('THIRD GUESS') 
x0 = [-2; 14]; 
[x_arr, i] = min_perf(delta, err, x0, Q, b, g); 
x_arr(end,:)
i 

% plot(x_arr(:,1), x_arr(:,2))

%% Problem 2 

clear; 

x = sym('x', [2 1]); 

% create performance index functions 
f = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2; 
f = matlabFunction(f); 
h = (x(1) + 0.5)^2 + (x(2) + 0.5)^2 - 0.25; 
h = matlabFunction(h); 

% initialize 
b    = 6;
a0   = 0.1; 
ak   = a0 / b; 
x0   = [1; 1]; 
xkp1 = x0; 
err  = 10e-4; 
i    = 0; 
h_err = h(xkp1(1), xkp1(2));

% iterate 
while h_err > err 
    
    i = i + 1; 
    
    xk   = xkp1; 
    ak   = b * ak; 

    phi  = @(x) f(x(1), x(2)) + 1/2 * ak * h(x(1), x(2))^2; 
    xkp1 = fminsearch(phi, xk); 
    
    h_err = h(xkp1(1), xkp1(2)); 
    
end 

%% Problem 3 

% initialize 
b    = 6;
a0   = 0.1; 
ak   = a0 / b; 
x0   = [1; 1]; 
xkp1 = x0; 
err  = 10e-4; 
i    = 0; 
h_err = h(xkp1(1), xkp1(2));
lmda0 = 10; 
lmdak = lmda0; 


% iterate 
while h_err > err 
    
    i = i + 1; 
    
    xk   = xkp1; 
    ak   = b * ak; 

    phi = @(x) f(x(1), x(2)) + lmdak * h(x(1), x(2)) + 1/2 * ak * h(x(1), x(2))^2; 
    xkp1 = fminsearch(phi, xk); 
    
    h_err = h(xkp1(1), xkp1(2)); 
    
end 


%% subfunctions 

function [x_arr, i] = min_perf(delta, err, x0, Q, b, g) 
% Minimize performance index 

% initialize 
i = 0; 
x_arr = x0'; 
xkp1  = x0; 

% start iterative process 
while delta > err
    
    % index loops 
    i = i + 1; 
    
    xk = xkp1; 
%     ak = ( 1/2 * xk' * Q * g(xk) + 1/2 * g(xk)' * Q * xk - g(xk)' ) * ... 
%         ( 1/2 * g(xk)' * Q * g(xk) + 1/2 * g(xk)' * Q * g(xk) )^-1; 
    ak = ( 1/2 * g(xk)' * Q * xk + 1/2 * xk' * Q * g(xk) - g(xk)' * b) * ... 
        ( g(xk)' * Q * g(xk) )^-1; 
    xkp1 = xk - ak * g(xk); 
    delta = norm(xkp1 - xk); 
    
    if isnan(delta) 
        disp('Contains NaNs') 
        break
    end 
    
    x_arr = [x_arr; xkp1']; 

end 

end 

