% Optimal control theory 
% HW 2
% Junette Hsin 

%% Problem 1 

% inputs Q and b 
Q = [2 0; 0 2]; 
b = [3; 5];

% g function 
g = @(x) Q*x - b; 

% error threshold 
e = 1e-3; 
delta = 1; 

% FIRST GUESS 
disp('FIRST GUESS')
x0 = [13; -2]; 
[x_arr, i] = min_perf(delta, e, x0, Q, b, g); 
x_arr(end,:)
i 

% SECOND GUESS 
disp('SECOND GUESS') 
x0 = [-10; 7]; 
[x_arr, i] = min_perf(delta, e, x0, Q, b, g); 
x_arr(end,:)
i 

% THIRD GUESS 
disp('THIRD GUESS') 
x0 = [-2; 14]; 
[x_arr, i] = min_perf(delta, e, x0, Q, b, g); 
x_arr(end,:)
i 


% plot(x_arr(:,1), x_arr(:,2))

%% subfunctions 

function [x_arr, i] = min_perf(delta, e, x0, Q, b, g) 

i = 0; 
x_arr = x0'; 
xkp1  = x0; 

% start iterative process 
while delta > e
    
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

