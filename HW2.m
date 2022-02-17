% Optimal control theory 
% HW 2
% Junette Hsin 

% inputs Q and b 
Q = [2 0; 0 2]; 
b = [3; 5];

% first guess 
x0 = [1; 2]; 
xkp1 = x0; 

% g function 
g = @(x) Q*x - b; 

% error threshold 
e = 1e-1; 
delta = 1; 

% start iterative process 
while delta > e
    
    xk = xkp1; 
    ak = ( 1/2 * xk' * Q * g(xk) + 1/2 * g(xk)' * Q * xk - g(xk)' ) * ... 
        ( 1/2 * g(xk)' * Q * g(xk) + 1/2 * g(xk)' * Q * g(xk) )^-1; 
    xkp1 = xk - ak * g(xk); 
    delta = norm(xkp1 - xk); 

end 

