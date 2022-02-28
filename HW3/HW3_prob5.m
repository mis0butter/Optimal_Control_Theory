
%% prob 5 

clear; 

syms x1 x2 x3 lmda
x = [x1; x2; x3]; 

% equations and performance index 
f = 2*x1 + 2*x2^2 + 2*x2*x3 + 4*x3^2; 
h = x1^2 + x2^2 + x3^2 - 1; 
l_eq = f + lmda * h; 

% partial derivatives 
grad_l_eq = gradient(l_eq); 
eq1 = grad_l_eq(1); 
eq2 = grad_l_eq(2); 
eq3 = grad_l_eq(3); 
eq4 = grad_l_eq(4); 

sol1 = solve(eq1, eq2, eq3, eq4); 

% compute hessians 
Hf = hessian(f); 
Hh = hessian(h); 
L = Hf + lmda * Hh; 
L_fn = matlabFunction(L); 

% find positive semidef L 
for i = 1:numel(sol1.lmda)
    
    sprintf('Soln = %d. x1 = %.3g, x2 = %.3g, x3 = %.3g, lmda = %.3g', ... 
        i, double(sol1.x1(i)), double(sol1.x2(i)), double(sol1.x3(i)), double(sol1.lmda(i))) 
    
    lmda_dbl = double(sol1.lmda(i)); 
    L_dbl(:,:,i) = L_fn(lmda_dbl); 
    eig_L = round(eig(L_dbl(:,:,i)), 8); 
    
    % pos, semi, neg, or in-definite 
    if all(eig_L >= 0)
        disp('Local min found!') 
        sprintf('eigvals = [ %.3g, %.3g, %.3g ]', eig_L(1), eig_L(2), eig_L(3)) 

    elseif all(eig_L <= 0)
        sprintf('Soln %d negative semidefinite', i)
        sprintf('eigvals = [ %.3g, %.3g, %.3g ]', eig_L(1), eig_L(2), eig_L(3))

    elseif all(eig_L < 0)
        sprintf('Soln %d negative definite', i)
        sprintf('eigvals = [ %.3g, %.3g, %.3g ]', eig_L(1), eig_L(2), eig_L(3))

    else
        sprintf('Soln %d indefinite', i)
        sprintf('eigvals = [ %.3g, %.3g, %.3g ]', eig_L(1), eig_L(2), eig_L(3))
        
        % check on tangent space 
        grad_h = gradient(h); 

        % get zLz 
        z = sym('z', [3 1]); 
        zLz = z.' * L * z; 
        zLz_fn = matlabFunction(zLz); 
        zLz_soln = subs(zLz, lmda, sol1.lmda(i)); 

        % substitute nullspace constraint 
        zLz_null_soln = subs(zLz_soln, z(3), -z(1)-z(2)) 
    
    end 
    
end 