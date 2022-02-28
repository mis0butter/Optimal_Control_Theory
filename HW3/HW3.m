x = sym('x', [2 1]); 

f = 1/2 * (x(1)^2 + 1/4 * x(2)^2); 

%% 

cvx_begin

    variable x1 
    variable x2 

    minimize( 2*x(1) + 2*x(2)^2 )
%     minimize( norm( Ustar ) )

    subject to
        x1 - x2 == 0; 
            
cvx_end

%% 

cvx_begin 

    variable x1 
    variable x2 
    variable x3 
    
    minimize ( 2*x1 + 2*x2^2 + 2*x2*x3 + 4*x3^2 )
    
    subject to 
    
        x1^2 + x2^2 + x3^2 - 1 == 0; 
        
cvx_end 
    
%% 

syms a b c d
eq1 = 64-2*a-b-6*d;
eq2 = 37-5*a-2*b;
eq3 = 66-7*b-2*c-2*d;
eq4 = 104-8*c-9*d;
sol = solve(eq1,eq2,eq3,eq4);

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
    
    lmda_dbl = double(sol1.lmda(i)); 
    L_dbl(:,:,i) = L_fn(lmda_dbl); 
    eig_L = eig(L_dbl(:,:,i)); 
    
    % pos, semi, neg, or in-definite 
    if all(eig_L > 0)
        
        disp('Local min found!') 
        sprintf('Soln = %d. x1 = %.3g, x2 = %.3g, x3 = %.3g, lmda = %.3g', ... 
            i, double(sol1.x1(i)), double(sol1.x2(i)), double(sol1.x3(i)), double(sol1.lmda(i))) 

    elseif all(eig_L < 0)
        sprintf('Soln %d negative definite', i)

    else
        sprintf('Soln %d indefinite', i)
        
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


%% prob 6 

syms x1 x2 l1 l2 

eq1 = 2 + 2*l1*(x1-3) + 2*l2*(x1-4); 
eq2 = 1 + 2*l1*x2 + 2*l2*x2; 
eq3 = (x1-3)^2 + x2^2 - 9; 
eq4 = (x1-4)^2 + x2^2 - 16; 

sol1 = solve(eq1, eq2, eq3, eq4); 

disp('4 systems of equations:')
if isempty(sol1.x1) 
    disp('Empty solution set') 
else 
    for i = 1:numel(sol1.x1)
        sprintf('Soln = %d. x1 = %.3g, x2 = %.3g, lmda1 = %.3g, lmda2 = %.3g', ... 
            i, double(sol1.x1(i)), double(sol1.x2(i)), double(sol1.l1(i)), double(sol1.l2(i))) 
    end 
end 

sol2 = solve(eq3, eq4); 

disp('2 systems of equations:')
if isempty(sol2.x1) 
    disp('Empty solution set') 
else 
    for i = 1:numel(sol2.x1)
        sprintf('Soln = %d. x1 = %.3g, x2 = %.3g', ... 
            i, double(sol2.x1(i)), double(sol2.x2(i))) 
    end 
end 








    