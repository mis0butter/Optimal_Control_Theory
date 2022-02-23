% Optimal control theory 
% HW 2
% Junette Hsin 

%% Problem 1 

clear; 

% inputs Q and b 
Q = [12 1; 1 2]; 
b = [3; 2];

% g function 
g = @(x) Q*x - b; 

% error threshold 
err   = 1e-3; 
delta = 1; 

% FIRST GUESS 
disp('FIRST GUESS')
x0 = [13; -2]; 
[x_arr, i] = min_perf(delta, err, x0, Q, b, g); 

    % plot 
    fname = 'P1 - Minimize Performance Index (Case 1)'; 
    plot_x1x2(fname, x_arr, i)

% SECOND GUESS 
disp('SECOND GUESS') 
x0 = [-10; 7]; 
[x_arr, i] = min_perf(delta, err, x0, Q, b, g); 

    % plot 
    fname = 'P1 - Minimize Performance Index (Case 2)'; 
    plot_x1x2(fname, x_arr, i)

% THIRD GUESS 
disp('THIRD GUESS') 
x0 = [-2; 14]; 
[x_arr, i] = min_perf(delta, err, x0, Q, b, g); 

    % plot 
    fname = 'P1 - Minimize Performance Index (Case 3)'; 
    plot_x1x2(fname, x_arr, i)

%% Problem 2 

clear; 

x = sym('x', [2 1]); 

% create performance index functions 
f = 100 * ( x(2) - x(1)^2 )^2 + ( 1 - x(1) )^2; 
f = matlabFunction(f); 
h = ( x(1) + 0.5 )^2 + ( x(2) + 0.5 )^2 - 0.25; 
h = matlabFunction(h); 

% initialize 
b    = 6;
a0   = 0.1; 
x0   = [1; 1]; 
err  = 10^-4; 

% 0 iteration 
j    = 0; 
akm1 = a0; 
xkm1 = x0; 
h_err = h(xkm1(1), xkm1(2));
x_arr = x0'; 

% iterate 
while h_err > err 
    
    % current index 
    j  = j + 1;  
    ak = b * akm1; 

    phi = @(x) f(x(1), x(2)) + 1/2 * ak * h(x(1), x(2))^2; 
    xk  = fminsearch(phi, xkm1); 
    
    % new penalty 
    h_err = norm(h(xk(1), xk(2))); 
    
    % save output 
    x_arr = [x_arr; xk']; 
    
    % set up next index 
    akm1 = ak; 
    xkm1 = xk; 
    
end 

% plot 
fname = 'P2 - Equality Constrained Optimization (Penalty)'; 
plot_x1x2(fname, x_arr, j, h_err)


%% Problem 3 

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
x0   = [1; 1]; 
err  = 10^-4; 
k    = 0; 
lmda0 = 10; 

% first iteration 
akm1 = a0; 
xkm1 = x0; 
lmdakm1 = lmda0; 
h_err = h(xkm1(1), xkm1(2));
x_arr = x0'; 

% iterate 
while h_err > err 
    
    % current index 
    k     = k + 1; 
    ak    = b * akm1; 
    lmdak = lmdakm1 + akm1 * h(xkm1(1), xkm1(2)); 

    phi = @(x) f(x(1), x(2)) + lmdak * h(x(1), x(2)) + 1/2 * ak * h(x(1), x(2))^2; 
    xk  = fminsearch(phi, xkm1); 
    
    % penalty 
    h_err = norm(h(xk(1), xk(2))); 
    
    % save output 
    x_arr = [x_arr; xk']; 
    
    % set up next index 
    akm1    = ak; 
    lmdakm1 = lmdak; 
    xkm1    = xk; 
    
end 

% plot 
fname = 'P3 - Equality Constrained Optimization (Lagrange Multiplier)'; 
plot_x1x2(fname, x_arr, k, h_err)

%% save as pdfs 

fig_h = findall(0,'Type','figure'); 

for i = 1:numel(fig_h)
    
    fname = fig_h(i).Name; 
    savePDF( fig_h(i), [fname '.pdf'], '.'); 
%     savePDF( '', ['' '.pdf'] ); 
    
end 

%  % Create filename 
%  fn = '.';  %in this example, we'll save to a temp directory.
%  
%  % Save first figure
%  export_fig(fn, '-pdf', figHandles(1))
%  
%  % Loop through figures 2:end
%  for i = 2:numel(figHandles)
%      export_fig(fn, '-pdf', figHandles(i), '-append')
%  end

%% subfunctions 

function [x_arr, i] = min_perf(delta, err, x0, Q, b, g) 
% Minimize performance index 

% initialize 
i     = 0; 
x_arr = x0'; 
xkm1  = x0; 

% start iterative process 
while delta > err
    
    % current index 
    i = i + 1; 
    
    % calc ak  
%     ak = ( 1/2 * g(xkm1)' * Q * xkm1 + 1/2 * xkm1' * Q * g(xkm1) - g(xkm1)' * b) * ... 
%         ( g(xkm1)' * Q * g(xkm1) )^-1;     
    ak = inv(g(xkm1)' * Q* g(xkm1)) * (xkm1' * Q * g(xkm1) - g(xkm1)' * b ); 
    
    % calc xk 
    xk = xkm1 - ak * g(xkm1); 
    delta = norm(xkm1 - xk); 
    
    if isnan(delta) 
        disp('Contains NaNs') 
        break
    end 
    
    % save output 
    x_arr = [x_arr; xk']; 
    
    % set up next index 
    xkm1 = xk;

end 

end 

% ------------------------------------------------------------------------ 

function plot_x1x2(fname, x_arr, k, h_err)

if ~exist('h_err', 'var') 
    h_err = NaN; 
end 
% plot x1-x2 plane 

figure('name', fname, 'position', [100 100 600 500])
    subplot(3,1,1:2) 
        hold on; grid on; 
        scatter(x_arr(1,1), x_arr(1,2), 40, 'linewidth', 2); 
        scatter(x_arr(2:end-1,1), x_arr(2:end-1,2), 10, 'filled')
        scatter(x_arr(end,1), x_arr(end,2), 40, 'linewidth', 2); 
        plot(x_arr(:,1), x_arr(:,2), '--k'); 
        bigger_lim; 

        xlabel('x1'); ylabel('x2'); 
        legend('start', 'middle', 'end', 'location', 'eastoutside')
        sgtitle( fname ); 
    subplot(3,1,3) 
        pos = get(gca, 'position'); 
%         text = {'test1'; 'test2'; 'test3'}; 
        
        text = { sprintf('Iterations = %d', k); 
            sprintf('h_{final} = %.4g', h_err); 
            sprintf( 'start = [%.4g , %.4g]', x_arr(1,1), x_arr(1,2) ); 
            sprintf( 'end = [%.4g , %.4g]', x_arr(end,1), x_arr(end,2) ) }; 
        
        annotation('textbox', pos, ...
          'String', text, ...
          'edgecolor', 'none');
        axis off 
        
end 

function bigger_lim
% Increase y-axis limits on plot by 30% on current axes 

    ylims       = get(gca, 'ylim'); 
    yrange      = ylims(2) - ylims(1);
    new_ylim    = [ ylims(1) - 0.15*yrange, ylims(2) + 0.15*yrange ]; 
    set(gca, 'ylim', new_ylim);
    
    xlims       = get(gca, 'xlim'); 
    xrange      = xlims(2) - xlims(1); 
    new_xlim    = [ xlims(1) - 0.15*xrange, xlims(2) + 0.15*xrange ]; 
    set(gca, 'xlim', new_xlim); 

end

