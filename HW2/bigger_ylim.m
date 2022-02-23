function bigger_ylim
% Increase y-axis limits on plot by 30% on current axes 

    ylimits     = get(gca, 'ylim'); 
    yrange      = ylimits(2) - ylimits(1);
    new_ylim    = [ ylimits(1) - 0.15*yrange, ylimits(2) + 0.15*yrange ]; 
    set(gca, 'ylim', new_ylim);

end