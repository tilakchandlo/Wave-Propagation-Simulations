% Function takes inputs for dimensions and number of iterations for the
% simulation.
function waterSim;
dim = 100;
iterations = 500;
% Just going to use simple values for the constants. Also since the spaces
% between each spot on the grid is 1, dx and dx are 1. The one we need to
% determine is dt. In this part of the project, mu, epsilon, dx, and dy are
% all 1 so for the sake of simplicity they will be omitted and then
% reincluded in the later parts of the project.
% mu = 1;
% ep = 1;
% dx = 1;
% dy = 1;
dt = .4;

% First step is to set up the ez, hx, hy grids using the input and also
% setting up the matricies to be used in the for loops for updating.
ez = zeros(dim);
hx = zeros(dim, dim-1);
hy = zeros(dim-1, dim);
nextEz = ez;
nextHx = hx;
nextHy = hy;

% This main loop will drive the update of the three components that we need
% to be updated.
for k = 1:iterations

    % This is the source of the electric field at the center of the grid.
    ez(floor(dim/2), floor(dim/2)) = cos(pi.*k/20.*dt);

    % This subloop updates the electric field for every point on the grid.
    for i = 2:dim-1
        for j = 2:dim-1
            nextEz(i,j) = ez(i,j) + dt.*((hy(i-1,j) - hy(i,j)) - (hx(i,j) - hx(i,j-1)));
        end
    end
    ez = nextEz;

    % This subloop updates the x component of the magnetic field for every
    % point on the grid.
    for i = 1:size(hx,1)
        for j = 1:size(hx,2)
            nextHx(i,j) = hx(i,j) - dt.*(ez(i,j+1) - ez(i,j));
        end
    end
    hx = nextHx;

    % This subloop updates the y component of the magnetic field for every
    % point on the grid.
    for i = 1:size(hy,1)
        for j = 1:size(hy,2)
            nextHy(i,j) = hy(i,j) - dt.*(ez(i+1,j) - ez(i,j));
        end
    end
    hy = nextHy;

    % These few lines just make the animation of Ez propagating in the
    % grid.
    figure(1);
    surf(ez);
    drawnow;
end

end