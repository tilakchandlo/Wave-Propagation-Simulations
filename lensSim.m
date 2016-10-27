% Function takes inputs for dimensions and number of iterations for the
% simulation.
function lensSim;
dim = 100;
iterations = 5000;
% Just going to use simple values for the constants. Also since the spaces
% between each spot on the grid is 1, dx and dx are 1. The one we need to
% determine is dt. In this part of the project, mu, epsilon, dx, and dy are
% all 1 so for the sake of simplicity they will be omitted and then
% reincluded in the later parts of the project.

% mu = 1;
% ep = 1; For areas outside of the lens the epsilon value will still be 1.
% dx = 1;
% dy = 1;
dt = .3;

% First step is to set up the ez, hx, hy grids using the input and also
% setting up the matricies to be used in the for loops for updating.
ez = zeros(dim);
hx = zeros(dim, dim-1);
hy = zeros(dim-1, dim);
nextEz = ez;
nextHx = hx;
nextHy = hy;

% Since we are making a lens in this part of the project we need a matrix
% of epsilon values that will be used to help focus the wave. Here we set
% up the matrix of epsilon values for the grid with special values of
% epsilon to simulate a lens.
ep = ones(dim);
n = 1:dim;
l1 = round(.05.*(n-round(dim/2)).^2 + round(dim.*(.4)));
l2 = -round(.05.*(n-round(dim/2)).^2 - round(dim.*(.6)));

for i = (dim*.35):(dim*.65)
    ep(i, l1(i):l2(i)) = 2.5;
end

% This main loop will drive the update of the three components that we need
% to be updated.
for k = 1:iterations

    % This is the source of the electric field at the center of the grid.
    ez(1:dim, 2) = cos(pi.*k/20.*dt);

    % This subloop updates the electric field for every point on the grid.
    for i = 2:dim-1
        for j = 2:dim-1
            nextEz(i,j) = ez(i,j) + (1./ep(i,j)).*dt.*((hy(i-1,j) - hy(i,j)) - (hx(i,j) - hx(i,j-1)));
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

    % The following four lines set the boundaries to 0 and simulate the
    % metallic walls surrounding the grid.
    ez(1,:) = 0;
    ez(:,1) = 0;
    ez(dim,:) = 0;
    ez(:,dim) = 0;

    % These few lines just make the animation of Ez propagating in the
    % grid.
   figure(1);
  % surf(ez);
   surf(ep);
   drawnow;
end

end