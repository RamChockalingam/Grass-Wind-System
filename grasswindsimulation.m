% Define the geometry of the metal plate
x = linspace(0, 1, 100);
y = linspace(0, 1, 100);
[X, Y] = meshgrid(x, y);
Z = sin(10*X).*sin(10*Y); % Grass-like structure
plate = struct('X', X, 'Y', Y, 'Z', Z);

% Define the geometry of the permanent magnet
R = 0.1; % Radius
H = 0.2; % Height
magnet = struct('R', R, 'H', H);

% Set up the simulation parameters
v = 10; % Air velocity
B = 1; % Magnetic field strength
theta = pi/4; % Magnetic field orientation
rho = 7800; % Density of metal plate
E = 2e11; % Young's modulus of metal plate

% Solve the FEA equations
tspan = [0 10];
sol = pdepe(0, @(x,t,u,DuDx) plate_motion(x,t,u,DuDx,v,B,theta,rho,E,plate,magnet), ...
            @(x) plate_initial(x,plate), @(xl,ul,xr,ur,t) plate_bc(xl,ul,xr,ur,t,plate));
u = sol(:,:,1);

% Calculate the energy produced
K = 0.5*rho*sum(sum(u.^2)); % Kinetic energy
E = 0.5*B^2*sum(sum(u.*curl(u))); % Electrical energy
energy = K + E;

% Visualize the results
figure;
surf(X, Y, u);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Metal Plate Motion');

figure;
quiver(X, Y, curl(u(:,:,1)), curl(u(:,:,2)));
xlabel('X');
ylabel('Y');
title('Magnetic Field');

fprintf('Total energy produced: %f\n', energy);

% Define the partial differential equations for FEA

function [c,f,s] = plate_motion(~,~,u,DuDx,v,B,theta,rho,E,plate,~)
    % Extract the displacement and velocity from the input vector
    ux = u(:,1);
    uy = u(:,2);
    duxdx = DuDx(:,1);
    duydy = DuDx(:,2);

    % Calculate the air drag force
    Fdragx = -0.5*rho*v*abs(v)*ux;
    Fdragy = -0.5*rho*v*abs(v)*uy;

    % Calculate the Lorentz force due to the magnetic field
    Fmagx = B*sin(theta)*duydy;
    Fmagy = -B*sin(theta)*duxdx;

    % Calculate the elastic force due to the deformation of the metal plate
    [Fx, Fy] = elastic_force(ux, uy, plate.X, plate.Y, plate.Z, E);

    % Calculate the total force on the plate
    Fx_total = Fdragx + Fmagx + Fx;
    Fy_total = Fdragy + Fmagy + Fy;

    % Calculate the acceleration of the plate
    ax = Fx_total./rho;
    ay = Fy_total./rho;

    % Return the PDE coefficients
    c = [1; 1];
    f = [duxdx; duydy];
    s = [ax; ay];
end

% Define the initial conditions for FEA
function u0 = plate_initial(~,plate)
    u0 = [plate.X(:), plate.Y(:)];
end

% Define the boundary conditions for FEA
function [pl, ql, pr, qr] = plate_bc(~,~,~,~,~,~)
    pl = [0; 0];
    ql = [1; 1];
    pr = [0; 0];
    qr = [1; 1];
end

% Define the function for calculating the elastic force
function [Fx, Fy] = elastic_force(ux, uy, ~, ~, Z, E)
    % Calculate the displacement gradients
    [dudx, dudy] = gradient(ux);
    [dvdx, dvdy] = gradient(uy);

    % Calculate the strain tensor
    exx = dudx;
    eyy = dvdy;
    exy = 0.5*(dudy + dvdx);

    % Calculate the stress tensor
    C = E/(1-2*0.3)*(eye(2)*0.3 + [1 -0.3; -0.3 1]);
    sigma_xx = C(1,1)*exx + C(1,2)*eyy;
    sigma_yy = C(2,1)*exx + C(2,2)*eyy;


    % Calculate the force vector
    Fx = -sigma_xx.*Z;
    Fy = -sigma_yy.*Z;
end 