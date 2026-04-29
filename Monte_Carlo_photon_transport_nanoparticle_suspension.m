% PHOTON TRANSPORT MONTE CARLO IN A NANOPARTICLE SUSPENSION
%
% This script simulates laser propagation through a nanoparticle suspension
% contained in a rectangular cuvette. Photons are launched from a defined
% beam profile and propagated through the solution until they are absorbed,
% scattered, or transmitted through one of the cuvette boundaries.
%
% The optical response of the suspension is described by the nanoparticle
% absorption cross-section, scattering cross-section, concentration, and
% scattering anisotropy factor g. The interaction length is sampled from
% the extinction coefficient, mu_ext = mu_abs + mu_sca. At each interaction
% event, the photon is either absorbed or scattered according to the
% relative probabilities mu_abs/mu_ext and mu_sca/mu_ext. Scattering angles
% are sampled from the Henyey-Greenstein phase function.
%
% The geometry is defined by the liquid volume of the cuvette and by a
% cubic voxel size. Absorption events are accumulated in a 3D array, Nabs,
% which is then converted into a volumetric heating rate q_abs.
%
% Each simulated photon represents an equal fraction of the optical power
% entering the cuvette, not one physical photon. Therefore, for an input
% power P_in and voxel volume dV, the heating rate is calculated as:
%
%     q_abs = P_in * Nabs / (N_inside * dV)
%
% where N_inside is the number of simulated photons that enter the liquid
% in the cuvette. The resulting q_abs is expressed in mW/cm^3.
%
% The implementation is inspired by the Monte Carlo light transport
% framework of Wang, Jacques, and Zheng, "MCML: Monte Carlo modeling of
% light transport in multi-layered tissues", Computer Methods and Programs
% in Biomedicine 47, 131-146 (1995), adapted here for a nanoparticle
% suspension in a finite cuvette geometry.

clc

% Define optical constants of materials
n_solvent = 1.333;         % Refractive index of solvent (water = 1.333)
n_cuvette = 1.44;      % Refractive index of cuvette (glass = 1.44)
n_air = 1;           % Refractive index of air

%% USER INPUT: optical parameters, geometry, grid, propagation, and laser

% Check whether all user-defined input parameters already exist
check_input = exist('N','var') && ...
    exist('g','var') && ...
    exist('Cabs','var') && ...
    exist('Csca','var') && ...
    exist('NPconc','var') && ...
    exist('x_length','var') && ...
    exist('y_length','var') && ...
    exist('z_length','var') && ...
    exist('voxel_size','var') && ...
    exist('x_beam','var') && ...
    exist('y_beam','var') && ...
    exist('z_beam','var') && ...
    exist('P_in','var') && ...
    exist('Gaussian_Beam','var') && ...
    exist('beam_radius','var');

if ~check_input
    
    prompt = { ...
        'Number of injected photons, N', ...
        'Scattering anisotropy factor, g (back = -1, isotropic = 0, forward = 1)', ...
        'Absorption cross-section, Cabs (cm^2)', ...
        'Scattering cross-section, Csca (cm^2)', ...
        'Nanoparticle concentration, NPconc (NPs/mL)', ...
        'x length of the liquid volume (cm)', ...
        'y length of the liquid volume (cm)', ...
        'z length of the liquid volume (cm)', ...
        'voxel side length (cm)', ...
        'Beam enters at x_beam (cm) (beam propagates along +x)', ...
        'Beam enters at y_beam (cm)', ...
        'Beam enters at z_beam (cm)', ...
        'Power entering the liquid (mW)', ...
        'Gaussian Beam (1) or Circular Beam (0)', ...
        'Gaussian beam 1/e^2 radius or circular beam radius (cm)'};
    
    dlg_title = 'Input: optical parameters, cuvette, grid, propagation, and laser beam';
    
    defaultans = { ...
        '1000000', ...
        '0', ...
        '6.80e-11', ...
        '9.47e-12', ...
        '2.97e10', ...
        '1.0', ...
        '0.2', ...
        '0.4', ...
        '0.005', ...
        '0', ...
        '0.1', ...
        '0.2', ...
        '400', ...
        '1', ...
        '0.1'};
    
    num_lines = [1 70];
    
    answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
    
    if isempty(answer)
        error('Input dialogue cancelled by user.')
    end
    
    N             = str2double(answer{1});
    g             = str2double(answer{2});
    Cabs          = str2double(answer{3});
    Csca          = str2double(answer{4});
    NPconc        = str2double(answer{5});
    x_length      = str2double(answer{6});
    y_length      = str2double(answer{7});
    z_length      = str2double(answer{8});
    voxel_size    = str2double(answer{9});
    x_beam        = str2double(answer{10});
    y_beam        = str2double(answer{11});
    z_beam        = str2double(answer{12});
    P_in          = str2double(answer{13});
    Gaussian_Beam = str2double(answer{14});
    beam_radius   = str2double(answer{15});
    
end

%% Derived geometry and grid variables

BC = [0, 0, 0, x_length, y_length, z_length];

dr = zeros(1,6);
ds = zeros(1,6);

nx = round(x_length/voxel_size);
ny = round(y_length/voxel_size);
nz = round(z_length/voxel_size);

dx = x_length/nx;
dy = y_length/ny;
dz = z_length/nz;

x = dx/2 : dx : x_length - dx/2;
y = dy/2 : dy : y_length - dy/2;
z = dz/2 : dz : z_length - dz/2;

%% Laser beam initialization

% Photons propagate along +x, axes perpendicular to propagation are y and z
u_i = [1, 0, 0];

r_i = zeros(N,3);

if Gaussian_Beam == 1
    sigma = beam_radius/2;
    r_i(:,1) = x_beam;                    % x: injection plane
    r_i(:,2) = y_beam + sigma*randn(N,1); % y
    r_i(:,3) = z_beam + sigma*randn(N,1); % z
else
    [y0, z0] = rand_circ(N, y_beam, z_beam, beam_radius);
    r_i(:,1) = x_beam;                    % x: injection plane
    r_i(:,2) = y0;                     % y
    r_i(:,3) = z0;                     % z
end

%% BASIC VALIDATIONS

if g < -1 || g > 1
    error('g must be between -1 and 1.')
end

if Cabs < 0 || Csca < 0 || NPconc < 0
    error('Cabs, Csca, and NPconc must be non-negative.')
end

if voxel_size <= 0
    error('voxel_size must be positive.')
end

if x_beam < 0 || x_beam > x_length || y_beam < 0 || y_beam > y_length || z_beam < 0 || z_beam > z_length
    error('Beam center must be inside or on the boundary of the cuvette.')
end

%% CALCULATIONS

% Extinction cross-section in cm^2
Cext = Csca + Cabs;

% Total internal reflection angle
a_TIR = asin(n_air/n_solvent);

% Scattering, absorption, and extinction coefficients in cm^(-1)
mu_sca = NPconc*Csca;
mu_abs = NPconc*Cabs;
mu_ext = NPconc*Cext;

% Nabs stores the number of absorption events in each voxel.
% Array order is Nabs(ix, iy, iz), corresponding to q_abs(x,y,z).
Nabs = zeros(nx, ny, nz);

% We create a vector to identify the exact boundaries where the photon hit
j = zeros(1,6);

% We create some variables to track the number of photons escaping the
% simulation if the are transmitted:
T_top =    0;
T_bottom = 0;
T_xplus =  0;
T_xminus = 0;
T_yplus =  0;
T_yminus = 0;

% Variables to study the code
% We create an index to track when a photon has been absorbed. For each
% absorption event, the variable i_abs increase one unit. That allows to
% creat a vector with ordered elements indexed by i_abs.
i_abs = 0;

% Preallocation
Event_X =   zeros(N,1);         % Event position in x
Event_Y =   zeros(N,1);         % Event position in y
Event_Z =   zeros(N,1);         % Event position in z
Event_sca = zeros(N,1);       % Track number of scattering events before absorption
Event_hit = zeros(N,1);       % Track number of photons that hit the BCs

%% THE MONTE CARLO METHOD STARTS HERE:

tic

N_skipped = 0;

for i = 1:N
    
    % Generate a random initial position for the photon with a spatial
    % distribution that follows the gaussian or circular beam profile.
    X = r_i(i,1); % random x coordinate (= 0)
    Y = r_i(i,2); % random y coordinate
    Z = r_i(i,3); % random z coordinate
    u = u_i;      % Initial direction of the photon packet
    
    % Discard launched photons outside the liquid volume:
    if X < 0 || X > x_length || ...
            Y < 0 || Y > y_length || ...
            Z < 0 || Z > z_length
        
        N_skipped = N_skipped + 1;
        continue
    end
    
    % Photon is active while W = 1; absorbed or transmitted photons terminate.
    W = 1;
    
    % Number of scattering events before absorption.
    Sca = 0;
    
    while W ~= 0
        
        % Record the old position
        X_m = X;
        Y_m = Y;
        Z_m = Z;
        
        % Sample distance to next interaction from exponential distribution
        % (Beer–Lambert law), mean free path = 1/mu_ext.
        e = rand;
        s = -log(e)/mu_ext;
        
        % Update position:
        X = X + u(1)*s;
        Y = Y + u(2)*s;
        Z = Z + u(3)*s;
        r = [X, Y, Z];
        
        % Check whether the photon crossed any boundary during this step.
        % If so, ds stores how far the photon travelled beyond each crossed
        % boundary. The largest ds corresponds to the first boundary reached
        % along the photon trajectory.
        for k = 1:3
            
            % Lower boundaries: x = 0, y = 0, z = 0
            dr(k) = BC(k) - r(k);
            
            % Upper boundaries: x = x_length, y = y_length, z = z_length
            dr(k + 3) = r(k) - BC(k + 3);
            
            % Crossing a lower boundary while moving in the negative
            % direction along axis k
            if u(k) < 0 && dr(k) > 0
                ds(k) = abs(dr(k)/u(k));
            else
                ds(k) = 0;
            end
            
            % Crossing an upper boundary while moving in the positive
            % direction along axis k
            if u(k) > 0 && dr(k + 3) > 0
                ds(k + 3) = abs(dr(k + 3)/u(k));
            else
                ds(k + 3) = 0;
            end
        end
        
        % Identify the first boundary encountered during the step.
        [M, I] = max(ds);
        
        % Boundary index vector:
        % j = [x_min, y_min, z_min, x_max, y_max, z_max]
        j = zeros(1,6);
        
        if M > 0
            j(I) = 1;
            Event_hit(i) = Event_hit(i) + 1;
        end
        
        % Handle interaction with cuvette boundaries.
        % If the photon step crosses a boundary, it is first moved back to the
        % intersection point. Fresnel reflection/transmission is then applied.
        % Transmitted photons leave the simulation; reflected photons continue.
        
        if j(1) == 1   % x_min: back reflection (water → glass)
            a_i = acos(abs(u(1)));
            
            X = 0;
            Y = Y_m + u(2)*(s - ds(I));
            Z = Z_m + u(3)*(s - ds(I));
            
            R = fresnel_R(n_solvent, n_cuvette, a_i);
            if rand <= R
                u(1) = -u(1);
            else
                T_xminus = T_xminus + W;
                break
            end
            
        elseif j(2) == 1   % y_min: front cuvette wall (water → glass)
            a_i = acos(abs(u(2)));
            
            X = X_m + u(1)*(s - ds(I));
            Y = 0;
            Z = Z_m + u(3)*(s - ds(I));
            
            R = fresnel_R(n_solvent, n_cuvette, a_i);
            if rand <= R
                u(2) = -u(2);
            else
                T_yminus = T_yminus + W;
                break
            end
            
        elseif j(3) == 1   % z_min: bottom of the cuvette (water → glass)
            a_i = acos(abs(u(3)));
            
            X = X_m + u(1)*(s - ds(I));
            Y = Y_m + u(2)*(s - ds(I));
            Z = 0;
            
            R = fresnel_R(n_solvent, n_cuvette, a_i);
            if rand <= R
                u(3) = -u(3);
            else
                T_bottom = T_bottom + W;
                break
            end
            
        elseif j(4) == 1   % x_max: forward transmission (water → glass)
            a_i = acos(abs(u(1)));
            
            X = x_length;
            Y = Y_m + u(2)*(s - ds(I));
            Z = Z_m + u(3)*(s - ds(I));
            
            R = fresnel_R(n_solvent, n_cuvette, a_i);
            if rand <= R
                u(1) = -u(1);
            else
                T_xplus = T_xplus + W;
                break
            end
            
        elseif j(5) == 1   % y_max: back cuvette wall (water → glass)
            a_i = acos(abs(u(2)));
            
            X = X_m + u(1)*(s - ds(I));
            Y = y_length;
            Z = Z_m + u(3)*(s - ds(I));
            
            R = fresnel_R(n_solvent, n_cuvette, a_i);
            if rand <= R
                u(2) = -u(2);
            else
                T_yplus = T_yplus + W;
                break
            end
            
        elseif j(6) == 1   % z_max: top interface (water → air)
            a_i = acos(abs(u(3)));
            
            X = X_m + u(1)*(s - ds(I));
            Y = Y_m + u(2)*(s - ds(I));
            Z = z_length;
            
            if a_i >= a_TIR   % total internal reflection
                u(3) = -u(3);
            else
                R = fresnel_R(n_solvent, n_air, a_i);
                if rand <= R
                    u(3) = -u(3);
                else
                    T_top = T_top + W;
                    break
                end
            end
        end
        
        if M == 0 % The photon has not crossed any boundary
            
            % Calculation of the absorption probability:
            Prob_abs = (mu_abs/mu_ext);
            
            if rand <= Prob_abs % The photon is absorbed and the loop is terminated
                
                i_abs = i_abs + 1;
                Event_sca(i_abs) = Sca;
                
                % Update the spatial distribution of the absorption:
                i_x = min(max(floor(X/dx) + 1, 1), nx);
                i_y = min(max(floor(Y/dy) + 1, 1), ny);
                i_z = min(max(floor(Z/dz) + 1, 1), nz);
                
                Nabs(i_x,i_y,i_z) = Nabs(i_x,i_y,i_z) + 1;
                
                Event_X(i_abs) = x(i_x);
                Event_Y(i_abs) = y(i_y);
                Event_Z(i_abs) = z(i_z);
                W = 0;
                
            else  % The photon is scattered
                Sca = Sca + 1;
                
                % Calculate the scattering and azimuthal angles:
                if abs(g) > 0
                    costheta = 1/(2*g)*(1 + g^2 - ((1 - g^2)/(1 - g +2*g*rand))^2);
                    sintheta = sqrt(1 - costheta^2);
                    phi = 2*pi*rand;
                    cosphi = cos(phi);
                    sinphi = sin(phi);
                else
                    costheta = (1 - 2*rand);
                    sintheta = sqrt(1 - costheta^2);
                    phi = 2*pi*rand;
                    cosphi = cos(phi);
                    sinphi = sin(phi);
                end
                
                % Calculate the new direction of propagation:
                
                ux = u(1);
                uy = u(2);
                uz = u(3);
                
                if abs(ux) > 0.9999
                    u_new = [ ...
                        sign(ux)*costheta, ...
                        sintheta*cosphi, ...
                        sintheta*sinphi];
                else
                    denom = sqrt(1 - ux^2);
                    
                    u_new = [ ...
                        ux*costheta - denom*sintheta*cosphi, ...
                        uy*costheta + sintheta*(ux*uy*cosphi - uz*sinphi)/denom, ...
                        uz*costheta + sintheta*(ux*uz*cosphi + uy*sinphi)/denom];
                end
                
                u = u_new / norm(u_new);
                
            end % End of absorption/scattering event
        end % End of step: no boundary crossing (interaction in bulk)
    end % End of photon trajectory (terminated by absorption or escape)
end % End of Monte Carlo loop over all photon packets

% Calculate the heat source density generated:
dV = dx*dy*dz;
Abs = sum(sum(sum(Nabs)));
N_inside = N - N_skipped;
q_abs = P_in * Nabs / (N_inside * dV);

% Check photon bookkeeping
Net = T_top + T_bottom + T_xplus + T_xminus + T_yplus + T_yminus + Abs + N_skipped - N;

% Calculate the effective optical density of the suspension
T_forward = T_xplus / N_inside;
OD_eff = -log10(T_forward);

fprintf('Skipped photons: %d / %d = %.2f %%\n', N_skipped, N, 100*N_skipped/N)
fprintf('Photon bookkeeping residual: %g\n', Net)
fprintf('Effective optical density: %g\n', OD_eff)

toc

%% Plots

[~, iy_beam] = min(abs(y - y_beam));
[~, iz_beam] = min(abs(z - z_beam));

figure(1); clf
set(gcf, 'DefaultAxesFontSize', 12)
set(gcf, 'DefaultTextFontSize', 12)
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
set(gcf, 'Color', 'w')

% Top left: stacked xy slices
subplot(2,2,1)
z_slices = linspace(0, z_length, 12);
hold on
[X,Y] = meshgrid(x,y);

for k = 1:length(z_slices)
    [~, iz] = min(abs(z - z_slices(k)));
    A = squeeze(q_abs(:,:,iz))';
    Z = z(iz) * ones(size(X));
    surf(X, Y, Z, A, ...
        'EdgeColor','none', ...
        'FaceAlpha',0.9)
end
colormap(hot)
cb = colorbar;
cb.Label.String = 'q_{abs} (mW/cm^3)';
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('z (cm)')
title('Heating rate q_{abs}')
axis tight
daspect([1 1 1])
grid on
box on
view(-35,25)
hold off

% Top right: xz slice through beam center y_beam
subplot(2,2,2)
surf(x,z,squeeze(q_abs(:,iy_beam,:))','EdgeColor','none')
colormap(hot)
cb = colorbar;
cb.Label.String = 'q_{abs}(x,z) (mW/cm^3)';
cb.Label.Interpreter = 'tex';
view(0,90)
xlabel('x (cm)')
ylabel('z (cm)')
title('q_{abs}(x,z) at y = y_{beam}')
axis equal tight

% Bottom left: xy slice through beam center z_beam
subplot(2,2,3)
surf(x,y,squeeze(q_abs(:,:,iz_beam))','EdgeColor','none')
colormap(hot)
cb = colorbar;
cb.Label.String = 'q_{abs}(x,y) (mW/cm^3)';
cb.Label.Interpreter = 'tex';
view(0,90)
xlabel('x (cm)')
ylabel('y (cm)')
title('q_{abs}(x,y) at z = z_{beam}')
axis equal tight

% Bottom right: histogram with analytical geometric distribution
subplot(2,2,4)
data = Event_sca(1:i_abs);
[counts, x_vals] = hist(data, unique(data));
bar(x_vals, counts, ...
    'FaceColor', [0.6 0.6 0.6], ...
    'EdgeColor', 'k')
hold on
p_abs = mu_abs / mu_ext; % absorption probability
p_sca = mu_sca / mu_ext; % scattering probability
k_vals = 0:max(x_vals);
expected_counts = i_abs * p_abs * p_sca.^k_vals; % geometric distribution probability
semilogy(k_vals, expected_counts, 'k--', 'LineWidth', 2)
set(gca,'YScale','log')
xlabel('Scattering events before absorption, S')
ylabel('Counts')
title('Photon transport statistics')
legend('Monte Carlo', 'Ideal geometric distribution = N_{abs}P_{abs}P_{sca}^S', 'Location', 'northeast')
grid on
hold off

%% Local functions

function [X,Y] = rand_circ(N,x,y,r)
% Generate N random points uniformly distributed in a circle
% of radius r centered at (x,y).

if nargin < 2
    x = 0; y = 0; r = 1;
end

Ns = round(1.28*N + 2.5*sqrt(N) + 100); % oversampling factor
X = rand(Ns,1)*(2*r) - r;
Y = rand(Ns,1)*(2*r) - r;

I = (X.^2 + Y.^2) <= r^2;

X = X(I);
Y = Y(I);

X = X(1:N) + x;
Y = Y(1:N) + y;
end

function R = fresnel_R(n1,n2,a_i)
% Calculate Fresnel reflection avoiding indeterminate solutions arising
% when photons propagate perpendicularly to the interface.

if abs(a_i) < 1e-12
    R = ((n1 - n2)/(n1 + n2))^2;
else
    a_t = asin(sin(a_i)*n1/n2);
    R = 0.5*((sin(a_i - a_t)/sin(a_i + a_t))^2 + ...
        (tan(a_i - a_t)/tan(a_i + a_t))^2);
end

end