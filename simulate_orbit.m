%% Single-planet smooth elliptical orbit (MATLAB / Octave)
clear; close all; clc;

% --- Constants ---
G      = 6.67430e-11;       % m^3 kg^-1 s^-2
M_sun  = 1.989e30;          % kg

% --- Orbital parameters (change 'e' to 0.0167 for Earth-like) ---
a = 1.496e11;   % semi-major axis (1 AU) [m]
e = 0.2;        % eccentricity (0.2 visible ellipse; set 0.0167 for Earth)

% initial position at perihelion and correct initial velocity (vis-viva)
r0 = [ a*(1 - e); 0; 0 ];                          % perihelion
v0y = sqrt( G*M_sun * ( 2 / norm(r0) - 1 / a ) ); % speed from vis-viva
v0 = [0; v0y; 0];                                  % perpendicular to r

Y0 = [r0; v0];   % state: [x;y;z; vx;vy;vz]

% orbital period (Kepler)
T_orbit = 2*pi*sqrt(a^3 / (G*M_sun));
fprintf('Orbital period = %.2f days\n', T_orbit/86400);

% --- Solve (let solver choose internal steps for accuracy) ---
opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
odefun = @(t,y) [ y(4:6); -G*M_sun * y(1:3) / norm(y(1:3))^3 ];
[t_sol, Y_sol] = ode45(odefun, [0 T_orbit], Y0, opts);

% --- Interpolate to uniform time grid for smooth animation ---
nFrames = 1200;                           % increase for smoother animation
t_anim = linspace(0, T_orbit, nFrames);
Y_anim = interp1(t_sol, Y_sol, t_anim);  % returns (nFrames x 6)

x = Y_anim(:,1); y = Y_anim(:,2); z = Y_anim(:,3);

% --- Plot setup ---
fig = figure('Color','k');
ax = axes('Parent',fig);
hold(ax,'on');
plot3(ax, 0,0,0, 'yo', 'MarkerSize', 20, 'MarkerFaceColor', 'y'); % Sun

h_planet = plot3(ax, x(1), y(1), z(1), 'o', 'MarkerSize', 8, ...
                 'MarkerFaceColor', 'r', 'MarkerEdgeColor','none');
h_trail  = plot3(ax, x(1), y(1), z(1), '-', 'LineWidth', 1.2, 'Color', [0.2 0.6 1]);

axis(ax,'equal');
maxR = max( sqrt(x.^2 + y.^2 + z.^2) );
axis(ax, [-1 1 -1 1 -1 1]*maxR*1.2);
grid(ax,'on');
view(ax,3);
xlabel(ax,'x [m]','Color','w'); ylabel(ax,'y [m]','Color','w'); zlabel(ax,'z [m]','Color','w');
set(ax,'Color',[0 0 0],'XColor','w','YColor','w','ZColor','w');

% --- Animation ---
frame_skip = 1;          % show every frame; increase to skip for speed
for k = 1:frame_skip:nFrames
    set(h_planet, 'XData', x(k), 'YData', y(k), 'ZData', z(k));
    set(h_trail,  'XData', x(1:k), 'YData', y(1:k), 'ZData', z(1:k));
    title(ax, sprintf('t = %.2f days', t_anim(k)/86400), 'Color','w');
    drawnow;            % in Octave use drawnow; in MATLAB you can use drawnow limitrate if desired
end
