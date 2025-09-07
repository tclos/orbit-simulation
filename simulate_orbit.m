clear; close all; clc;

% --- Constants ---
G      = 6.67430e-11;       % m^3 kg^-1 s^-2
M_sun  = 1.989e30;          % kg

% --- Planet parameters ---
planets = {
    'Mercury', 0.6e11, 0.15, [0.8 0.8 0.8];
    'Venus',   0.9e11, 0.05, [1.0 0.7 0.2];
    'Earth',   1.0e11, 0.10, [0.2 0.6 1.0];
    'Mars',    1.2e11, 0.20, [1.0 0.4 0.3];
    'Jupiter', 1.5e11, 0.05, [0.9 0.8 0.5];
};

nFrames = 1500;          % Number of animation frames
frame_skip = 2;         % Skipping some frames for performance
fig = figure('Color', 'k');
ax = axes('Parent', fig);
hold(ax, 'on');
plot3(ax, 0, 0, 0, 'yo', 'MarkerSize', 20, 'MarkerFaceColor', 'y'); % Sun

maxR_all = 0;
planet_handles = [];
trail_handles  = [];
orbit_data     = {};

% Precompute planetary orbits
for p = 1:size(planets, 1)
    name = planets{p, 1};
    a    = planets{p, 2};
    e    = planets{p, 3};
    col  = planets{p, 4};

    % Random initial angle (true anomaly)
    theta0 = rand() * 2 * pi;

    % Random inclination (between 0 and 90 degrees)
    inclination = rand() * pi / 4;  % Random inclination between 0 and 45 degrees

    % Distance from Sun at this true anomaly (using the vis-viva equation)
    r_mag = a * (1 - e^2) / (1 + e * cos(theta0));

    % Position in orbital plane (accounting for inclination)
    r0 = [r_mag * cos(theta0) * cos(inclination);
          r_mag * sin(theta0) * cos(inclination);
          r_mag * sin(inclination)];

    % Orbital speed from vis-viva equation
    v_mag = sqrt(G * M_sun * (2 / r_mag - 1 / a));

    % Velocity direction (perpendicular to position vector)
    v_dir = [-sin(theta0); cos(theta0); 0];
    v0 = v_mag * v_dir;

    % Combine into initial state vector
    Y0 = [r0; v0];

    % Orbital period
    T_orbit = 2 * pi * sqrt(a^3 / (G * M_sun));

    % Solve orbital equation using ode45 (Octave uses the same function)
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
    odefun = @(t, y) [y(4:6); -G * M_sun * y(1:3) / norm(y(1:3))^3];
    [t_sol, Y_sol] = ode45(odefun, [0 T_orbit], Y0, opts);

    % Interpolate for animation
    t_anim = linspace(0, T_orbit, nFrames);
    Y_anim = interp1(t_sol, Y_sol, t_anim);
    x = Y_anim(:, 1);
    y = Y_anim(:, 2);
    z = Y_anim(:, 3);

    orbit_data{p} = struct('x', x, 'y', y, 'z', z, 'color', col, 'name', name);
    maxR_all = max(maxR_all, max(sqrt(x.^2 + y.^2 + z.^2)));

    % Plot initial position and trails
    planet_handles(p) = plot3(ax, x(1), y(1), z(1), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'none');
    trail_handles(p)  = plot3(ax, x(1), y(1), z(1), '-', 'LineWidth', 1.2, ...
        'Color', col);
end

% Set axis limits and grid
axis(ax, 'equal');
axis(ax, [-1 1 -1 1 -1 1] * maxR_all * 1.3);
grid(ax, 'on');
view(ax, 3);
xlabel(ax, 'x [m]', 'Color', 'w');
ylabel(ax, 'y [m]', 'Color', 'w');
zlabel(ax, 'z [m]', 'Color', 'w');
set(ax, 'Color', [0 0 0], 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');

% --- Animation ---
for k = 1:frame_skip:nFrames
    for p = 1:size(planets, 1)
        % Update planet position
        set(planet_handles(p), 'XData', orbit_data{p}.x(k), ...
                               'YData', orbit_data{p}.y(k), ...
                               'ZData', orbit_data{p}.z(k));

        % Update trail
        set(trail_handles(p), 'XData', orbit_data{p}.x(1:k), ...
                              'YData', orbit_data{p}.y(1:k), ...
                              'ZData', orbit_data{p}.z(1:k));
    end
    % Update title for current frame
    title(ax, sprintf('Frame %d / %d', k, nFrames), 'Color', 'w');
    drawnow;
end
s
