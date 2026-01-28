function projectile_with_drag_optimization()
% PROJECTILE_WITH_DRAG_OPTIMIZATION Simulates projectile motion with drag
% and finds the optimal launch angle for maximum horizontal distance.
%
% This script models a projectile launched from a specified height with
% atmospheric drag proportional to velocity squared. It sweeps through
% launch angles to determine the optimal angle for maximum range.
%
% Drag force model: F_drag = k * v^2
% where k is the drag coefficient and v is velocity magnitude
% The simulation is parameterized by k/m (drag coefficient / mass)

    %% Simulation Parameters
    % Initial conditions
    h0 = 0;              % Initial height (m)
    v0 = 100;              % Initial velocity magnitude (m/s)
    k_over_m = 0.0057;     % Drag coefficient over mass (1/m) - Golf Ball approx
    
    % Physical constants
    g = 9.81;             % Gravitational acceleration (m/s^2)
    
    % Angle sweep parameters
    angles = 15:0.01:60;   % Launch angles to test (degrees)
    num_angles = length(angles);
    
    %% Run simulations for each angle
    ranges = zeros(1, num_angles);
    
    fprintf('Simulating projectile trajectories with drag...\n');
    fprintf('Initial height: %.1f m\n', h0);
    fprintf('Initial velocity: %.1f m/s\n', v0);
    fprintf('k/m (drag parameter): %.4f 1/m\n\n', k_over_m);
    
    for i = 1:num_angles
        theta = angles(i);
        ranges(i) = simulate_trajectory(h0, v0, theta, k_over_m, g);
    end
    
    %% Find optimal angle
    [max_range, max_idx] = max(ranges);
    optimal_angle = angles(max_idx);
    
    fprintf('Results:\n');
    fprintf('========\n');
    fprintf('Optimal launch angle: %.2f degrees\n', optimal_angle);
    fprintf('Maximum range: %.2f m\n\n', max_range);
    
    %% Plot results
    figure('Position', [100, 100, 1200, 500]);
    
    % Subplot 1: Range vs Launch Angle
    subplot(1, 2, 1);
    plot(angles, ranges, 'b-', 'LineWidth', 2);
    hold on;
    plot(optimal_angle, max_range, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('Launch Angle (degrees)', 'FontSize', 12);
    ylabel('Horizontal Range (m)', 'FontSize', 12);
    title('Range vs Launch Angle (With Drag)', 'FontSize', 14);
    legend('Range', sprintf('Optimal: %.1f°', optimal_angle), 'Location', 'best');
    
    % Subplot 2: Optimal trajectory
    subplot(1, 2, 2);
    [~, x_opt, y_opt] = simulate_trajectory(h0, v0, optimal_angle, k_over_m, g);
    plot(x_opt, y_opt, 'b-', 'LineWidth', 2);
    hold on;
    plot(0, h0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(max_range, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('Horizontal Distance (m)', 'FontSize', 12);
    ylabel('Height (m)', 'FontSize', 12);
    title(sprintf('Optimal Trajectory (%.1f°)', optimal_angle), 'FontSize', 14);
    legend('Trajectory', 'Launch Point', 'Landing Point', 'Location', 'best');
    axis equal;
    xlim([0, max_range * 1.1]);
    ylim([0, max(y_opt) * 1.2]);
    
    %% Compare with no-drag case
    fprintf('Comparison with no-drag scenario:\n');
    fprintf('==================================\n');
    
    % No-drag optimal angle (theoretical: 45° from ground level, adjusted for height)
    range_no_drag = simulate_trajectory(h0, v0, optimal_angle, 0, g);
    fprintf('Range at %.1f° with drag: %.2f m\n', optimal_angle, max_range);
    fprintf('Range at %.1f° without drag: %.2f m\n', optimal_angle, range_no_drag);
    fprintf('Range reduction due to drag: %.2f m (%.1f%%)\n', ...
            range_no_drag - max_range, 100 * (range_no_drag - max_range) / range_no_drag);
    
    % Find no-drag optimal
    ranges_no_drag = zeros(1, num_angles);
    for i = 1:num_angles
        ranges_no_drag(i) = simulate_trajectory(h0, v0, angles(i), 0, g);
    end
    [max_range_no_drag, max_idx_no_drag] = max(ranges_no_drag);
    optimal_angle_no_drag = angles(max_idx_no_drag);
    
    fprintf('\nOptimal angle without drag: %.1f°\n', optimal_angle_no_drag);
    fprintf('Maximum range without drag: %.2f m\n', max_range_no_drag);
    fprintf('Angle shift due to drag: %.1f°\n', optimal_angle_no_drag - optimal_angle);
end

function [range, x_traj, y_traj] = simulate_trajectory(h0, v0, theta_deg, k_over_m, g)
% SIMULATE_TRAJECTORY Simulates a single projectile trajectory with drag
%
% Inputs:
%   h0        - Initial height (m)
%   v0        - Initial velocity magnitude (m/s)
%   theta_deg - Launch angle (degrees)
%   k_over_m  - Drag coefficient over mass (1/m)
%   g         - Gravitational acceleration (m/s^2)
%
% Outputs:
%   range     - Horizontal distance traveled (m)
%   x_traj    - Horizontal position array (m)
%   y_traj    - Vertical position array (m)

    % Convert angle to radians
    theta = deg2rad(theta_deg);
    
    % Initial conditions
    vx0 = v0 * cos(theta);
    vy0 = v0 * sin(theta);
    
    % State vector: [x, y, vx, vy]
    initial_state = [0, h0, vx0, vy0];
    
    % Time span for integration (generous upper bound)
    tspan = [0, 20];
    
    % ODE options - stop when projectile hits ground
    options = odeset('Events', @ground_event, 'RelTol', 1e-6, 'AbsTol', 1e-8);
    
    % Solve ODE
    [~, state] = ode45(@(t, s) projectile_ode(t, s, k_over_m, g), tspan, initial_state, options);
    
    % Extract trajectory
    x_traj = state(:, 1);
    y_traj = state(:, 2);
    
    % Range is the final x position
    range = x_traj(end);
end

function dstate = projectile_ode(~, state, k_over_m, g)
% PROJECTILE_ODE Defines the ODEs for projectile motion with drag
%
% State vector: [x, y, vx, vy]
% Drag force: F_drag = -k * v^2 * (v_hat)
%            = -k * v * v_vector

    % Extract state variables
    vx = state(3);
    vy = state(4);
    
    % Velocity magnitude
    v = sqrt(vx^2 + vy^2);
    
    % Drag acceleration components
    % a_drag = -(k/m) * v * v_vector
    ax_drag = -k_over_m * v * vx;
    ay_drag = -k_over_m * v * vy;
    
    % Total acceleration
    ax = ax_drag;
    ay = ay_drag - g;
    
    % State derivatives
    dstate = [vx; vy; ax; ay];
end

function [value, isterminal, direction] = ground_event(~, state)
% GROUND_EVENT Event function to detect when projectile hits the ground
    value = state(2);      % y-position
    isterminal = 1;        % Stop integration
    direction = -1;        % Only detect when decreasing (falling)
end
