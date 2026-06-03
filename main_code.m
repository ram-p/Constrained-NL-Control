% This code implements the algorithm in the paper "Ignore Drift, Embrace
% Simplicity: Constrained Nonlinear Control through Driftless
% Approximation" on https://arxiv.org/abs/2509.06188.
% 
% The results of this code exactly duplicate Example 2 in the paper, for
% the forced Van der Pol oscillator. Modifying the code for a particular
% user-defined system is simple using the functions f and g. The main
% design parameter is t(2) = t^*, as described in the paper.

close all
clear
clc
format compact

T = 30;             % Final time to cut simulation
x0 = [3; 0.1];      % Initial state
xtg = [0; 0];       % Target state

% Horizon of first piece: t(1) = t_0, t(2) = t^*.
t(1) = 0;
t(2) = 7;

% Input matrix used for driftless approximation
B = [0; 1];
n = size(B, 1);     % State dimension
m = size(B, 2);     % Input dimension

% Constructing the sequence of time instants until the horizon T.
i = 2;
while t(i) < T
    t(i+1) = t(i) + t(2)/i;
    i = i+1;
end
t(end) = T;
pieces = length(t)-1;

% Initializing trajectory arrays
xtraj = [];
ttraj = [];
utraj = [];

u = zeros(m, pieces);           % Initializing input
xpoint = zeros(n, pieces+1);    % Initializing state reached at each point in sequence
xpoint(:,1) = x0;

% Method
for i = 1:pieces
    delta = t(i+1)-t(i);
    tspan = [t(i) t(i+1)];
    xBtilde = pinv(B)*(xtg - xpoint(:,i));
    u(:,i) = (1/delta)*xBtilde;             % Control applied over tspan
    
    % Solving system dynamics in sys_eg for the control above, over tspan.
    [t_hor,x] = ode45(@(t,x) sys_eg(t,x,u(:,i)), tspan, xpoint(:,i));

    xpoint(:,i+1) = x(end,:)';

    % Updating trajectory arrays
    xtraj = horzcat(xtraj, x');
    ttraj = vertcat(ttraj, t_hor);
    u_piece = repmat(u(:,i), 1, length(t_hor));
    utraj = horzcat(utraj, u_piece);
end

% Plots
figure(1)
set(gcf, 'DefaultLineLineWidth', 2.5)
plot(ttraj, xtraj);
grid on
box on
xlabel('Time')
ylabel('State trajectory')
set(gca, 'FontSize', 25, 'FontName', 'Times New Roman')
legend('$x$', '$y$', 'Interpreter', 'latex')

figure(2)
set(gcf, 'DefaultLineLineWidth', 2.5)
plot(ttraj, utraj);
grid on
box on
xlabel('Time')
ylabel('Input trajectory')
set(gca, 'FontSize', 25, 'FontName', 'Times New Roman')
