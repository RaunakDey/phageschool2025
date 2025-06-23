% Define parameters
pars.phi = 6.7e-10;
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 5e6;
pars.beta = 50;
pars.eta = 0.5;

% Disease-Free Equilibrium (DFE) of S
S_DFE = pars.K * (1 - pars.d / pars.b);

% Components of R0
prod = pars.beta;
proba_lyse = (pars.phi * S_DFE) / (pars.phi * S_DFE + pars.m);
proba_inf = pars.eta / (pars.eta + pars.d);

% Compute R0
R0 = prod * proba_lyse * proba_inf;

% Display
fprintf('R0 = %.4f\n', R0);




% Time
t0 = 0;
tf = 2000;
T = linspace(t0, tf, 500);  % same as Python's linspace

% Parameters
pars.phi = 6.7e-10;
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 5e6;
pars.beta = 50;
pars.eta = 0.5;

% Disease-Free Equilibrium
S_DFE = pars.K * (1 - pars.d / pars.b);

% Initial Conditions [S, I, V]
y0 = [S_DFE; 0; 100];

% Solve the system
[T_out, Y] = ode45(@(t, y) SIV_ode(t, y, pars), T, y0);
S = Y(:, 1);
I = Y(:, 2);
V = Y(:, 3);

% Plotting
figure;
semilogy(T_out, S + I, 'k', 'LineWidth', 2); hold on;
semilogy(T_out, S, 'g', 'LineWidth', 2);
semilogy(T_out, V, 'r', 'LineWidth', 2);
scatter(0, S_DFE, 80, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'red');

legend('S + I', 'S', 'V', 'DFE', 'Location', 'best', 'Box', 'on');
xlabel('Time (h)');
ylabel('Density (cells or viruses/mL)');
title('SIV model dynamics from DFE');
grid on;
set(gca,'FontSize',24);
saveas(gcf, 'siv_dynamics.png');



% Time
t0 = 0;
tf = 2000;
T = linspace(t0, tf, 500);

% Parameters (same as before, but beta changed)
pars.phi = 6.7e-10;
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 5e6;
pars.beta = 10; % changed here
pars.eta = 0.5;

% Disease-Free Equilibrium
S_DFE = pars.K * (1 - pars.d / pars.b);

% Initial Conditions [S, I, V]
y0 = [S_DFE; 0; 100];

% Integrate system
[T_out, Y] = ode45(@(t, y) SIV_ode(t, y, pars), T, y0);
S = Y(:, 1);
I = Y(:, 2);
V = Y(:, 3);

% Plotting
figure;
semilogy(T_out, S + I, 'k', 'LineWidth', 2); hold on;
semilogy(T_out, S, 'g', 'LineWidth', 2);
semilogy(T_out, V, 'r', 'LineWidth', 2);
scatter(0, S_DFE, 80, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'red');

legend('S + I', 'S', 'V', 'DFE', 'Location', 'best', 'Box', 'on');
xlabel('Time (h)');
ylabel('Density (cells or viruses/mL)');
title('SIV Model Dynamics (\beta = 10)');
set(gca,'FontSize',24);
grid on;
saveas(gcf, 'siv_beta_10.png');


%%%%%%%


% Parameters
pars.phi = 6.7e-10;
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 5e6;
pars.eta = 0.5;

% Disease-Free Equilibrium
S_DFE = pars.K * (1 - pars.d / pars.b);

% Critical burst size beta_crit = numerator / denominator
num = (pars.eta + pars.d) * (pars.phi * S_DFE + pars.m);
den = pars.eta * pars.phi * S_DFE;
beta_crit = num / den;

% Sweep over burst sizes
Brange = linspace(1, 100, 20);
R0_vals = zeros(size(Brange));

figure; hold on;

% Loop over burst sizes
for i = 1:length(Brange)
    beta = Brange(i);
    pars.beta = beta;

    % Components of R0
    prod = beta;
    proba_lyse = (pars.phi * S_DFE) / (pars.phi * S_DFE + pars.m);
    proba_inf = pars.eta / (pars.eta + pars.d);
    R0 = prod * proba_lyse * proba_inf;

    % Store and plot
    R0_vals(i) = R0;
    c = [0.75, 0.75, 0.75] * i / length(Brange);
    scatter(beta, R0, 50, 'filled', 'MarkerFaceColor', c);
end

% Highlight beta_crit
scatter(beta_crit, 1, 100, 'r', 'filled', 'DisplayName', '\beta_{crit}');
yline(1, 'k--', 'DisplayName', 'R_0 = 1');

xlabel('Burst Size (\beta)');
ylabel('R_0');
title('Invasion Threshold for Varying Burst Size');
xlim([0, 105]);
ylim([0, max(R0_vals) + 0.5]);
%legend('Location', 'southeast');
grid on;
saveas(gcf, 'siv_beta_variation.png');



%%%%%%%%%%%%%%%%%%%
% Time setup
t0 = 0;
tf = 2000;
T = linspace(t0, tf, 500);

% Parameters
pars.phi = 6.7e-10;
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 5e6;
pars.beta = 50;
pars.eta = 0.1; % Updated eta

% Disease-Free Equilibrium
S_DFE = pars.K * (1 - pars.d / pars.b);

% Initial Conditions: [S, I, V]
y0 = [S_DFE; 0; 100];

% Integrate the ODE
[T_out, Y] = ode45(@(t, y) SIV_ode(t, y, pars), T, y0);
S = Y(:,1);
I = Y(:,2);
V = Y(:,3);

% Plotting
figure;
semilogy(T_out, S + I, 'k', 'LineWidth', 2); hold on;
semilogy(T_out, S, 'g', 'LineWidth', 2);
semilogy(T_out, V, 'r', 'LineWidth', 2);
scatter(0, S_DFE, 80, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'red');

legend('S + I', 'S', 'V', 'DFE', 'Location', 'best', 'Box', 'on');
xlabel('Time (h)');
ylabel('Density (cells or viruses/mL)');
title('SIV Model Dynamics (\eta = 0.1 h^{-1})');
grid on;
set(gca,'FontSize',24);
saveas(gcf, 'siv_eta_0pt1.png');




%%%%%%%%%

% Parameters
pars.phi = 6.7e-10;
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 5e6;
pars.beta = 50;

% Disease-Free Equilibrium
S_DFE = pars.K * (1 - pars.d / pars.b);

% Numerator and Denominator
num = pars.d * (pars.phi * S_DFE + pars.m);
den = pars.beta * pars.phi * S_DFE - pars.phi * S_DFE - pars.m;

% Critical lysis rate
eta_crit = num / den;

% Display
fprintf('Critical lysis rate (eta_crit) = %.4f h^{-1}\n', eta_crit);



%%%%%%%

% Time setup
t0 = 0;
tf = 5000;
dt = 0.1;
T = t0:dt:tf;

% Fixed parameters
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 5e6;
pars.beta = 50;
pars.eta = 0.5;

% DFE for plotting reference
S_DFE = pars.K * (1 - pars.d / pars.b);

% List of phi values and corresponding initial conditions
phi_list = [1e-12, 1e-9, 5e-9, 1e-8];
y0_list = {[1000; 0; 100], [1000; 0; 100], [100; 0; 1e6], [100; 0; 1e6]};
titles = {'Low \phi', 'Medium \phi', 'Medium \phi (high V_0)', 'High \phi'};

% Create subplot figure
figure;
for i = 1:4
    pars.phi = phi_list(i);
    y0 = y0_list{i};

    % Integrate SIV model
    [~, Y] = ode45(@(t, y) SIV_ode(t, y, pars), T, y0);
    S = Y(:, 1);
    V = Y(:, 3);

    % Plot in subplot
    subplot(2, 2, i);
    loglog(S, V, 'LineWidth', 2); hold on;
    scatter(y0(1), y0(3), 50, 'red', 'filled', 'DisplayName', 'Initial');
    scatter(S(end), V(end), 50, 'green', 'filled', 'DisplayName', 'Final');
    if i == 1
        scatter(S_DFE, 1e-20, 50, 'blue', 'filled', 'DisplayName', 'DFE');
    end
    xlabel('Microbe S');
    ylabel('Virus V');
    title(titles{i});
    legend('Location', 'best');
    grid on;
end

sgtitle('SIV Model: Phase Plane for Different \phi');

saveas(gcf, 'siv_phi_variation.png');