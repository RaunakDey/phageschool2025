
% Time settings
t0 = 0;      % Initial time (hours)
tf = 24;     % Final time (hours)
T = linspace(t0, tf, 250);  % Time vector

% Initial conditions
y0 = [1e4; 0];  % [S0; V0] (V0 is 0 by default)

% Parameters
pars.phi = 6.7e-10;
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 7.5e7;
pars.beta = 50;

% Solve ODE
[T_out, Y] = ode45(@(t,y) SV_ode(t, y, pars), T, y0);

% Extract S and compute equilibrium value
S = Y(:,1);
Seq_num = S(end);

% Plotting
figure;
plot(T_out, S, 'k', 'LineWidth', 2); hold on;
scatter(T_out(end), Seq_num, 80, 'red', 'filled');
legend('S(t)', 'S^* numerically', 'Location', 'best', 'Box', 'on');
xlabel('Time (h)');
ylabel('S cells/mL');
title('SV numerical integration');
set(gca,'FontSize',20);
grid on;

% Save figure
saveas(gcf, 'SV_model_plot.png');



% Analytical equilibrium
Seq_ana = pars.K * (1 - pars.d / pars.b);

% Plotting
figure;
plot(T_out, S, 'k', 'LineWidth', 2); hold on;
scatter(T_out(end), Seq_num, 80, 'red', 'filled');
plot(T_out, Seq_ana * ones(size(T_out)), 'g--', 'LineWidth', 2);
legend('S(t)', 'S^* numerically', 'S^* analytically', 'Location', 'best', 'Box', 'on');
xlabel('Time (h)');
ylabel('S cells/mL');
title('SV numerical integration with equilibrium');
set(gca,'FontSize',20);
grid on;

% Save figure
saveas(gcf, 'SV_model_equilibrium_plot.png');


% Analytical equilibrium
Seq_ana = pars.K * (1 - pars.d / pars.b);

% Compute R0
R0 = (pars.beta * pars.phi * Seq_ana) / (pars.phi * Seq_ana + pars.m);
fprintf('R0 = %.4f\n', R0);



% Updated initial conditions (start from DFE + virus)
y0 = [Seq_ana; 100];

% Solve ODE
[T_out, Y] = ode45(@(t,y) SV_ode(t, y, pars), T, y0);

% Extract solutions
S = Y(:,1);
V = Y(:,2);

% Plot semilog
figure;
semilogy(T_out, S, 'k', 'LineWidth', 2); hold on;
semilogy(T_out, V, 'r', 'LineWidth', 2);
scatter(0, Seq_ana, 80, 'red', 'filled', 'DisplayName', 'DFE');
legend('S(t)', 'V(t)', 'DFE', 'Location', 'best', 'Box', 'on');
xlabel('Time (h)');
ylabel('Density (cells or viruses/mL)');
title('SV numerical integration from DFE');
set(gca,'FontSize',24);
grid on;

% Save figure
saveas(gcf, 'SV_model_DFE_plot.png');





% Take log of virus counts
logV = log(V);

% Select exponential growth phase (e.g., T(2:15))
idx = 2:15;  % MATLAB is 1-indexed
x_reg = T_out(idx)';
y_reg = logV(idx);

% Fit linear regression model
coeffs = polyfit(x_reg, y_reg, 1);
r_num = coeffs(1);  % Slope of the line (exponential rate)
y_fit = polyval(coeffs, x_reg);

% Calculate R^2
SS_res = sum((y_reg - y_fit).^2);
SS_tot = sum((y_reg - mean(y_reg)).^2);
r_sq = 1 - SS_res / SS_tot;

% Display numerical results
fprintf('Numerical estimate:\n');
fprintf('r_num = %.5f\n', r_num);
fprintf('R^2 = %.4f\n', r_sq);

% Plot log(V)
figure;
scatter(T_out, logV, 'filled'); hold on;
scatter(x_reg, y_reg, 60, 'r', 'filled');
plot(x_reg, y_fit, 'b-', 'LineWidth', 2);
xlabel('Time (h)');
ylabel('log(V)');
title('Log(V) and exponential phase linear fit');
legend('log(V)', 'Exponential phase', 'Linear fit', 'Location', 'best');
grid on;
set(gca,'FontSize',24);
saveas(gcf, 'growth_estimate.png');




% Time parameters
t0 = 0;
tf = 5000;
dt = 0.1;
T = t0:dt:tf;

% Shared parameters
pars.m = 1/24;
pars.d = 1/4;
pars.b = 0.95;
pars.K = 7.5e7;
pars.beta = 50;

% Different phi values
phi_values = [1e-12, 1e-9, 1e-5];
titles = {'Low \phi', 'Medium \phi', 'High \phi'};

figure;

for i = 1:3
phi = phi_values(i);
pars.phi = phi;

% Initial condition & equilibrium calculations
if phi == 1e-12
    y0 = [1000; 100];
    Seq_ana = pars.K * (1 - pars.d / pars.b);
    dy = ode45(@(t, y) SV_ode(t, y, pars), T, y0);
    Y = deval(dy, T);
    S = Y(1, :);
    V = Y(2, :);
else
    Seq = pars.m / ((pars.beta - 1) * pars.phi);
    Veq = (1 / pars.phi) * (pars.b * (1 - Seq / pars.K) - pars.d);
    if phi == 1e-5
        y0 = [Seq + 250; Veq + 500];
    else
        y0 = [1000; 100];
    end
    dy = ode45(@(t, y) SV_ode(t, y, pars), T, y0);
    Y = deval(dy, T);
    S = Y(1, :);
    V = Y(2, :);
end

% Plot in subplot
subplot(1, 3, i)
loglog(S, V, 'LineWidth', 2); hold on;
scatter(y0(1), y0(2), 60, 'r', 'filled', 'DisplayName', 'Initial condition');
scatter(S(end), V(end), 60, 'g', 'filled', 'DisplayName', 'Final condition');

if phi == 1e-12
    scatter(Seq_ana, 1e-20, 60, 'b', 'filled', 'DisplayName', 'DFE');
    xline(Seq_ana, 'k--', 'DisplayName', 'S* (DFE)');
else
    scatter(Seq, Veq, 60, 'k', 'filled', 'DisplayName', 'Endemic Eq.');
    xline(Seq, 'k--', 'DisplayName', 'S*');
    yline(Veq, 'k--', 'DisplayName', 'V*');
end

xlabel('Microbe S');
ylabel('Virus V');
title(titles{i});
legend('Location', 'best');
grid on;
set(gca,'YScale','log');
set(gca,'XScale','log');
end

sgtitle('Phase Plane (log-log) for Different \phi Values');
saveas(gcf, 'sv_phase_plot.png');