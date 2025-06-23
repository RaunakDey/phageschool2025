function dydt = SIV_ode(t, y, pars)
    % Unpack state variables
    S = y(1);
    I = y(2);
    V = y(3);

    % Unpack parameters
    phi = pars.phi;
    m = pars.m;
    d = pars.d;
    b = pars.b;
    K = pars.K;
    beta = pars.beta;
    eta = pars.eta;

    % Define ODEs
    dS = b * S * (1 - S / K) - phi * S * V - d * S;
    dI = phi * S * V - eta * I - d * I;
    dV = beta * eta * I - phi * S * V - m * V;

    dydt = [dS; dI; dV];
end