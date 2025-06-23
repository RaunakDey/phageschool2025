
function dydt = SV_ode(t, y, pars)
    % Standard SV model ODEs
    S = y(1);
    V = y(2);
    phi = pars.phi;
    m = pars.m;
    d = pars.d;
    b = pars.b;
    K = pars.K;
    beta = pars.beta;

    dS = b*S*(1 - S/K) - phi*S*V - d*S;
    dV = beta*phi*S*V - phi*S*V - m*V;

    dydt = [dS; dV];
end

