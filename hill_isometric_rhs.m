function dFdt = hill_isometric_rhs(t, F, x, x_star, kse, kpe, b)
    % Activation function A(t) — simple pulse at 100 ms
    if t >= 0.1 && t <= 0.11
        A = 1;
    else
        A = 0;
    end

    % Force-length scaling factor: 1 when x = x_star
    s = exp(-40 * (x - x_star)^2);

    % ODE from the Hill model
    dFdt = (kse / b) * ( kpe * (x - x_star) - (1 + kpe / kse) * F + A * s );
end
