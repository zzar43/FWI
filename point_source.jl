
# ====================
# This file includes some point sources for personal usage.
# ====================

function source_spike(Nt, spike_loc)
    s = zeros(Nt);
    s[spike_loc] = 1;
    return s;
end

function source_sin(Nt, fre)
    x = linspace(0,1-1/Nt,Nt);
    s = sin.(2*pi*fre*x);
    return s;
end

function source_ricker(Nt, dt, center_point, sigma)
    t = linspace(0,(Nt-1)*dt,Nt);
    s = (ones(Nt)-((t-center_point*dt)/sigma).^2) .* exp.(-(t-center_point*dt).^2/(2*sigma^2));
    s = s / maximum(s);
    return s;
end
