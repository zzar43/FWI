# Source Function

function source_ricker(center_fre, center_time, t)
    x = (1 - 2*pi^2*center_fre^2*(t-center_time).^2) .* exp.(-pi^2*center_fre^2*(t-center_time).^2);
end;
function source_gaus(center_fre, center_time, t)
    x =  exp.(-pi^2*center_fre^2*(t-center_time).^2);
end;

source_vec0 = 100*source_ricker(25, 0.05, t);
source_vec = zeros(source_num, Nt);
for i = 1:source_num
    source_vec[i,:] = source_vec0;
end

# Test
# using PyPlot

# sample_fre = 500; # Hertz
# dt = 1/sample_fre;
# Nt = 500;
# t = linspace(0,(Nt-1)*dt,Nt);

# plot(t,source_func)
# gcf()
# clf()

# fre = sample_fre * linspace(0,1-1/Nt,Nt);
# plot(fre,abs(fft(source_func)))
# xlim([20,40])
# gcf()
# clf()