# Setup the model parameters.

# Grid dimension
dx = 5; # step size in meter
Nx, Nz = 201, 201; # grid number

# True velocity model
vel_true = ones(Nx,Nz);
vel_true = 2000 * vel_true;
# vel_true[:,30:60] = 2200;
vel_true[:,100:end] = 2200;
# vel_true[95:105,100:105] = 2200;

# Initial velocity model
vel_init = 2000 * ones(Nx,Nz);

# Time
sample_fre = 1000; # Hertz
dt = 1/sample_fre;
Nt = 1000;
t = linspace(0,(Nt-1)*dt,Nt);

# Source
# Source function
function source_ricker(center_fre, center_time, t)
    x = (1 - 2*pi^2*center_fre^2*(t-center_time).^2) .* exp.(-pi^2*center_fre^2*(t-center_time).^2);
    return x;
end;
function source_sine(center_fre, t)
    x = sin.(2*pi*center_fre.*t);
    return x;
end

source_num = 1;
source_coor = zeros(Int,source_num,2);
for i = 1:source_num
    source_coor[i,1] = 50;
    source_coor[i,2] = 50;
end
# for i = 6:10
#     source_coor[i,1] = 16*(i-5);
#     source_coor[i,2] = 100;
# end

source_vec0 = 100*source_ricker(25, 0.05, t);
# source_vec0 = 100*source_sine(50, t);
source_vec = zeros(source_num, Nt);
for i = 1:source_num
    source_vec[i,:] = source_vec0;
end

# Receiver
receiver_num = Nx;
receiver_coor = zeros(Int,receiver_num,2);
for i = 1:receiver_num
    receiver_coor[i,1] = i;
    receiver_coor[i,2] = 1;
end

# PML
pml_len = 20;
pml_alpha = 300;

# Display the model
# using PyPlot;
# pcolor(vel_true)
# scatter(source_coor[:,2], source_coor[:,1])
# scatter(receiver_coor[:,2], receiver_coor[:,1], alpha=0.1)
# gcf()
# clf()