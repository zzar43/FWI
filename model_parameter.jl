# Setup the model parameters.
# include("draw_graph.jl");

# Grid dimension
h = 5; # step size in meter
Nx, Ny = 101, 101; # grid number

# True velocity model
vel_true = 2000*ones(Nx,Ny);
for i = 1:Nx
    vel_true[i,:] = linspace(1800,2200,Ny);
end
vel_true[:,40:45] = 1850;
vel_true[:,70:end] = 2200;
# vel_true = 2000 * vel_true;
# vel_true[:,30:60] = 2200;
# vel_true[:,41:end] = 2200;
# vel_true[95:105,100:105] = 2200;
# using ImageFiltering
# vel_init = imfilter(vel_true, Kernel.gaussian(10));

# Initial velocity model
vel_init = 2000*ones(Nx,Ny);
for i = 1:Nx
    vel_init[i,:] = linspace(1800,2200,Ny);
end

# Time
sample_fre = 1000; # Hertz
dt = 1/sample_fre;
Nt = 1000;
t = linspace(0,(Nt-1)*dt,Nt);

# Source
function source_ricker(center_fre, center_time, t)
    x = (1 - 2*pi^2*center_fre^2*(t-center_time).^2) .* exp.(-pi^2*center_fre^2*(t-center_time).^2);
    return x;
end;
function source_sine(center_fre, t)
    x = sin.(2*pi*center_fre.*t);
    return x;
end
source_num = 6;
source_coor = zeros(Int,source_num,2);
for i = 1:source_num
    source_coor[i,1] = 20*(i-1)+1;
    source_coor[i,2] = 5;
end
source_vec0 = 100*source_ricker(25, 0.05, t);
# source_vec0 = 100*source_sine(40, t);
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

# Display model
println("Source number: ", source_num)
println("Receiver number: ", receiver_num)
# draw_model(vel_true, vel_init, receiver_coor,source_coor);