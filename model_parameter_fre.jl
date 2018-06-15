# Setup the model parameters.
# include("draw_graph.jl");
# ===================================================
# Head function
struct acquisition_fre
    # space and frequency
    Nx::Int
    Ny::Int
    h::Float32
    # time
    Nt::Int
    dt
    t
    # frequency
    frequency::Array{Float32}
    fre_num::Int
    # source
    source_num::Int
    # receiver
    receiver_num::Int
    projection_op
    # PML
    pml_len::Int
    pml_alpha::Float32
end

function source_ricker(center_fre, center_time, t)
    x = (1 - 2*pi^2*center_fre^2*(t-center_time).^2) .* exp.(-pi^2*center_fre^2*(t-center_time).^2);
    return x;
end

function build_source_vec(acq_fre,source_coor,center_fre,center_time,fre_position,ricker=true)
    if ricker == true
        source_time = source_ricker(center_fre, center_time, acq_fre.t);
        source_fre = fft(source_time);
        source_func = source_fre[fre_position];
    else
        source_func = ones(length(fre_position));
    end

    source_vec = zeros(Complex64,acq_fre.Nx,acq_fre.Ny,acq_fre.fre_num,acq_fre.source_num);
    # the ith source
    # for ind_source = 1:acq_fre.source_num
    #     for ind_fre = 1:acq_fre.fre_num
    #         source_vec[source_coor[ind_source,1],source_coor[ind_source,2],ind_fre,ind_source] = source_func[ind_fre];
    #     end
    # end
    return source_vec
end

function build_proj_op(Nx,Ny,receiver_coor,receiver_num)
    R = spzeros(Int64,Nx*Ny,Nx*Ny);
    receiver_ind = receiver_coor[:,1] + (receiver_coor[:,2]-1)*Nx;
    for i = 1:receiver_num
        R[receiver_ind[i],receiver_ind[i]] = 1;
    end
    return R
end
# ===================================================
# Space
Nx = 101;
Ny = 101;
h = 1/100;
vel_true = 2ones(Float32,Nx,Ny);
vel_true[46:55,46:55] = 2.5;
# vel_true[:,51:end] = 1.25;
vel_init = 2ones(Float32,Nx,Ny);
# using ImageFiltering
# vel_init = imfilter(vel_true, Kernel.gaussian(15));

# ===================================================
# Time
sample_fre = 1000.0; # Hertz
dt = 1/sample_fre;
Nt = 1000;
t = linspace(0,(Nt-1)*dt,Nt);
fre = sample_fre * linspace(0,1-1/Nt,Nt);
fre_position = 1:35;
frequency = fre[fre_position];
fre_num = length(frequency);
println("Frequency: ", frequency)

# ===================================================
# Source
source_num = 1;
source_coor = zeros(Int,source_num,2);
for i = 1:source_num
    source_coor[i,1] = 51;
    source_coor[i,2] = 51;
end
println("Source number: ", source_num)

# ===================================================
# Receiver
receiver_num = Nx;
receiver_coor = zeros(Int,receiver_num,2);
for i = 1:receiver_num
    receiver_coor[i,1] = i;
    receiver_coor[i,2] = 1;
end
println("Receiver number: ", receiver_num)

# Projection operator
R = spzeros(Int64,Nx*Ny,Nx*Ny);
receiver_ind = receiver_coor[:,1] + (receiver_coor[:,2]-1)*Nx;
for i = 1:receiver_num
    R[receiver_ind[i],receiver_ind[i]] = 1;
end

# PML
pml_len = 30;
pml_alpha = 1;

# Display model
# draw_model(vel_true, vel_init, receiver_coor,source_coor);

# Make acquisition
acq_fre = acquisition_fre(Nx,Ny,h,Nt,dt,t,frequency,fre_num,source_num,receiver_num,R,pml_len,pml_alpha)

# Source vec
source_vec = build_source_vec(acq_fre,source_coor,15,0.1,fre_position,true);
