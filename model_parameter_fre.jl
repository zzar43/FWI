# Setup the model parameters.
# include("draw_graph.jl");

# ===================================================
# Head function
struct acquisition_fre
    # space and frequency
    Nx::Int64
    Ny::Int64
    h::Float32
    # time
    Nt::Int64
    dt
    t
    # frequency
    frequency::Array{Float32}
    fre_num::Int64
    # source
    source_num::Int64
    source_coor
    # receiver
    receiver_num::Int64
    receiver_coor
    projection_op
    projection_op_pml
    # PML
    pml_len::Int64
    pml_alpha::Float32
    Nx_pml::Int64
    Ny_pml::Int64
end

function source_ricker(center_fre, center_time, t)
    x = (1 - 2*pi^2*center_fre^2*(t-center_time).^2) .* exp.(-pi^2*center_fre^2*(t-center_time).^2);
    return x;
end

function build_source_multi(center_fre,center_time,t,fre_position,source_num,fre_num,ricker=true)
    if ricker == true
        source_time = source_ricker(center_fre, center_time, t);
        source_fre = fft(source_time);
        source_func = source_fre[fre_position];
    else
        source_func = ones(length(fre_position));
    end

    # Outer is source index, middle is frequency index
    source_multi = zeros(Complex64,Nx*Ny,fre_num,source_num);
    for ind_source = 1:source_num
        for ind_fre = 1:fre_num
            source_mat = zeros(Complex64,Nx,Ny);
            source_mat[source_coor[ind_source,1], source_coor[ind_source,2]] = source_func[ind_fre];
            source_mat = reshape(source_mat,Nx*Ny,1);
            source_multi[:,ind_fre,ind_source] = source_mat;
        end
    end

    return source_multi
end

function build_proj_op(Nx,Ny,receiver_coor,receiver_num)
    R = spzeros(Int64,Nx*Ny,Nx*Ny);
    receiver_ind = receiver_coor[:,1] + (receiver_coor[:,2]-1)*Nx;
    for i = 1:receiver_num
        R[receiver_ind[i],receiver_ind[i]] = 1;
    end
    return R
end

function build_proj_op_pml(Nx,Ny,receiver_coor,receiver_num,pml_len)
    # This is for build the adjoint source during the adjoint method
    Nx_pml = Nx + 2pml_len - 2;
    Ny_pml = Ny + 2pml_len - 2;
    R = spzeros(Int64,Nx_pml*Ny_pml,Nx_pml*Ny_pml);
    receiver_coor += (pml_len-1);
    receiver_ind = receiver_coor[:,1] + (receiver_coor[:,2]-1)*Nx_pml;
    for i = 1:receiver_num
        R[receiver_ind[i],receiver_ind[i]] = 1;
    end
    return R
end

function extend_vel(vel, acq_fre)
    # return to the vector version extended velocity
    pml_len = acq_fre.pml_len;
    Nx_pml = acq_fre.Nx + 2pml_len;
    Ny_pml = acq_fre.Ny + 2pml_len;

    vel_ex = zeros(Nx_pml,Ny_pml);
    vel_ex[pml_len+1:end-pml_len,pml_len+1:end-pml_len] = vel;
    for i = 1:pml_len
        vel_ex[i,:] = vel_ex[pml_len+1,:];
        vel_ex[end-i+1,:] = vel_ex[end-pml_len,:];
        vel_ex[:,i] = vel_ex[:,pml_len+1];
        vel_ex[:,end-i+1] = vel_ex[:,end-pml_len];
    end
    vel_ex_vec = reshape(vel_ex,Nx_pml*Ny_pml,1);
    return vel_ex_vec
end

# ===================================================
# Space
# Nx = 101;
# Ny = 101;
# h = 1/100;
# vel_true = 2ones(Float32,Nx,Ny);
# # vel_true[46:55,46:55] = 2.5;
# vel_true[:,33:66] = 2.5;
# vel_true[:,67:end] = 3;
# vel_init = 2ones(Float32,Nx,Ny);

# Read Marmousi
using MAT;
vars = matread("marmousi_dz10.mat");
vel_true = vars["vel"]; vel_true = vel_true.';
vel_true = convert(Array{Float32,2},vel_true)
Nx, Ny = size(vel_true);
h = 10;
using ImageFiltering
vel_init = imfilter(vel_true, Kernel.gaussian(15));
matshow(vel_true')
savefig("vel_true.png")


# PML
pml_len = 50;
pml_alpha = 1;
Nx_pml = Nx + 2pml_len;
Ny_pml = Ny + 2pml_len;

# ===================================================
# Time
sample_fre = 1000.0; # Hertz
dt = 1/sample_fre;
Nt = 1000;
t = linspace(0,(Nt-1)*dt,Nt);
fre = sample_fre * linspace(0,1-1/Nt,Nt);
fre_position = 10:10;
frequency = fre[fre_position];
fre_num = length(frequency);
println("Frequency: ", frequency)

# ===================================================
# Source
source_num = 1;
source_coor = zeros(Int,source_num,2);
for i = 1:source_num
    source_coor[i,1] = 500;
    source_coor[i,2] = 1;
end
# for i = 7:12
#     source_coor[i,1] = 1 + 20*(i-7);
#     source_coor[i,2] = 101;
# end
println("Source number: ", source_num)
source_multi = build_source_multi(15,0.1,t,fre_position,source_num,fre_num,true);
source_vec = zeros(Nx,Ny);
source_vec[20,20] = 1;
source_vec = reshape(source_vec,Nx*Ny,1);

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
R = build_proj_op(Nx,Ny,receiver_coor,receiver_num);
R_pml = build_proj_op_pml(Nx,Ny,receiver_coor,receiver_num,pml_len);

# Display model
# draw_model(vel_true, vel_init, receiver_coor,source_coor);

# Make acquisition
acq_fre = acquisition_fre(Nx,Ny,h,Nt,dt,t,frequency,fre_num,source_num,source_coor,receiver_num,receiver_coor,R,R_pml,pml_len,pml_alpha,Nx_pml,Ny_pml);
