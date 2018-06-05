include("helmholtz_pml.jl";)
include("draw_graph.jl")

using IterativeSolvers;


Nx = 101;
Ny = 111;
h = 1/200;
fre = 22.0;
omega = fre*2*pi;
vel_true = ones(Complex128,Nx,Ny);
vel_true[:,85:end] = 1.3;
# vel_true[95:105,81:91] = 0.5;
vel_init = ones(Complex128,Nx,Ny);
source_coor = [51 1];
source_func = [100.0+0im];
# source_func = zeros(Nx,Ny);
# source_func[source_coor[:,1], source_coor[:,2]] = -1e5;
# Receiver
receiver_num = Nx;
receiver_coor = zeros(Int,receiver_num,2);
for i = 1:receiver_num
    receiver_coor[i,1] = i;
    receiver_coor[i,2] = 1;
end

pml_len = 30;
pml_alpha = 2.0;

@time u_true, r_true = helmholtz_pml1(vel_true,Nx,Ny,omega,h,source_coor,source_func,receiver_coor,pml_len,pml_alpha);

@time u1, r1 = helmholtz_pml1(vel_init,Nx,Ny,omega,h,source_coor,source_func,receiver_coor,pml_len,pml_alpha);

r_back = conj(r_true - r1);

@time u2, r2 = helmholtz_pml1(vel_init,Nx,Ny,omega,h,source_coor,r_back,receiver_coor,pml_len,pml_alpha);

# dJ = u1 + u2;
# draw_real(dJ)
# @time u_true, r_true = helmholtz_pml(vel_true::Array{Complex{Float64},2},Nx::Int64,Ny::Int64,omega::Float64,h::Float64,source_coor::Array{Int64,2},source_func::Array{Complex128,1},receiver_coor::Array{Int64,2},pml_len::Int64,pml_alpha::Float64);
# @time u1, r1 = helmholtz_pml(vel_init::Array{Complex{Float64},2},Nx::Int64,Ny::Int64,omega::Float64,h::Float64,source_coor::Array{Int64,2},source_func::Array{Complex128,1},receiver_coor::Array{Int64,2},pml_len::Int64,pml_alpha::Float64);
# # Backward
# r_back = conj(r_true-r1);
# @time u2, r2 = helmholtz_pml(vel_init::Array{Complex{Float64},2},Nx::Int64,Ny::Int64,omega::Float64,h::Float64,receiver_coor::Array{Int64,2},r_true::Array{Complex128,1},source_coor::Array{Int64,2},pml_len::Int64,pml_alpha::Float64);

# F(G(m)-d);
# dJ = 2pi .* u1 .* omega^2 .* u2;
# dJ = dJ ./ maximum(abs.(dJ));

# draw_model(vel_true, receiver_coor, source_coor)
# draw_real(u1);
# draw_real(dJ)

# matshow(real(dJ'), cmap="gray", clim=[-0.5,0.5])