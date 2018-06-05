
# This file is for check the wave solver and sensitivity in time domain.
# Before run this file, model_parameter.jl should be checked.
# There should only be one source and one parameter.
# The vel_true and vel_init may be different.
# For the source function, both ricker wavelet and sine wave can be used.

include("model_parameter.jl");
include("2d_wave_solver.jl");

# Make data
ind_source = 1
@time true_wavefield, received_data = wave_solver_2d_pml(vel_true,Nx,Nz,h,Nt,dt,pml_len,pml_alpha,source_coor[ind_source,:]',source_vec[ind_source,:]',receiver_coor);
draw_real(true_wavefield[:,:,300])

# Forward and backward
@time forward_wavefield, received_data_forward = wave_solver_2d_pml(vel_init,Nx,Nz,h,Nt,dt,pml_len,pml_alpha,source_coor[ind_source,:]',source_vec[ind_source,:]',receiver_coor);
draw_real(forward_wavefield[:,:,300])

# adjoint source
adjoint_source = flipdim(received_data_forward-received_data, 2);
@time backward_wavefield, received_data1 = wave_solver_2d_pml(vel_init,Nx,Nz,h,Nt,dt,pml_len,pml_alpha,receiver_coor,adjoint_source,source_coor[ind_source,:]');
backward_wavefield = flipdim(backward_wavefield,3);
draw_real(backward_wavefield[:,:,700])

forward_wavefield_tt = zeros(Nx, Nz, Nt);
forward_wavefield_tt[:,:,2:Nt-1] = (forward_wavefield[:,:,1:Nt-2] - 2*forward_wavefield[:,:,2:Nt-1] + forward_wavefield[:,:,3:Nt])/(dt^2);
S = forward_wavefield_tt .* backward_wavefield;
S = sum(S,3);
S = S[:,:,1];
S = -2 ./ (vel_init).^3 .* S;
S = S ./ maximum(abs.(S));

draw_real(S)
