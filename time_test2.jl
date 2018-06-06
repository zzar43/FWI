# This file is for multiple source and multiple receivers.
# Before run this file, model_parameter.jl should be checked.

include("model_parameter.jl");
include("2d_wave_solver.jl");
include("FWI_time.jl")

# Make data
true_wavefield = zeros(Nx,Ny,Nt,source_num);
received_data = zeros(receiver_num,Nt,source_num);
for ind_source = 1:source_num
    @time true_wavefield[:,:,:,ind_source], received_data[:,:,ind_source] = wave_solver_2d_pml(vel_true,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor[ind_source,:]',source_vec[ind_source,:]',receiver_coor);
    println("Source: ", ind_source, " done.")
    # draw_real(true_wavefield[:,:,200,ind_source])
    # draw_real(received_data[:,:,ind_source])
end
# Clear ram
true_wavefield = [];

println("Start source loop.")
# S = zeros(Nx,Ny);
# # Forward and backward
# for ind_source = 1:source_num
#     println("Source: ", ind_source, " begin. Start forward.")
#     @time forward_wavefield, received_data_forward = wave_solver_2d_pml(vel_init,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor[ind_source,:]',source_vec[ind_source,:]',receiver_coor);
#     # draw_real(forward_wavefield[:,:,300])
#     # draw_real(received_data_forward)

#     # adjoint source
#     println("Start backward.")
#     adjoint_source = flipdim(received_data_forward-received_data[:,:,ind_source], 2);
#     draw_real(adjoint_source)
#     @time backward_wavefield, received_data1 = wave_solver_2d_pml(vel_init,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,receiver_coor,adjoint_source,source_coor[ind_source,:]');
#     backward_wavefield = flipdim(backward_wavefield,3);
#     # draw_real(backward_wavefield[:,:,700])

#     println("Build sensitivity")
#     forward_wavefield_tt = zeros(Nx, Ny, Nt);
#     forward_wavefield_tt[:,:,2:Nt-1] = (forward_wavefield[:,:,1:Nt-2] - 2*forward_wavefield[:,:,2:Nt-1] + forward_wavefield[:,:,3:Nt])/(dt^2);
#     S0 = forward_wavefield_tt .* backward_wavefield;
#     S0 = sum(S0,3);
#     S0 = S0[:,:,1];
#     S0 = -2 ./ (vel_init).^3 .* S0;
#     S0 = S0 ./ maximum(abs.(S0));
#     S = S + S0;
#     println("Source: ", ind_source, " done.")
# end
# S = S./maximum(S);

S = source_loop(vel_init, Nx, Ny, h, Nt, dt, pml_len, pml_alpha, source_coor, source_vec, receiver_coor, received_data);

draw_real(S)