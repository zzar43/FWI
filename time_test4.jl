# This file is for the steepest gradient method and line search algorithm.

using NPZ;
@everywhere include("model_parameter.jl");
@everywhere include("2d_wave_solver.jl");
@everywhere include("FWI_time_parallel.jl")

# Make data
println("Start make received data.")
@time true_wavefield, received_data = make_data(vel_true,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);
true_wavefield = [];
println("Received data has been made.")
# Source Loop
println("Start main loop.")
iter_time = 3;
for iter_main = 1:iter_time
    println("Main iter time: ", iter_main)
    @time S, received_data_forward = source_loop(vel_init, Nx, Ny, h, Nt, dt, pml_len, pml_alpha, source_coor, source_vec, receiver_coor, received_data);

    alpha0 = [100; 50; 20];
    @time alpha = line_search(vel_init, S, alpha0[iter_main], vel_true,received_data_forward, received_data, Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor)

    if alpha == 0
        break;
    end

    vel_init = vel_init - alpha * S;
    # draw_real(vel_init)
    # draw_real(S)
end
npzwrite("vel.npy", vel_init);