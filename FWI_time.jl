# This file contains some functions needed in FWI algorithm.
include("AcousticWaveSolver.jl");

function make_data(vel, Nx, Ny, h, Nt, dt, source_vec, source_coor, receiver_coor, pml_alpha, pml_len)
    # println("Making data.")
    # wavefield = SharedArray{Float32}(Nx,Ny,Nt,source_num);
    # recorded_data = SharedArray{Float32}(receiver_num,Nt,source_num);
    wavefield = zeros(Float32,Nx,Ny,Nt,source_num);
    recorded_data = zeros(Float32,receiver_num,Nt,source_num);
    for ind_source = 1:source_num
        source_vec0 = source_vec[ind_source,:]';
        source_coor0 = source_coor[ind_source,:]';
        wavefield[:,:,:,ind_source], recorded_data[:,:,ind_source] = AcousticWaveSolver2d_PML(vel, Nx, Ny, h, Nt, dt, source_vec0, source_coor0, receiver_coor, pml_alpha, pml_len, true);
        
        # println("   Source ", ind_source, " done.")
    end
    # println("Data has been made.")
    return wavefield, recorded_data;
end

function compute_gradient(vel_init, Nx, Ny, h, Nt, dt, source_vec, source_coor, source_num, receiver_coor, pml_alpha, pml_len, recorded_data_true)
    println("Computing gradient.")
    # gradient = SharedArray{Float32}(Nx,Ny,source_num);
    gradient = zeros(Float32,Nx,Ny,source_num);

    # Forward
    wavefield_forward, recorded_data_forward  = make_data(vel_init, Nx, Ny, h, Nt, dt, source_vec, source_coor, receiver_coor, pml_alpha, pml_len);
    local forward_wavefield_tt = zeros(Float32, Nx, Ny, Nt, source_num);
    forward_wavefield_tt[:,:,2:end-1,:] = (wavefield_forward[:,:,3:end,:] - 2*wavefield_forward[:,:,2:end-1,:] + wavefield_forward[:,:,1:end-2,:])/(dt^2);

    # source loop
    for ind_source = 1:source_num
        source_vec0 = source_vec[ind_source,:]';
        source_coor0 = source_coor[ind_source,:]';

        # Backward
        adjoint_source = flipdim(recorded_data_true[:,:,ind_source] - recorded_data_forward[:,:,ind_source],2);
        wavefield_back, recorded_data_back = AcousticWaveSolver2d_PML(vel_init, Nx, Ny, h, Nt, dt, adjoint_source, receiver_coor, source_coor0, pml_alpha, pml_len, false);
        wavefield_back = flipdim(wavefield_back,3);

        # Gradient
        gradient0 = - forward_wavefield_tt[:,:,:,ind_source] .* wavefield_back;
        gradient0 = sum(gradient0,3);
        gradient0 = gradient0[:,:,1];
        gradient[:,:,ind_source] = gradient0;

        println("   Source ", ind_source, " done.")
    end
    gradient = sum(gradient,3);
    gradient = gradient[:,:,1];
    gradient = gradient ./ maximum(gradient);

    # misfit func
    J = 1/2 * sum((recorded_data_forward - recorded_data_true).^2);

    println("Finish gradient.")
    return J, gradient;
end

function source_loop(vel_init, Nx, Ny, h, Nt, dt, pml_len, pml_alpha, source_coor, source_vec, receiver_coor, received_data)
    S = zeros(Nx,Ny);
    # Forward and backward
    println("Source loop started. Computing sensitivity.")
    forward_wavefield, received_data_forward = make_data(vel_init,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor)
    for ind_source = 1:source_num
        # adjoint source
        # println("Start backward.") 
        adjoint_source = received_data_forward[:,:,ind_source]-received_data[:,:,ind_source];
        adjoint_source = adjoint_source[:,:,1];
        adjoint_source = -flipdim(adjoint_source, 2);
        # draw_real(adjoint_source)
        backward_wavefield, received_data1 = wave_solver_2d_pml(vel_init,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,receiver_coor,adjoint_source,source_coor[ind_source,:]');
        backward_wavefield = flipdim(backward_wavefield,3);
        # draw_real(backward_wavefield[:,:,700])

        # println("Build sensitivity")
        forward_wavefield_tt = zeros(Nx, Ny, Nt);
        forward_wavefield_tt[:,:,2:Nt-1] = (forward_wavefield[:,:,1:Nt-2,ind_source] - 2*forward_wavefield[:,:,2:Nt-1,ind_source] + forward_wavefield[:,:,3:Nt,ind_source])/(dt^2);
        S0 = forward_wavefield_tt .* backward_wavefield;
        S0 = sum(S0,3);
        S0 = S0[:,:,1];
        S0 = -2 ./ (vel_init).^3 .* S0;
        # S0 = S0 ./ maximum(abs.(S0));
        S = S + S0;
        println("   Source: ", ind_source, " done.")
    end
    S = S./maximum(S);
    println("Finish source loop.")
    return S, received_data_forward;
end

# That is not a backtracking line search method. Just some poor algorithm which works...
function line_search(vel, recorded_data_true, vel_new, Nx, Ny, h, Nt, dt, source_vec, source_coor, receiver_coor, pml_alpha, pml_len)
    println("Start line search")
    N_alpha = 10;
    Max_alpha = 200;
    alpha_vec = linspace(Max_alpha,0,N_alpha);
    J_vec = zeros(N_alpha)
    for i = 1:N_alpha-1
        alpha = alpha_vec[i];
        vel_new = vel - alpha * gradient;
        J_vec[i] = misfit_func(recorded_data_true, vel_new, Nx, Ny, h, Nt, dt, source_vec, source_coor, receiver_coor, pml_alpha, pml_len);
    end
    println("alpha: ", alpha_vec)
    println("J: ", J_vec)
    alpha = J_vec[indmin(alpha)]
    return alpha
end

# function backtracking_line_search(vel0, gradient,c,tau,iter_max,alpha,received_data_forward, received_data, Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor)
#     println("Start backtracking line search")
#     J0 = sum((received_data_forward - received_data).^2);
#     gradient_vec = reshape(gradient,Nx*Ny,1);
#     p = -1*gradient_vec;
#     m = p'*gradient_vec;
#     t = -c*m;
#     alpha = alpha / tau;

#     for i = 1:iter_max
#         alpha = tau * alpha;
#         vel_new = vel0 - alpha * gradient;
#         wavefield_new, received_data_new = make_data(vel_new,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);
#         wavefield_new = [];
#         J = sum((received_data_forward - received_data_new).^2);
#         println("J0 - J is ", J0-J, " αt is ", alpha*t[1])
#         println("Iter time: ", i, " alpha: ", alpha, " J0: ", J0, " J: ", J)
#         if (J0-J) >= (alpha*t[1])
#             println("Line search done.")
#             break;
#         end
#         println("Line search done.")
#     end
#     return alpha
# end