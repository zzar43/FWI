# This file contains some functions needed in FWI algorithm.

function make_data(vel_true,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor)
    true_wavefield = SharedArray{Float64}(Nx,Ny,Nt,source_num);
    received_data = SharedArray{Float64}(receiver_num,Nt,source_num);
    # println("Start to build the received data.")
    @sync @parallel for ind_source = 1:source_num
        true_wavefield[:,:,:,ind_source], received_data[:,:,ind_source] = wave_solver_2d_pml(vel_true,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor[ind_source,:]',source_vec[ind_source,:]',receiver_coor);
        # println("Source: ", ind_source, " done.")
        # draw_real(true_wavefield[:,:,200,ind_source])
        # draw_real(received_data[:,:,ind_source])
    end
    # println("Received data is built.")
    return true_wavefield, received_data
end


function source_loop(vel_init, Nx, Ny, h, Nt, dt, pml_len, pml_alpha, source_coor, source_vec, receiver_coor, received_data)
    S = SharedArray{Float64}(Nx,Ny,source_num);
    # Forward and backward
    println("Source loop started. Computing sensitivity.")
    forward_wavefield, received_data_forward = make_data(vel_init,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor)
    @sync @parallel for ind_source = 1:source_num
        # forward_wavefield, received_data_forward = wave_solver_2d_pml(vel_init,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor[ind_source,:]',source_vec[ind_source,:]',receiver_coor);
        # draw_real(forward_wavefield[:,:,300])
        # draw_real(received_data_forward)

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
        S0 = S0 ./ maximum(abs.(S0));
        S[:,:,ind_source] = S0;
        println("   Source: ", ind_source, " done.")
    end
    S = sum(S,3);
    S = S[:,:,1];
    S = S./maximum(S);
    println("Finish source loop.")
    return S, received_data_forward;
end

# That is not a backtracking line search method. Just some poor algorithm which works...
function line_search(c_init, S, alpha, vel_true,received_data_forward, received_data, Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor)
    println("Start line search")
    iter_max = 7;
    alpha_vec = SharedArray{Float64}(iter_max,1);
    J = SharedArray{Float64}(iter_max,1);
    alpha_vec[iter_max] = 0;
    alpha_vec[1] = alpha;
    for i = 2:iter_max-1
        alpha_vec[i] = alpha_vec[i-1] / 2;
    end
    J[iter_max] = sum((received_data_forward - received_data).^2);
    @sync @parallel for i = 1:iter_max-1
        alpha = alpha_vec[i]
        wavefield_new, received_data_new = make_data(vel_init-alpha*S,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);
        wavefield_new = [];
        J[i] = sum((received_data_new - received_data).^2);
        # println("alpha = ", alpha_vec[i], " J = ", J[i])
    end
    println("Linear search result alpha: ", alpha_vec)
    println("Misfit function value: ", J)
    alpha = alpha_vec[indmin(J)];
    println("alpha is: ", alpha)
    return alpha
end

# function backtracking_line_search(c_init, S, alpha, vel_true,received_data_forward, received_data, Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor)
#     println("Start line search")
#     J_old = sum((received_data_forward - received_data).^2);
#     wavefield_new, received_data_new = make_data(vel_init-alpha*S,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);
#     wavefield_new = [];
#     J_new = sum((received_data_new - received_data).^2);

#     if J_old <= J_new
#         alpha = 0;
#         println("alpha = ", alpha, " J_new = ", J_new, " J_old = ", J_old)
#     else
#         iter_line = 1;
#         iter_max = 6;
#         println("alpha = ", alpha, " J_new = ", J_new, " J_old = ", J_old)
#         while (iter_line <= iter_max) && (J_new < J_old)
#             alpha = alpha / 2;
#             J_old = J_new;
#             wavefield_new, received_data_new = make_data(vel_init-alpha*S,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);
#             wavefield_new = [];
#             J_new = sum((received_data_new - received_data).^2);
#             iter_line += 1;
#             println("alpha = ", alpha, " J_new = ", J_new, " J_old = ", J_old)
#         end
#     end
#     alpha = 2*alpha;
#     println("Linear search result alpha = ", alpha)
#     return alpha
# end



    
#     while (iter_line <= iter_max) && (J_new <= J_old)
#         if iter_line > 1
#             alpha = alpha / 2;
#             J_old = J_new;
#         end
#         wavefield_new, received_data_new = make_data(vel_init-alpha*S,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);
#         wavefield_new = [];
#         J_new = sum((received_data_new - received_data).^2);
#         iter_line += 1;
#         println("alpha = ", alpha, " J_new = ", J_new, " J_old = ", J_old)
#     end
#     for iter_line = 1:iter_max
#         wavefield_new, received_data_new = make_data(vel_init-alpha*S,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);
#         wavefield_new = [];
#         J_new = sum((received_data_new - received_data).^2);
#     end
#     return alpha*2
# end