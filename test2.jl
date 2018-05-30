# This file is for test forward modelling

# using PyPlot
using NPZ
@everywhere include("model_parameter.jl");
@everywhere include("2d_wave_solver.jl");
npzwrite("vel_true.npy", vel_true);
npzwrite("source_coor.npy", source_coor);
npzwrite("receiver_coor.npy", receiver_coor);
println("Model saved.")

# plot model
# imshow(vel_true'); colorbar()
# scatter(source_coor[:,1], source_coor[:,2])
# scatter(receiver_coor[:,1], receiver_coor[:,2])
# savefig("model.pdf", format="pdf")

# Make data
println("==============================")
println("Generating received data.")
tic()
r_u = SharedArray{Float64}(receiver_num,Nt,source_num);
@sync @parallel for i = 1:source_num
    wavefield, r_u[:,:,i] = wave_solver_2d_pml(vel_true,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,source_coor[i,:]',source_vec[i,:]',receiver_coor);
    # r_u[:,:,i] = r_u_temp;
    println("Source ", i, " data generated.")
end
toc()
println("Generate received data complete.")
println("==============================\n")

# Main Loop
iter_time = 1;
println("==============================")
println("Start main loop. Iter time is: ", iter_time, ".\n")
for iter_main = 1:iter_time
    println("Iter: ", iter_main, ".")

        dJ = SharedArray{Float64}(Nx,Nz,source_num);
    @sync @parallel for ind_source = 1:source_num
        # solve u_0
        s_u0, r_u0 = wave_solver_2d_pml(vel_init,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,source_coor[ind_source,:]',source_vec[ind_source,:]',receiver_coor);
        # adjoint source
        adjoint_source = flipdim(r_u0-r_u[:,:,ind_source],2);
        # adjoint_source = flipdim(r_u0,2);
        s_u1, r_u1 = wave_solver_2d_pml(vel_init,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,receiver_coor,adjoint_source,source_coor[ind_source,:]');
        s_u1 = flipdim(s_u1,3);

        # adjoint kernal
        # second order derivative of u0
        s_u0_tt = zeros(Nx, Nz, Nt);
        for i = 2:(Nt-1)
            s_u0_tt[:,:,i] = (s_u0[:,:,i-1] - 2s_u0[:,:,i] + s_u0[:,:,i+1])/(dt^2);
        end
        dJ_temp = s_u1 .* s_u0_tt;
        dJ_temp = sum(dJ_temp,3);
        dJ_temp = dJ_temp[:,:,1];
        dJ_temp = -2 ./ (vel_init).^3 .* dJ_temp;
        
        dJ[:,:,source_num] = dJ_temp;

        println("Sensitivity on source ", ind_source, " done.")
    end
    dJ = sum(dJ,3);
    dJ = dJ ./ maximum(abs.(dJ));
    dJ = dJ[:,:,1];
    npzwrite("dJ.npy", dJ)
    println("dJ saved.")

    # Choose a step size
    # use the 1st source to test
    test_source = 1;
    alpha = linspace(1,10,5);
    alpha = alpha.^2;
    val_alpha = SharedArray{Float64}(length(alpha));
    r_u_alpha = SharedArray{Float64}(receiver_num,Nt,source_num);
    @sync @parallel for i = 1:length(alpha)
        wavefield, r_u_alpha[:,:,test_source] = wave_solver_2d_pml(vel_init+alpha[i]*dJ,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,source_coor[test_source,:]',source_vec[test_source,:]',receiver_coor);
        val_alpha[i] = norm(r_u[:,:,test_source] - r_u_alpha[:,:,test_source]);
    end
    println("val_alpha: ", val_alpha, "\n")
    vel_init = vel_init+alpha[indmin(val_alpha)]*dJ;
end
println("==============================")
npzwrite("vel.npy", vel_init)
println("data saved.")
println("All Done!")

