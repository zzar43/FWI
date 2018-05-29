# This file is for test forward modelling

using PyPlot
include("model_parameter.jl");
include("2d_wave_solver.jl");

# plot model
imshow(vel_true'); colorbar()
scatter(source_coor[:,1], source_coor[:,2])
scatter(receiver_coor[:,1], receiver_coor[:,2])
savefig("model.pdf", format="pdf")

# Make data
r_u = zeros(receiver_num,Nt,source_num);
for i = 1:source_num
    wavefield, r_u[:,:,i] = wave_solver_2d_pml(vel_true,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,source_coor[i,:]',source_vec[i,:]',receiver_coor);
    print("\nSource ", i, " data generated.")
end

# Main Loop
iter_time = 3
for iter_main = 1:iter_time

        dJ0 = zeros(Nx,Nz);
    for ind_source = 1:source_num
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
        dJ = s_u1 .* s_u0_tt;
        dJ = sum(dJ,3);
        dJ = dJ[:,:,1];
        dJ = -2 ./ (vel_init).^3 .* dJ;
        dJ = dJ ./ maximum(abs.(dJ));

        dJ0 = dJ0 + dJ;

        print("\nSource ", ind_source, " done.")
    end

    matshow(dJ0', clim=[-1,1], cmap="gray"); colorbar()
    scatter(source_coor[:,1], source_coor[:,2])
    scatter(receiver_coor[:,1], receiver_coor[:,2])
    savefig("dJ.png", format="png")


    # Choose a step size
    # alpha = [10,20,30,40,50];
    alpha = linspace(1,50,5);
    val_alpha = zeros(length(alpha));
    r_u_alpha = zeros(receiver_num,Nt,source_num);
    for i = 1:length(alpha)
        wavefield, r_u_alpha[:,:,1] = wave_solver_2d_pml(vel_init+alpha[i]*dJ0,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,source_coor[1,:]',source_vec[1,:]',receiver_coor);
        val_alpha[i] = norm(r_u[:,:,1] - r_u_alpha[:,:,1]);
        print("\n", i)
    end
    print("\nval_alpha: ", val_alpha)
    vel_init = vel_init+alpha[indmin(val_alpha)]*dJ0;
    # vel_init = vel_init - 10*dJ0;
    matshow(vel_init', cmap="gray"); colorbar()
    scatter(source_coor[:,1], source_coor[:,2])
    scatter(receiver_coor[:,1], receiver_coor[:,2])
    savefig("vel_init.png", format="png")
    
    print("\nMain iter: ", iter_main);
end




# v = maximum(abs.(s_u1)) * 0.2;
# matshow(s_u1[:,:,300]', cmap="gray", clim=[-v,v])



# forward plot
# matshow(s_u[:,:,600]', cmap="gray", clim=[-50,50])
# matshow(r_u', cmap="gray", clim=[-10,10])

# matshow(s_u0[:,:,500]', clim=[-5,5], cmap="gray")

# matshow(adjoint_source', cmap="gray", clim=[-5,5])

# backward plot
# matshow(s_u1[:,:,1]', cmap="gray")
# scatter(source_coor[:,1], source_coor[:,2])
# scatter(receiver_coor[:,1], receiver_coor[:,2])
# plot(r_u1')


