using PyPlot
include("model_parameter.jl");
include("2d_wave_solver.jl");

imshow(vel_true'); colorbar()
scatter(source_coor[:,1], source_coor[:,2])
scatter(receiver_coor[:,1], receiver_coor[:,2],alpha=0.1)

# Make data
r_u = zeros(receiver_num,Nt,source_num);
for i = 1:source_num
    wavefield, r_u[:,:,i] = wave_solver_2d_pml(vel_true,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,source_coor[i,:]',source_vec[i,:]',receiver_coor);
    wavefield = 0;
    print("\nSource ", i, " data generated.")
end

dJ0 = zeros(Nx,Nz);
ind_source = 1;
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

matshow(dJ')

matshow(r_u1)