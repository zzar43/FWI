
# ====================
# This file includes two 2 dimensional acoustic wave equation solvers.
# The wave equation is in second order and we are using explicit finite difference method for fast computation.
# The boundary condition for the first solver is u = 0 on the boundary. And the second solver using Perfectly Matched Layers (PML) to simulate free boundary.
# Initial conditions for both solvers are u = 0 and u_t = 0 when t = 0.
# In the PML condition, we set the coefficients in linear relation to the grid. More sophisicated condition may be used in the future.
# ====================

function wave_solver_2d_dirichlet(c,Nx,Ny,h,Nt,dt,source_coor,source_func,receiver_coor)

    u0 = zeros(Nx,Ny);
    u1 = zeros(Nx,Ny);
    u2 = zeros(Nx,Ny);
    snaps_u = zeros(Nx,Ny,Nt);
    receiver_num = size(receiver_coor,1);
    received_data = zeros(receiver_num,Nt);

    lambda = (c*dt/h).^2 .* ones(Nx,Ny);

    source = zeros(Nx,Ny,Nt);
    for i = 1:size(source_coor, 1)
        source[source_coor[i,1], source_coor[i,2], :] = (1/dt)^2*source_func[i,:];
    end

    # Main loop
    for iter_t = 1:Nt
        coef_1 = lambda[2:end-1,2:end-1].*(u1[1:end-2,2:end-1] - 2*u1[2:end-1,2:end-1] + u1[3:end,2:end-1]);
        coef_2 = lambda[2:end-1,2:end-1].*(u1[2:end-1,1:end-2] - 2*u1[2:end-1,2:end-1] + u1[2:end-1,3:end]);
        u2[2:end-1,2:end-1] = 2*u1[2:end-1,2:end-1] - u0[2:end-1,2:end-1] + coef_1 + coef_2 + dt^2*source[2:end-1,2:end-1,iter_t];

        u0[:] = u1; u1[:] = u2;
        snaps_u[:,:,iter_t] = u2;
    end

    for i = 1:receiver_num
        received_data[i,:] = snaps_u[receiver_coor[i,1],receiver_coor[i,2],:];
    end
    return u2, snaps_u, received_data;
end

function wave_solver_2d_pml(c,Nx,Ny,h,Nt,dt,pml_len,pml_alpha,source_coor,source_func,receiver_coor)

    # ====================
    # Build PML Area
    # ====================

    # pml coef with linear relation
    pml_value = linspace(0,pml_alpha,pml_len);
    # pml_value = 1./(pml_alpha-pml_value+1) * 100;
    # pml_value = pml_value.^2;

    sigma_x = zeros(Nx+2*pml_len,Ny+2*pml_len);
    for i = 1:pml_len
        sigma_x[pml_len+1-i,:] = pml_value[i];
        sigma_x[pml_len+Nx+i,:] = pml_value[i];
    end

    sigma_y = zeros(Nx+2*pml_len,Ny+2*pml_len);
    for i = 1:pml_len
        sigma_y[:,pml_len+1-i] = pml_value[i];
        sigma_y[:,pml_len+Ny+i] = pml_value[i];
    end

    # ====================
    # Extend velocity
    # ====================
    c_ex = zeros(Nx+2*pml_len,Ny+2*pml_len);
    for i in 1:pml_len
        c_ex[i,pml_len+1:pml_len+Ny] = c[1,:];
        c_ex[pml_len+Nx+i , pml_len+1 : pml_len+Ny] = c[end,:];
        c_ex[pml_len+1:pml_len+Nx,i] = c[:,1];
        c_ex[pml_len+1 : pml_len+Nx , pml_len+Ny+i] = c[:,end];
    end
    c_ex[1:pml_len,1:pml_len] = c[1,1];
    c_ex[1:pml_len,pml_len+Ny+1:end] = c[1,end];
    c_ex[pml_len+Nx+1:end,1:pml_len] = c[end,1];
    c_ex[pml_len+Nx+1:end,pml_len+Ny+1:end] = c[end,end];
    c_ex[pml_len+1:pml_len+Nx, pml_len+1:pml_len+Ny] = c;

    # ====================
    # Initialize
    # ====================
    u0 = zeros(Nx+2*pml_len,Ny+2*pml_len);
    u1 = zeros(Nx+2*pml_len,Ny+2*pml_len);
    u2 = zeros(Nx+2*pml_len,Ny+2*pml_len);
    snaps_u = zeros(Nx,Ny,Nt);

    vx2 = zeros(Nx+2*pml_len,Ny+2*pml_len);
    vx1 = zeros(Nx+2*pml_len,Ny+2*pml_len);
    vy2 = zeros(Nx+2*pml_len,Ny+2*pml_len);
    vy1 = zeros(Nx+2*pml_len,Ny+2*pml_len);

    receiver_num = size(receiver_coor,1);
    received_data = zeros(receiver_num,Nt);

    A = ones(Nx+2*pml_len,Ny+2*pml_len) ./ c_ex.^2;
    B = (sigma_x + sigma_y) ./ c_ex.^2;
    C = ones(Nx+2*pml_len,Ny+2*pml_len) ./ c_ex.^2;

    source = zeros(Nx+2*pml_len,Ny+2*pml_len,Nt);
    # change source coordinate
    source_coor_ex = source_coor + pml_len*ones(Int,size(source_coor,1),size(source_coor,2));
    # here we changed the maximum value of source by (1/dx)^2 just for convenience
    for i = 1:size(source_coor, 1)
        source[source_coor_ex[i,1], source_coor_ex[i,2], :] = (1/dt)^2*source_func[i,:];
    end

    # ====================
    # Main loop
    # ====================
    for iter_t = 1:Nt

        # original equation
        coef_1 = 2*u1[2:end-1,2:end-1] - u0[2:end-1,2:end-1] - (dt^2.*B[2:end-1,2:end-1])./(A[2:end-1,2:end-1]*dt).*(u1[2:end-1,2:end-1]-u0[2:end-1,2:end-1]) - dt^2.*C[2:end-1,2:end-1]./A[2:end-1,2:end-1].*u1[2:end-1,2:end-1];

        coef_2 = dt^2./(A[2:end-1,2:end-1].*(2h)).*(vx1[3:end,2:end-1] - vx1[1:end-2,2:end-1] + vy1[2:end-1,3:end] - vy1[2:end-1,1:end-2]);

        coef_3 = dt^2./(A[2:end-1,2:end-1]*h^2).*(u1[3:end,2:end-1] - 2*u1[2:end-1,2:end-1] + u1[1:end-2,2:end-1] + u1[2:end-1,3:end] - 2*u1[2:end-1,2:end-1] + u1[2:end-1,1:end-2]);

        u2[2:end-1,2:end-1] = coef_1 + coef_2 + coef_3 + dt^2*source[2:end-1,2:end-1,iter_t];

        # auxiliary equation
        vx2[2:end-1,2:end-1] = vx1[2:end-1,2:end-1] - dt.*sigma_x[2:end-1,2:end-1].*vx1[2:end-1,2:end-1] - dt/(2h).*(sigma_x[2:end-1,2:end-1]-sigma_y[2:end-1,2:end-1]).*(u1[3:end,2:end-1]-u1[1:end-2,2:end-1]);

        vy2[2:end-1,2:end-1] = vy1[2:end-1,2:end-1] - dt.*sigma_y[2:end-1,2:end-1].*vy1[2:end-1,2:end-1] - dt/(2h).*(sigma_y[2:end-1,2:end-1]-sigma_x[2:end-1,2:end-1]).*(u1[2:end-1,3:end]-u1[2:end-1,1:end-2]);

        # time update
        vx1[:] = vx2; vy1[:] = vy2;
        u0[:] = u1; u1[:] = u2;

        # record time domain wavefield
        snaps_u[:,:,iter_t] = u2[pml_len+1:pml_len+Nx,pml_len+1:pml_len+Ny];
    end

    # ====================
    # Record
    # ====================
    for i = 1:receiver_num
        received_data[i,:] = snaps_u[receiver_coor[i,1],receiver_coor[i,2],:];
    end

    # ====================
    # Output
    # ====================
    # return u2, snaps_u, received_data;
    return snaps_u, received_data;

end
