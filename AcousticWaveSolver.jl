# This is a new version for time domain acoustic wave solver in 2 dimensional.
# For both Dirichlet boundary condition and pml version boundary condition.

function AcousticWaveSolver2d(vel, Nx, Ny, h, Nt, dt, source_vec, source_coor, receiver_coor)
    u0 = zeros(Nx,Ny);
    u1 = zeros(Nx,Ny);
    u2 = zeros(Nx,Ny);
    receiver_num = size(receiver_coor,1);
    recorded_data = zeros(receiver_num,Nt);
    wavefield = zeros(Nx,Ny,Nt);
    # Prepare coefficients
    lambda = (vel*dt/h).^2;
    coef = zeros(Nx,Ny);
    # Source and receiver index in vector form
    source_coor_vec = source_coor[:,1] + (source_coor[:,2]-1)*Nx;
    receiver_coor_vec = receiver_coor[:,1] + (receiver_coor[:,2]-1)*Nx;
    # Time loop
    for iter_t in 1:Nt
        coef[2:end-1,2:end-1] = lambda[2:end-1,2:end-1] .* (u1[1:end-2,2:end-1] + u1[3:end,2:end-1] + u1[2:end-1,1:end-2] + u1[2:end-1,3:end] - 4*u1[2:end-1,2:end-1]);
        u2[:] = 2*u1 - u0 + coef;
        # Should be += ?????
        u2[source_coor_vec] += vel[source_coor_vec].^2 * dt.*2 .* (-1) .* source_vec[iter_t];
        # update
        u0[:] = u1; u1[:] = u2;
        # record
        wavefield[:,:,iter_t] = u2;
        recorded_data[:,iter_t] = u2[receiver_coor_vec];
    end
    
    return wavefield, recorded_data;
end

function AcousticWaveSolver2d_PML(vel, Nx, Ny, h, Nt, dt, source_vec, source_coor, receiver_coor, pml_alpha, pml_len, forward=true)
    # Extend Model
    local Nx_pml = Nx + 2pml_len;
    local Ny_pml = Ny + 2pml_len;
    local source_coor_ex = source_coor + pml_len;
    local receiver_coor_ex = receiver_coor + pml_len;
    # Extend velocity
    local vel_ex = zeros(Float32,Nx_pml,Ny_pml);
    vel_ex[pml_len+1:end-pml_len,pml_len+1:end-pml_len] = vel;
    for i = 1:pml_len
        vel_ex[i,:] = vel_ex[pml_len+1,:];
        vel_ex[end-i+1,:] = vel_ex[end-pml_len,:];
        vel_ex[:,i] = vel_ex[:,pml_len+1];
        vel_ex[:,end-i+1] = vel_ex[:,end-pml_len];
    end
    # Vectorization source and receiver coordinates
    local source_coor_ex_vec = source_coor_ex[:,1] + (source_coor_ex[:,2]-1)*Nx_pml;
    local receiver_coor_ex_vec = receiver_coor_ex[:,1] + (receiver_coor_ex[:,2]-1)*Nx_pml;

    # PML coefficients
    local pml_value = linspace(0,pml_alpha,pml_len);
    local sigma_x = zeros(Nx_pml,Ny_pml);
    local sigma_y = zeros(Nx_pml,Ny_pml);
    for i = 1:pml_len
        sigma_x[pml_len+1-i,:] = pml_value[i];
        sigma_y[:,pml_len+1-i] = pml_value[i];
        sigma_x[end-pml_len+i,:] = pml_value[i];
        sigma_y[:,end-pml_len+i] = pml_value[i];
    end

    # Coefficients
    local u0 = zeros(Float32,Nx_pml,Ny_pml);
    local u1 = zeros(Float32,Nx_pml,Ny_pml);
    local u2 = zeros(Float32,Nx_pml,Ny_pml);
    local ae_x0 = zeros(Float32,Nx_pml,Ny_pml);
    local ae_x1 = zeros(Float32,Nx_pml,Ny_pml);
    local ae_y0 = zeros(Float32,Nx_pml,Ny_pml);
    local ae_y1 = zeros(Float32,Nx_pml,Ny_pml);
    local part1 = zeros(Float32,Nx_pml,Ny_pml);
    local part2 = zeros(Float32,Nx_pml,Ny_pml);
    local part3 = zeros(Float32,Nx_pml,Ny_pml);
    local part4 = zeros(Float32,Nx_pml,Ny_pml);
    local A = dt*(sigma_x+sigma_y);
    local B = dt.^2 .* sigma_x .* sigma_y;
    local lambda1 = vel_ex.^2 * dt^2 / h^2;
    local lambda2 = vel_ex.^2 * dt^2 / (2h);

    local receiver_num = size(receiver_coor,1);
    local recorded_data = zeros(Float32,receiver_num,Nt);
    local wavefield = zeros(Float32,Nx_pml,Ny_pml,Nt);

    # Time loop
    for iter_t in 1:Nt
        part1[2:end-1,2:end-1] = lambda1[2:end-1,2:end-1] .* (u1[3:end,2:end-1] - 2u1[2:end-1,2:end-1] + u1[1:end-2,2:end-1]);
        part2[2:end-1,2:end-1] = lambda1[2:end-1,2:end-1] .* (u1[2:end-1,3:end] - 2u1[2:end-1,2:end-1] + u1[2:end-1,1:end-2]);
        part3[2:end-1,2:end-1] = lambda2[2:end-1,2:end-1] .* (ae_x1[3:end,2:end-1] - ae_x1[1:end-2,2:end-1]);
        part4[2:end-1,2:end-1] = lambda2[2:end-1,2:end-1] .* (ae_y1[2:end-1,3:end] - ae_y1[2:end-1,1:end-2]);

        u2 = 2u1 - u0 - A.*(u1-u0) - B.*u1 + part1 + part2 + part3 + part4;
        # Source point ????? += ?????
        if forward == true
            u2[source_coor_ex_vec] += vel_ex[source_coor_ex_vec].^2 .* dt^2 .* (-1) .* source_vec[:,iter_t];
        else
            u2[source_coor_ex_vec] = vel_ex[source_coor_ex_vec].^2 .* dt^2 .* (-1) .* source_vec[:,iter_t];
        end

        # auxiliary equation
        ae_x1[2:end-1,2:end-1] = ae_x0[2:end-1,2:end-1] - dt.*sigma_x[2:end-1,2:end-1].*ae_x0[2:end-1,2:end-1] - dt/(2h).*(sigma_x[2:end-1,2:end-1]-sigma_y[2:end-1,2:end-1]).*(u1[3:end,2:end-1]-u1[1:end-2,2:end-1]);
        ae_y1[2:end-1,2:end-1] = ae_y0[2:end-1,2:end-1] - dt.*sigma_y[2:end-1,2:end-1].*ae_y0[2:end-1,2:end-1] - dt/(2h).*(sigma_y[2:end-1,2:end-1]-sigma_x[2:end-1,2:end-1]).*(u1[2:end-1,3:end]-u1[2:end-1,1:end-2]);

        u0[:] = u1; u1[:] = u2;
        ae_x0[:] = ae_x1; ae_y0[:] = ae_y1;

        # record
        wavefield[:,:,iter_t] = u2;
        recorded_data[:,iter_t] = u2[receiver_coor_ex_vec];
    end

    wavefield = wavefield[pml_len+1:end-pml_len,pml_len+1:end-pml_len,:];

    return wavefield, recorded_data;
end