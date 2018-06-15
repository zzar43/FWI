# helmholtz_pml1 is the simplified pml version. With (1+i\beta)k^2 u + \Delta u = -f.
function make_diff_operator(h,omega,vel,beta,Nx,Ny)
    coef = (1 + im*beta) .* (h^2*omega.^2) ./ (vel.^2);
    coef = coef - 4;
    A = spzeros(Complex64, (Nx-2)*(Ny-2), (Nx-2)*(Ny-2));
    # A = zeros((Nx-2)*(Ny-2),(Nx-2)*(Ny-2));
    # Center area
    for i = 2:Nx-3
        for j = 2:Ny-3
            ind_row = (j-1)*(Nx-2)+i;
            A[ind_row,ind_row] = coef[i+1,j+1];
            A[ind_row,ind_row-1] = 1;
            A[ind_row,ind_row+1] = 1;
            A[ind_row,ind_row-Nx+2] = 1;
            A[ind_row,ind_row+Nx-2] = 1;
        end
    end
    # Top
    i = 1;
    for j = 2:Ny-3
        ind_row = (j-1)*(Nx-2)+i;
        A[ind_row,ind_row] = coef[i+1,j+1];
        # A[ind_row,ind_row-1] = 1;
        A[ind_row,ind_row+1] = 1;
        A[ind_row,ind_row-Nx+2] = 1;
        A[ind_row,ind_row+Nx-2] = 1;
    end
    # Bottom
    i = Nx-2;
    for j = 2:Ny-3
        ind_row = (j-1)*(Nx-2)+i;
        A[ind_row,ind_row] = coef[i+1,j+1];
        A[ind_row,ind_row-1] = 1;
        # A[ind_row,ind_row+1] = 1;
        A[ind_row,ind_row-Nx+2] = 1;
        A[ind_row,ind_row+Nx-2] = 1;
    end
    # Left
    j = 1;
    for i = 2:Nx-3
        ind_row = (j-1)*(Nx-2)+i;
        A[ind_row,ind_row] = coef[i+1,j+1];
        A[ind_row,ind_row-1] = 1;
        A[ind_row,ind_row+1] = 1;
        # A[ind_row,ind_row-Nx+2] = 1;
        A[ind_row,ind_row+Nx-2] = 1;
    end
    # Right
    j = Ny-2;
    for i = 2:Nx-3
        ind_row = (j-1)*(Nx-2)+i;
        A[ind_row,ind_row] = coef[i+1,j+1];
        A[ind_row,ind_row-1] = 1;
        A[ind_row,ind_row+1] = 1;
        A[ind_row,ind_row-Nx+2] = 1;
        # A[ind_row,ind_row+Nx-2] = 1;
    end
    # Corner
    i = 1; j = 1;
    ind_row = (j-1)*(Nx-2)+i;
    A[ind_row,ind_row] = coef[i+1,j+1];
    # A[ind_row,ind_row-1] = 1;
    A[ind_row,ind_row+1] = 1;
    # A[ind_row,ind_row-Nx+2] = 1;
    A[ind_row,ind_row+Nx-2] = 1;
    i = Nx-2; j = 1;
    ind_row = (j-1)*(Nx-2)+i;
    A[ind_row,ind_row] = coef[i+1,j+1];
    A[ind_row,ind_row-1] = 1;
    # A[ind_row,ind_row+1] = 1;
    # A[ind_row,ind_row-Nx+2] = 1;
    A[ind_row,ind_row+Nx-2] = 1;
    i = 1; j = Ny-2;
    ind_row = (j-1)*(Nx-2)+i;
    A[ind_row,ind_row] = coef[i+1,j+1];
    # A[ind_row,ind_row-1] = 1;
    A[ind_row,ind_row+1] = 1;
    A[ind_row,ind_row-Nx+2] = 1;
    # A[ind_row,ind_row+Nx-2] = 1;
    i = Nx-2; j = Ny-2;
    ind_row = (j-1)*(Nx-2)+i;
    A[ind_row,ind_row] = coef[i+1,j+1];
    A[ind_row,ind_row-1] = 1;
    # A[ind_row,ind_row+1] = 1;
    A[ind_row,ind_row-Nx+2] = 1;
    # A[ind_row,ind_row+Nx-2] = 1;
    return A;
end;

function acoustic_helmholtz_solver(vel,Nx,Ny,h,fre,source_vec,pml_len,pml_alpha,return_vec=true)
    Nx_pml = Nx + 2pml_len;
    Ny_pml = Ny + 2pml_len;
    omega = fre * 2 * pi;
    pml_value = linspace(0,pml_alpha,pml_len);

    # (1+i\beta)k^2 u + \Delta u = -f.
    beta = zeros(Nx_pml,Ny_pml);
    for i = 1:pml_len
        beta[pml_len+1-i,:] = pml_value[i];
        beta[end-pml_len+i,:] = pml_value[i];
        beta[:,pml_len+1-i] = pml_value[i];
        beta[:,end-pml_len+i] = pml_value[i];
    end

    vel_ex = zeros(Nx_pml,Ny_pml);
    vel_ex[pml_len+1:end-pml_len,pml_len+1:end-pml_len] = vel;
    for i = 1:pml_len
        vel_ex[i,:] = vel_ex[pml_len+1,:];
        vel_ex[end-i+1,:] = vel_ex[end-pml_len,:];
        vel_ex[:,i] = vel_ex[:,pml_len+1];
        vel_ex[:,end-i+1] = vel_ex[:,end-pml_len];
    end

    # Source term
    source = zeros(Complex64,Nx_pml-2,Ny_pml-2);
    source[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny] = reshape(source_vec,Nx,Ny);
    source_vec = reshape(source, (Nx_pml-2)*(Ny_pml-2), 1);
    # source_vec = -1*conj(source_vec);
    source_vec = -1source_vec;

    A = make_diff_operator(h,omega,vel_ex,beta,Nx_pml,Ny_pml);

    u_vec = A\source_vec;
    u = reshape(u_vec,Nx_pml-2,Ny_pml-2);
    u = u[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny];
    if return_vec == true
        u = reshape(u, Nx*Ny, 1);
    end
    return u;
end


function acoustic_helmholtz_solver_pml(vel,Nx,Ny,omega,h,source_vec,pml_len,pml_alpha,return_vec=true)
    Nx0 = Nx + 2pml_len;
    Ny0 = Ny + 2pml_len;
    # pml_alpha = pml_alpha * omega.^2 / (4*pi^2*20*maximum(abs.(vel))) * (100*h);
    # println("pml alpha: ", pml_alpha)
    # pml_alpha = 0.1;
    pml_value = linspace(0,pml_alpha,pml_len);

    # Extend velocity
    vel_ex = zeros(Complex64,Nx0,Ny0);
    for i in 1:pml_len
        vel_ex[i,pml_len+1:pml_len+Ny] = vel[1,:];
        vel_ex[pml_len+Nx+i , pml_len+1 : pml_len+Ny] = vel[end,:];
        vel_ex[pml_len+1:pml_len+Nx,i] = vel[:,1];
        vel_ex[pml_len+1 : pml_len+Nx , pml_len+Ny+i] = vel[:,end];
    end
    vel_ex[1:pml_len,1:pml_len] = vel[1,1];
    vel_ex[1:pml_len,pml_len+Ny+1:end] = vel[1,end];
    vel_ex[pml_len+Nx+1:end,1:pml_len] = vel[end,1];
    vel_ex[pml_len+Nx+1:end,pml_len+Ny+1:end] = vel[end,end];
    vel_ex[pml_len+1:pml_len+Nx, pml_len+1:pml_len+Ny] = vel;

    # PML Coef
    sigma_x = zeros(Complex64,Nx0,Ny0);
    sigma_y = zeros(Complex64,Nx0,Ny0);
    for i = 1:pml_len
        sigma_x[pml_len+1-i,:] = pml_value[i];
        sigma_x[pml_len+Nx+i,:] = pml_value[i];
        sigma_y[:,pml_len+1-i] = pml_value[i];
        sigma_y[:,pml_len+Ny+i] = pml_value[i];
    end

    # Sx = 1 ./ (ones(Complex64, Nx0, Ny0) - im .* sigma_x / omega);
    # Sy = 1 ./ (ones(Complex64, Nx0, Ny0) - im .* sigma_y / omega);
    Sx = 1 ./ (ones(Complex64, Nx0, Ny0) - im .* sigma_x / omega .* vel_ex);
    Sy = 1 ./ (ones(Complex64, Nx0, Ny0) - im .* sigma_y / omega .* vel_ex);


    coef = (omega./vel_ex).^2 - (2/h^2)*(Sx + Sy);
    coef_x = zeros(Complex64,Nx0,Ny0);
    coef_y = zeros(Complex64,Nx0,Ny0);
    coef_x[2:end-1,2:end-1] = (Sx[2:end-1,2:end-1] .* (Sx[3:end,2:end-1]-Sx[1:end-2,2:end-1])) ./ (4*h.^2);
    coef_x[1,:] = coef_x[2,:]; coef_x[end,:] = coef_x[end-1,:];
    coef_x[:,1] = coef_x[:,2]; coef_x[:,end] = coef_x[:,end-1];
    coef_y[2:end-1,2:end-1] = (Sy[2:end-1,2:end-1] .* (Sy[2:end-1,3:end]-Sy[2:end-1,1:end-2])) ./ (4*h.^2);
    coef_y[1,:] = coef_y[2,:]; coef_y[end,:] = coef_y[end-1,:];
    coef_y[:,1] = coef_y[:,2]; coef_y[:,end] = coef_y[:,end-1];
    coef_x1 = coef_x + Sx.^2/h^2;
    coef_x2 = -1*coef_x + Sx.^2/h^2;
    coef_y1 = coef_y + Sy.^2/h^2;
    coef_y2 = -1*coef_y + Sy.^2/h^2;

    # build differential matrix
    A = spzeros(Complex64,(Nx0-2)*(Ny0-2),(Nx0-2)*(Ny0-2));
    for i = 2:Nx0-3
        for j = 2:Ny0-3
            ind_row = (j-1)*(Nx0-2) + i;
            A[ind_row,ind_row] = coef[i+1,j+1];
            A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
            A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
            A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
            A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];
        end
    end
    # Top
    i = 1;
    for j = 2:Ny0-3
        ind_row = (j-1)*(Nx0-2) + i;
        A[ind_row,ind_row] = coef[i+1,j+1];
        # A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
        A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
        A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
        A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];
    end
    # Bottom
    i = Nx0-2;
    for j = 2:Ny0-3
        ind_row = (j-1)*(Nx0-2) + i;
        A[ind_row,ind_row] = coef[i+1,j+1];
        A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
        # A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
        A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
        A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];
    end
    # Left
    j = 1;
    for i = 2:Nx0-3
        ind_row = (j-1)*(Nx0-2) + i;
        A[ind_row,ind_row] = coef[i+1,j+1];
        A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
        A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
        # A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
        A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];
    end
    # Right
    j = Ny0-2;
    for i = 2:Nx0-3
        ind_row = (j-1)*(Nx0-2) + i;
        A[ind_row,ind_row] = coef[i+1,j+1];
        A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
        A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
        A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
        # A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];
    end
    # Top Left
    i = 1; j = 1;
    ind_row = (j-1)*(Nx0-2) + i;
    A[ind_row,ind_row] = coef[i+1,j+1];
    # A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
    A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
    # A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
    A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];
    # Top Right
    i = 1; j = Ny0-2;
    ind_row = (j-1)*(Nx0-2) + i;
    A[ind_row,ind_row] = coef[i+1,j+1];
    # A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
    A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
    A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
    # A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];
    # Bottom Left
    i = Nx0-2; j = 1;
    ind_row = (j-1)*(Nx0-2) + i;
    A[ind_row,ind_row] = coef[i+1,j+1];
    A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
    # A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
    # A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
    A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];
    # Bottom Right
    i = Nx0-2; j = Ny0-2;
    ind_row = (j-1)*(Nx0-2) + i;
    A[ind_row,ind_row] = coef[i+1,j+1];
    A[ind_row,ind_row-1] = coef_x2[i+1,j+1];
    # A[ind_row,ind_row+1] = coef_x1[i+1,j+1];
    A[ind_row,ind_row-Nx0+2] = coef_y2[i+1,j+1];
    # A[ind_row,ind_row+Nx0-2] = coef_y1[i+1,j+1];

    # Source Term
    source = zeros(Complex64,Nx0-2,Ny0-2);
    # source_coor += pml_len;
    # source_ind = source_coor[:,1]-1 + (source_coor[:,2]-2)*(Nx0-2);
    # source[source_ind] = -1*source_func;
    # source_vec = reshape(source, (Nx0-2)*(Ny0-2), 1);
    source[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny] = reshape(source_vec,Nx,Ny);
    source_vec = reshape(source, (Nx0-2)*(Ny0-2), 1);

    u_vec = A\source_vec;
    # Iterative method
    # u_vec = bicgstabl(A,source_vec);

    u = reshape(u_vec,Nx0-2,Ny0-2);
    u = u[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny];
    if return_vec == true
        u = reshape(u, Nx*Ny, 1);
    end
    # received_data = u[receiver_coor[:,1]+(receiver_coor[:,2]-1)*Nx];
    return u
end
