

function make_diff_operator(h,omega,vel,Nx,Ny)
    coef = (h.^2 .* omega.^2) ./ (vel.^2);
    coef = coef - 4;
    A = spzeros(Complex128, (Nx-2)*(Ny-2), (Nx-2)*(Ny-2));
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

function make_source_vec(source_coor, source_func, h, Nx, Ny)
    source = zeros(Complex128,Nx-2,Ny-2);
    source_ind = source_coor[:,1]-1 + (source_coor[:,2]-2)*(Nx-2);
    source[source_ind] = -1*h^2*source_func;
    source_vec = reshape(source, (Nx-2)*(Ny-2), 1);
    return source_vec
end

function wave_solver_2d_helmholtz(Nx, Ny, h, omega, vel, source_coor, source_func)

    A = make_diff_operator(h,omega,vel,Nx,Ny);
    source_vec = make_source_vec(source_coor, source_func, h, Nx, Ny);

    u = zeros(Complex128,Nx,Ny);
    u_vec = A\source_vec;
    u[2:Nx-1,2:Ny-1] = reshape(u_vec,Nx-2,Ny-2);

    return u;
end