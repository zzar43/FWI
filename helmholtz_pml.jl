
function make_diff_operator(h,omega,vel,beta,Nx,Ny)
    coef = (1 + im*beta) .* (h^2*omega.^2) ./ (vel.^2);
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

function helmholtz_pml_1(c,Nx,Ny,h,omega,source_coor,source_func,pml_len,pml_alpha)
    Nx0 = Nx + 2pml_len;
    Ny0 = Ny + 2pml_len;
    pml_value = linspace(0,pml_alpha,pml_len);

    beta = zeros(Nx0,Ny0);
    for i = 1:pml_len
        beta[pml_len+1-i,:] = pml_value[i];
        beta[end-pml_len+i,:] = pml_value[i];
        beta[:,pml_len+1-i] = pml_value[i];
        beta[:,end-pml_len+i] = pml_value[i];
    end
    # beta = zeros(Nx,Ny);
    # for i = 1:pml_len
    #     beta[pml_len+1-i,:] = pml_value[i];
    #     beta[end-pml_len+i,:] = pml_value[i];
    #     beta[:,pml_len+1-i] = pml_value[i];
    #     beta[:,end-pml_len+i] = pml_value[i];
    # end

    c_ex = zeros(Nx0,Ny0);
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

    A = make_diff_operator(h,omega,c_ex,beta,Nx0,Ny0);
    source_vec = make_source_vec(source_coor, source_func, h, Nx0, Ny0);

    u = zeros(Complex128,Nx,Ny);
    u_vec = A\source_vec;
    u_vec = reshape(u_vec,Nx0-2,Ny0-2);
    u = u_vec[pml_len-1:end-pml_len+1,pml_len-1:end-pml_len+1];
    return u;
end

function helmholtz_pml_2(c,Nx,Ny,omega,h,source_coor,source_func,pml_len,pml_alpha)
    Nx0 = Nx + 2pml_len;
    Ny0 = Ny + 2pml_len;
    pml_value = linspace(0,pml_alpha,pml_len);

    # Extend velocity
    c_ex = zeros(Nx0,Ny0);
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


    # PML Coef
    sigma_x = zeros(Nx0,Ny0);
    sigma_y = zeros(Nx0,Ny0);
    for i = 1:pml_len
        sigma_x[pml_len+1-i,:] = pml_value[i];
        sigma_x[pml_len+Nx+i,:] = pml_value[i];
        sigma_y[:,pml_len+1-i] = pml_value[i];
        sigma_y[:,pml_len+Ny+i] = pml_value[i];
    end

    Sx = zeros(Complex128, Nx0, Ny0);
    Sy = zeros(Complex128, Nx0, Ny0);
    Sx = 1 ./ (ones(Complex128, Nx0, Ny0) - im .* sigma_x / omega);
    Sy = 1 ./ (ones(Complex128, Nx0, Ny0) - im .* sigma_y / omega);

    coef = (omega./c_ex).^2 - (2/h^2)*(Sx + Sy);
    coef_x = zeros(Complex128,Nx0,Ny0);
    coef_y = zeros(Complex128,Nx0,Ny0);
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
    A = spzeros(Complex128,(Nx0-2)*(Ny0-2),(Nx0-2)*(Ny0-2));
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
    source_vec = make_source_vec(source_coor, source_func, h, Nx0, Ny0);
    source_vec = source_vec / h^2; # Just for the coef on source right.
    u_vec = A\source_vec;
    u = reshape(u_vec,Nx0-2,Ny0-2);
    u = u[pml_len-1:end-pml_len+1,pml_len-1:end-pml_len+1];
    return u;
end 