include("draw_graph.jl")

Nx = 101;
Ny = 101;
h = 1/100;
vel_true = ones(Nx,Ny);
# vel_true[:,50:end] = 1.1;
# vel_true[:,80:end] = 1.2;
using ImageFiltering
vel_init = imfilter(vel_true, Kernel.gaussian(30));

# vel_init = ones(Nx,Ny);

fre = [8];
fre_len = length(fre);
# omega = fre * 2pi;

pml_len = 20;
pml_alpha = 2;

Nx_pml = Nx + 2pml_len;
Ny_pml = Ny + 2pml_len;

# Source and receiver
source_num = 1;
source_coor = zeros(Int,source_num,2);
for i = 1:source_num
    source_coor[i,1] = 20;
    source_coor[i,2] = 1;
end
source_func = 100*ones(source_num,1);
# source_func = [100];

receiver_num = 1;
receiver_coor = zeros(Int,receiver_num,2);
for i = 1:receiver_num
    receiver_coor[i,1] = 50;
    receiver_coor[i,2] = 50;
end

draw_model(vel_true,receiver_coor,source_coor)
draw_model(vel_init,receiver_coor,source_coor)

# Line Search
alpha = linspace(0.01,0.5,5);
alpha = alpha.^2;

function extend_acquisition(c,c_init,source_num,source_func,source_coor,receiver_coor,Nx_pml,Ny_pml,pml_len)
    # Extend velocity
    c_ex = zeros(Nx_pml,Ny_pml);
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

    # Extend velocity
    c_init_ex = zeros(Nx_pml,Ny_pml);
    for i in 1:pml_len
        c_init_ex[i,pml_len+1:pml_len+Ny] = c_init[1,:];
        c_init_ex[pml_len+Nx+i , pml_len+1 : pml_len+Ny] = c_init[end,:];
        c_init_ex[pml_len+1:pml_len+Nx,i] = c_init[:,1];
        c_init_ex[pml_len+1 : pml_len+Nx , pml_len+Ny+i] = c_init[:,end];
    end
    c_init_ex[1:pml_len,1:pml_len] = c_init[1,1];
    c_init_ex[1:pml_len,pml_len+Ny+1:end] = c_init[1,end];
    c_init_ex[pml_len+Nx+1:end,1:pml_len] = c_init[end,1];
    c_init_ex[pml_len+Nx+1:end,pml_len+Ny+1:end] = c_init[end,end];
    c_init_ex[pml_len+1:pml_len+Nx, pml_len+1:pml_len+Ny] = c_init;

    # Coor
    source_coor += pml_len;
    receiver_coor += pml_len;
    source_coor_vec = source_coor[:,1] + (source_coor[:,2]-1)*Nx_pml;
    receiver_coor_vec = receiver_coor[:,1] + (receiver_coor[:,2]-1)*Nx_pml

    source_func_vec = zeros(Complex128,Nx_pml*Ny_pml,source_num);
    for i = 1:source_num
        source_func_vec[source_coor_vec[i],i] = -1*source_func[i];
    end
    return c_ex, c_init_ex, source_coor_vec, receiver_coor_vec, source_func_vec;
end

function make_differential_operator(c_ex,Nx_pml,Ny_pml,pml_len,pml_alpha,omega)
    # PML Coef
    pml_alpha =  pml_alpha * omega.^2 / (4*pi^2*20*maximum(abs.(c_ex))) * (100*h);
    pml_value = linspace(0,pml_alpha,pml_len);

    sigma_x = zeros(Nx_pml,Ny_pml);
    sigma_y = zeros(Nx_pml,Ny_pml);
    for i = 1:pml_len
        sigma_x[pml_len+1-i,:] = pml_value[i];
        sigma_x[pml_len+Nx+i,:] = pml_value[i];
        sigma_y[:,pml_len+1-i] = pml_value[i];
        sigma_y[:,pml_len+Ny+i] = pml_value[i];
    end

    Sx = zeros(Complex128, Nx_pml, Ny_pml);
    Sy = zeros(Complex128, Nx_pml, Ny_pml);
    Sx = 1 ./ (ones(Complex128, Nx_pml, Ny_pml) - im .* sigma_x / omega);
    Sy = 1 ./ (ones(Complex128, Nx_pml, Ny_pml) - im .* sigma_y / omega);

    coef = (omega./c_ex).^2 - (2/h^2)*(Sx + Sy);
    coef_x = zeros(Complex128,Nx_pml,Ny_pml);
    coef_y = zeros(Complex128,Nx_pml,Ny_pml);
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
    A = spzeros(Complex128,Nx_pml*Ny_pml,Nx_pml*Ny_pml);
    for i = 2:Nx_pml-1
        for j = 2:Ny_pml-1
            ind_row = (j-1)*Nx_pml + i;
            A[ind_row,ind_row] = coef[i,j];
            A[ind_row,ind_row-1] = coef_x2[i,j];
            A[ind_row,ind_row+1] = coef_x1[i,j];
            A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
            A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];
        end
    end
    # Top
    i = 1;
    for j = 2:Ny_pml-1
        ind_row = (j-1)*Nx_pml + i;
        A[ind_row,ind_row] = coef[i,j];
        # A[ind_row,ind_row-1] = coef_x2[i,j];
        A[ind_row,ind_row+1] = coef_x1[i,j];
        A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
        A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];
    end
    # Bottom
    i = Nx_pml;
    for j = 2:Ny_pml-1
        ind_row = (j-1)*Nx_pml + i;
        A[ind_row,ind_row] = coef[i,j];
        A[ind_row,ind_row-1] = coef_x2[i,j];
        # A[ind_row,ind_row+1] = coef_x1[i,j];
        A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
        A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];
    end
    # Left
    j = 1;
    for i = 2:Nx_pml-1
        ind_row = (j-1)*Nx_pml + i;
        A[ind_row,ind_row] = coef[i,j];
        A[ind_row,ind_row-1] = coef_x2[i,j];
        A[ind_row,ind_row+1] = coef_x1[i,j];
        # A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
        A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];
    end
    # Right
    j = Ny_pml;
    for i = 2:Nx_pml-1
        ind_row = (j-1)*Nx_pml + i;
        A[ind_row,ind_row] = coef[i,j];
        A[ind_row,ind_row-1] = coef_x2[i,j];
        A[ind_row,ind_row+1] = coef_x1[i,j];
        A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
        # A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];
    end
    # Top Left
    i = 1; j = 1;
    ind_row = (j-1)*Nx_pml + i;
    A[ind_row,ind_row] = coef[i,j];
    # A[ind_row,ind_row-1] = coef_x2[i,j];
    A[ind_row,ind_row+1] = coef_x1[i,j];
    # A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
    A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];
    # Top Right
    i = 1; j = Ny_pml;
    ind_row = (j-1)*Nx_pml + i;
    A[ind_row,ind_row] = coef[i,j];
    # A[ind_row,ind_row-1] = coef_x2[i,j];
    A[ind_row,ind_row+1] = coef_x1[i,j];
    A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
    # A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];
    # Bottom Left
    i = Nx_pml; j = 1;
    ind_row = (j-1)*Nx_pml + i;
    A[ind_row,ind_row] = coef[i,j];
    A[ind_row,ind_row-1] = coef_x2[i,j];
    # A[ind_row,ind_row+1] = coef_x1[i,j];
    # A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
    A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];
    # Bottom Right
    i = Nx_pml; j = Ny_pml;
    ind_row = (j-1)*Nx_pml + i;
    A[ind_row,ind_row] = coef[i,j];
    A[ind_row,ind_row-1] = coef_x2[i,j];
    # A[ind_row,ind_row+1] = coef_x1[i,j];
    A[ind_row,ind_row-Nx_pml] = coef_y2[i,j];
    # A[ind_row,ind_row+Nx_pml] = coef_y1[i,j];

    return A;
end

function vec2mat(u_vec, Nx_pml, Ny_pml, pml_len)
    u = reshape(u_vec,Nx_pml,Ny_pml);
    u = u[pml_len+1:end-pml_len,pml_len+1:end-pml_len];
    return u
end

function locate_receiver(u_vec,receiver_coor_vec)
    r_vec = zeros(Complex128,Nx_pml*Ny_pml,1);
    r_vec[receiver_coor_vec] = u_vec[receiver_coor_vec];
    return r_vec
end

function single_frequency_loop(vel_init_ex, Nx_pml, Ny_pml, pml_alpha, pml_len, omega, source_num, source_func_vec, receiver_num, receiver_coor_vec, r_true_vec)

    # Build forward
    A_1 = make_differential_operator(vel_init_ex,Nx_pml,Ny_pml,pml_len,pml_alpha,omega);
    u_1_vec = zeros(Complex128,Nx_pml*Ny_pml,source_num);
    r_1_vec = zeros(Complex128,Nx_pml*Ny_pml,source_num);
    for ind_source = 1:source_num
        u_1_vec[:,ind_source] = A_1\source_func_vec[:,ind_source];
        r_1_vec[:,ind_source] = locate_receiver(u_1_vec[:,ind_source],receiver_coor_vec);
    end
    # draw_real(vec2mat(u_1_vec[:,1], Nx_pml, Ny_pml, pml_len))
    # draw_real(vec2mat(r_1_vec[:,1], Nx_pml, Ny_pml, pml_len))
    println("   Forward is built.")

    # Build backward
    source_back_vec = -conj(r_true_vec - r_1_vec);
    u_2_vec = zeros(Complex128,Nx_pml*Ny_pml,source_num);
    r_2_vec = zeros(Complex128,Nx_pml*Ny_pml,source_num);
    for ind_source = 1:source_num
        u_2_vec[:,ind_source] = A_1\source_back_vec[:,ind_source];
        # r_2_vec[:,ind_source] = locate_receiver(u_2_vec[:,ind_source],receiver_coor_vec);
    end
    # draw_real(vec2mat(u_2_vec[:,1], Nx_pml, Ny_pml, pml_len))
    println("   Back is built.")

    # Sensitivity

    S_vec = zeros(Nx_pml*Ny_pml,1);
    for ind_source = 1:source_num
        S_vec += omega^2 .* u_1_vec[:,ind_source] .* u_2_vec[:,ind_source];
    end
    S_vec = real(S_vec);
    S_vec = S_vec ./ maximum(S_vec);
    
    println("   Sensitivity is built.")

    return S_vec, r_1_vec, u_2_vec, source_back_vec
end

function line_search(alpha, vel_init_ex, r_true_vec, S, omega)
    search_value = zeros(length(alpha),1);
    for i = 1:length(alpha)
        alpha0 = -alpha[i];
        vel_init_ex1 = zeros(Complex128, Nx_pml, Ny_pml);
        vel_init_ex1[pml_len+1:end-pml_len,pml_len+1:end-pml_len] += alpha0 * S;

        A_search = make_differential_operator(vel_true_ex,Nx_pml,Ny_pml,pml_len,pml_alpha,omega);
        ind_source = 1;
        u_search_vec = A_search\source_func_vec[:,ind_source];
        r_search_vec = locate_receiver(u_search_vec,receiver_coor_vec);

        search_value[i] = norm(r_search_vec-r_true_vec[:,ind_source])
    end
    return alpha[findmin(search_value)[2]]
end

S_vec = zeros(Complex128,Nx_pml*Ny_pml,1);

# Make acquisition
vel_true_ex, vel_init_ex, source_coor_vec, receiver_coor_vec, source_func_vec = extend_acquisition(vel_true,vel_init,source_num,source_func,source_coor,receiver_coor,Nx_pml,Ny_pml,pml_len);
println("Acquisition made.")

for ind_fre = 1:fre_len
    omega = 2 * pi * fre[ind_fre];
    println("Omega is: ", omega)

    # Build true solution
    A_true = make_differential_operator(vel_true_ex,Nx_pml,Ny_pml,pml_len,pml_alpha,omega);
    u_true_vec = zeros(Complex128,Nx_pml*Ny_pml,source_num);
    u_true = zeros(Complex128,Nx_pml,Ny_pml,source_num)
    r_true_vec = zeros(Complex128,Nx_pml*Ny_pml,source_num);
    for ind_source = 1:source_num
        u_true_vec[:,ind_source] = A_true\source_func_vec[:,ind_source];
        r_true_vec[:,ind_source] = locate_receiver(u_true_vec[:,ind_source],receiver_coor_vec);
    end
    println("Solution built.")
    draw_real(vec2mat(u_true_vec[:,1], Nx_pml, Ny_pml, pml_len))

    S_vec, r1, u2, source_back_vec = single_frequency_loop(vel_init_ex, Nx_pml, Ny_pml, pml_alpha, pml_len, omega, source_num, source_func_vec, receiver_num, receiver_coor_vec, r_true_vec);
    S = real(vec2mat(S_vec, Nx_pml, Ny_pml, pml_len));
    draw_real(S)
    println("Frequency ", fre[ind_fre], " Hz has been done.")

    step_size = line_search(alpha, vel_init_ex, r_true_vec, S, omega);
    println("Line search step size: ", step_size)
    vel_init_ex[pml_len+1:end-pml_len,pml_len+1:end-pml_len] += step_size * S;
end



# u2 = vec2mat(u2, Nx_pml, Ny_pml, pml_len)
# r_back = vec2mat(source_back_vec, Nx_pml, Ny_pml, pml_len)
# S = vec2mat(S_vec, Nx_pml, Ny_pml, pml_len)
vel = vec2mat(vel_init_ex, Nx_pml, Ny_pml, pml_len)
draw_real(vel)

# draw_real(vec2mat(u2[:,4], Nx_pml, Ny_pml, pml_len))
# u_back_vec = A_true\conj(r_true_vec);
# u_back = vec2mat(u_back_vec, Nx_pml, Ny_pml, pml_len);

# dJ_vec = omega.^2 .* u_vec_true .* u_back_vec;
# dJ = vec2mat(dJ_vec, Nx_pml, Ny_pml, pml_len)
# draw_real(dJ)