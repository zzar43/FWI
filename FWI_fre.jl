include("scalar_helmholtz_solver.jl")

# Compute gradient by adjoint method.
function compute_gradient(vel, recorded_data, source_multi, acq_fre)
    Nx_pml = acq_fre.Nx_pml;
    Ny_pml = acq_fre.Ny_pml;
    pml_len = acq_fre.pml_len;
    Nx = acq_fre.Nx;
    Ny = acq_fre.Ny;
    frequency = acq_fre.frequency;
    omega = frequency * 2 * pi;
    fre_num = acq_fre.fre_num;
    source_num = acq_fre.source_num;

    # Initialize
    gradient = zeros(Float32, Nx*Ny, fre_num, source_num);
    # Extend area
    beta, vel_ex = extend_area(vel, acq_fre);
    # Source term
    # Size: Nx_pml-2 * Ny_pml-2
    source_vec = change_source(source_multi, acq_fre);

    # Main loop
    for ind_fre = 1:fre_num
        A = make_diff_operator(h,omega[ind_fre],vel_ex,beta,Nx_pml,Ny_pml);
        F = lufact(A);
        for ind_source = 1:source_num
            source = source_vec[:,ind_fre,ind_source];
            # Forward
            u_forward_vec = F\source;
            u_forward = reshape(u_forward_vec,Nx_pml-2,Ny_pml-2);
            u_forward = u_forward[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny];
            # Adjoint source
            r_forward_vec = acq_fre.projection_op * reshape(u_forward,Nx*Ny,1);
            source_adjoint = conj(r_forward_vec - recorded_data[:,ind_fre,ind_source]);
            source_adjoint0 = zeros(Complex64,Nx_pml-2,Ny_pml-2);
            source_adjoint0[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny] = reshape(source_adjoint,Nx,Ny);
            source_adjoint = reshape(source_adjoint0, (Nx_pml-2)*(Ny_pml-2), 1);
            source_adjoint = -source_adjoint;
            # Backward
            u_back_vec = F\source_adjoint;
            u_back = reshape(u_back_vec,Nx_pml-2,Ny_pml-2);
            u_back = u_back[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny];
            # Gradient
            grad = real(-omega[ind_fre].^2 ./ (vel.^3) .* u_forward .* u_back);
            gradient[:,ind_fre,ind_source] = reshape(grad, Nx*Ny, 1);

            println("Frequency ", frequency[ind_fre], "Hz source ", ind_source, " complete.")
        end
    end
    return gradient
end

# function make_data(vel,Nx,Ny,h,frequency,source_vec,pml_len,pml_alpha,R)
#     size_source = size(source_vec);
#     fre_num = size_source[3];
#     source_num = size_source[4];
#
#     wavefield = zeros(Complex64,Nx,Ny,fre_num,source_num);
#     recorded_vec = zeros(Complex64,Nx*Ny,1,fre_num,source_num);
#
#     for ind_source = 1:source_num
#         for ind_fre = 1:fre_num
#             fre = frequency[ind_fre];
#             source = source_vec[:,:,ind_fre,ind_source];
#             u_vec = AcousticHelmholtzSolver(vel_true,Nx,Ny,h,fre,source,pml_len,pml_alpha,true);
#             u = reshape(u_vec, Nx, Ny);
#             r_vec = R * u_vec;
#             wavefield[:,:,ind_fre,ind_source] = u;
#             recorded_vec[:,:,ind_fre,ind_source] = r_vec;
#             println("Source: ",ind_source, ", frequency: ", fre, " done.")
#         end
#     end
#     return wavefield, recorded_vec
# end


# function compute_gradient(recorded_true_vec,vel,Nx,Ny,h,frequency,source_vec,pml_len,pml_alpha,R)
#     size_source = size(source_vec);
#     fre_num = size_source[3];
#     source_num = size_source[4];
#
#     gradient = zeros(Complex64,Nx,Ny,fre_num,source_num);
#
#     Nx_pml = Nx + 2pml_len;
#     Ny_pml = Ny + 2pml_len;
#     pml_value = linspace(0,pml_alpha,pml_len);
#
#     beta = zeros(Nx_pml,Ny_pml);
#     for i = 1:pml_len
#         beta[pml_len+1-i,:] = pml_value[i];
#         beta[end-pml_len+i,:] = pml_value[i];
#         beta[:,pml_len+1-i] = pml_value[i];
#         beta[:,end-pml_len+i] = pml_value[i];
#     end
#
#     vel_ex = zeros(Nx_pml,Ny_pml);
#     vel_ex[pml_len+1:end-pml_len,pml_len+1:end-pml_len] = vel;
#     for i = 1:pml_len
#         vel_ex[i,:] = vel_ex[pml_len+1,:];
#         vel_ex[end-i+1,:] = vel_ex[end-pml_len,:];
#         vel_ex[:,i] = vel_ex[:,pml_len+1];
#         vel_ex[:,end-i+1] = vel_ex[:,end-pml_len];
#     end
#
#     for ind_source = 1:source_num
#         # Source term
#         source = zeros(Complex64,Nx_pml-2,Ny_pml-2);
#         source_vec0 = source_vec[:,:,ind_fre,ind_source];
#         source[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny] = reshape(source_vec0,Nx,Ny);
#         source_vec0 = reshape(source, (Nx_pml-2)*(Ny_pml-2), 1);
#         # frequency loop
#         for ind_fre = 1:fre_num
#             fre = frequency[ind_fre];
#             omega = fre * 2 * pi;
#             A = make_diff_operator(h,omega,vel_ex,beta,Nx_pml,Ny_pml);
#             # Forward
#             u_forward_vec = A\source_vec0;
#             u_forward = reshape(u_forward_vec,Nx_pml-2,Ny_pml-2);
#             u_forward = u_forward[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny];
#             u_forward_vec = reshape(u_forward, Nx*Ny, 1);
#             r_forward_vec = R * u_forward_vec;
#             # Adjoint
#             # has a -1 ???
#             source_adjoint_vec = -1 * conj(recorded_true_vec[:,:,ind_fre]-r_forward_vec);
#             source_adjoint = zeros(Complex64,Nx_pml-2,Ny_pml-2);
#             source_adjoint[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny] = reshape(source_adjoint_vec,Nx,Ny);
#             source_adjoint_vec = reshape(source_adjoint, (Nx_pml-2)*(Ny_pml-2), 1);
#             # Backward
#             u_back_vec = A\source_adjoint_vec;
#             u_back = reshape(u_back_vec,Nx_pml-2,Ny_pml-2);
#             u_back = u_back[pml_len:pml_len-1+Nx,pml_len:pml_len-1+Ny];
#             # Gradient
#             gradient[:,:,ind_fre,ind_source] = -2 * omega^2 ./ (vel.^3) .* u_forward .* u_back;
#             println("Source: ",ind_source, ", frequency: ", fre, " done.")
#         end
#     end
#     return gradient;
# end

# function compute_gradient(recorded_true_vec,vel,Nx,Ny,h,frequency,source_vec,pml_len,pml_alpha,R)
#
#     size_source = size(source_vec);
#     fre_num = size_source[3];
#     source_num = size_source[4];
#
#     gradient = zeros(Complex64,Nx,Ny,fre_num,source_num);
#
#     for ind_source = 1:source_num
#         for ind_fre = 1:fre_num
#             source_vec0 = source_vec[:,:,ind_fre,ind_source];
#             fre = frequency[ind_fre];
#             omega = 2*pi*fre;
#             # Forward
#             u_forward_vec = AcousticHelmholtzSolver(vel,Nx,Ny,h,fre,source_vec0,pml_len,pml_alpha,true);
#             r_forward_vec = R * u_forward_vec;
#             u_forward = reshape(u_forward_vec, Nx, Ny);
#             # Backward
#             source_adjoint = -conj(recorded_true_vec[:,:,ind_fre]-r_forward_vec);
#             u_back_vec = AcousticHelmholtzSolver(vel,Nx,Ny,h,fre,source_adjoint,pml_len,pml_alpha,true);
#             u_back = reshape(u_back_vec, Nx, Ny);
#             # gradient
#             gradient[:,:,ind_fre,ind_source] = -2omega.^2 ./(vel.^3) .* u_forward .* u_back;
#             println("Source: ",ind_source, ", frequency: ", fre, " done.")
#         end
#     end
#
#     return gradient
# end


# function make_data(vel_true,Nx,Ny,fre,h,source_num,source_func,pml_len,pml_alpha,R)
#     fre_num = length(fre);
#     u_true_vec = Array{Complex128}(Nx*Ny,source_num,fre_num);
#     u_true = Array{Complex128}(Nx,Ny,source_num,fre_num);
#     r_true_vec = Array{Complex128}(Nx*Ny,source_num,fre_num);
#     for ind_fre = 1:fre_num
#         omega = 2*pi*fre[ind_fre];
#         for ind_source = 1:source_num
#             # Build the source
#             source_vec = zeros(Nx,Ny);
#             source_vec[source_coor[ind_source,1], source_coor[ind_source,2]] =source_func[ind_source];
#             source_vec = reshape(source_vec,Nx*Ny,1);
#
#             u_true_vec[:,ind_source,ind_fre] = helmholtz_solver_2d(vel_true,Nx,Ny,omega,h,source_vec,pml_len,pml_alpha,true);
#             u_true[:,:,ind_source,ind_fre] = reshape(u_true_vec[:,ind_source], Nx, Ny);
#             r_true_vec[:,ind_source,ind_fre] = R * u_true_vec[:,ind_source];
#         end
#     end
#
#     if fre_num == 1
#         u_true_vec = u_true_vec[:,:,1];
#         u_true = u_true[:,:,:,1];
#         r_true_vec = r_true_vec[:,:,1];
#     end
#     # return u_true_vec, u_true, r_true_vec
#     return u_true, r_true_vec
# end

# function source_loop(vel_init,Nx,Ny,omega,h,source_coor,source_func,pml_len,pml_alpha,r_true_vec)
#     println("Computing gradient.")
#     gradient = Array{Complex128}(Nx,Ny,source_num);
#     for ind_source = 1:source_num
#         # Build the source
#         source_vec = zeros(Nx,Ny);
#         source_vec[source_coor[ind_source,1], source_coor[ind_source,2]] = source_func[ind_source];
#         source_vec = reshape(source_vec,Nx*Ny,1);
#
#         # Forward
#         u_forward_vec = helmholtz_solver_2d(vel_init,Nx,Ny,omega,h,source_vec,pml_len,pml_alpha,true);
#         u_forward = reshape(u_forward_vec, Nx, Ny);
#         r_forward_vec = R * u_forward_vec;
#
#         # Adjoint source
#         source_adjoint = conj(r_true_vec[:,ind_source]-r_forward_vec);
#
#         # Backward
#         u_back_vec = helmholtz_solver_2d(vel_init,Nx,Ny,omega,h,source_adjoint,pml_len,pml_alpha,true);
#         u_back = reshape(u_back_vec, Nx, Ny);
#
#         # Sensitivity
#         gradient[:,:,ind_source] = -2 * omega.^2 .* (1./vel_init).^3 .* u_forward .* u_back;
#         println("   Source: ", ind_source, " complete.")
#     end
#     gradient = real(gradient);
#     gradient = sum(gradient,3); gradient = gradient[:,:,1];
#     gradient = gradient ./ maximum(gradient);
#
#     println("Gradient Done.")
#     return gradient
# end
#
# function objective_func(vel,Nx,Ny,omega,h,source_num,source_func,pml_len,pml_alpha,R,r_true_vec)
#     # println("start objective func")
#     u_vec = Array{Complex128}(Nx*Ny,source_num);
#     r_vec = Array{Complex128}(Nx*Ny,source_num);
#     for ind_source = 1:source_num
#         # Build the source
#         source_vec = zeros(Nx,Ny);
#         source_vec[source_coor[ind_source,1], source_coor[ind_source,2]] =source_func[ind_source];
#         source_vec = reshape(source_vec,Nx*Ny,1);
#
#         u_vec[:,ind_source] = helmholtz_solver_2d(vel,Nx,Ny,omega,h,source_vec,pml_len,pml_alpha,true);
#         r_vec[:,ind_source] = R * u_vec[:,ind_source];
#     end
#     r_true_vec = real(r_true_vec);
#     r_vec = real(r_vec);
#     J = 1/2*sum((r_true_vec-r_vec).^2);
#     return J;
# end
#
# function backtracking_line_search(vel0,Nx,Ny,omega,h,source_num,source_func,pml_len,pml_alpha,R,r_true_vec,gradient,c,tau,alpha,iter_max)
#     println("Start line search")
#     gradient_vec = reshape(gradient, Nx*Ny, 1);
#     J0 = objective_func(vel0,Nx,Ny,omega,h,source_num,source_func,pml_len,pml_alpha,R,r_true_vec);
#     p = -1*gradient_vec;
#     m = p'*gradient_vec;
#     t = -c*m;
#     alpha = alpha / tau;
#     for i = 1:iter_max
#         alpha = tau * alpha;
#         vel_new = vel0 - alpha * gradient;
#         J = objective_func(vel_new,Nx,Ny,omega,h,source_num,source_func,pml_len,pml_alpha,R,r_true_vec);
#         println("J0 - J is ", J0-J, " αt is ", alpha*t)
#         println("Iter time: ", i, " alpha: ", alpha, " J0: ", J0, " J: ", J)
#         if (J0-J) >= (alpha*t[1])
#             println("Line search done.")
#             break;
#         end
#         println("Line search done.")
#     end
#     return alpha;
# end
