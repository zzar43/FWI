# using PyPlot
include("model_parameter.jl");
include("2d_wave_solver.jl");

# imshow(vel_true'); colorbar()
# scatter(source_coor[:,1], source_coor[:,2])
# scatter(receiver_coor[:,1], receiver_coor[:,2],alpha=0.1)

# Make data

@time u2, snaps_u,rr = wave_solver_2d_dirichlet(vel_true,Nx,Nz,dx,Nt,dt,source_coor,source_vec,receiver_coor);

@code_warntype wave_solver_2d_dirichlet(vel_true,Nx,Nz,dx,Nt,dt,source_coor,source_vec,receiver_coor);

@time wavefield1, r1= wave_solver_2d_pml1(vel_true,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);

@time wavefield2, r2 = wave_solver_2d_pml(vel_true,Nx,Nz,dx,Nt,dt,pml_len,pml_alpha,source_coor,source_vec,receiver_coor);

matshow(r1);
matshow(r2)
matshow(r1-r2)
maximum(abs.(r1-r2))

t = 120;
matshow(wavefield1[:,:,t])
matshow(wavefield2[:,:,t])
#matshow(wavefield1[1:20,1:20,t] - wavefield2[1:20,1:20,t])

 function pisum(w,n)
     u_sum = 0;
     for vi = 1:w
         u_sum = 0;
         for vj = 1:n
             u_sum += 1.0/(vj*vj);
         end
     end
     u_sum
 end