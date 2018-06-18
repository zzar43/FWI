include("model_parameter_fre.jl");
include("scalar_helmholtz_solver.jl");
include("FWI_fre.jl");
# using NPZ;
using PyPlot;
matshow(real(vel_init')); colorbar()

@time wavefield_true, recorded_data_true = scalar_helmholtz_solver(vel_true, source_multi, acq_fre, "all");
npzwrite("data_wavefield_true.npy", wavefield_true)
npzwrite("data_recorded_true.npy", recorded_data_true)
