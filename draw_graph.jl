using PyPlot
function draw_abs(mat)
    matshow(abs.(mat'));
    colorbar();
end
function draw_real(mat)
    matshow(real(mat'),cmap="gray");
    colorbar();
end
function draw_imag(mat)
    matshow(imag(mat'));
    colorbar();
end

function draw_model(mat,receiver_coor,source_coor)
    matshow(real(mat)',cmap="gray");
    colorbar()
    scatter(receiver_coor[:,1], receiver_coor[:,2],alpha = 0.2);
    scatter(source_coor[:,1], source_coor[:,2]);
end