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

function draw_model(vel_true,vel_init,receiver_coor,source_coor)
    fig = figure("Model Setting")
    receiver_coor -= 1;
    source_coor -= 1;
    ax1 = subplot(121)
    ax1[:matshow](real(vel_true)',cmap="gray");
    ax1[:scatter](receiver_coor[:,1], receiver_coor[:,2],alpha = 0.2);
    ax1[:scatter](source_coor[:,1], source_coor[:,2]);
    ax2 = subplot(122)
    ax2[:matshow](real(vel_init)',cmap="gray")
    ax2[:scatter](receiver_coor[:,1], receiver_coor[:,2],alpha = 0.2);
    ax2[:scatter](source_coor[:,1], source_coor[:,2]);
    tight_layout()
end