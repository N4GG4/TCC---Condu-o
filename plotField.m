function plotField(mesh, T)

    imagesc(mesh.x, mesh.y, T);
    set(gca,'YDir','normal');     % y crescente para cima
    axis equal tight;
    colorbar;

    xlabel('x (m)');
    ylabel('y (m)');
    title('Distribuição de temperatura (FVM – pixel a pixel)');

end

