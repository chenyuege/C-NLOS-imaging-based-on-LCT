figure
for i=1:1:292
    subplot(131)
    imagesc(tic_x,tic_y,squeeze(vol(i,:,:)));
    title('Front view');
    set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
    set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
    xlabel('x (m)');
    ylabel('y (m)');
    colormap('gray');
    axis square;
    pause(0.1)
end