figure
v = 1-vol;
h = vol3d('cdata',v,'texture','2D');
view(3); 
% % Update view since 'texture' = '2D'
vol3d(h);  colormap('gray')
alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')