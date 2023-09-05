function fig = set_fig_size(fig,fig_sz)
%
% Set figure size based on current screen size
%
sz = get(0, 'ScreenSize'); 
fig_pos = [fig_sz(1)*sz(3),fig_sz(2)*sz(4),fig_sz(3)*sz(3),fig_sz(4)*sz(4)];
set(fig,'Position',fig_pos); 
% set(gcf,'Color','w'); % make figure background white 
