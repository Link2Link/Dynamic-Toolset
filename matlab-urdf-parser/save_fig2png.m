function save_fig2png(fig,png_name,VERBOSE)
%
% Save figure to png file
%

% set(fig,'PaperPositionMode','auto')
% print (fig , '-dpng', figname) ;

% fprintf('Press any key to save.\n');
% pause;

% fig.InvertHardcopy = 'off';
saveas(fig, png_name)
if VERBOSE
    fprintf('%s SAVED. \n', png_name);
end