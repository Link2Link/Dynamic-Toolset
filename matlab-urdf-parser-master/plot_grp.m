function fig = plot_grp(grp,opt)
%
% Plot GRP mean 2std
%

% Parse inputs
n_sample = getfield_safe(opt,'n_sample',10);
seed = getfield_safe(opt,'seed',0);
names = getfield_safe(opt,'names','');
title_str = getfield_safe(opt,'title_str','');

if seed % fix random seed
    rng(seed);
end

fig = figure(); set_fig_size(fig,[0.1,0.4,0.8,0.4]);
subaxes(fig,1,1,1);
hold on; grid on;
colors = hsv(grp.xdim);

% GRP varaince with shaded araes
for d_idx = 1:grp.xdim
    h_fill = fill([grp.t_test;flipud(grp.t_test)],...
        [grp.mu(:,d_idx)-2*grp.std;flipud(grp.mu(:,d_idx))+2*grp.std],...
        colors(d_idx,:),'LineStyle','none'); % grp 2std
    set(h_fill,'FaceAlpha',0.2);
end


% GRP sampled paths
sampled_path_list = sample_grp(grp,n_sample);
for i = 1:n_sample
    for d_idx = 1:grp.xdim
        plot(grp.t_test,sampled_path_list{i}(:,d_idx),...
            '-','Color',0.8*colors(d_idx,:),'linewidth',0.5);
    end
end

% Data points (circles when #data<100, otherwise line)
for d_idx = 1:grp.xdim
    if grp.n_anchor > 100
        % If #data is bigger than 100, do not plot
        plot(grp.t_anchor,grp.x_anchor(:,d_idx),'-',...
            'MarkerSize',12,'linewidth',2,'Color','k'); % plot data
    else
        plot(grp.t_anchor,grp.x_anchor(:,d_idx),'ko',...
            'MarkerSize',12,'linewidth',2,'MarkerFaceColor',colors(d_idx,:)); % plot data
    end
end

% GRP mu
for d_idx = 1:grp.xdim
    hs(d_idx) = plot(grp.t_test,grp.mu(:,d_idx),'-',...
        'linewidth',3,'color',colors(d_idx,:)); % plot GRP mu
    strs{d_idx} = sprintf('%d',d_idx);
end

set(gcf,'color','w');
if isempty(names)
    legend(hs,strs,'fontsize',15,'interpreter','none');
else
    legend(hs,names','fontsize',15,'interpreter','none');
end

% Title
if title_str
    title(title_str,'fontsize',25,'interpreter','none');
end