function fig = plot_node(node,opt)
%
% Plot node information
%
persistent h
if isempty(h)
    h.cnt = 1;
else
    h.cnt = h.cnt + 1;
end

% Parse inputs
fig_idx = getfield_safe(opt,'fig_idx',1);
fig_size = getfield_safe(opt,'fig_size',[0.1,0.4,0.3,0.3]);
title_str = getfield_safe(opt,'title_str','Title string');


nodes4plot = zeros(1,length(node));
for i = 1:length(node)
    if ~isempty(node(i).parent)
        nodes4plot(i) = node(i).parent;
    end
end
p = nodes4plot;
[x,y,height]=treelayout(p);
f = find(p~=0);
pp = p(f);
X = [x(f); x(pp); NaN(size(f))];
Y = [y(f); y(pp); NaN(size(f))];
X = X(:);
Y = Y(:);
n = length(p);

fig = figure(fig_idx); subaxes(fig,1,1,1,0.05,0.15);
hold on;
set_fig_size(fig,fig_size);
set(gcf,'Color','w');
if n < 500
    plot (X, Y, 'k-','MarkerSize',15,'LineWidth',2,'MarkerFaceColor','w');
    n = length(x);
    % Color based on depth
    if isfield(node,'depth')
        colors = jet(max([node.depth])+1);
        for i = 1:n
            ith_depth = node(i).depth;
            plot (x(i),y(i),'ko',...
                'MarkerSize',14,'LineWidth',2,'MarkerFaceColor',colors(ith_depth,:));
            str = sprintf('   %02d (%s)[d=%d]'...
                ,i,node(i).name,ith_depth);
            text(x(i),y(i),str,'FontSize',13,'Interpreter','none');
        end
    else
        for i = 1:n
            plot (x(i),y(i),'ko',...
                'MarkerSize',14,'LineWidth',2,'MarkerFaceColor','g');
            str = sprintf('   %02d (%s)'...
                ,i,node(i).name);
            text(x(i),y(i),str,'FontSize',13,'Interpreter','none');
        end
    end
else
    plot (X, Y, 'k-');
end
xlabel(['height = ' int2str(height)],'FontSize',15);
grid on;
axis([0.1 1 0 1]);
axis off;


% Title
title(title_str,'FontSize',20);
