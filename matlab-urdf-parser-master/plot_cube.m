function handler = plot_cube(t,xyz_min,xyz_len,R,opt)
%
% Plot 3-d cube
%
persistent h

% Parse options
color = getfield_safe(opt,'color','b');
alpha = getfield_safe(opt,'alpha',0.5);
sub_idx = getfield_safe(opt,'sub_idx',1);
line_style = getfield_safe(opt,'line_style','none');
line_color = getfield_safe(opt,'line_color','k');

% Make enough handlers at the first
if isempty(h)
    for i = 1:100, h{i}.first_flag = true; end
end

vertex_matrix = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
faces_matrix = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];

vertex_matrix = vertex_matrix.*xyz_len;
vertex_matrix = vertex_matrix + xyz_min; % do basic translation 
vertex_matrix = vertex_matrix*R; % do rotation
vertex_matrix = vertex_matrix + t; % do final translation 

if h{sub_idx}.first_flag
    h{sub_idx}.first_flag = false;    
    h{sub_idx}.patch = patch('Vertices',vertex_matrix,...
        'Faces',faces_matrix,...
        'FaceColor',color,'FaceAlpha',alpha,...
        'LineStyle',line_style,'EdgeColor',line_color,...
        'EdgeAlpha',0.8,'LineWidth',1);
else
    h{sub_idx}.patch.Vertices = vertex_matrix;
end
handler = h{sub_idx}.patch;
