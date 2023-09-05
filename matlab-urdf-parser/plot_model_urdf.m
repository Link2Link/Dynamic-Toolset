function fig = plot_model_urdf(model,opt)
%
% Plot kinematic chain
%
% opt = struct('ball_r',0.05,'coord_r',0.05,'subfig_idx',1,...
%     'view_info',[82,7],'fig_size',[0.1,0.3,0.3,0.5],...
%     'title_str','Initialize Theo','title_fs',25,...
%     'SHOW_TEXT',true,'PLOT_LINKS',true,'PLOT_JOINTS',true,'PLOT_COORDINATES',true...
%     ,'PLOT_ROTATIONAL_AXES',true);
%
persistent h

% Configurations
ball_r = getfield_safe(opt,'ball_r',0.05);
coord_r = getfield_safe(opt,'coord_r',0.05);
fig_idx = getfield_safe(opt,'fig_idx',1);
subfig_idx = getfield_safe(opt,'subfig_idx',1);
view_info = getfield_safe(opt,'view_info',[82,7]);
fig_size = getfield_safe(opt,'fig_size',[0.1,0.4,0.3,0.5]);
title_str = getfield_safe(opt,'title_str','plot_model_urdf');
title_fs = getfield_safe(opt,'title_fs',20);
axis_info = getfield_safe(opt,'axis_info',[-1.0,+1.0,-1.0,+1.0,0.5,2.0]);

SHOW_TEXT = getfield_safe(opt,'SHOW_TEXT',true);
PLOT_LINKS = getfield_safe(opt,'PLOT_LINKS',true);
PLOT_JOINTS = getfield_safe(opt,'PLOT_JOINTS',true);
PLOT_COORDINATES = getfield_safe(opt,'PLOT_COORDINATES',true);
PLOT_ROTATIONAL_AXES = getfield_safe(opt,'PLOT_ROTATIONAL_AXES',true);
PLOT_MESH = getfield_safe(opt,'PLOT_MESH',true);
PLOT_CUBE = getfield_safe(opt,'PLOT_CUBE',true);

% Make enough handlers at the first
if isempty(h)
    for i = 1:10
        for j = 1:5
            h{i,j}.first_flag = true;
        end
    end
end

% Plot using update handler
if h{fig_idx,subfig_idx}.first_flag % first flag
    % Initialize figure
    h{fig_idx,subfig_idx}.fig = figure(fig_idx); hold on; grid on; axis equal;
    view(view_info(1),view_info(2));
    xlabel('X','fontsize',15); ylabel('Y','fontsize',15); zlabel('Z','fontsize',15);
    if ~isempty(axis_info)
        axis(axis_info);
    end
    set_fig_size(h{fig_idx,subfig_idx}.fig,fig_size);
    set(gcf,'color','w');
    
    % ============================= Make handler =============================
    h{fig_idx,subfig_idx}.first_flag = false;
    % Plot all the links
    if PLOT_LINKS
        idx_base = get_matched_idx_in_node_name(model.node,model.base_name);
        for i = 1:model.n_joint % for each joint
            for j = model.node(i).childs % for each child of current joint
                if (i == idx_base) % If it stems from the origin, thin dashed line
                    h{fig_idx,subfig_idx}.links{i,j} = plot3([model.node(i).T(1),model.node(j).T(1)],...
                        [model.node(i).T(2),model.node(j).T(2)],...
                        [model.node(i).T(3),model.node(j).T(3)],...
                        '--','LineWidth',1,'Color',0.4*[1,1,1]);
                else
                    % Think lines
                    switch model.joint_names{j}(1)
                        case 'r' % right part links: red
                            line_color = [0.8,0.3,0.3];
                        case 'l' % left part links: blue
                            line_color = [0.3,0.3,0.8];
                        otherwise % others: black
                            line_color = [0.3,0.3,0.3];
                    end
                    h{fig_idx,subfig_idx}.links{i,j} = plot3([model.node(i).T(1),model.node(j).T(1)],...
                        [model.node(i).T(2),model.node(j).T(2)],...
                        [model.node(i).T(3),model.node(j).T(3)],...
                        '-','LineWidth',3,'Color',line_color);
                end
            end
        end
    end % if PLOT_LINKS
    
    % Plot joints with spheres
    if PLOT_JOINTS
        colorsDepth = jet(max([model.node.depth]));
        for i = 1:model.n_joint
            if isequal(model.joint_types{i},'revolute')
                % Plot only joints that move.
                pos = model.node(i).T;
                r = ball_r;
                col = colorsDepth(model.node(i).depth,:);
                [x,y,z] = ellipsoid(pos(1),pos(2),pos(3),r,r,r,15);
                h{fig_idx,subfig_idx}.joints{i} = surf(x,y,z,'EdgeColor','none','FaceColor',col,'FaceAlpha',0.2);
            end
            if SHOW_TEXT
                h{fig_idx,subfig_idx}.texts{i} = text(model.node(i).T(1),model.node(i).T(2),model.node(i).T(3),...
                    sprintf(' [%d]%s',i,model.joint_names{i}),'fontsize',11,'interpreter','none');
            end
        end
    end % if PLOT_JOINTS
    
    % Plot coordinates
    if PLOT_COORDINATES
        for i = 1:model.n_joint
            if isequal(model.joint_types{i},'revolute'), colorful = true;
            else, colorful = false;
            end
            p = model.node(i).T; R = model.node(i).R; len = coord_r; lw = 2;
            ex = R(:,1); ey = R(:,2); ez = R(:,3);
            if colorful
                h{fig_idx,subfig_idx}.x_coords{i} = quiver3(p(1),p(2),p(3),...
                    coord_r*ex(1),coord_r*ex(2),coord_r*ex(3),...
                    'Color','r','lineWidth',lw,'MaxHeadSize',1.5,'MarkerSize',10);
                h{fig_idx,subfig_idx}.y_coords{i} = quiver3(p(1),p(2),p(3),...
                    coord_r*ey(1),coord_r*ey(2),coord_r*ey(3),...
                    'Color','g','lineWidth',lw,'MaxHeadSize',1.5,'MarkerSize',10);
                h{fig_idx,subfig_idx}.z_coords{i} = quiver3(p(1),p(2),p(3),...
                    coord_r*ez(1),coord_r*ez(2),coord_r*ez(3),...
                    'Color','b','lineWidth',lw,'MaxHeadSize',1.5,'MarkerSize',10);
            else
                h{fig_idx,subfig_idx}.x_coords{i} = plot3([p(1),p(1)+len*ex(1)],[p(2),p(2)+len*ex(2)],[p(3),p(3)+len*ex(3)],...
                    '-','LineWidth',lw,'color',0.5*[1,1,1]);
                h{fig_idx,subfig_idx}.y_coords{i} = plot3([p(1),p(1)+len*ey(1)],[p(2),p(2)+len*ey(2)],[p(3),p(3)+len*ey(3)],...
                    '-','LineWidth',lw,'color',0.5*[1,1,1]);
                h{fig_idx,subfig_idx}.z_coords{i} = plot3([p(1),p(1)+len*ez(1)],[p(2),p(2)+len*ez(2)],[p(3),p(3)+len*ez(3)],...
                    '-','LineWidth',lw,'color',0.5*[1,1,1]);
            end
        end
    end % if PLOT_COORDINATES
    
    % Plot rotation axis
    if PLOT_ROTATIONAL_AXES
        for i = 1:model.n_joint % for all joints
            if isequal(model.joint_types{i},'revolute')
                p = model.node(i).T; R = model.node(i).R;
                ex = R(:,1); ey = R(:,2); ez = R(:,3);
                curr_axis = model.joint_axes{i};
                nz_idx = find(curr_axis~=0); % find the rotational axis
                len = coord_r*1.2; lw = 2;
                switch nz_idx
                    case 1
                        h{fig_idx,subfig_idx}.axes{i} = plot3([p(1),p(1)+len*ex(1)],[p(2),p(2)+len*ex(2)],[p(3),p(3)+len*ex(3)],...
                            'r-','LineWidth',lw);
                    case 2
                        h{fig_idx,subfig_idx}.axes{i} = plot3([p(1),p(1)+len*ey(1)],[p(2),p(2)+len*ey(2)],[p(3),p(3)+len*ey(3)],...
                            'g-','LineWidth',lw);
                    case 3
                        h{fig_idx,subfig_idx}.axes{i} = plot3([p(1),p(1)+len*ez(1)],[p(2),p(2)+len*ez(2)],[p(3),p(3)+len*ez(3)],...
                            'b-','LineWidth',lw);
                end
            end
        end % for i = 1:model.n_joint % for all joints
    end % if PLOT_ROTATIONAL_AXES
    
    % Plot mesh information
    if PLOT_MESH || PLOT_CUBE % mesh from stl files
        camlight('r'); 
        material('dull');
        cube_colors = jet(model.n_joint);
        for i = 1:model.n_joint % for all joint
            % get child link of the current joint
            curr_child_link_name = model.child_links{i};
            curr_child_link_idx = ...
                get_matched_idx_list_in_cell(model.link_names,curr_child_link_name);
            % get the stl of child link (mesh)
            curr_link_filename = model.link_filenames{curr_child_link_idx};
            if curr_link_filename % if stl exists, load and plot
                % get item
                fv = model.link_fvs{curr_child_link_idx};
                p_link = model.link_ps{curr_child_link_idx};
                R_link = model.link_Rs{curr_child_link_idx};
                v_min = model.link_vmins{curr_child_link_idx};
                v_len = model.link_vlens{curr_child_link_idx};
                % plot mesh
                if PLOT_MESH
                    mesh_color = [0.5 0.55 0.6]; % color of the body links
                    mesh_alpha = 0.4; % transparency
                    h{fig_idx,subfig_idx}.patch{i} = ...
                        patch(fv,'FaceColor', mesh_color, ...
                        'EdgeColor', 'none','FaceLighting','gouraud',...
                        'AmbientStrength', 0.15, 'FaceAlpha', mesh_alpha);
                end
                % plot bounding box
                if PLOT_CUBE
                    opt = struct('color',cube_colors(i,:),'alpha',0.1,'sub_idx',i,...
                        'line_style','-','line_color',cube_colors(i,:));
                    plot_cube(p_link',v_min,v_len,R_link',opt);
                end
            end % if curr_link_filename % if stl exists, load and plot
            
        end % for i = 1:model.n_joint % for all joint
    end % if PLOT_MESH || PLOT_CUBE % mesh from stl files
    
    % Title
    if title_str
        h{fig_idx,subfig_idx}.title = title(title_str,'fontsize',title_fs,'interpreter','none');
    end
    
else % not for the first time
    % ============================= Update handler =============================
    
    % Plot all the links
    if PLOT_LINKS
        for i = 1:model.n_joint % for each joint
            for j = model.node(i).childs % for each child of current joint
                h{fig_idx,subfig_idx}.links{i,j}.XData = [model.node(i).T(1),model.node(j).T(1)];
                h{fig_idx,subfig_idx}.links{i,j}.YData = [model.node(i).T(2),model.node(j).T(2)];
                h{fig_idx,subfig_idx}.links{i,j}.ZData = [model.node(i).T(3),model.node(j).T(3)];
            end
        end
    end % if PLOT_LINKS
    
    % Plot joints with spheres
    if PLOT_JOINTS
        for i = 1:model.n_joint
            if isequal(model.joint_types{i},'revolute')
                % Plot only joints that move.
                pos = model.node(i).T;
                r = ball_r;
                [x,y,z] = ellipsoid(pos(1),pos(2),pos(3),r,r,r,15);
                h{fig_idx,subfig_idx}.joints{i}.XData = x;
                h{fig_idx,subfig_idx}.joints{i}.YData = y;
                h{fig_idx,subfig_idx}.joints{i}.ZData = z;
            end
            if SHOW_TEXT
                h{fig_idx,subfig_idx}.texts{i}.Position = model.node(i).T;
                % h{fig_idx,subfig_idx}.texts{i}.String = sprintf(' [%d]%s',i,model.joint_names{i});
            end
        end
    end % if PLOT_JOINTS
    
    % Plot coordinates
    if PLOT_COORDINATES
        for i = 1:model.n_joint
            if isequal(model.joint_types{i},'revolute'), colorful = true;
            else, colorful = false;
            end
            p = model.node(i).T; R = model.node(i).R; len = coord_r; lw = 2;
            ex = R(:,1); ey = R(:,2); ez = R(:,3);
            
            if colorful
                h{fig_idx,subfig_idx}.x_coords{i}.XData = p(1);
                h{fig_idx,subfig_idx}.x_coords{i}.YData = p(2);
                h{fig_idx,subfig_idx}.x_coords{i}.ZData = p(3);
                h{fig_idx,subfig_idx}.x_coords{i}.UData = len*ex(1);
                h{fig_idx,subfig_idx}.x_coords{i}.VData = len*ex(2);
                h{fig_idx,subfig_idx}.x_coords{i}.WData = len*ex(3);
                
                h{fig_idx,subfig_idx}.y_coords{i}.XData = p(1);
                h{fig_idx,subfig_idx}.y_coords{i}.YData = p(2);
                h{fig_idx,subfig_idx}.y_coords{i}.ZData = p(3);
                h{fig_idx,subfig_idx}.y_coords{i}.UData = len*ey(1);
                h{fig_idx,subfig_idx}.y_coords{i}.VData = len*ey(2);
                h{fig_idx,subfig_idx}.y_coords{i}.WData = len*ey(3);
                
                h{fig_idx,subfig_idx}.z_coords{i}.XData = p(1);
                h{fig_idx,subfig_idx}.z_coords{i}.YData = p(2);
                h{fig_idx,subfig_idx}.z_coords{i}.ZData = p(3);
                h{fig_idx,subfig_idx}.z_coords{i}.UData = len*ez(1);
                h{fig_idx,subfig_idx}.z_coords{i}.VData = len*ez(2);
                h{fig_idx,subfig_idx}.z_coords{i}.WData = len*ez(3);
                
            else
                h{fig_idx,subfig_idx}.x_coords{i}.XData = [p(1),p(1)+len*ex(1)];
                h{fig_idx,subfig_idx}.x_coords{i}.YData = [p(2),p(2)+len*ex(2)];
                h{fig_idx,subfig_idx}.x_coords{i}.ZData = [p(3),p(3)+len*ex(3)];
                
                h{fig_idx,subfig_idx}.y_coords{i}.XData = [p(1),p(1)+len*ey(1)];
                h{fig_idx,subfig_idx}.y_coords{i}.YData = [p(2),p(2)+len*ey(2)];
                h{fig_idx,subfig_idx}.y_coords{i}.ZData = [p(3),p(3)+len*ey(3)];
                
                h{fig_idx,subfig_idx}.z_coords{i}.XData = [p(1),p(1)+len*ez(1)];
                h{fig_idx,subfig_idx}.z_coords{i}.YData = [p(2),p(2)+len*ez(2)];
                h{fig_idx,subfig_idx}.z_coords{i}.ZData = [p(3),p(3)+len*ez(3)];
            end
        end
    end % if PLOT_COORDINATES
    
    % Plot rotation axis
    if PLOT_ROTATIONAL_AXES
        for i = 1:model.n_joint % for all joints
            if isequal(model.joint_types{i},'revolute') % onlt for revolute joionts
                p = model.node(i).T; R = model.node(i).R;
                ex = R(:,1); ey = R(:,2); ez = R(:,3);
                curr_axis = model.joint_axes{i};
                nz_idx = find(curr_axis~=0); % find the rotational axis
                len = coord_r*1.2;
                switch nz_idx
                    case 1
                        h{fig_idx,subfig_idx}.axes{i}.XData = [p(1),p(1)+len*ex(1)];
                        h{fig_idx,subfig_idx}.axes{i}.YData = [p(2),p(2)+len*ex(2)];
                        h{fig_idx,subfig_idx}.axes{i}.ZData = [p(3),p(3)+len*ex(3)];
                    case 2
                        h{fig_idx,subfig_idx}.axes{i}.XData = [p(1),p(1)+len*ey(1)];
                        h{fig_idx,subfig_idx}.axes{i}.YData = [p(2),p(2)+len*ey(2)];
                        h{fig_idx,subfig_idx}.axes{i}.ZData = [p(3),p(3)+len*ey(3)];
                    case 3
                        h{fig_idx,subfig_idx}.axes{i}.XData = [p(1),p(1)+len*ez(1)];
                        h{fig_idx,subfig_idx}.axes{i}.YData = [p(2),p(2)+len*ez(2)];
                        h{fig_idx,subfig_idx}.axes{i}.ZData = [p(3),p(3)+len*ez(3)];
                end
            end % if isequal(model.joint_types{i},'revolute') % onlt for revolute joionts
        end % for i = 1:model.n_joint % for all joints
    end % if PLOT_ROTATIONAL_AXES
    
    % Update mesh information
    if PLOT_MESH || PLOT_CUBE % mesh from stl files
        cube_colors = hsv(model.n_joint);
        for i = 1:model.n_joint % for all joint
            % get child link of the current joint
            curr_child_link_name = model.child_links{i};
            curr_child_link_idx = ...
                get_matched_idx_list_in_cell(model.link_names,curr_child_link_name);
            % get the stl of child link (mesh)
            curr_link_filename = model.link_filenames{curr_child_link_idx};
            if curr_link_filename % if stl exists, load and plot
                % get item
                fv = model.link_fvs{curr_child_link_idx};
                p_link = model.link_ps{curr_child_link_idx};
                R_link = model.link_Rs{curr_child_link_idx};
                v_min = model.link_vmins{curr_child_link_idx};
                v_len = model.link_vlens{curr_child_link_idx};
                % plot mesh
                if PLOT_MESH
                    h{fig_idx,subfig_idx}.patch{i}.Vertices = fv.vertices;
                end
                % plot bounding box
                if PLOT_CUBE
                    opt = struct('color',cube_colors(i,:),'alpha',0.1,'sub_idx',i,...
                        'line_style','-','line_color',cube_colors(i,:));
                    plot_cube(p_link',v_min,v_len,R_link',opt);
                end
            end % if curr_link_filename % if stl exists, load and plot
            
        end % for i = 1:model.n_joint % for all joint
    end % if PLOT_MESH || PLOT_CUBE % mesh from stl files
    
    if title_str
        h{fig_idx,subfig_idx}.title.String = title_str;
    end
end % if h{fig_idx,subfig_idx}.first_flag % first flag

rotate3d on

fig = h{fig_idx,subfig_idx}.fig;