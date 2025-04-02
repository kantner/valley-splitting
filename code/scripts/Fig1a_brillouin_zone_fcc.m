

%% parameters

scaling_factor=3;

%% lists
corner_list = [
    1.5, 0, -3;
    0, 1.5, -3;
    -1.5, 0, 3;
    0, -1.5, 3;
    -1.5, 3, 0;
    0, 3, -1.5;
    1.5, -3, 0;
    0, -3, 1.5;
    3, -1.5, 0;
    3, 0, -1.5;
    -3, 1.5, 0;
    -3, 0, 1.5;
    0, 1.5, 3;
    1.5, 0, 3;
    0, 3, 1.5;
    1.5, 3, 0;
    3, 0, 1.5;
    3, 1.5, 0;
    0, -1.5, -3;
    -1.5, 0, -3;
    0, -3, -1.5;
    -1.5, -3, 0;
    -3, 0, -1.5;
    -3, -1.5, 0
]/scaling_factor;

edge_list = {...
    [0, 1.5, -3 ; 1.5, 0, -3], ...
    [0, -1.5, 3 ; -1.5, 0, 3], ...
    [0, 3, -1.5 ; 0, 1.5, -3], ...
    [0, 3, -1.5 ; -1.5, 3, 0], ...
    [0, -3, 1.5 ; 0, -1.5, 3], ...
    [0, -3, 1.5 ; 1.5, -3, 0], ...
    [3, -1.5, 0 ; 1.5, -3, 0], ...
    [3, 0, -1.5 ; 1.5, 0, -3], ...
    [3, 0, -1.5 ; 3, -1.5, 0], ...
    [-3, 1.5, 0 ; -1.5, 3, 0], ...
    [-3, 0, 1.5 ; -1.5, 0, 3], ...
    [-3, 0, 1.5 ; -3, 1.5, 0], ...
    [0, 1.5, 3 ; -1.5, 0, 3], ...
    [1.5, 0, 3 ; 0, -1.5, 3], ...
    [1.5, 0, 3 ; 0, 1.5, 3], ...
    [0, 3, 1.5 ; -1.5, 3, 0], ...
    [0, 3, 1.5 ; 0, 1.5, 3], ...
    [1.5, 3, 0 ; 0, 3, -1.5], ...
    [1.5, 3, 0 ; 0, 3, 1.5], ...
    [3, 0, 1.5 ; 3, -1.5, 0], ...
    [3, 0, 1.5 ; 1.5, 0, 3], ...
    [3, 1.5, 0 ; 3, 0, -1.5], ...
    [3, 1.5, 0 ; 1.5, 3, 0], ...
    [3, 1.5, 0 ; 3, 0, 1.5], ...
    [0, -1.5, -3 ; 1.5, 0, -3], ...
    [-1.5, 0, -3 ; 0, 1.5, -3], ...
    [-1.5, 0, -3 ; 0, -1.5, -3], ...
    [0, -3, -1.5 ; 1.5, -3, 0], ...
    [0, -3, -1.5 ; 0, -1.5, -3], ...
    [-1.5, -3, 0 ; 0, -3, 1.5], ...
    [-1.5, -3, 0 ; 0, -3, -1.5], ...
    [-3, 0, -1.5 ; -3, 1.5, 0], ...
    [-3, 0, -1.5 ; -1.5, 0, -3], ...
    [-3, -1.5, 0 ; -3, 0, 1.5], ...
    [-3, -1.5, 0 ; -1.5, -3, 0], ...
    [-3, -1.5, 0 ; -3, 0, -1.5], ...   
};


surface_idx_list = {...
    [3,4,14,13], ...
    [11,12,24,23], ...
    [1,2,20,19], ...
    [16,15,5,6], ...
    [17,18,10,9], ...
    [21,22,8,7], ...
    [5,15,13,3,12,11], ...
    [24,12,3,4,8,22], ...
    [8,4,14,17,9,7], ...
    [18,17,14,13,15,16], ...
    [2,6,5,11,23,20], ...
    [20,23,24,22,21,19], ...
    [19,21,7,9,10,1], ...
    [1,10,18,16,6,2], ...
};

%%
%{
for i = 1:size(corner_list,1)
    p = corner_list(i,:);
    plot3(p(1),p(2),p(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    label = string(i);
    text(p(1),p(2),p(3), label, 'FontSize', 17, 'Color', 'blue', 'HorizontalAlignment', 'left');
end
%}

figure(1); clf;
hold on;
for i = 1:size(edge_list,2)
    p = edge_list{i};
    p = p/scaling_factor;
    p1 = p(1,:);
    p2 = p(2,:);
    
    %plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], 'k-', 'LineWidth', 1);
end

delete(findall(gcf,'Type','light'))
%light('Position', [-5, -10, 5], 'Style', 'infinite'); % Add an infinite light source
h=light;

set(h,'Color',[1 0.9 0.9])
set(h,'Position',[-0.5 -1 1])
set(h,'Style','infinite')


%lighting phong
%camlight headlight
lighting gouraud; % Use Gouraud shading for smooth lighting effects


% Adjust material properties
%material shiny; % Make the faces shiny



for i = 1:size(surface_idx_list,2)

    obj=patch('Vertices', corner_list, ...
          'Faces', surface_idx_list{i}, ...
          'FaceColor', [1 1 1]*0.8, ...
          'FaceAlpha', 0.4, ... 
          'EdgeColor', [1 1 1]*0.5, ...
          'EdgeAlpha', 0, ...
          'LineWidth',1);


    %obj.FaceLighting = 'gouraud';
    %{
obj.AmbientStrength = 0.3;
obj.DiffuseStrength = 1.0;
obj.SpecularStrength = 0.9;
obj.SpecularExponent = 25;
obj.SpecularColorReflectance=0.5;
obj.BackFaceLighting = 'reverselit';
    %}
    %
    set(obj,'AmbientStrength',0.5)
    set(obj,'DiffuseStrength',1.0)
    set(obj,'SpecularStrength',0.5)
    set(obj,'SpecularExponent',25)
    set(obj,'SpecularColorReflectance',0.5)
    set(obj,'BackFaceLighting','unlit');
    %}
    %set(obj,'BackFaceLighting','unlit');    
end
%material shiny

% Adjust plot appearance
xlabel('X');
ylabel('Y');
zlabel('Z');
%grid off;
axis equal;
axis off
%view([-68 10]);
%view([-47 15]);
view(3)


% Add lighting
%light('Position', [5, 5, 5], 'Style', 'infinite'); % Add an infinite light source







%ellipsoids
color_x = [0 0.2 1];
color_y = [0 0.2 1];
color_z = [1 0 0.1];

k0 = 0.78;
ellip_wid = 0.10;
ellip_len = 0.2;
ellip_res = 50;




%axes

Xmax = max(abs(corner_list(:,1)));
Ymax = max(abs(corner_list(:,2)));
Zmax = max(abs(corner_list(:,3)));
axisRange = 1.5;%*max([Xmax, Ymax, Zmax]);

arrowSize = 1; % Adjust arrow head size as needed
lineWidth = 1;

%quiver3(0,0,0,axisRange,0,0,'k','LineWidth',lineWidth,'MaxHeadSize',arrowSize);
%line([-axisRange 0],[0 0],[0 0],'Color','k','LineWidth',lineWidth);

%quiver3(0,0,0,0,axisRange,0,'k','LineWidth',lineWidth,'MaxHeadSize',arrowSize);
%line([0 0],[-axisRange 0],[0 0],'Color','k','LineWidth',lineWidth);

%quiver3(0,0,0,0,0,axisRange,'k','LineWidth',lineWidth,'MaxHeadSize',arrowSize);
%line([0 0],[0 0],[-axisRange 0],'Color','k','LineWidth',lineWidth);

%text(1,0,0,'k_x','FontSize',12,'Color','k','HorizontalAlignment','left','VerticalAlignment','middle');
%text(0,1,0,'k_y','FontSize',12,'Color','k','HorizontalAlignment','left','VerticalAlignment','bottom');
%text(0,0,1,'k_z','FontSize',12,'Color',k','HorizontalAlignment','left','VerticalAlignment','bottom');
col_arrow = [1 1 1]*0.4;
arrow_length = 1.55;
arrowHandle =arrow3D([0,0,0], [0,0,arrow_length] ,col_arrow, 0.9);
arrowHandle =arrow3D([0,0,0], [0,0,-arrow_length] , col_arrow, 0.9);
arrowHandle =arrow3D([0,0,0], [0,arrow_length,0] , col_arrow, 0.9);
arrowHandle =arrow3D([0,0,0], [0,-arrow_length,0] , col_arrow, 0.9);
arrowHandle =arrow3D([0,0,0], [arrow_length,0,0] , col_arrow, 0.9);
arrowHandle =arrow3D([0,0,0], [-arrow_length,0,0] , col_arrow, 0.9);

face_alpha = 1;
[x_e, y_e, z_e] = ellipsoid(k0,0,0,ellip_len,ellip_wid,ellip_wid,ellip_res);
surf(x_e, y_e, z_e, 'FaceColor', color_x, 'FaceAlpha', face_alpha, 'EdgeColor', 'none');
[x_e, y_e, z_e] = ellipsoid(-k0,0,0,ellip_len,ellip_wid,ellip_wid,ellip_res);
surf(x_e, y_e, z_e, 'FaceColor', color_x, 'FaceAlpha', face_alpha, 'EdgeColor', 'none');

[x_e, y_e, z_e] = ellipsoid(0,k0,0,ellip_wid,ellip_len,ellip_wid,ellip_res);
surf(x_e, y_e, z_e, 'FaceColor', color_y, 'FaceAlpha', face_alpha, 'EdgeColor', 'none');
[x_e, y_e, z_e] = ellipsoid(0,-k0,0,ellip_wid,ellip_len,ellip_wid,ellip_res);
surf(x_e, y_e, z_e, 'FaceColor', color_y, 'FaceAlpha', face_alpha, 'EdgeColor', 'none');


[x_e, y_e, z_e] = ellipsoid(0,0,k0,ellip_wid,ellip_wid,ellip_len,ellip_res);
surf(x_e, y_e, z_e, 'FaceColor', [color_z], 'FaceAlpha', face_alpha, 'EdgeColor', 'none');
[x_e, y_e, z_e] = ellipsoid(0,0,-k0,ellip_wid,ellip_wid,ellip_len,ellip_res);
surf(x_e, y_e, z_e, 'FaceColor', [color_z], 'FaceAlpha', face_alpha, 'EdgeColor', 'none');


%plot3([0, h.Position(1)],[0 h.Position(2)],[0 h.Position(3)],'r-s')
%%

drawnow


exportgraphics(gcf,'Fig1a_brillouin_zone_fcc.png','ContentType','image','Resolution',1200,'BackgroundColor','white')
%exportgraphics(gcf,'Fig1a_brillouin_zone_fcc.pdf','ContentType','vector')


