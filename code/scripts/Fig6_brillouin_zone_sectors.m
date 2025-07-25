

% high-symmetry points for Brillouin zone segment for [0 0 k0] in units of 2pi/a0
  X      = [0 0 1]';
  Gamma  = [0 0 0]';

  L{1} = [+0.5 +0.5 +0.5]';
  L{2} = [-0.5 +0.5 +0.5]';
  L{3} = [-0.5 -0.5 +0.5]';
  L{4} = [+0.5 -0.5 +0.5]';

  W{1}     = [+0.5, 0, 1]';
  W{2}     = [0, +0.5, 1]';
  W{3}     = [-0.5, 0, 1]';  
  W{4}     = [0, -0.5, 1]';

  K{1} = [0.75 0 0.75]';
  K{2} = [0 0.75 0.75]';
  K{3} = [-0.75 0 0.75]';
  K{4} = [0 -0.75 0.75]';  

% top facet: four W points  
  facet{1} = [W{1}, W{2}, W{3}, W{4}];

% upper side facets: L, two K and two W points
  facet{2} = [L{1}, K{1}, W{1}, W{2}, K{2}];
  facet{3} = [L{2}, K{2}, W{2}, W{3}, K{3}];
  facet{4} = [L{3}, K{3}, W{3}, W{4}, K{4}];
  facet{5} = [L{4}, K{4}, W{4}, W{1}, K{1}];
  
% lower side facets: Gamma, two L and one K point
  facet{6} = [Gamma, L{1}, K{2}, L{2}];
  facet{7} = [Gamma, L{2}, K{3}, L{3}];
  facet{8} = [Gamma, L{3}, K{4}, L{4}];
  facet{9} = [Gamma, L{4}, K{1}, L{1}];
  

% unit vectors
  ex = [1 0 0]';
  ey = [0 1 0]';
  ez = [0 0 1]';

% rotation matrix
  Rz = @(angle) [cos(angle), -sin(angle), 0; sin(angle), cos(angle), 0; 0 0 1];
  Ry = @(angle) [cos(angle), 0, -sin(angle); 0, 1, 0; sin(angle), 0, cos(angle)];
  Rx = @(angle) [1 0 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];

% angle list
  angles{1} = [ 0 0 0];
  angles{2} = [+1 0 0]*pi/2;
  angles{3} = [-1 0 0]*pi/2;
  angles{4} = [0 +1 0]*pi/2;
  angles{5} = [0 -1 0]*pi/2;
  angles{6} = [2 0 0]*pi/2;

% shift vector list  
  shift{1} = [0 0 +1]';
  shift{2} = [0 -1 0]';
  shift{3} = [0 +1 0]';
  shift{4} = [-1 0 0]';
  shift{5} = [+1 0 0]';
  shift{6} = [0 0 -1]';


% color list
  color{1} = [0.8 0.1 0.3]*1.1;
  color{2} = [1 1 1]*0.8;
  color{3} = color{2};
  color{4} = color{2};
  color{5} = color{2};
  color{6} = color{1};



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%    
% plot Brillouin zone sectors
  figure(61); clf; hold all;

% axes
  view([123 9]);
  axis equal
  axis off

% sweep over all sectors and facets  
  for s = 1 : 6

    % set rotation
      R = Rx(angles{s}(1)) * Ry(angles{s}(2)) * Rz(angles{s}(3));
      
   % set shift
     shift_scale = 0.1; % in units of 2pi/a0
    % plot facets  
      for i = 1 : length(facet)      
        pts = shift_scale * shift{s} + R * facet{i};
        p = patch(ex'*pts, ey'*pts, ez'*pts,'k-','LineWidth',1);
        set(p,'FaceColor',color{s},'FaceAlpha',1.0,'EdgeColor',color{s}*0.5,'EdgeAlpha',1.0);
        %set(p,'EdgeColor','none')
      end
  end


% draw axes
  col_arrow    = [1 1 1]*0.4;
  arrow_length = 3.;
  tipLength    = 0.2;
  tipWidth     = 0.06;
  stemWidth    = 0.007;
  %arrowHandle =arrow3D([0,0,0], [0,0,arrow_length] ,col_arrow, 0.9);
  %arrowHandle =arrow3D([0,0,0], [0,arrow_length,0] ,col_arrow, 0.9);
  %arrowHandle =arrow3D([0.5,0,0], [arrow_length,0,0] ,col_arrow, 0.9);
  mArrow3([0 0 -0*arrow_length/2], [0 0 arrow_length/2+tipLength], 'color',col_arrow,'stemWidth',stemWidth,'tipWidth',tipWidth,'tipLength',tipLength)
  mArrow3([0 -0*arrow_length/2 0], [0 arrow_length/2+1.5*tipLength 0], 'color',col_arrow,'stemWidth',stemWidth,'tipWidth',tipWidth,'tipLength',tipLength)
  mArrow3([-0*arrow_length/2 0 0], [arrow_length/2+4*tipLength 0 0], 'color',col_arrow,'stemWidth',stemWidth,'tipWidth',tipWidth,'tipLength',tipLength)



% set light  
  material dull
  %lighting gouraud
  camlight headlight




  drawnow

% save  
  exportgraphics(gcf,'Fig6a_brillouin_zone_sectors.png','ContentType','image','Resolution',1200,'BackgroundColor','white')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%    
% plot Brillouin zone sectors
  figure(62); clf; hold all;

% axes
  view([123 9]);
  axis equal
  axis off

% sweep over all sectors and facets  
  for s = [1,6]

    % set rotation
      R = Rx(angles{s}(1)) * Ry(angles{s}(2)) * Rz(angles{s}(3));

    % set shift
      shift_scale = 0.0; % in units of 2pi/a0
      

    % plot valley vector
      k0_scale = 0.8;
      plot3(k0_scale * ex'*shift{s},k0_scale *  ey'*shift{s},k0_scale *  ez'*shift{s}, 'k.','MarkerSize',25)

    % plot facets  
      for i = 1 : length(facet)      
        pts = shift_scale * shift{s} + R * facet{i};
        p = patch(ex'*pts, ey'*pts, ez'*pts,'k-','LineWidth',1);
        set(p,'FaceColor',color{s},'FaceAlpha',0.5,'EdgeColor',color{s}*0.5,'EdgeAlpha',1.0);
        %set(p,'EdgeColor','none')      
      end
  end


% draw axes
  arrow_length = 2.6;
  %arrowHandle =arrow3D([0,0,-arrow_length], [0,0,2*arrow_length] ,col_arrow, 0.9);
  %arrowHandle =arrow3D([0,-arrow_length,0], [0,2*arrow_length,0] ,col_arrow, 0.9);
  %arrowHandle =arrow3D([-arrow_length,0,0], [2*arrow_length,0,0] ,col_arrow, 0.9);
  mArrow3([0 0 -arrow_length/2], [0 0 arrow_length/2+tipLength], 'color',col_arrow,'stemWidth',stemWidth,'tipWidth',tipWidth,'tipLength',tipLength)
  mArrow3([0 -arrow_length/2 0], [0 arrow_length/2+tipLength 0], 'color',col_arrow,'stemWidth',stemWidth,'tipWidth',tipWidth,'tipLength',tipLength)
  mArrow3([-arrow_length/2 0 0], [arrow_length/2+tipLength 0 0], 'color',col_arrow,'stemWidth',stemWidth,'tipWidth',tipWidth,'tipLength',tipLength)

% set light  
  material dull
  %lighting gouraud
  camlight headlight




  drawnow

% save  
  exportgraphics(gcf,'Fig6b_brillouin_zone_sectors.png','ContentType','image','Resolution',1200,'BackgroundColor','white')

    
    
