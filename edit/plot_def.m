function plot_def(displ_blown_up, femodel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite Element Reliability using Matlab, FERUM, Version 3.0       %
%                                                                   %
% This program is free software; you can redistribute it and/or     %
% modify it under the terms of the GNU General Public License       %
% as published by the Free Software Foundation; either version 2    %
% of the License, or (at your option) any later version.            %
%                                                                   %
% This program is distributed in the hope that it will be useful,   %
% but WITHOUT ANY WARRANTY; without even the implied warranty of    %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     %
% GNU General Public License for more details.                      %
%                                                                   %
% A copy of the GNU General Public License is found in the file     %
% <gpl.txt> following this collection of program files.             %
%                                                                   %
% Developed under the sponsorship of the Pacific                    %
% Earthquake Engineering (PEER) Center.                             %
%                                                                   %
% For more information, visit: http://www.ce.berkeley.edu/~haukaas  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Extract finite element model data
node       = femodel.node;
el         = femodel.el;
ndfpernode = femodel.ndf;
nnode      = max(size(node));
elsize     = size(el);
nel        = elsize(1);

% FIND OUTER PLOT LIMITS
xmin = 0;
xmax = 0;
ymin = 0;
ymax = 0;

for i = 1 : nnode
   if node(i,1) < xmin
      xmin = node(i,1);
   elseif node(i,1) > xmax
      xmax = node(i,1);
   elseif node(i,2) < ymin
      ymin = node(i,2);
   elseif node(i,2) > ymax
      ymax = node(i,2);
   end
end

xmin = xmin - (xmax-xmin) * 0.2; 
ymin = ymin - (ymax-ymin) * 0.2; 
xmax = xmax + (xmax-xmin) * 0.2; 
ymax = ymax + (ymax-ymin) * 0.2; 


% FIND A SCALING FACTOR FOR THE DISPLACEMENT
x_disp_max = 0;
y_disp_max = 0;
for i = 1 : nnode
   dofx = (i-1)*ndfpernode+1;
   dofy = (i-1)*ndfpernode+2;
   if abs( displ_blown_up(dofx) ) > x_disp_max
      x_disp_max = abs( displ_blown_up(dofx) );
   end
   if abs( displ_blown_up(dofy) ) > y_disp_max
      y_disp_max = abs (displ_blown_up(dofy) );
   end
end

charact_dim=sqrt((xmax-xmin)^2+(ymax-ymin)^2); 
if x_disp_max > y_disp_max
   scale_factor = (0.05 * charact_dim) / x_disp_max;
else
   scale_factor = (0.05 * charact_dim) / y_disp_max;
end

% ADD MAX. DISPLACEMENT AROUND PLOT
xmin = xmin - x_disp_max*scale_factor;
ymin = ymin - y_disp_max*scale_factor;
xmax = xmax + x_disp_max*scale_factor;
ymax = ymax + y_disp_max*scale_factor;


% CHECK IF STRUCTURE IS 1-D
if (ymax-ymin) == 0
   ymin_new = ymin - y_disp_max*scale_factor*1.2;
   ymax = ymin + y_disp_max*scale_factor*1.2;
   ymin = ymin_new;
elseif (xmax - xmin) == 0
   xmin_new = xmin - x_disp_max*scale_factor*1.2;
   xmax = xmin + x_disp_max*scale_factor*1.2;
   xmin = xmin_new;
end


figure(1)
clf;
% axis([xmin,xmax,ymin,ymax]) 
title('Deformed Shape of Structure')

% COORDINATES OF POINTS THAT ARE TO BE PLOTTED
for i = 1 : nel
   if ( el(i,1)==1 | el(i,1)==2 )
      node1 = el(i,2);
      node2 = el(i,3);
      xcoord(1) = node(node1,1) + displ_blown_up((node1-1)*ndfpernode+1) * scale_factor;
      ycoord(1) = node(node1,2) + displ_blown_up((node1-1)*ndfpernode+2) * scale_factor;
      xcoord(2) = node(node2,1) + displ_blown_up((node2-1)*ndfpernode+1) * scale_factor;
      ycoord(2) = node(node2,2) + displ_blown_up((node2-1)*ndfpernode+2) * scale_factor;
      hold on;
      plot(xcoord,ycoord,'b-square');
      
   elseif ( el(i,1)==3 )
      node1 = el(i,2);
      node2 = el(i,3);
      node3 = el(i,4);
      node4 = el(i,5);
      xcoord(1) = node(node1,1) + displ_blown_up((node1-1)*ndfpernode+1) * scale_factor;
      ycoord(1) = node(node1,2) + displ_blown_up((node1-1)*ndfpernode+2) * scale_factor;
      xcoord(2) = node(node2,1) + displ_blown_up((node2-1)*ndfpernode+1) * scale_factor;
      ycoord(2) = node(node2,2) + displ_blown_up((node2-1)*ndfpernode+2) * scale_factor;
      hold on;
      plot(xcoord,ycoord,'b-square');
      xcoord(1) = node(node3,1) + displ_blown_up((node3-1)*ndfpernode+1) * scale_factor;
      ycoord(1) = node(node3,2) + displ_blown_up((node3-1)*ndfpernode+2) * scale_factor;
      xcoord(2) = node(node2,1) + displ_blown_up((node2-1)*ndfpernode+1) * scale_factor;
      ycoord(2) = node(node2,2) + displ_blown_up((node2-1)*ndfpernode+2) * scale_factor;
      hold on;
      plot(xcoord,ycoord,'b-square');
      xcoord(1) = node(node3,1) + displ_blown_up((node3-1)*ndfpernode+1) * scale_factor;
      ycoord(1) = node(node3,2) + displ_blown_up((node3-1)*ndfpernode+2) * scale_factor;
      xcoord(2) = node(node4,1) + displ_blown_up((node4-1)*ndfpernode+1) * scale_factor;
      ycoord(2) = node(node4,2) + displ_blown_up((node4-1)*ndfpernode+2) * scale_factor;
      hold on;
      plot(xcoord,ycoord,'b-square');
      xcoord(1) = node(node1,1) + displ_blown_up((node1-1)*ndfpernode+1) * scale_factor;
      ycoord(1) = node(node1,2) + displ_blown_up((node1-1)*ndfpernode+2) * scale_factor;
      xcoord(2) = node(node4,1) + displ_blown_up((node4-1)*ndfpernode+1) * scale_factor;
      ycoord(2) = node(node4,2) + displ_blown_up((node4-1)*ndfpernode+2) * scale_factor;
      hold on;
      plot(xcoord,ycoord,'b-square');
      
   end
   
   if i == nel
      hold off;
   end
   
end

