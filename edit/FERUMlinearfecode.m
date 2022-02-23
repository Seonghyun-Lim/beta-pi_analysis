function [g,dg] = FERUMlinearfecode(lsf,x,grad_flag,gfundata,femodel,randomfield)

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


% Extract data from passed data structures
ndfpernode   = femodel.ndf;
node         = femodel.node;
el           = femodel.el;
loading      = femodel.loading;
nodal_spring = femodel.nodal_spring;
fixed_dof    = femodel.fixed_dof;
id           = femodel.id;
resp         = gfundata(lsf).resp;
lim          = gfundata(lsf).lim;
if length(gfundata(lsf).scenario) == 1
    component = gfundata(lsf).scenario;
end

% Number of random variables 'nrv'
nrv = length(x);

% Number of nodes 'nnode'
ndim = size(node);
nnode = ndim(1);

% Number of elements 'nel'
eldim = size(el);
nel = eldim(1);

% Total number of degrees of freedom
ndof = nnode*ndfpernode;

% Number of loads
load_dim = size(loading);
nload = load_dim(1);

% Check if nodal springs are given, and in case how many
spring_dim = size(nodal_spring);
one_if_no_springs = spring_dim(2);
nspring = spring_dim(1);

% Position random field and compute its Jacobian (if specified)
% or: Use 'id' array to position random variables 'x' in the FE model
if length(randomfield.mesh) ~= 1
   [femodel,rf_jac]=positionrandomfield(x,femodel,randomfield,'no ');
   one_if_no_rf = 0;
else
   one_if_no_rf = 1;
   rf_jac = 1;
   [femodel]=position(x,femodel);
end

% Extract updated data from passed data structures
el           = femodel.el;
loading      = femodel.loading;
nodal_spring = femodel.nodal_spring;

% Establish matrix for introducing boundary conditions 
bound = create_bound(fixed_dof,ndfpernode,nnode);

% Assemble stiffness matrix 
K = zeros( ndof );
for i = 1 : nel
   if ( el(i,1)==1 ) % Elastic truss
      node1 = el(i,2);
      node2 = el(i,3);
      nodelist    = [ node1 node2 ];
      nodalcoords = [ node(node1,1) node(node1,2) node(node2,1) node(node2,2) ];
      paramlist   = [ el(i,4) el(i,5) ];
      conn = create_conn(nodelist,ndfpernode,nnode);
      k = ke_truss(nodalcoords,paramlist);
   elseif ( el(i,1)==2 ) % Elastic beam
      node1 = el(i,2);
      node2 = el(i,3);
      nodelist    = [ node1 node2 ];
      nodalcoords = [ node(node1,1) node(node1,2) node(node2,1) node(node2,2) ];
      paramlist   = [ el(i,4) el(i,5) el(i,6) ];
      conn = create_conn(nodelist,ndfpernode,nnode);
      k = ke_beam(nodalcoords,paramlist);
   elseif ( el(i,1)==3 ) % Elastic quad4
      node1 = el(i,2);
      node2 = el(i,3);
      node3 = el(i,4);
      node4 = el(i,5);
      nodelist    = [ node1 node2 node3 node4 ];
      nodalcoords = [ node(node1,1) node(node1,2) node(node2,1) node(node2,2) ...
            node(node3,1) node(node3,2) node(node4,1) node(node4,2) ];
      paramlist   = [ el(i,6) el(i,7) el(i,8) el(i,9) ];
      conn = create_conn(nodelist,ndfpernode,nnode);
      k = ke_quad4(nodalcoords,paramlist);
   else
      disp('ERROR: Unknown linear elastic element type');
   end
   K = K + conn' * k * conn;
end

% Add nodal spring stiffnesses
if one_if_no_springs ~= 1
   for i = 1 : nspring
      the_dof = nodal_spring(i,1)*ndfpernode - ndfpernode + nodal_spring(i,3);
      K(the_dof,the_dof) = K(the_dof,the_dof) + nodal_spring(i,2); 
   end
end

K = bound' * K * bound;

% Assemble load vector 
F = zeros( ndof , 1 );
for i = 1 : nload
   F( loading(i,1)*ndfpernode-(ndfpernode-abs(loading(i,3))))=sign(loading(i,3))*loading(i,2);
end
F = bound' * F;

% Compute displacement response 
iK = inv(K);
displ = iK * F;

% Evaluate limit-state function
displ_blown_up = bound * displ;

resp_lx = el(component,2)*ndfpernode - ndfpernode + 1;
resp_ly = el(component,2)*ndfpernode - ndfpernode + 2;
resp_rx = el(component,3)*ndfpernode - ndfpernode + 1;
resp_ry = el(component,3)*ndfpernode - ndfpernode + 2;

disp_x = displ_blown_up(resp_rx)-displ_blown_up(resp_lx);
disp_y = displ_blown_up(resp_ry)-displ_blown_up(resp_ly);
C = node(el(component,3),:)-node(el(component,2),:);
Le = norm(C);
T = C/Le;
G = el(component,4)*el(component,5)/Le;
Tj = G*T;
F = sum(Tj.*[disp_x disp_y]);
St = F/el(component,5);


g = lim - abs(St);

% Plot deformed shape
plot_def(displ_blown_up,femodel);


% SENSITIVITY COMPUTATIONS (compute dgdz = dgdu*dudz, where dudz=inv(K)*(dFdz-dKdz*u) )

if grad_flag == 'yes'; % If the user wants gradients to be computed by DDM
   
   % First compute vector dgdu (for now assuming 'displacementlimit' gfun)
   dgdu = zeros(ndof,1);
   dgdu(resp_dof) = -1;  
   dgdu = bound' * dgdu; 
   
   % Using the Adjoint Method we then find the 'lambda' vector
   lambda_vector = iK * dgdu;    
   
   % (# of gradients should be equal to # of r.v.'s from the id-array)
   dim_id = size(id);
   nid = dim_id(1);
   nrv_from_id = max( id(:,1) );
   
   for i = 1 : nrv_from_id   % (For each of these loops, one gradient is computed)
      
      % i : count which gradient is computed now
      % j : count where we are in the id-array   
      
      dK = zeros( ndof );
      dF = zeros( ndof , 1 );
      
      for j = 1 : nid % (Step through the id-array, and pick up contributions)
         
         if id(j,1) == i  % (If r.v. # in this id-array line equals the current gradient #)
            
            parameter = id(j,2);
            
            if parameter==1 % Nodal load
               load_node = loading( (id(j,3)) , 1 );
               load_dir  = abs( loading( (id(j,3)) , 3 ) );
               if load_dir == 1
                  the_dof = (load_node-1)*ndfpernode + 1;
               elseif load_dir == 2j
                  the_dof = (load_node-1)*ndfpernode + 2;
               end
               if sign( loading( (id(j,3)) , 3 ) ) == 1
                  dF(the_dof) = loading( (id(j,3)) , 4 );
               else
                  dF(the_dof) = -loading( (id(j,3)) , 4 );
               end
               
            elseif parameter==5 % Nodal spring
               the_dof = nodal_spring(( id(j,3) ),1)*ndfpernode - ndfpernode + nodal_spring(( id(j,3) ),3);
               dK(the_dof,the_dof) = 1;
               
            else
               elno = id(i,3);
               elementtype = el(elno,1);
               if ( elementtype==1 ) % Elastic truss
                  node1 = el(elno,2);
                  node2 = el(elno,3);
                  nodelist    = [ node1 node2 ];
                  nodalcoords = [ node(node1,1) node(node1,2) node(node2,1) node(node2,2) ];
                  paramlist   = [ el(elno,4) el(elno,5) ];
                  conn = create_conn(nodelist,ndfpernode,nnode);
                  dk = dke_truss(nodalcoords,paramlist,parameter);
               elseif ( elementtype==2 ) % Elastic beam
                  node1 = el(elno,2);
                  node2 = el(elno,3);
                  nodelist    = [ node1 node2 ];
                  nodalcoords = [ node(node1,1) node(node1,2) node(node2,1) node(node2,2) ];
                  paramlist   = [ el(elno,4) el(elno,5) el(elno,6) ];
                  conn = create_conn(nodelist,ndfpernode,nnode);
                  dk = dke_beam(nodalcoords,paramlist,parameter);
               elseif ( elementtype==3 ) % Elastic quad4

                  node1 = el(elno,2);
                  node2 = el(elno,3);
                  node3 = el(elno,4);
                  node4 = el(elno,5);
                  nodelist    = [ node1 node2 node3 node4 ];
                  nodalcoords = [ node(node1,1) node(node1,2) node(node2,1) node(node2,2) ...
                        node(node3,1) node(node3,2) node(node4,1) node(node4,2) ];
                  paramlist   = [ el(elno,6) el(elno,7) el(elno,8) el(elno,9) ];
                  conn = create_conn(nodelist,ndfpernode,nnode);
                  dk = dke_quad4(nodalcoords,paramlist,parameter);
               else
                  disp('ERROR: Unknown element type');
               end
               dK = dK + conn' * dk * conn;
            end
         end
      end
      
      dF = bound' * dF;
      dK = bound' * dK * bound;
      
      % Find gradient vector by multiplying by the 'lambda' vector
      dg(i) = lambda_vector' * ( dF - (dK * displ) );
      
   end
      
   % Mult. by the Jacobian of the random field transformation 
   if one_if_no_rf ~= 1
      dg = dg * rf_jac;
   end
   
else % Assign dummy value if grad_flag==0 (gradients not computed by DDM) 
   dg = 0; 
end
