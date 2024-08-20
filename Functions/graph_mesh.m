function graph_mesh(x,mesh,Gypsilab)
% Plots the values of x defined over mesh and displays the edges of the 
% mesh.
% See documentation of Gypsilab.
%
% Arguments:
% x ('double'): Values to be plot over mesh.
% mesh ('msh'): Mesh used to define x.
%               See documentation of Gypsilab.
% Gypsilab ('logical'): 1 means that x is displayed using Gypsilab function
%                       graph. 0 means that x is displayed using
%                       triangulation.
%
% Returns:
% Plot the results on a figure.

T = triangulation(mesh.elt,mesh.vtx(:,1),mesh.vtx(:,2),x);

if Gypsilab
    Vh = fem(mesh, 'P1');
    
    graph(Vh, x);
    hold on
    trisurf(T);
    hold off
else
    trisurf(T);
end
end