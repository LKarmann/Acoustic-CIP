function graph_mesh(x,mesh,Gypsilab)

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