function y = funFn(X,wh,mesh,submeshb,Proposal)
% Calculates the normal derivative of wh on the boundary using discrete
% derivative. Errors might be observed on vertices of the square.
% See Appendix D.2.
%
% Arguments:
% X (Nx3 'double'): Position for a calculation of an integral.
%                   See documentation of Gypsilab.
% wh ('double'): Laplace transform of the solution of a forward problem.
% mesh ('msh'): Mesh on which wh has been calculated.
%               See documentation of Gypsilab.
% submeshb ('msh'): Mesh of the boundary of the sub-open set.
%                   See documentation of Gypsilab.
% Proposal ('int'): Proposal for a calculation at the vertices of the
%                   square (see Appendix D.2).
%
% Returns:
% y ('double'): Value which is used for the integration.
%               See documentation of Gypsilab.


N = sqrt(mesh.size(1)/2);                    % Parameter of the mesh

fn = zeros([size(submeshb,1) 1]);

dx = mesh.vtx(N+2,1) - mesh.vtx(1,1);
dy = mesh.vtx(2,2) - mesh.vtx(1,2);

for k = 1:size(submeshb,1)
    x = submeshb.vtx(k,1);
    y = submeshb.vtx(k,2);

    [~,l] = min((mesh.vtx(:,1)-x).^2 + (mesh.vtx(:,2)-y).^2);

    if x == 0 && y > 0 && y < 1
        fn(k) = (wh(l-N-1) - wh(l+N+1))/(2*dx);
    elseif x == 1 && y > 0 && y < 1
        fn(k) = (wh(l+N+1) - wh(l-N-1))/(2*dx);
    elseif y == 0 && x > 0 && x < 1
        fn(k) = (wh(l-1) - wh(l+1))/(2*dy);
    elseif y == 1 && x > 0 && x < 1
        fn(k) = (wh(l+1) - wh(l-1))/(2*dy);

    % Proposal 1: Mean value discretisation
    elseif x == 0 && y == 0 && Proposal == 1
        fn(k) = ((wh(l-N-1) - wh(l+N+1))/(2*dx) + (wh(l-1) - wh(l+1))/(2*dy))/sqrt(2);
    elseif x == 0 && y == 1 && Proposal == 1
        fn(k) = ((wh(l-N-1) - wh(l+N+1))/(2*dx) + (wh(l+1) - wh(l-1))/(2*dy))/sqrt(2);
    elseif x == 1 && y == 0 && Proposal == 1
        fn(k) = ((wh(l+N+1) - wh(l-N-1))/(2*dx) + (wh(l-1) - wh(l+1))/(2*dy))/sqrt(2);
    elseif x == 1 && y == 1 && Proposal == 1
        fn(k) = ((wh(l+N+1) - wh(l-N-1))/(2*dx) + (wh(l+1) - wh(l-1))/(2*dy))/sqrt(2);

    % Proposal 2: Diagonal discretisation
    elseif x == 0 && y == 0 && Proposal == 2
        fn(k) = (wh(l-N-2) - wh(l+N+2))/sqrt(dx^2 + dy^2);
    elseif x == 0 && y == 1 && Proposal == 2
        fn(k) = (wh(l-N) - wh(l+N))/sqrt(dx^2 + dy^2);
    elseif x == 1 && y == 0 && Proposal == 2
        fn(k) = (wh(l+N) - wh(l-N))/sqrt(dx^2 + dy^2);
    elseif x == 1 && y == 1 && Proposal == 2
        fn(k) = (wh(l+N+2) - wh(l-N-2))/sqrt(dx^2 + dy^2);
    end

end


v = submeshb.vtx;


% Proposal 3: Use the nearest node
if Proposal == 3
    v(v(:,1) == 0 & v(:,2) == 0, 1:2) = [2,2];
    v(v(:,1) == 0 & v(:,2) == 1, 1:2) = [2,2];
    v(v(:,1) == 1 & v(:,2) == 0, 1:2) = [2,2];
    v(v(:,1) == 1 & v(:,2) == 1, 1:2) = [2,2];
end


y = zeros([size(X,1) 1]);

for k = 1:size(y,1)
    [~, I] = min((v(:,1)-X(k,1)).^2 + (v(:,2)-X(k,2)).^2);
    y(k) = fn(I);
end

end