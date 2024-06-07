function u = undamped_newmark(M,K,B,u0,u1,dt,Nt,delta,theta)
% Solves the Newmark numerical integration when C = 0.
%
% Arguments:
% M (NxN 'double'): Mass matrix of the Newmark method.
% K (NxN 'double'): Stiffness matrix of the Newmark method.
% B ('function_handle'): Right-side term of the Newmark method. Must be a
%                        function of time with vector values of size Nx1.
% u0 (Nx1 'double'): Initial condition of the ODE.
% u1 (Nx1 'double'): Initial derivative condition of the ODE.
% dt ('scalar'): Time step.
% Nt ('integer'): Number of time steps.
% delta ('scalar'): Parameter of the Newmark method.
%                   Should be between 1/2 and 1.
% theta ('scalar'): Parameter of the Newmark method.
%                   Should be between 0 and 1/2.
%
% Returns:
% u (Nx(Nt+1) 'double'): Vectors of the sequence solution of the undamped
%                        Newmark numerical integration. Each column
%                        corresponds at the same time value k*dt for
%                        k = 0:(Nt+1)*dt.


u = zeros([size(u0,1), Nt+1]);


u(:,1) = u0;
u(:,2) = u0 + dt*u1;


A = M + theta*dt^2*K;

dA = decomposition(A);



for k = 3:Nt+1

    b = dt^2*(theta*B((k-1)*dt) + (1/2 + delta - 2*theta)*B((k-2)*dt) + (1/2 - delta + theta)*B((k-3)*dt)) ...
        + M*(2*u(:,k-1) - u(:,k-2)) + dt^2*K*(-(1/2 + delta - 2*theta)*u(:,k-1) - (1/2 - delta + theta)*u(:,k-2));


    u(:,k) = dA\b;

end
end