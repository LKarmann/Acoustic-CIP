function [u, u_dot, u_ddot] = damped_newmark(M,C,K,B,u0,u1,dt,Nt,delta,theta,Progression)
% Solves the Newmark numerical integration. Returns bad results if C=0.
% See undamped_newmark for better calculation when C=0.
%
% Arguments:
% M (NxN 'double'): Mass matrix of the Newmark method.
% C ('function_handle'): Damping matrix of the Newmark method. Can be
%                        time-dependent.
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
% Progression ('logical'): 1 means that the progression step is displayed.
%
% Returns:
% u (Nx(Nt+1) 'double'): Vectors of the sequence solution of the damped
%                        Newmark numerical integration. Each column
%                        corresponds at the same time value k*dt for
%                        k = 0:(Nt+1)*dt.
% u_dot (Nx(Nt+1) 'double'): Vectors of the sequence solution of the damped
%                            Newmark numerical integration.
% u_ddot (Nx(Nt+1) 'double'): Vectors of the sequence solution of the damped
%                             Newmark numerical integration.


u = zeros([size(u0,1), Nt+1]);
u_dot = zeros([size(u0,1), Nt+1]);
u_ddot = zeros([size(u0,1), Nt+1]);



u(:,1) = u0;
u_dot(:,1) = u1;
u_ddot(:,1) = M\(B(0)-K*u0-C(0)*u1);

if theta == 0

    for k = 2:Nt+1

        u(:,k) = u(:,k-1) + dt*u_dot(:,k-1) + dt^2/2*u_ddot(:,k-1);
    
        u_ddot(:,k) = (M+delta*dt*C((k-1)*dt))\(B((k-1)*dt)-K*u(:,k)-C((k-1)*dt)*u_dot(:,k-1)-(1-delta)*dt*C((k-1)*dt)*u_ddot(:,k-1));
        u_dot(:,k) = u_dot(:,k-1) + dt*(delta*u_ddot(:,k)+(1-delta)*u_ddot(:,k-1));

        if Progression && mod(k,10)==0
            disp(k)
        end


    end


else
    
    
    for k = 2:Nt+1
    
        b = theta*dt^2*B((k-1)*dt) + (M+delta*dt*C((k-1)*dt))*u(:,k-1)...
            + dt*(M+(delta-theta)*dt*C((k-1)*dt))*u_dot(:,k-1)...
            + dt^2/2*((1-2*theta)*M+(delta-2*theta)*dt*C((k-1)*dt))*u_ddot(:,k-1);


        A = M + delta*dt*C((k-1)*dt) + theta*dt^2*K;
    
    
        u(:,k) = A\b;
    
        u_dot(:,k) = (1-delta/theta)*u_dot(:,k-1) + delta/theta*(u(:,k)-u(:,k-1))/dt...
                    + dt/2*(2-delta/theta)*u_ddot(:,k-1);
        u_ddot(:,k) = (u(:,k)-u(:,k-1))/(theta*dt^2) - u_dot(:,k-1)/(theta*dt) ...
                    - (1-2*theta)/(2*theta)*u_ddot(:,k-1);
    

        if Progression && mod(k,10)==0
            disp(k)
        end
    
    end
end
end