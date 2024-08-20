function [u, u_dot, u_ddot] = damped_newmark(M,C1,C2,K,B,u0,u1,dt,Nt,t_1,delta,theta)
% Solves the Newmark numerical integration assuming that C =! 0 and C takes
% two values, C1 for t <= t_1 and C1 + C2 for t > t_1. 
% Returns bad results if C1 = 0. Used for the resolution of Problems (2) 
% and (3).
% See Appendix C.
%
% Arguments:
% M (NxN 'double'): Mass matrix of the Newmark method.
% C1 (NxN 'double'): Damping matrix of the Newmark method for t <= t_1.
% C2 (NxN 'double'): Additional damping matrix of the Newmark method for 
%                    t > t_1.
% K (NxN 'double'): Stiffness matrix of the Newmark method.
% B ('function_handle'): Right-side term of the Newmark method. Must be a
%                        function of time with vector values of size Nx1.
% u0 (Nx1 'double'): Initial condition of the ODE.
% u1 (Nx1 'double'): Initial derivative condition of the ODE.
% dt ('scalar'): Time step for the time discretisation.
% Nt ('integer'): Number of time steps. Corresponds to Nt+1 time steps or
%                 T = Nt x dt.
% t_1 ('scalar'): Time instant for changing the damping matrix.
% delta ('scalar'): Parameter for the Newmark numerical integration.
%                   Must be between 1/2 and 1.
% theta ('scalar'): Parameter for the Newmark numerical integration.
%                   Must be between 0 and 1/2.
%
% Returns:
% u (Nx(Nt+1) 'double'): Vectors of the sequence u, solution of the damped
%                        Newmark numerical integration. Each column
%                        corresponds at the same time value k*dt for
%                        0 <= k <= Nt+1.
% u_dot (Nx(Nt+1) 'double'): Vectors of the sequence u_dot, solution of the
%                            damped Newmark numerical integration.
% u_ddot (Nx(Nt+1) 'double'): Vectors of the sequence u_double_dot solution
%                             of the damped Newmark numerical integration.


u = zeros([size(u0,1), Nt+1]);
u_dot = zeros([size(u0,1), Nt+1]);
u_ddot = zeros([size(u0,1), Nt+1]);


u(:,1) = u0;
u_dot(:,1) = u1;
u_ddot(:,1) = M\(B(0)-K*u0-C1*u1);

if theta == 0

    for k = 2:Nt+1

        u(:,k) = u(:,k-1) + dt*u_dot(:,k-1) + dt^2/2*u_ddot(:,k-1);

        if (k-1)*dt <= t_1
            u_ddot(:,k) = (M+delta*dt*C1)\(B((k-1)*dt)-K*u(:,k)-C1*u_dot(:,k-1)-(1-delta)*dt*C1*u_ddot(:,k-1));
        else
            u_ddot(:,k) = (M+delta*dt*(C1+C2))\(B((k-1)*dt)-K*u(:,k)-(C1+C2)*u_dot(:,k-1)-(1-delta)*dt*(C1+C2)*u_ddot(:,k-1));
        end

        u_dot(:,k) = u_dot(:,k-1) + dt*(delta*u_ddot(:,k)+(1-delta)*u_ddot(:,k-1));

    end

else

    A1 = M + delta*dt*C1 + theta*dt^2*K;
    A2 = M + delta*dt*(C1+C2) + theta*dt^2*K;

    dA1 = decomposition(A1);
    dA2 = decomposition(A2);
     
    for k = 2:Nt+1

        if (k-1)*dt <= t_1
    
            b = theta*dt^2*B((k-1)*dt) + (M+delta*dt*C1)*u(:,k-1)...
                + dt*(M+(delta-theta)*dt*C1)*u_dot(:,k-1)...
                + dt^2/2*((1-2*theta)*M+(delta-2*theta)*dt*C1)*u_ddot(:,k-1);
        
        
            u(:,k) = dA1\b;
        
            u_dot(:,k) = (1-delta/theta)*u_dot(:,k-1) + delta/theta*(u(:,k)-u(:,k-1))/dt...
                        + dt/2*(2-delta/theta)*u_ddot(:,k-1);
            u_ddot(:,k) = (u(:,k)-u(:,k-1))/(theta*dt^2) - u_dot(:,k-1)/(theta*dt) ...
                        - (1-2*theta)/(2*theta)*u_ddot(:,k-1);

        else

            b = theta*dt^2*B((k-1)*dt) + (M+delta*dt*(C1+C2))*u(:,k-1)...
                + dt*(M+(delta-theta)*dt*(C1+C2))*u_dot(:,k-1)...
                + dt^2/2*((1-2*theta)*M+(delta-2*theta)*dt*(C1+C2))*u_ddot(:,k-1);
        
        
            u(:,k) = dA2\b;
        
            u_dot(:,k) = (1-delta/theta)*u_dot(:,k-1) + delta/theta*(u(:,k)-u(:,k-1))/dt...
                        + dt/2*(2-delta/theta)*u_ddot(:,k-1);
            u_ddot(:,k) = (u(:,k)-u(:,k-1))/(theta*dt^2) - u_dot(:,k-1)/(theta*dt) ...
                        - (1-2*theta)/(2*theta)*u_ddot(:,k-1);

        end
    
    end
end
end