function w = discrete_laplace_transform(u,t,s,Quadrature)
% Calculates an approximation of the Laplace transform of u known only on
% discrete values of t. Various quadratures are available.
%
% Arguments:
% u (NxM 'double'): Array of the discrete values of u(x,t) for various t.
%                   Each column correspond to a time value and each row to
%                   a x value.
% t ('double'): Array of the corresponding time values. Must have a size
%               Mx1 or 1xM.
% s ('scalar'): Pseudo-frequency at which the Laplace transform is
%               calculated.
% Quadrature ('string'): Quadrature method used to calculate the integral.
%                        Available values are "LeftRectangle",
%                        "RightRectangle" and "Trapezoidal".
%
% Returns:
% w (Nx1 'double'): Laplace transform of u for each x value.


w = zeros([size(u,1), 1]);

for k = 1:length(t)-1

    if Quadrature == "LeftRectangle"

        a = u(:,k) * exp(-s*t(k));

    elseif Quadrature == "RightRectangle"

        a = u(:,k+1) * exp(-s*t(k+1));

    elseif Quadrature == "Trapezoidal"

        a = (u(:,k) * exp(-s*t(k)) + u(:,k+1) * exp(-s*t(k+1)))/2;

    end

    w = w + a*(t(k+1)-t(k));

end
end