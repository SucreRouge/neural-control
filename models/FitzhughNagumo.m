classdef FitzhughNagumo
    % A class representing a Fitzhugh-Nagumo type dynamical system.
    %
    % We write the Fitzhugh-Nagumo equations as
    % 
    %     dx/dt = f(x) ; i.e.
    %
    %     dv/dt = v - v^3/3 - w + I_ext
    %     dw/dt = (v + a - bw)/tau
    
    properties
        A
        B
        Tau
        Iext
    end
    
    methods
        function res = modelfun(x, y)
            % Encodes the RHS of the Fitzhugh-Nagumo equation. In other
            % words, this function encodes f(x) in the ODE:
            %    dx/dt = f(x)

            dx = x - (x^3)/3 - y + this.Iext;
            dy = (x + this.A - this.B * y) / this.Tau;
            res = [dx, dy];
        end

        function res = d_modelfun_mtx(x, y)
            % Returns the matrix Df(x, y).

            res = zeros(2);
            res(1,1) = x^2;
            res(1,2) = 0;
            res(2,1) = 1;
            res(2,2) = -fnparams.B;
        end
    end
end

