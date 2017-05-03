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
            %
            %    dx/dt = f(x)

            dx = x - (x^3)/3 - y + this.Iext;
            dy = (x + this.A - this.B * y) / this.Tau;
            res = [dx, dy];
        end
        
        function res = modelfunT(t, y, T)
            % This function encodes the RHS of the ODE solved by @solve_bvp
            % and is used in the process of finding limit cycle solutions.

            dx = T*(y(1) - (y(1)^3)/3 - y(2) + this.Iext);
            dy = T*((y(1) + this.A - this.B * y(2)) / this.Tau);
            res = [dx, dy];
        end
        
        function res = bdryT(ya, yb, T, fnparams)
            % Companion function to modelfunT. This function encodes the
            % boundary conditions and phase condition for the @solve_bvp.
            % These conditions are needed to derive the period of a 
            % limit-cycle.

            res = [ ya(1) - yb(1)
                    ya(2) - yb(2)
                    T*((ya(1) + fnparams.A - fnparams.B * ya(2)) / fnparams.Tau) ];
        end
        
        function g = guessFn(t)
            % Encodes a circle. Used as an initial guess for the 
            % bvp solver called in @solve_bvp. 

            g = [ sin(2*pi*t)
                  cos(2*pi*t) ];
        end

        function res = d_modelfun_mtx(x, y)
            % Returns the matrix Df(x, y).

            res = zeros(2);
            res(1,1) = x^2;
            res(1,2) = 0;
            res(2,1) = 1;
            res(2,2) = -this.B;
        end
        
        function lim_cycle = solve_bvp(fnparams)
            % In order to find limit-cycle/periodic solutions to our ODE
            % we need to switch to the equivalent system
            %
            % 1) dx/ds = T * f(x) 
            %
            % 2) x(0) = x(T)
            % 
            % 3) *some phase condition here*
            %
            % This function solves the BVP above. 
            
            solinit = bvpinit(linspace(0, 1, 5), @guessFn, 2*pi);
            sol = bvp4c(@(x0, y0, T0) modelfunT(x0, y0, T0, fnparams), @(a0, b0, T0) bdryT(a0, b0, T0, fnparams), solinit);

            lim_cycle = LimitCycle;
            lim_cycle.T = sol.parameters;
            lim_cycle.t = lim_cycle.T * sol.x;
            lim_cycle.x = sol.y(1,:);
            lim_cycle.y = sol.y(2,:);
        end
    end
end

