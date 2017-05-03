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
        function res = modelfun(this, x, y)
            % Encodes the RHS of the Fitzhugh-Nagumo equation. In other
            % words, this function encodes f(x) in the ODE:
            %
            %    dx/dt = f(x)

            dx = x - (x^3)/3 - y + this.Iext;
            dy = (x + this.A - this.B * y) / this.Tau;
            res = [dx, dy];
        end
        
        function res = modelfunT(this, ~, y, T)
            % This function encodes the RHS of the ODE solved by @solve_bvp
            % and is used in the process of finding limit cycle solutions.

            dx = T*(y(1) - (y(1)^3)/3 - y(2) + this.Iext);
            dy = T*((y(1) + this.A - this.B * y(2)) / this.Tau);
            res = [dx, dy];
        end
        
        function res = bdryT(this, ya, yb, T)
            % Companion function to modelfunT. This function encodes the
            % boundary conditions and phase condition for the @solve_bvp.
            % These conditions are needed to derive the period of a 
            % limit-cycle.

            res = [ ya(1) - yb(1)
                    ya(2) - yb(2)
                    T*((ya(1) + this.A - this.B * ya(2)) / this.Tau) ];
        end
        
        function g = guessFn(~, t)
            % Encodes a circle. Used as an initial guess for the 
            % bvp solver called in @solve_bvp. 

            g = [ sin(2*pi*t)
                  cos(2*pi*t) ];
        end

        function res = d_modelfun_mtx(this, x, ~)
            % Returns the matrix Df(x, y).

            res = zeros(2);
            res(1,1) = x^2;
            res(1,2) = 0;
            res(2,1) = 1;
            res(2,2) = -this.B;
        end
        
        function lim_cycle = solve_bvp(this)
            % In order to find limit-cycle/periodic solutions to our ODE
            % we need to switch to the equivalent system
            %
            % 1) dx/ds = T * f(x) 
            %
            % 2) x(0) = x(T)
            % 
            % 3) *some phase condition here*
            %
            % This function solves the BVP above. Note that the limit 
            % cycle solution returned is already rescaled to have period
            % T. 
            
            solinit = bvpinit(linspace(0, 1, 5), @(x,y) this.guessFn(x, y), 2*pi);
            sol = bvp4c(@(t, y, T) this.modelfunT(t, y, T), @(ya, yb, T) this.bdryT(ya, yb, T), solinit);

            lim_cycle = LimitCycle;
            lim_cycle.T = sol.parameters;
            lim_cycle.t = lim_cycle.T * sol.x;
            lim_cycle.x = sol.y(1,:);
            lim_cycle.y = sol.y(2,:);
        end
        
        function floquet_matrix(this)
            % Comments

            lim_cycle         = this.solve_bvp();
            transients        = zeros(2, 2, 1);
            transients(:,:,1) = eye(2);

            for i=1:(length(lim_cycle.t) - 1)
                A = this.d_modelfun_mtx(lim_cycle.x(i), lim_cycle.y(i));
                [solt, solv] = fundamental_solution(A, transients(:,:,i), lim_cycle.t(i), lim_cycle.t(i+1));
                transients(:,:,i+1) = solv(length(solt));
            end

            transients
        end
        
        function plot_vector_field(this)
            % Comments

            [x, y] =  meshgrid(-3:0.2:3,-2:0.2:4);
            [dx, dy] = meshvec(x, y, @(x0, y0) this.modelfun(x0, y0));

            figure
            quiver(x, y, dx, dy)
        end
        
        function plot_limit_cycle(this)
            % Comments

            lim_cycle = this.solve_bvp();

            figure
            plot_vector_field(fnparams)
            hold on
            plot(lim_cycle.x, lim_cycle.y)
            hold off
        end
    end
end

