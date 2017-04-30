function fitzhugh_nagumo(task, fnparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if task == TaskManager.VECTOR_FIELD
    plot_vector_field(fnparams)
elseif task == TaskManager.LIMIT_CYCLE
    plot_limit_cycle(fnparams)
elseif task == TaskManager.FUN_SOL
    floquet_matrix(fnparams)
end
end

function res = modelfun(x, y, fnparams)
% Comments

dx = x - (x^3)/3 - y + fnparams.Iext;
dy = (x + fnparams.A - fnparams.B * y) / fnparams.Tau;
res = [dx, dy];
end

function res = d_modelfun_mtx(x, y, fnparams)
% Comments

res = zeros(2);
res(1,1) = x^2;
res(1,2) = 0;
res(2,1) = 1;
res(2,2) = -fnparams.B;
end

function plot_vector_field(fnparams)
% Comments

[x, y] =  meshgrid(-3:0.2:3,-2:0.2:4);
[dx, dy] = meshvec(x, y, @(x0, y0) modelfun(x0, y0, fnparams));

figure
quiver(x, y, dx, dy)
end

function res = modelfunT(t, y, T, fnparams)
% Comments

dx = T*(y(1) - (y(1)^3)/3 - y(2) + fnparams.Iext);
dy = T*((y(1) + fnparams.A - fnparams.B * y(2)) / fnparams.Tau);
res = [dx, dy];
end

function res = bdryT(ya, yb, T, fnparams)
% Comments

res = [ ya(1) - yb(1)
        ya(2) - yb(2)
        T*((ya(1) + fnparams.A - fnparams.B * ya(2)) / fnparams.Tau) ];
end

function g = guessFn(t)
% Comments

g = [ sin(2*pi*t)
      cos(2*pi*t) ];
end

function lim_cycle = solve_bvp(fnparams)
% Comments

solinit = bvpinit(linspace(0, 1, 5), @guessFn, 2*pi);
sol = bvp4c(@(x0, y0, T0) modelfunT(x0, y0, T0, fnparams), @(a0, b0, T0) bdryT(a0, b0, T0, fnparams), solinit);

lim_cycle = LimitCycle;
lim_cycle.T = sol.parameters;
lim_cycle.t = lim_cycle.T * sol.x;
lim_cycle.x = sol.y(1,:);
lim_cycle.y = sol.y(2,:);
end

function plot_limit_cycle(fnparams)
% Comments

lim_cycle = solve_bvp(fnparams);

figure
plot_vector_field(fnparams)
hold on
plot(lim_cycle.x, lim_cycle.y)
hold off
end

function floquet_matrix(fnparams)
% Comments

lim_cycle = solve_bvp(fnparams);
transients = zeros(2, 2, 1);
transients(:,:,1) = eye(2);

for i=1:(length(lim_cycle.t) - 1)
    A   = d_modelfun_mtx(lim_cycle.x(i), lim_cycle.y(i), fnparams);
    [solt, solv] = fundamental_solution(A, transients(:,:,i), lim_cycle.t(i), lim_cycle.t(i+1));
    transients(:,:,i+1) = solv(length(solt));
end

transients
end

function dydx = odefun(t, X, A)
% Comments

dydx = reshape(A*reshape(X, [2, 2]), [4,1]);
end

function [t,v]  = fundamental_solution(A, init, t0, t1)
% Comments

[t, v] = ode45(@(t, X) odefun(t, X, A), [t0, t1], init);
end