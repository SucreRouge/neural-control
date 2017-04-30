function fitzhugh_nagumo(task, fnparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if task == TaskManager.VECTOR_FIELD
    plot_vector_field(fnparams)
elseif task == TaskManager.LIMIT_CYCLE
    plot_limit_cycle(fnparams)
end
end

function res = modelfun(x, y, fnparams)
% Comments

dx = x - (x^3)/3 - y + fnparams.Iext;
dy = (x + fnparams.A - fnparams.B * y) / fnparams.Tau;
res = [dx, dy];
end

function plot_vector_field(fnparams)
% Comments

[x, y] =  meshgrid(-5:0.2:5,-5:0.2:5);
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

function sol = solve_bvp(fnparams)
% Comments

solinit = bvpinit(linspace(0, 1, 5), @guessFn, 2*pi);
sol = bvp4c(@(x0, y0, T0) modelfunT(x0, y0, T0, fnparams), @(a0, b0, T0) bdryT(a0, b0, T0, fnparams), solinit);
end

function plot_limit_cycle(fnparams)
% Comments

sol = solve_bvp(fnparams);
sol.y(1,:)

figure
plot_vector_field(fnparams)
hold on
plot(sol.y(1,:), sol.y(2,:))
hold off
end