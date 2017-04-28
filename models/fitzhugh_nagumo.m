function fitzhugh_nagumo(a, b, tau, iext)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fnparams = FitzhughNagumoParams;
fnparams.A    = a;
fnparams.B    = b;
fnparams.Tau  = tau;
fnparams.Iext = iext;

plot_vector_field(fnparams)
end

function res = modelfun(x, y, fnparams)
% Comments

  dx = x - (x^3)/3 - y + fnparams.Iext;
  dy = (x + fnparams.A - fnparams.B * y) / fnparams.Tau;
  res = [dx, dy];
end

function plot_vector_field(fnparams)
% Comments

[x, y] =  meshgrid(0:0.2:2,0:0.2:2);
[dx, dy] = meshvec(x, y, @(x0, y0) modelfun(x0, y0, fnparams));

figure
quiver(x, y, dx, dy)
end