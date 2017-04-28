function [u, v] = meshvec(x, y, f)
% Returns [ ..., (f_1(x(i), y(i)), f_2(x(i), y(i)), ... ]

szx = size(x)
u = zeros(szx)
v = zeros(szx)

for i=1:numel(x)
   vec = f(x(i), y(i));
   u(i) = vec(1);
   v(i) = vec(2);
   
u
end
