% Based on the python function matplotlib.mlab.slopes(x, y)

%% Inputs
% x
% y

%% Outputs
% yp = y'(x) (an estimate)

%%
function [yp] = slope(x, y)

yp=zeros(length(y), 1);

dx= x(2:end) - x(1:(end-1));
dy = y(2:end) - y(1:(end-1));

dydx = dy./dx;

yp(2:(end-1)) = (dydx(1:(end-1)) .* dx(2:end) + dydx(2:end) .* dx(1:(end-1))) ./ (dx(2:end) + dx(1:(end-1)));

yp(1) = 2 * dy(1)/dx(1) - yp(2);

%2*(last y difference)/(last x difference) - the slope of the previous 
yp(end) = 2 * dy(end)/dx(end) - yp(end-1);

end