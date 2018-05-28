function [ xl xr ] = cover( keeper, attacker )
% Computes the uncovered parts

xk = keeper(1);
yk = keeper(2);
theta = keeper(3);

xa = attacker(1);
ya = attacker(2);

xkl = xk - d * cos(theta);
ykl = yk - d * sin(theta);

xkr = xk + d * cos(theta);
ykr = yk + d * sin(theta);

xcl = (ya*xkl - xa*ykl)/(ya-ykl);
xcr = (xkr*ya-xa*ykr)/(ya - ykr);

xl = xcl + 7.32/2;
xr = 7.32/2 - xcr;
end
