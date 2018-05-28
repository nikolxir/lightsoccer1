function [ J ] = costfunction( keeper, attacker )
% Computes cost function

[xl, xr] = cover(keeper, attacker);
J = xl^2 + xr^2;
end