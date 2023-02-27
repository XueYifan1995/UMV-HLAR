%  m,u are vectors
%add control inputs to the GP inputs
function [tileM, tileS] = jointDistribution(m, s, u)
    dimS = size(s, 1);
    dimU = size(u,1);
    tileM = [m; u];
    tileS = zeros(dimS+dimU,dimS+dimU);
    tileS(1:dimS, 1:dimS) = s;
end