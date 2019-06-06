function [fval]=obj_func(params,MW,zOil1,WaxPercentOil1,WOil1)
a1=params(1);
a2=params(2);
a3=params(3);
a4=params(4);
a5=params(5);
R = 10.732;
T = WaxPercentOil1(:,1);
Tf = 374.5+0.02617.*MW-20172./MW;
dL25 = 0.8155+0.6273e-4.*MW-13.06./MW;
V = MW./dL25;
nL = 0.5;
nS = 0.5;
NL = T*0;
NS = NL;
WP = NL; % weight percent
K = [zeros(8,1);ones(24,1)];
xL = [zOil1(1:8);zeros(24,1)];
xS = zeros(32,1);
for i=1:numel(T)
    err = inf;
    %DCp = linspace(20,85,24)';
    deltaLPed = 7.41 + a1.*(log(7:30)-log(7))';
    deltaSPed = 8.50 + a2.*(log(7:30)-log(7))';
    DH = a3.*0.1426.*MW.*Tf; % cal/mole
    DCp = a4.*MW+a5.*MW.*T(i);
    while err > 0.1
        f = @(nS)(1-sum(zOil1./(1+nS*(K-1))));
        options = optimoptions('fsolve','Display','off');
        nS = fsolve(f,0.5,options)
        nL = 1 - nS;
        xL = zOil1./(1+K.*nS);
        xS = zOil1./(1+1./K.*nL);
        PhiL = xL(9:32).*V./sum(xL(9:32).*V);
        PhiS = xS(9:32).*V./sum(xS(9:32).*V);
        deltaAvgL = sum(PhiL.*deltaLPed);
        deltaAvgS = sum(PhiS.*deltaSPed);
        gammaL = exp(V.*(deltaAvgL-deltaLPed).^2./(R.*T(i)));
        gammaS = exp(V.*(deltaAvgS-deltaSPed).^2./(R.*T(i)));
        KC = K(9:32);
        %K(9:32) = gammaL./gammaS.*exp(DH./(R.*T(i)).*(1-T(i)./Tf)-a4.*MW./R.*(Tf./T(i)-1-log(Tf./T(i))));
        K(9:32) = gammaL./gammaS.*exp(DH./(R.*T(i)).*(1-T(i)./Tf)+DCp./R.*(1-Tf./T(i)+log(Tf./T(i))));
        %K(9:32) = gammaL./gammaS.*exp(DH./(R.*T(i)).*(1-T(i)./Tf)-a4.*MW./R.*(Tf./T(i)-1-log(Tf./T(i)))-a5.*MW./(2.*R).*(Tf.^2./T(i)+T(i)-2.*Tf.^2));
        err = sum((KC-K(9:32)).^2./(KC.*K(9:32)));
%         [DH./(R.*T(i)).*(1-T(i)./Tf), DH./(R.*T(i)).*(1-T(i)./Tf)-a4.*MW./R.*(Tf./T(i)-1-log(Tf./T(i)))-a5.*MW./(2.*R).*(Tf.^2./T(i)+T(i)-2.*Tf.^2)]
%         pause
    end
    WP(i) = sum(xS(9:end).*nS.*MW)./WOil1;
    NS(i) = nS;
    NL(i) = nL;
end
fval=sum((WP-WaxPercentOil1(:,2)).^2);
end