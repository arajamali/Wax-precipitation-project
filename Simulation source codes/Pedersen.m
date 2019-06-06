clear all
close all
clc

R = 10.732;
%R = 1.99;

Components = char({'C1';'C2';'C3';'iC4';'C4';'iC5';'C5';'C6';'C7';'C8';'C9';'C10';'C11';'C12';'C13';'C14';'C15';'C16';'C17';'C18';'C19';'C20';'C21';'C22';'C23';'C24';'C25';'C26';'C27';'C28';'C29';'C30+'});
zOil1 = [1.139;0.507;0.481;0.563;0.634;1.113;0.515;2.003;5.478;8.756;7.222;5.414;4.323;4.547;5.289;4.720;4.445;3.559;3.642;3.104;2.717;2.597;1.936;2.039;1.661;1.616;1.421;1.233;1.426;1.343;1.300;13.234]/100;
MW = [90.9;105.0;117.7;132.0;148.0;159.0;172.0;185.0;197.0;209.0;227.0;243.0;254.0;262.0;281.0;293.0;307.0;320.0;333.0;346.0;361.0;374.0;381.0;624.0]; % for the 24 components from C7 to C30+
MWC1_C6 = [16.043;30.070;44.097;58.123;58.123;72.15;72.15;86.177];

WOil1 = sum(zOil1.*[MW;MWC1_C6]);% sample weight, one mole
WaxPercentOil1=[233.3,6.3;238.7,4.4;243.1,3.3;248.5,1.9;253.6,2.0;258.2,1.2;263.4,1.1;268.1,1.2;273.5,1.3;278.6,1.0;283.2,0.5;288.7,0.7;293.4,0.6;298.3,0.7;303.5,0.3];

% for the 24 components from C7 to C30 and from table 11.12 of the book
%DH = [9564;12053;15238;repmat(19401,2,1);24357;repmat(29196,2,1);repmat(36420,2,1);repmat(43335,2,1);repmat(49128,2,1);repmat(58487,4,1);repmat(73756,6,1)];
%Tf = [166.9;188.8;211.0;repmat(233,2,1);253.4;repmat(268.3,2,1);repmat(285.5,2,1);repmat(298,2,1);repmat(306.4,2,1);repmat(317.1,4,1);repmat(330.1,6,1)]
T = (325:-5:235)';
Tf = 374.5+0.02617.*MW-20172./MW;

%[a1,a2,a3,a4,a5]=deal(0.5914,5.763,0.5148,0.3033,-4.635e-4);
[a1,a2,a3,a4,a5]=deal(0.1,5.763,0.26,0.36,-10e-4); % I had to slightly change
% some of the fitting parameters to get a better match
deltaLPed = 7.41 + a1.*(log(7:30)-log(7))';
deltaSPed = 8.50 + a2.*(log(7:30)-log(7))';
deltaLWon = [7.41;7.53;7.63;7.71;7.78;7.83;7.88;7.92;7.96;7.99;8.02;8.05;8.07;8.09;8.11;8.13;8.15;8.17;8.18;8.20;8.21;8.22;8.24;8.25;];%8.26;8.27;8.28;8.29;8.30;8.31;8.32;8.33;8.34;8.35];
deltaSWon = [8.5;8.78;9;9.17;9.32;9.44;9.55;9.64;9.72;9.79;9.86;9.92;9.97;10;10.1;10.1;10.1;10.2;10.2;10.3;10.3;10.3;10.3;10.4;];%10.4;10.4;10.4;10.4;10.5;10.5;10.5;10.5;10.5;10.6];
deltaLPed = deltaLWon;
deltaSPed = deltaSWon;
% deltaLPed = [7.41;7.49;7.56;7.62;7.68;7.13;7.78;7.82;7.86;7.90;7.93;7.91;8;8.03;8.06;8.09;8.11;8.14;8.16;8.19;8.21;8.23;8.25;8.27;];%8.29;8.31;8.33;8.34;8.36;8.38;8.39;8.41;8.43;8.44];
% deltaSPed = [8.5;9.27;9.95;10.6;11.1;11.6;12.1;12.5;12.9;13.3;13.6;13.9;14.3;14.6;14.8;15.1;15.4;15.6;15.8;16.1;16.3;16.5;16.7;16.9;];%17.1;17.3;17.4;17.6;17.8;17.9;18.1;18.2;18.4;18.5];
DH = a3.*0.1426.*MW.*Tf; % cal/mole
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
    DCp = a4.*MW+a5.*MW.*T(i);
    while err > 0.1
        f = @(nS)(1-sum(zOil1./(1+nS*(K-1))));
        options = optimoptions('fsolve','Display','off');
        nS = fsolve(f,0.5,options);
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
        %K(9:32) = gammaL./gammaS.*exp(DH./(R.*T(i)).*(1-T(i)./Tf));
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
% calculation of weight percent
scatter(WaxPercentOil1(:,1),WaxPercentOil1(:,2),'+');
hold on
plot(T,WP*100)

