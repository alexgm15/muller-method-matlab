clear; clc;

% p=[1 -5 -1 2 12 -1 -4]; % se da polinomul
% p=[1 2 3 4 5 6 7 8];
%p=randn(1,7);
p=@(x) x^5-x^4-x^3-x^2-x-1;
p=[1 -1 -1 -1 -1 -1];

z=[-5, 0.1, 10]; % se dau primele 3 aproximari pentru radacinile lui p
itmax=100; % setam numarul maxim de iteratii

[r]=Radacini_Muller(p,z,itmax); % se apeleaza functia de calcul a radacinilor polinomului p si afisarea grafica a acestora
disp('Radacinile polinomului sunt: '); disp(r); % se afiseaza vectorul care contine radacinile polinomului

disp('Testare prin functia roots() din Matlab: '); disp(roots(p)'); % testam corectitudinea cu ajutorul functiei de gasire a radacinilor din Matlab