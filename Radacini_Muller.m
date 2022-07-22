function [r] = Radacini_Muller(p,z,itmax)

pc=p; % facem o copie a polinoumlui in pc care o sa ne foloseasca la afisarea graficului
n=length(p)-1; % apelam functia length() din Matlab pentru determinarea numarului de elemente din p si il scadem cu 1 pentru a afla gradul polinomului

r=zeros(1,n); % initializam vectorul r cu 0, iar in el vom retine radacinile polinomului p
i=1; % initializam i cu 1 pentru a gasi prima radacina

while i<n % pana cand gradul polinomului p este diferit de 0, adica pana cand se gasesc toate radacinile
    
    if i~=1 % dupa ce gasim prima radacina
        p=deconv(p,poly(r(i-1))); % reducem polinomul prin impartire la x-r, unde r este ultima radacina calculata
    end
    
    cz=z; % facem o copie a vectorului z care contine aproximarile pentru radacini
    it=itmax; % facem o copie a variabilei itmax pentru a pastra acelasi numar de iteratii la fiecare calcul al radacinilor
    
    p1=polyval(p,cz(1)); % p1 este valoarea polinomului p in punctul in care se gaseste prima aproximare a radacinilor
    p2=polyval(p,cz(2)); % p2 este valoarea polinomului p in punctul in care se gaseste a doua aproximare a radacinilor
    p3=polyval(p,cz(3)); % p3 este valoarea polinomului p in punctul in care se gaseste a treia aproximare a radacinilor
    zr=cz(1); % initializam zr cu valoarea uneia dintre cele 3 aproximari pentru a intra in bucla

    % aplicam algoritmul lui Muller pana cand norma lui p(zr) este foarte aproape de 0, adica zr sa fie o radacina pt polinomul p
    while norm(polyval(p,zr)) > 1e-5

    % se testeaza daca s-a atins numarul maxim de iteratii
        if it==0
            disp('S-a atins numarul maxim de iteratii fara a se gasi o radacina exacta pentru functia f');
            return;
        end

    % se aplica formulele din metoda lui Muller
        q=(cz(3)-cz(2))/(cz(2)-cz(1));
        A=q*p3-q*(1+q)*p2+q^2*p1;
        B=(2*q+1)*p3-(1+q)^2*p2+q^2*p1;
        C=(1+q)*p3;

    % se verifica daca valoarea radicalului este mai mica decat 0, iar in acest caz o inmultim cu (-1)
        rad=sqrt(B^2-4*A*C);
        if rad<0
            rad=-rad;
        end

    % se afla urmatoarea aproximare a unei radacini pentru polinomul p
        zr=cz(3)-(cz(3)-cz(2))*(2*C/(B+rad));

    % se modifica valoarile vectorului z cu ultimele doua valori din vector si valoarea calculata mai sus a lui zr
        cz(1)=cz(2);
        cz(2)=cz(3);
        cz(3)=zr;

        p1=polyval(p,cz(1)); % p1 primeste valoarea polinomului p in noul punct cz(1)
        p2=polyval(p,cz(2)); % p2 primeste valoarea polinomului p in noul punct cz(2)
        p3=polyval(p,cz(3)); % p3 primeste valoarea polinomului p in noul punct cz(3)

        it=it-1; % scadem numarul de iteratii cu 1
    end
    
r(i)=zr; % retinem valoarea radacinii actuale in vectorul r in pozitia i
i=i+1; % incrementam variabila i pentru a trece la urmatoarea radacina
end

p=deconv(p,poly(r(i-1)));
r(i)=-p(2)/p(1);

% cautam limita minima si limita maxima a graficului
min=real(r(1)); % min e partea reala a primei radacini
max=min; % pe max il egal cu min, adica tot partea reala a primei radacini
i=2; % incepem cu al doilea element din vectorul de radacini pentru ca prima radacina a fost folosita la initializarea variabilelor min si max
while i<=n % pana ajungem la sfarsitul vectorului
    if min>r(i)
        min=real(r(i)); % gasim partea reala pentru cea mai mica radacina 
    elseif max<r(i)
        max=real(r(i)); % gasim partea reala pentru cea mai mare radacina
    end
    i=i+1; % incrementam i
end

% alegem o marja suficienta pentru a afisa graficul astfel incat sa se observe radacinile din punctele extreme
min=min-2;
max=max+2;

% se va face afisarea grafica a polinomului si a radacinilor acestuia
hold on
x=min:0.001:max;
plot(x,polyval(pc,x)); grid on; % se afiseaza grafic polinomul p

plot(r,'xr','LineWidth',2,'MarkerSize',10); % se afiseaza grafic radacinile polinomului p
hold off

end