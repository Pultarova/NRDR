% metoda siti pro rovnici -a Laplace u + pu = f na (0,L1)x(0,L2)
% okrajova podminka u=g na hranici
function abc
L1 = 3; % delka oblasti ve smeru x1
L2 = 2; % delka oblasti ve smeru x2
N1 = 18; % pocet intervalu ve smeru x1
N2 = 17; % pocet intervalu ve smeru x2
h1 = L1/N1;
h2 = L2/N2;
NN = (N1+1)*(N2+1); % celkovy pocet uzlu
A = zeros(NN); % matice soustavy
B = zeros(NN,1); % prava strana
uzel_hra = zeros(NN,1); % 1 pro hranicni uzel, jinak 0
upom = zeros(NN,1); % ulozim sem zname okrajove hodnoty
% jdu po uzlech:
kde = 0; % poradove cislo aktualniho uzlu
for k2 = 1:N2+1
    for k1 = 1:N1+1
        kde = kde+1; % poradove cislo aktualniho uzlu
        if (k1==1 || k1==N1+1 || k2==1 || k2==N2+1) 
            uzel_hra(kde) = 1; % je to hranicni uzel
            upom(kde) = Fceg((k1-1)*h1,(k2-1)*h2);
        else
            x1 = (k1-1)*h1; % souradnice aktualniho uzlu
            x2 = (k2-1)*h2; % souradnice aktualniho uzlu
            fa = Fcea(x1,x2); 
            A(kde,kde) = 2*fa*(1/h1^2+1/h2^2)+Fcep(x1,x2);
            A(kde,[kde-1,kde+1]) = -fa/h1^2;
            A(kde,[kde-N1-1,kde+N1+1]) = -fa/h2^2;
            B(kde) = Fcef(x1,x2);
        end;
    end;
end;
for k = NN:-1:1 % vymazani radku a sloupcu pro hranicni uzly a uprava B
    if (uzel_hra(k)==1)
        B = B - A(:,k)*upom(k);
        A(k,:) = []; A(:,k) = []; B(k) = [];
    end;
end;
% spy(A)
u = A\B; % vyreseni soustavy lin. rovnic
u2 = reshape(u,N1-1,N2-1); % preusporadani do tvaru oblasti - vnitrni uzly
u3 = reshape(upom,N1+1,N2+1); % preusporadani - vsechny uzly
u3(2:N1,2:N2) = u2; % doplneni hodnot vnitrnich uzlu

x1 = linspace(0,L1,N1+1);
x2 = linspace(0,L2,N2+1);
[xx1,xx2] = meshgrid(x1,x2);
surf(xx1,xx2,u3'); % vykresleni reseni na skutecne siti bodu

function z = Fcea(x1,x2) % funkce a
z = 2+x2;
function z = Fcep(x1,x2) % funkce p
z = 1;
function z = Fcef(x1,x2) % prava strana f
z = 1+3*x2-2*x1;
function z = Fceg(x1,x2) % okrajova podminka
z = (1+x1)/50;
