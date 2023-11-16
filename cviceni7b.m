% metoda siti pro rovnici -au''+pu=f a OP u=g
function abc
L1 = 3; % delka oblasti ve smeru x1
L2 = 2; % delka oblasti ve smeru x2
N1 = 8; % pocet intervalu ve smeru x1
N2 = 7; % pocet intervalu ve smeru x2
h1 = L1/N1;
h2 = L2/N2;
NN = (N1+1)*(N2+1); % celkovy pocet uzlu
A = zeros(NN); % matice soustavy
B = zeros(NN,1); % prava strana
uzel_hra = zeros(NN,1); % 1 pro hranicni uzel, jinak 0
% jdu po uzlech:
kde = 0; % poradove cislo uzlu
for k2 = 1:N2+1
    for k1 = 1:N1+1
        kde = kde+1; % poradove cislo aktualniho uzlu
        if (k1==1 || k1==N1+1 || k2==1 || k2==N2+1) 
            uzel_hra(kde) = 1; % je to hranicni uzel
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
% pocitame zatim s tim, ze g=0 !!!
for k = NN:-1:1 % vymazani radku a sloupcu pro hranicni uzly
    if (uzel_hra(k)==1)
        A(k,:) = []; A(:,k) = []; B(k) = [];
    end;
end;
% spy(A)
u = A\B; % vyreseni soustavy
u2 = reshape(u,N1-1,N2-1); % % preusporadani do tvaru oblasti
surf(u2); % vykresleni hodnot na vnitrnich uzlech

function z = Fcea(x1,x2)
z = 1+3*x1;
function z = Fcep(x1,x2)
z = 1;
function z = Fcef(x1,x2)
z = 1;
function z = Fceg(x1,x2)
z = 0;