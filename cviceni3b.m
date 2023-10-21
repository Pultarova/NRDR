% resim y''-3y'+2y = x*cos(x), y(0)=1, y(2)=3
% metoda siti - zatim jen rychle napsany jednoduchy kod - bude doplnen
N = 70; % pocet intervalu
h = 2/N; % delka kroku
x = linspace(0,2,N+1); % sit bodu x_k

% sestavime matici ulohy a vektor prave strany:
A = zeros(N+1); % matice
B = zeros(N+1,1); % prava strana
for r = 2:N
    A(r,r-1:r+1) = [1/h^2+3/2/h,   -2/h^2+2,   1/h^2-3/2/h];
    B(r) = x(r)*cos(x(r));
end;
B = B - A(:,1)*1; % prevedeni na pravou stranu
B = B - A(:,N+1)*3; % prevedeni na pravou stranu

% vymazani prvnich a poslednich sloupcu a radku:
A(1,:) = [];
A(end,:) = [];
A(:,1) = [];
A(:,end) = [];
B(1) = [];
B(end) = [];

y = A\B
plot(x,[1,y',3])

