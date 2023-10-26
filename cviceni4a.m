% resim y''+py'+qy = f, y(0)=..., y(L)=... nebo y'(L)=...
% metoda siti 
function abc
L = 2; % delka intervalu
OP2_H = 1; % 0 pro derivaci v prave OP, 1 pro hodnotu v prave OP
a1 = 1; % OP vlevo
a2 = 1; % OP vpravo
a2d = -4; % OP vpravo - derivace
N = 70; % pocet intervalu

h = L/N; % delka kroku
x = linspace(0,L,N+1); % sit bodu x_k

% sestavime matici ulohy a vektor prave strany:
A = zeros(N+1); % matice
B = zeros(N+1,1); % prava strana
for r = 2:N
    A(r,r-1:r+1) = [1/h^2-FceP(x(r))/2/h,   -2/h^2+FceQ(x(r)),   1/h^2+FceP(x(r))/2/h];
    B(r) = FceF(x(r));
end;
B = B - A(:,1)*a1; % prevedeni na pravou stranu

if (OP2_H==1) % pro Dirichletovu podminku vravo
    B = B - A(:,N+1)*a2; % prevedeni na pravou stranu
    % vymazani prvnich a poslednich sloupcu a radku:
    A(1,:) = [];
    A(end,:) = [];
    A(:,1) = [];
    A(:,end) = [];
    B(1) = [];
    B(end) = [];
    y = A\B
    plot(x,[a1,y',a2]);
else % pro Neumannovu podminku vravo
    A(N+1,N:N+1) = [2/h^2,  -2/h^2+FceQ(x(N+1))];
    B(N+1) = -a2d*2/h - FceP(x(N+1))*a2d + FceF(x(N+1));
    A(1,:) = [];    
    A(:,1) = [];    
    B(1) = [];    
    y = A\B
    plot(x,[a1,y']);
end;

function y = FceP(x)
y = 2+x;
function y = FceQ(x)
y = 6+x^2;
function y = FceF(x)
y = 1/(x+1);







