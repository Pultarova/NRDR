% resim y''+py'+qy = f, y(0)=a1, y(L)=a2 nebo y'(L)=a2d
% metoda siti 
function abc
L = 2; % delka intervalu
OP2_H = 0; % 0 pro derivaci v prave OP, 1 pro hodnotu v prave OP
a1 = 5; % OP vlevo
a2 = 1; % OP vpravo
a2d = 0; % OP vpravo - derivace

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

cla; hold on;

if (OP2_H==1) % pro Dirichletovu podminku vpravo
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
else % pro Neumannovu podminku vpravo
    A(N+1,N:N+1) = [2/h^2,  -2/h^2+FceQ(x(N+1))];
    B(N+1) = -a2d*2/h - FceP(x(N+1))*a2d + FceF(x(N+1));
    A(1,:) = [];    
    A(:,1) = [];    
    B(1) = [];    
    y = A\B
    plot(x,[a1,y']);
end;

% kolokacni metoda s polynomem 4. stupne:
Akol = zeros(5);
Bkol = zeros(5,1);
Akol(1,1) = 1;
Bkol(1) = a1;
xx = [0,1/2,1,3/2,2];
for r = 2:4
    z = xx(r);
    p0 = [1,z,z^2,z^3,z^4];
    p1 = [0,1,2*z,3*z^2,4*z^3];
    p2 = [0,0,2,6*z,12*z^2];
    Akol(r,:) = p2+p1*FceP(z)+p0*FceQ(z);
    Bkol(r) = FceF(z);
end;
z = L;
if (OP2_H==1) % pro Dirichletovu podminku vpravo
    Akol(5,:) = [1,z,z^2,z^3,z^4];
    Bkol(5) = a2;
else % pro Neumannovu podminku vpravo
    Akol(5,:) = [0,1,2*z,3*z^2,4*z^3];
    Bkol(5) = a2d;
end;
Akol
c = Akol\Bkol
plot(x,c(1)+c(2)*x+c(3)*x.^2+c(4)*x.^3+c(5)*x.^4,'r')

legend('metoda siti','kolokacni metoda');


function y = FceP(x)
y = 2+x;
function y = FceQ(x)
y = 6+x^2;
function y = FceF(x)
y = 1/(x+1);






