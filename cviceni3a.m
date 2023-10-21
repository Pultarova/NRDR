% resim y''=f(x,y,y'), y(0)=1, y(5*pi/2)=3 na (0,5*pi/2)
% Eulerova metoda: N kroku
% metoda strelby
function abc
clc % vymaze command window
T = 5/2*pi; % nastaveni delky casoveho intervalu
y_poc = 1; % nastaveni pocatecni hodnoty y
% y_der_poc = -1; % nastaveni pocatecni derivace y
N = 1000; % nastaveni poctu kroku
h = T/N; % vypocet delky 1 kroku
x = linspace(0,T,N+1);
tol = 1e-3; % tolerance pro y(T)

% metoda puleni intervalu:
a = -1  % y(end)<3
b = 4  % y(end)>3
for k = 1:100
    c = (a+b)/2
    y = Resic(T,y_poc,c,N);
    yn = y(end)
    if (yn<3) a=c; end;
    if (yn>3) b=c; end;
    if (abs(yn-3)<tol) break; end;
end;

cla;
hold on;
plot(x,y,'b') % priblizne reseni


function z = F(x,y1,y2) % funkce popisujici difer. rovnici
z = -y1;

function  yy  = Resic(T,y_poc,y_der_poc,N) % Eulerova metoda
h = T/N;
y1 = zeros(1,N+1);
y2 = zeros(1,N+1);
x = linspace(0,T,N+1);
y1(1) = y_poc;
y2(1) = y_der_poc;
for k = 1:N
    y1(k+1) = y1(k) + h*y2(k); 
    y2(k+1) = y2(k) + h*F(x(k),y1(k),y2(k));
end;
yy = y1; % vraci cele reseni



