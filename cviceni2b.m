% resim y''=f(x,y,y'), y(0)=a, y'(0)=b na (0,T)
% Eulerova metoda: N kroku
function abc
clc % vymaze command window
T = 5; % nastaveni delky casoveho intervalu
y_poc = 1; % nastaveni pocatecni hodnoty y
y_der_poc = -1; % nastaveni pocatecni derivace y
N = 100; % nastaveni poctu kroku
h = T/N; % vypocet delky 1 kroku

% pocitame y na (0,T):
y1 = zeros(1,N+1);
y2 = zeros(1,N+1);
x = linspace(0,T,N+1);
y1(1) = y_poc;
y2(1) = y_der_poc;
for k = 1:N
    y1(k+1) = y1(k) + h*y2(k); 
    y2(k+1) = y2(k) + h*F(x(k),y1(k),y2(k));
end;
cla;
hold on;
plot(x,y1,'b') % priblizne reseni
plot(x,cos(x)-sin(x),'r') % presne reseni

function z = F(x,y1,y2)
z = -y1;