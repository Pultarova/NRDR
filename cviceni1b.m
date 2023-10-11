% resim y'=f(x,y), y(0)=y0 na (0,T)
% Eulerova metoda: N kroku
function abc
clc % vymaze command window
T = 5; % nastaveni delky casoveho intervalu
y_poc = 1; % nastaveni pocatecni podminky
N = 200; % nastaveni poctu kroku
h = T/N; % vypocet delky 1 kroku

% pocitame jen y(T):
y0 = y_poc;
x0 = 0;
for k = 1:N
    %y1 = y0 + h*(-3)*y0;
    % y1 = y0 + h*y0*sin(x0*y0); %% <<---
    y1 = y0 + h*F(x0,y0);
    y0 = y1;
    x0 = x0 + h;
end;
y1 % vysleddne priblizne y(T)
% y_presne = exp(-3*T); % presne y(T) pro urcitou zadanou funkci f
% pocitame y na (0,T):
y = zeros(1,N+1);
x = linspace(0,5,N+1);
y(1) = y_poc;
for k = 1:N
    % y(k+1) = y(k) + h*y(k)*sin(x(k)*y(k)); %% funcki f zadavame zde
    y(k+1) = y(k) + h*F(x(k),y(k)); %% funkci f definujeme jako matlabovskou funkci F a volame
end;
cla;
hold on;
plot(x,y,'b')
% set(gca,'YScale','log')
% plot(x,exp(-3*x),'r')

function z = F(x,y)
z = y*sin(x*y);
