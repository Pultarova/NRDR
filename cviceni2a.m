% resim y'=f(x,y), y(0)=y0 na (0,T)
% Eulerova metoda: N kroku
% Runge-Kutta 4. radu
function abc
clc % vymaze command window
T = 5; % nastaveni delky casoveho intervalu
y_poc = 1; % nastaveni pocatecni podminky
N = 10; % nastaveni poctu kroku
h = T/N; % vypocet delky 1 kroku

% pocitame jen y(T):
y0 = y_poc;
yRK0 = y_poc;
x0 = 0;
for k = 1:N    
    y1 = y0 + h*F(x0,y0); % Eulerova metoda
    y0 = y1;    

    k1 = F(x0,yRK0);
    k2 = F(x0+h/2,yRK0+h*k1/2);
    k3 = F(x0+h/2,yRK0+h*k2/2);
    k4 = F(x0+h,yRK0+h*k3);
    yRK1 = yRK0 + h/6*(k1+2*k2+2*k3+k4); % Rungova-Kuttova metoda
    yRK0 = yRK1;

    x0 = x0 + h;
end;
y1 % vysledne priblizne y(T) - m Euler
yRK1 % vysledne priblizne y(T) - m RK
y_presne = exp(-3*T) % presne y(T) pro urcitou zadanou funkci f
% return

% pocitame y na (0,T):
y = zeros(1,N+1);
yRK = zeros(1,N+1);
x = linspace(0,5,N+1);
y(1) = y_poc;
yRK(1) = y_poc;
for k = 1:N
    % Eulerova metoda:
    y(k+1) = y(k) + h*F(x(k),y(k)); %% funkci f definujeme jako matlabovskou funkci F a volame
    % Rungova-Kuttova metoda:
    x0 = x(k);
    yRK0 = yRK(k);
    k1 = F(x0,yRK0);
    k2 = F(x0+h/2,yRK0+h*k1/2);
    k3 = F(x0+h/2,yRK0+h*k2/2);
    k4 = F(x0+h,yRK0+h*k3);
    yRK(k+1) = yRK(k) + h/6*(k1+2*k2+2*k3+k4);
end;
cla;
hold on;
plot(x,y,'b') % kresleni E. m.
plot(x,yRK,'k') % kresleni RK. m.
plot(x,exp(-3*x),'r--') % kresleni presneho reseni
% set(gca,'YScale','log')

function z = F(x,y)
% z = y*sin(x*y);
z = -3*y;

