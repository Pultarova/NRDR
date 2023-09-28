% Eulerova metoda pro y'=f(x,y), y(0)=y0;

function EM
x0 = 1;
y0 = 1.001;
T = 50; % celkovy casovy interval (0,T)
N = 200; % pocet casovych kroku
h = T/N;
for j = 1:N
    y1 = y0 + h*f(x0,y0);
    y0 = y1;
    x0 = x0+h;
end;
y1
PresneReseni = 0.1*exp(3*T); % nahodou zname presne reseni


function z = f(x,y)
z = 3*y;
z = (y-1);



