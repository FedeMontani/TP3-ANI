Volumen = [0; 66; 481; 948];
Cota = [74; 76; 78; 80];

#Interpolacion por Lagrange

produ_denominadores_delta  = [];
for i = 1:4
    denominadores_delta = [];
    k = 0;
    for j = 1:4
            if j != i
                    k = k+1;
                      denominadores_delta(k) = (Volumen(i) - Volumen(j));
            endif
    endfor
    produ_denominadores_delta(i) = prod(denominadores_delta(1,:));
endfor

function L = l(x, Volumen,produ_denominadores_delta, Cota)  #polinomio interpolante
    delta0 = (x-Volumen(2)).*(x-Volumen(3)).*(x-Volumen(4))/produ_denominadores_delta(1);
    delta1 = (x-Volumen(1)).*(x-Volumen(3)).*(x-Volumen(4))/produ_denominadores_delta(2);
    delta2 = (x-Volumen(1)).*(x-Volumen(2)).*(x-Volumen(4))/produ_denominadores_delta(3);
    delta3 = (x-Volumen(1)).*(x-Volumen(2)).*(x-Volumen(3))/produ_denominadores_delta(4);
    L = delta0.*Cota(1) + delta1.*Cota(2) + delta2.*Cota(3) + delta3.*Cota(4);
endfunction

x = 0:25:1000;
plot(x,l(x, Volumen,produ_denominadores_delta, Cota))
hold on;
scatter(Volumen, Cota,100, 'filled')
hold off;
xlabel('Volumen [Hm^3]','fontsize',14)
ylabel('Cota [m]','fontsize',14)
title('Polinomio Interpolador de Lagrange','fontsize',14,'color','blue')
grid

#Ajuste cuadratico

fi0 = [1; 1; 1; 1];
fi1 = [0; 66; 481; 948];
fi2 = [0; 66^2; 481^2; 948^2];

A = [fi0'*fi0 fi1'*fi0 fi2'*fi0; fi0'*fi1 fi1'*fi1 fi2'*fi1; fi0'*fi2 fi1'*fi2 fi2'*fi2]; #matriz de sist. de ecuaciones normales
b = [Cota'*fi0; Cota'*fi1; Cota'*fi2]; #vector independiente del sist. de ecuaciones normales

c = inv(A)*b; #vector de coeficientes de la funcion de ajuste

function ajuste = a(x,c)
    ajuste = c(1) + c(2).*x + c(3).*x.^2;
endfunction

x = 0:25:1000;
plot(x,a(x,c))
hold on;
scatter(Volumen, Cota,100, 'filled')
hold off;
xlabel('Volumen [Hm^3]','fontsize',14)
ylabel('Cota [m]','fontsize',14)
title('Ajuste cuadratico','fontsize',14,'color','blue')
grid
