// codigo do Scilab:
clear()


// Definição dos parâmetros iniciais da analogia:
m1 = 556.454;
m2 = 121.554;
m0 = 15408.992;
k1 = 192842.44;
k2 = 71003.48;
omega1 = 18.616;
omega2 = 24.169;
c1 = 6.289;
c2 = 1.784;
h1 = -0.1185548;
h2 = 0.03291344;
h0 = 1.4996E-6;
l = 2.6;
I0 = 647.31;
c_d = 144.524;
k = 5400;
g = 9.8;

// Definição dos coeficientes da matriz A:
If = m0 * h0^2 + I0 + m1 * h1^2 + m2 * h2^2;
a61 = m1 * g / If;
a62 = m2 * g / If;
a63 = (m1 * h1 + m2 * h2 - m0 * h0 - k * l^2/2) / If;
a66 = -c_d * l^2 / (2 * If);

// Definição das matrizes do espaço de estados:
A = [0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 1; -k1/m1, 0, g, -c1/m1, 0, 0; 0, -k2/m2, g, 0, -c2/m2, 0; a61, a62, a63, 0, 0, a66];
B = [0; 0; 0; 0; 0; 1/If];
C = [0, 0, 1, 0, 0, 0];

// Encontrar os autovalores da matriz A (polos do sistema):
auto = spec(A);

// Plotagem dos polos no plano complexo:
scf(1);
for i = 1:length(auto),
    plot(real(auto(i)), imag(auto(i)), 'marker', 'o');
end
set(gca(), 'grid', [1 1]);
a = gca();
a.x_location = 'origin';
a.y_location = 'origin';
xlabel('Polos no plano complexo');

// Determinação da função de transferência do sistema a partir do espaço de estados:
s1 = syslin('c', A, B, C);
s = poly([0 1], 's', 'coeff');
I = eye(A);
h = C * ((s * I - A)^(-1)) * B;
func_transferencia = syslin('c', h);

// Obtenção do numerador e do denominador:
num = h.num;
den = h.den;

// Tabela de Routh-Horwitz:
estab = routh_t(den);

// Plotagem do gráfico de BODE:
scf(2);
bode(func_transferencia, 0.001, 1E5);
title('Resposta no domínio da frequência');

scf(3);
bode(func_transferencia, 0.1, 10);
title('Resposta no domínio da frequência');
// Definição do vetor do tempo:
t0 = 0;
tf = 40;
delta_t = 0.01;
t = t0:delta_t:tf;

// Condições iniciais:
x0 = [0; 0; 0; 0; 0; 0];

// Resolução do sistema pela matriz de transição (código desenvolvido na lista F):
function [x] = integrador3(A, B, C, u, t, x0)
    n = length(x0);
    delta_t = t(2) - t(1);
    ordem_de_expansao = 15;

    function [f] = phi(delta_t)
        f = eye(n, n);
        for k = 1:ordem_de_expansao - 1
            f = f + (A^k) * (delta_t^k) / factorial(k);
        end
    end

    function [g] = gama(delta_t)
        g = delta_t * eye(n, n);
        for k = 1:ordem_de_expansao - 1
            g = g + (A^k) * (delta_t^(k+1)) / factorial(k+1);
        end
    end

    y = zeros(n, length(t));
    x(:, 1) = x0;
    y(:, 1) = C * x(:, 1);

    for i = 2:length(t)
        x(:, i) = phi(delta_t) * x(:, i-1) + gama(delta_t) * B * u(:, i-1);
        y(:, i) = C * x(:, i);
    end
end

// Entradas do sistema na forma de impulso:
u = zeros(1, length(t));
u(1, 1) = 1;
z = integrador3(A, B, C, u, t, x0);

scf(1);
plot2d(t, z(3, :), 2);
xlabel('Tempo (s)');
ylabel('Posição psi do tanque (rad)');
title('Resposta do sistema a um impulso unitário');

// Impulso real:
u = zeros(1, length(t));
u(1, 1) = (2 * 1.24 * 16087.5 * (60 / 3.6)^2) / 50;
z = integrador3(A, B, C, u, t, x0);

scf(2);
plot2d(t, z(3, :), 2);
xlabel('Tempo (s)');
ylabel('Posição psi do tanque (rad)');
title('Resposta do sistema a um impulso de uma curva de R=50m');

scf(3);
plot2d(t, z(1, :), 2);
xlabel('Tempo (s)');
ylabel('Posição x1 da massa do primeiro modo de vibração (m)');
title('Resposta do sistema a um impulso de uma curva de R=50m');

// Entrada do sistema na forma de degrau:
u = zeros(1, length(t));
for i = 1:length(t)/2
    u(1, i) = 1;
end
z = integrador3(A, B, C, u, t, x0);

scf(4);
plot2d(t, z(3, :), 2);
xlabel('Tempo (s)');
ylabel('Ângulo psi do tanque (rad)');
title('Resposta do sistema a um degrau unitário agindo por 20s');

// Degrau real:
u = zeros(1, length(t));
for i = 1:length(t)/4
    u(1, i) = (2 * 1.24 * 16087.5 * (60 / 3.6)^2) / 50;
end
z = integrador3(A, B, C, u, t, x0);

scf(5);
plot2d(t, z(3, :), 2);
xlabel('Tempo (s)');
ylabel('Ângulo psi do tanque (rad)');
title('Resposta do sistema a um degrau em uma curva de R=50m');

// Entrada do sistema na forma senoidal:
u = zeros(1, length(t));
f = 1; // Frequência de 1Hz
omega = 2 * %pi * f;
u(1, :) = sin(omega * t);
z = integrador3(A, B, C, u, t, x0);

scf(6);
z(3, 1) = 1E-6;
plot(t, z(3, :), 2);
xlabel('Tempo (s)');
ylabel('Ângulo psi do tanque (rad)');
title('Resposta a uma entrada senoidal com frequência de 1Hz');
