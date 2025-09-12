function main()
clear all
close all
open_figure(1);
axis_plot(6.0,0.5);

%intervallo di definizione delle curve di bezier
a = 0;
b = 1;
np = 50;

%% interpolazione dei punti per la prima S
Q1 = [5,5;1.9,4.3;1.9,3.1;3.1,1.9;3.1,0.7;0,0]; %punti da interpolare
param = 0;                                      %parametrizzazione uniforme
bezQ1 = curv2_bezier_interp(Q1, a, b, param);   %interpolazione
curv2_bezier_plot(bezQ1, np, 'r-');             %disegno della curva 1

m=6;                        %numero di punti, grado è m-1 
t=linspace(0,1,m);          %genera i punti equispaziati su metà curva
Q=decast_val(bezQ1,t);      %valuta la curva iniziale nei punti equispaziati tramite de Casteljau
point_plot(Q,'ko');         %disegna i punti in cui è campionata la curva
param=0;                    %parametrizzazione uniforme

%crea una curva per interpolazione dei punti trovati
bezD = curv2_bezier_interp(Q, a, b, param); 
curv2_bezier_plot(bezD, np, 'r-',2);            %la disegna

%% Simmetria della curva rispetto alla retta verticale x = 5
bezP1 = bezD;                             %crea una copia della curva Bezier originale su cui applicare la trasformazione
bezP1.cp(:,1) = -bezP1.cp(:,1) + 2 * 5;   %riflette i punti di controllo rispetto a x = 5

%% Rotazione della curva riflessa di pi/4 rispetto al punto [5, -5]
b = [5, -5];               %punto intorno a cui ruotare
T = get_mat_trasl(b);      %matrice di traslazione per portare al punto b
R = get_mat2_rot(pi/4);    %matrice di rotazione di pi/4
T_inv = get_mat_trasl(-b); %matrice di traslazione inversa 
M = T * R * T_inv;                        %combinazione delle trasformazioni
bezP1.cp = point_trans(bezP1.cp, M);      %applica la trasformazione ai punti di controllo
curv2_ppbezier_plot(bezP1, np, 'b-',2);   %traccia la curva blu

%% definisco un petalo della figura 
%creo un'unica curva tramite intersezione
[A, t1, t2] = curv2_intersect(bezP1, bezD); %calcola i punti di intersezione tra le curve di Bezier
%point_plot(A', 'ro', 1, 'r');              %disegna i punti di intersezione

%suddivisione della curva in due sottocurve in corrispondenza di un parametro t di intersezione
[sx, dx] = ppbezier_subdiv(bezP1, t1(1));
%curv2_ppbezier_plot(sx, np, 'm-', 2);      %traccia la parte della curva che mi interessa utilizzare

[sx1, dx1] = ppbezier_subdiv(bezD, t2(1));
%curv2_ppbezier_plot(sx1, np, 'm-', 2); 

bezA = curv2_ppbezier_join(sx1, sx, 1e0);   %join delle curve ottenute tramite la subdiv
%curv2_ppbezier_plot(bezA,np,'m-',2);       %fa il plot del join
%% figura due, disegno completo 
open_figure(2);
%inizializzo una curva che in un primo momento è copia di quella trovata
%con il join, la traslo e ruoto utilizzando la matrice di trasformazione
%precedente, unisco i tratti generati nel for nella curva a tratti "perimetro"
perimetro = bezA;  
for i = 1: 7
    bezA.cp = point_trans(bezA.cp, M);
    perimetro = curv2_ppbezier_join(perimetro, bezA, 1.0e-2);
end

%disegno e coloro la curva a tratti chiusa appena ottenuta
Pxy = curv2_ppbezier_plot(perimetro, -np);      %-np così che non disegni il tratto ma lo calcoli solo
point_fill(Pxy, 'b', 'r-',2 );

%% calcolo dell'area
%utilizzo la function per il calcolo dell'area, la function utilizza
%integral di matlab che fa uso del metodo adattivo per l'integrazione
val = curv2_ppbezier_area(perimetro);
fprintf('area della regione blu: %e\n',val);
%per ottenere un valore dell'area minore è necessario scalare la figura

%% figura due, disegno completo
set(gca, 'color', [1.0, 1.0, 0.6]) %sfondo

%% Circonferenza interna
%Punti della curva da interpolare campionati da una circonferenza
%a parametri equispaziati (parametrizzazione uniforme)
a=0; 
b=2*pi;
tpar=linspace(a,b,m);                       %trova m punti equispaziati nell'intervallo
[xp,yp,xp1,yp1]=cp2_circle(tpar, 1, -1);    %utilizzo la function per formare la circonferenza
C=[xp',yp'];                                %lista di punti da interpolare
C1=[xp1',yp1'];                             %lista vettori tangenti nei punti da interpolare

%chiama la funzione per interpolare valori e derivate, restituisce curva
%cubica a tratti
circle = curv2_ppbezierCC1_interp_der(C,C1,tpar);
%ritorna una curva di bezier a tratti

%% circonferernza con raggio maggiore
%scala la circonferenza, applico la trasformazione geometrica per modificare le dimensioni della curva 
s = 5;                                          %raggio della circonferenza
S = get_mat_scale([s, s]);
circle.cp = point_trans(circle.cp, S);          %applico la matrice di scala ai punti
Cxy=curv2_ppbezier_plot(circle,-np);            %-np così che non disegni il tratto ma lo calcoli solo
point_fill(Cxy, 'r');                           %disegno e coloro

end

%% funzioni
function [x, y, xp, yp] = cp2_circle(t, x0, y0)
    % espressione parametrica della curva circonferenza con traslazione
    x = cos(t) + x0;  % Traslazione sull'asse x, per modificare origine cerchio
    y = sin(t) + y0;  % Traslazione sull'asse y
    xp = -sin(t);
    yp = cos(t);
end