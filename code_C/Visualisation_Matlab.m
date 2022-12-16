
%% Visualisation potentiel & champ electrique Matlab
% Importer les fichiers txt suivants ( -> Numeric Matrix )
% Vrelax pour potentiel
% Ex Ey pour champ electrique 
% E pour amplitude (Attention: fixer Delimited, Column delimiters: space)

clf;
V = Vrelax;

[Ny,Nx] = size(V);
x = (1:Nx);
y = (1:Ny);

% affichage potentiel electrique
figure(1)
contour_range_V = -101:0.5:101;
contour(x,y,V,contour_range_V,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',13);
xlabel('Axe Y (m)','fontsize',13);
ylabel('Axe X (m)','fontsize',13);
title('Potentiels V(x,y) en V','fontsize',13);
axis equal;
h1=gca;
set(h1,'fontsize',13);
fh1 = figure(1); 
set(fh1, 'color', 'white')

% affichage ligne de champ electrique
figure(2)
contour(x,y,E,'linewidth',0.5);
hold on, quiver(x,y,Ey,Ex,2) %quiver(x,y,Ey,Ex,2) pour charges uniquement
title('Lignes de champ electrique, E (x,y) en V/m','fontsize',13);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',13);
xlabel('Axe Y (m)','fontsize',13);
ylabel('Axe X (m)','fontsize',13);
axis equal;
h2=gca;
set(h2,'fontsize',13);
fh3 = figure(2); 
set(fh3, 'color', 'white')


