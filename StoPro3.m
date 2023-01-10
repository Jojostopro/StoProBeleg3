clc
clear
format long
close all

%Daten einlesen
data=load('aufgabe12w.txt');

t=data(:,2);
zr1=data(:,3);

clear data

%% Aufgabe 1 - Filter
zr1_36=zr1;
zr1_36(1:18)=[];
zr1_36(end-17:end)=[];
t32=t;
t32(1:16)=[];
t32(end-15:end)=[];
t36=t;
t36(1:18)=[];
t36(end-17:end)=[];

%Tiefpassfilter
m1=32;

for j = (1+m1/2) : (size(zr1,1)-m1/2)
 
  summe=0;
  for k=0:m1
      
    summe(:,1)=summe(:,1)+nchoosek(m1,k)*zr1(j+k-(m1/2),1);
  end
    tief32_raw(j-m1/2,1)=(1/(2^m1))*summe(:,1); %Formel 88, Kapitel 5.3.1
end
tief32_plot=tief32_raw;

%Hochpassfilter
m2=36;

for j = (1+m2/2) : (size(zr1,1)-m2/2)
  summe=0;
  for k=0:m2
    summe(:,1)=summe(:,1)+nchoosek(m2,k)*zr1(j+k-(m2/2),1);
  end
    tief36_raw(j-m2/2,1)=(1/(2^m2))*summe(:,1);
end
tief36=tief36_raw;

hoch36=zr1_36-tief36; %Formel Kapitel 5.3.3



%Bandpassfilter
tief32=tief32_raw;
tief32(1)=[];
tief32(2)=[];
tief32(end-1)=[];
tief32(end)=[];
band36=tief32-tief36;


% plots
figure('Name','Aufgabe 1')

subplot(4,1,1)
plot(t,zr1)
title('Ursprungsdatenreihe')
xlabel('Zeit in h')

subplot(4,1,2)
plot(t32,tief32_plot);
title('Tiefpassfilter m=32')
xlabel('Zeit in h')

subplot(4,1,3)
plot(t36,hoch36)
title('Hochpassfilter m=36')
xlabel('Zeit in h')

subplot(4,1,4)
plot(t36,band36)
title('Bandpassfilter m1=32 m2=36')
xlabel('Zeit in h')

%% Aufgabe 2
% Um wie viele Werte wird bei Anwendung des jeweiligen Filters die resultierende Zeitreihe 
% gegenueber der ursprünglichen Zeitreihe am Anfang und am Ende gekürzt? 

%% Aufgabe 3 Durchlasscharakteristik

deltaT=t(2)-t(1);

omegaG=pi/deltaT; 
%omegaG_grad = omegaG*180/pi;
T=1000*deltaT;
deltaOmega=2*pi/T;
frequenzachse=(0:deltaOmega:omegaG)'; %Einheit [rad]
frequenzachse_plot=frequenzachse.*180/pi;%Einheit [°]

G_tief(:,1)=cos((pi/2)*(frequenzachse/omegaG)).^m1; %sh. Seite 62

G_hoch(:,1)=1-cos((pi/2)*(frequenzachse/omegaG)).^m2; %Johnathan: .^m1

G_band(:,1)= G_tief(:,1) + G_hoch - 1; %Ausgeschrieben: =cos((pi/2)*(frequenzachse/omegaG)).^m1-cos((pi/2)*(frequenzachse/omegaG)).^m2;

%plots
figure('Name','Aufgabe 3')
subplot(3,1,1)
plot(frequenzachse_plot,G_tief)
title('Durchlasscharakteristik Tiefpassfilter')
ylabel('G(\omega)')
xlabel('Frequenz in °/h')

subplot(3,1,2)
plot(frequenzachse_plot,G_hoch)
title('Durchlasscharakteristik Hochpassfilter')
ylabel('G(\omega)')
xlabel('Frequenz in °/h')

subplot(3,1,3)
plot(frequenzachse_plot,G_band)
title('Durchlasscharakteristik Bandpassfilter')
ylabel('G(\omega)')
xlabel('Frequenz in °/h')

%% Aufgabe 4 Maximum Bandpassfilter

[max_band, pos_band] = max(G_band);
omegaMax_band = frequenzachse(pos_band)*180/pi; %Einheit [°]

fprintf('Maximum Bandpassfilter: omegaMax = %2.4f GMax = %2.4f\n',omegaMax_band,max_band);

verhaeltnis=omegaMax_band/(omegaG*180/pi);
fprintf('Verhaeltniss: omega_Max/omega_G = %2.5f\n',verhaeltnis)
%1/max_band;
band36_skaliert=band36/max_band;

G_band_skaliert=G_band/max_band;

figure('Name','Aufgabe 4')
subplot(2,1,1)
plot(t36,band36_skaliert)
title('Bandpassfilter skaliert')
xlabel('Zeit in h')

subplot(2,1,2)
plot(frequenzachse_plot,G_band_skaliert)
title('Durchlasscharakteristik Bandpassfilter skaliert') %So dass G(v) maximal = 1 ist
ylabel('G(\omega)')
xlabel('Frequenz in °/h')

%% Aufgabe 5 Amplitudenspektrum (Beleg 1)
[c_0, phi_0, w_0]=amplitudenspektrum(t, zr1); %Ausgangszeitreihe %Hier stimmt was nicht
[c_a, phi_a, w]=amplitudenspektrum(t36, tief36); %Mit Tiefpassfilter
[c_b, phi_b, w]=amplitudenspektrum(t36, hoch36); %Mit Hochpassfilter %Hier stimmt was nicht
[c_c, phi_c, w]=amplitudenspektrum(t36, band36); %Mit Bandpassfilter

figure('Name','Aufgabe 5')
plot(w_0(1:500), c_0,  w(1:482),c_a, w(1:482), c_b, w(1:482), c_c)
title ('Amplitudenspektrum')
xlabel('\omega [°/h]')
ylabel('c_n [m]')
legend('Ausgangszeitreihe','Mit Tiefpassfilter','Mit Hochpassfilter','Mit Bandpassfilter')
