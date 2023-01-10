function [c, phi, w]=amplitudenspektrum(zeit, messwerte)
    trend=polyfit(zeit,messwerte,1);
    werte=messwerte-trend(1)*zeit-trend(2); %um Trend und Offset reduzierte Zeitreihe

% a_n und b_n
    I=length(werte);
%a
    a_0= 1/I*(sum(werte)); 

    for n=1:(I/2-1)         
        a_sum(n)=werte(1);
        for i=2:I
           a_sum(n)=a_sum(n)+cos(2*pi*n*(i-1)/I)*werte(i);
        end
        a(n)=a_sum(n)*2/I; 
    end
    
    a500=werte(1);         
    for i=2:I
       a500=a500+(-1)^(i-1)*werte(i);
    end
    a_500=a500/I;

    aTotal=zeros((I/2+1),1); aTotal(1)=a_0; aTotal(2:I/2)=a; aTotal(I/2+1)=a_500;        
    
%b
    for n=1:(I/2-1)           
        b_sum(n)=0;
        for i=1:I
           b_sum(n)=b_sum(n)+sin(2*pi*n*(i-1)/I)*werte(i);
        end
        b(n)=b_sum(n)*2/I;   
    end   

    bTotal=zeros(I/2+1,1); bTotal(2:I/2)=b; 
    
% Amplituden Spektren 
c=zeros(I/2,1); phi=zeros(I/2,1);
    for n=1:I/2
        c(n)= sqrt(aTotal(n)^2+bTotal(n)^2);
        phi(n)=atan2(bTotal(n),aTotal(n))*180/pi; %in grad
    end

 % w, zur Darstellung
delta_t= zeit(2)-zeit(1);
V_G=1/(2*delta_t);
W_G=pi/delta_t;
T=length(zeit)*delta_t;
delta_W=(2*pi)/T;

%Frequenzachse
for i=1:I/2+1 
    wRad(i)=(i-1)*delta_W;
end
w = wRad*180/(pi); %in grad
end

