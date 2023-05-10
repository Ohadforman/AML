clc,clear
%Author: Hassan Najjar
%%%%%%%%%%%%%%%%%%%%% PART A:  BUTTERWORTH FILTER DESIGN %%%%%%%%%%%%%%%%%%%%
%equal ripple low-pass filter prototypes (g0=1, wc=1, N=1:10, 0,5 dB & 3.0
%dB ripple)
%we use the Pi network prototype (Pozar P. 448)

TYPE=input('FILTER TYPE:\n    OPTION (A) LOW-PASS\n    OPTION (B) HIGH-PASS\n    OPTION (C) BAND-PASS\n    OPTION (D) BAND-STOP\n FILTER TYPE OPTION: A,B,C,OR D (in single quotes) = ');
RIPPLE=input('RIPPLE: \n    OPTION (1) 0.5 dB\n    OPTION (2) 3 dB \n    OPTION (3) MAXIMALLY FLAT\n RIPPLE OPTION: 1,2,OR 3 (without quotes) = ');
N=input('FILTER ORDER: (INTEGER N 1 TO 10)\n N = ');
if ( (TYPE=='A') || (TYPE=='B') )
    fc=input('CUTOFF FREQUENCY (GHz)\n fc = ');
    fc=fc*10^9;
elseif ( (TYPE=='C') || (TYPE=='D') )
     f1=input('FIRST EDGE FREQUENCY (GHz)\n f1 = ');
     f2=input('SECOND EDGE FREQUENCY (GHz)\n f2 = ');
     f1=f1*10^9;
     f2=f2*10^9;    %convert to Hz
end
Z0=input('INPUT LINE CHARACTERISTIC IMPEDANCE (Ohms)\n Z0 = ');
d=input('MICROSTRIP HEIGHT (mm)\n d = ');
d=d*10^-3;
er=input('SUBSTRATE DIELECTRIC CONSTANT\n er = ');
%----------------------------------------0.5 dB RIPPLE-------------------------------------%
g(:,:,1) = [0.6986 1.0000 0      0      0      0      0      0      0      0      0;  ...
            1.4029 0.7071 1.9841 0      0      0      0      0      0      0      0;  ...
            1.5963 1.0967 1.5963 1.0000 0      0      0      0      0      0      0; ...
            1.6703 1.1926 2.3661 0.8419 1.9841 0      0      0      0      0      0; ...
            1.7058 1.2296 2.5408 1.2296 1.7058 1.0000 0      0      0      0      0; ...
            1.7254 1.2479 2.6064 1.3137 2.4758 0.8696 1.9841 0      0      0      0; ...
            1.7372 1.2583 2.6381 1.3444 2.6381 1.2483 1.7372 1.0000 0      0      0; ...
            1.7451 1.2647 2.6564 1.3590 2.6964 1.3389 2.5093 0.8796 1.9841 0      0; ...
            1.7504 1.2690 2.6678 1.3673 2.7239 1.3673 2.6678 1.2690 1.7504 1.0000 0; ...
            1.7543 1.2721 2.6754 1.3725 2.7392 1.3806 2.7231 1.3485 2.5239 0.8842 1.9841];
%----------------------------------------3 dB RIPPLE---------------------------------------%
g(:,:,2) = [1.9953 1.0000 0      0      0      0      0      0      0      0      0; ...
            3.1013 0.5339 5.8095 0      0      0      0      0      0      0      0; ...
            3.3487 0.7117 3.3487 1.0000 0      0      0      0      0      0      0; ...
            3.4389 0.7483 4.3471 0.5920 5.8095 0      0      0      0      0      0; ...
            3.4817 0.7618 4.5381 0.7618 3.4817 1.0000 0      0      0      0      0; ...
            3.5045 0.7685 4.6061 0.7929 4.4641 0.6933 5.8095 0      0      0      0; ...
            3.5182 0.7723 4.6386 0.8039 4.6386 0.7723 3.5182 1.0000 0      0      0; ...
            3.5277 0.7745 4.6575 0.8089 4.6990 0.8018 4.4990 0.6073 5.8095 0      0; ...
            3.5340 0.7760 4.6692 0.8118 4.7272 0.8118 4.6692 0.7760 3.5340 1.0000 0; ...
            3.5384 0.7771 4.6768 0.8136 4.7425 0.8164 4.7260 0.8051 4.5142 0.6091 4.8095];
%----------------------------------------MAXIMALLY FLAT-----------------------------------%
g(:,:,3) = [2.0000 1.0000 0      0      0      0      0      0      0      0      0; ...
            1.4142 1.4142 1.0000 0      0      0      0      0      0      0      0; ... 
            1.0000 2.0000 1.0000 1.0000 0      0      0      0      0      0      0; ... 
            0.7654 1.8478 1.8478 0.7654 1.0000 0      0      0      0      0      0; ...
            0.6180 1.6180 2.0000 1.6180 0.6180 1.0000 0      0      0      0      0; ...
            0.5176 1.4142 1.9318 1.9318 1.4142 0.5176 1.0000 0      0      0      0; ...
            0.4450 1.2470 1.8019 2.0000 1.8019 1.2470 0.4450 1.0000 0      0      0; ... 
            0.3902 1.1111 1.6629 1.9615 1.9615 1.6629 1.1111 0.3902 1.0000 0      0; ...
            0.3473 1.0000 1.5321 1.8794 2.0000 1.8794 1.5321 1.0000 0.3473 1.0000 0; ...
            0.3129 0.9080 1.4142 1.7820 1.9754 1.9754 1.7820 1.4142 0.9080 0.3129 1.0000];
%----------------------------------------LOWPASS PROTOTYPE-------------------------------%
RS = 1;                      %normalized values
for k=1:N
    if ( mod(k,2)==0 )
        L(k) = g(N,k,RIPPLE);
    else 
        C(k) = g(N,k,RIPPLE);
    end
RL = g(N,N+1,RIPPLE);
end
                     %Impedance and frequency scaling
RS = Z0;
RL = Z0*RL;
C = C/Z0;
L = L*Z0;
if ( TYPE == 'A' )
    wc = 2*pi*fc;
    L = L/wc;
    C = C/(wc);
end
%--------------------------------------HIGHPASS TRANSFORMATION-----------------------------%
if ( TYPE == 'B' )
    wc = 2*pi*fc;
    C = 1/(wc*L);
    L = 1/(wc*C);
end
%--------------------------------------BANDPASS TRANSFORMATION-----------------------------%
if ( TYPE == 'C' )
     f0=sqrt(f1*f2);    %center frequency (geometric mean)
     w0=2*pi*f0;
     DELTA = (f2 - f1)/f0;
     for k=1:2:N
        Lp(k) = DELTA/(w0*C(k));  
        Cp(k) = C(k)/(DELTA*w0);  %shunt elements
     end
     for l=2:2:N
        Lp(l) = L(l)/(DELTA*w0);
        Cp(l) = DELTA/(w0*L(l));  %series elements
     end
    L = Lp;
    C = Cp;
end
%--------------------------------------BANDSTOP TRANSFORMATION-----------------------------%
if ( TYPE == 'D' )
     f0=sqrt(f1*f2);    %center frequency (geometric mean)
     w0=2*pi*f0;
     DELTA = (f2 - f1)/f0;
     for k=1:2:N
        Lp(k) = 1/(w0*DELTA*C(k));  
        Cp(k) = DELTA*C(k)/w0;  %shunt elements
     end
     for l=2:2:N
        Lp(l) = DELTA*L(l)/w0;
        Cp(l) = 1/(w0*DELTA*L(l));  %series elements
     end
    L = Lp;
    C = Cp;
end
%%%%%%%%%%%%%%%%%%%%%% PART B: MICROSTRIP IMPLEMENTATION %%%%%%%%%%%%%%%%%%%%%%
if ( (TYPE == 'A') || TYPE == 'B' ) %for lowpass/highpass implement stepped-impedance microstrip
    disp('FOR STEPPED-IMPEDANCE IMPLEMENTATION');
    Zh = input('HIGHEST PRACTICAL LINE IMPEDANCE (Ohms):\n Zh = ');
    Zl = input('LOWEST PRACTICAL LINE IMPEDANCE (Ohms):\n Zl = ');
    for k=1:N
        if (mod(k,2) == 0)
            bl(k) = g(N,k,RIPPLE)*Z0/Zh;
        else
            bl(k) = g(N,k,RIPPLE)*Zl/Z0;
        end
    end
    c=299792458;                                  %speed of light
    lambda0=c/fc;                                 %vacuum wavelength
    k0=2*pi/lambda0;                             % ===   wavenumber
    %find the width of the low, and high impedance strips in order to
    %calculate the lenghts of each
    Al = Zl/60*sqrt((er+1)/2)+(er-1)/(er+1)*(0.23+0.11/er);
    Ah = Zh/60*sqrt((er+1)/2)+(er-1)/(er+1)*(0.23+0.11/er);    
    Bl = 377*pi/(2*Zl*sqrt(er));
    Bh = 377*pi/(2*Zh*sqrt(er));
    Wratiodl(1) = 8*exp(Al)./(exp(2*Al)-2);      % w/d < 2
    Wratiodl(2) = 2/pi*(Bl-1-log(2*Bl-1)+(er-1)/(2*er).*(log(Bl-1)+0.39-0.61/er)); %w/d>2
    Wratiodh(1) = 8*exp(Ah)./(exp(2*Ah)-2);      % w/d < 2
    Wratiodh(2) = 2/pi*(Bh-1-log(2*Bh-1)+(er-1)/(2*er).*(log(Bh-1)+0.39-0.61/er)); %w/d>2
    if (Wratiodl(1)<2)
        Wl=d*Wratiodl(1);
    elseif (Wratiodl(2)>2)
        Wl=d*Wratiodl(2);
    end
    if (Wratiodh(1)<2)
        Wh=d*Wratiodh(1);
    elseif (Wratiodh(2)>2)
        Wh=d*Wratiodh(2);
    end
    epsel=(er+1)/2+(er-1)/2*1./(sqrt(1+12*d./Wl));   %effective dielctric constant in microstrip
    epseh=(er+1)/2+(er-1)/2*1./(sqrt(1+12*d./Wh));   %effective dielctric constant in microstrip
    betal=k0*sqrt(epsel);                            %propagation constant in microstrip
    betah=k0*sqrt(epseh);                            %propagation constant in microstrip
    ll = bl/betal;
    lh = bl/betah;
    disp('NUMBER LOW, AND HIGH MICROSTRIP LINES:'), N
    disp('LENGTHS OF LOW, AND HIGH MICROSTRIP LINES(mm):');
    disp('Llow = '), ll*10^3
    disp('Lhigh = '), lh*10^3
    disp('WIDTHS OF LOW, AND HIGH MICROSTRIP LINES (mm):');
    disp('Wlow = '), Wl*10^3
    disp('Whigh = '), Wh*10^3
    ZL = RL;
    disp('THE LOAD RL MUST BE (Ohms)'), RL
    
elseif ( (TYPE == 'C') || (TYPE == 'D') ) %for bandpass implement coupled line strips
    Z0J(1) = sqrt(pi*DELTA/(2*g(N,1,RIPPLE)));
    for k=2:N+1
        Z0J(k) = sqrt(pi*DELTA/(2*g(N,k,RIPPLE)*g(N,k-1,RIPPLE)));
    end
    Z0e=Z0*(1+Z0J+Z0J.^2);
    Z0o=Z0*(1-Z0J+Z0J.^2);

%%find width of microstrips
%find the width of the microstrip
A = Z0e/60*sqrt((er+1)/2)+(er-1)/(er+1)*(0.23+0.11/er);
B = 377*pi./(2*Z0e*sqrt(er));

Wratiod(1,:) = 8*exp(A)./(exp(2*A)-2);      % w/d < 2
Wratiod(2,:) = 2/pi*(B-1-log(2*B-1)+(er-1)/(2*er).*(log(B-1)+0.39-0.61/er)); %w/d>2

for k=1:length(Wratiod)
    if (Wratiod(1,k)<2)
        W(1,k)=d*Wratiod(1,k);
    elseif (Wratiod(2,k)>2)
        W(2,k)=d*Wratiod(2,k);
    end
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('LUMPED ELEMENT VALUES'), L, C
disp('FOR BANDPASS IMPLEMENTED IN MICROSTRIP, USE COUPLED LINES WITH EVEN, AND ODD')
disp('CHARACTERISTIC IMPEDANCES (Ohms)'), Z0e, Z0o
disp('EACH LINE HAVING A WIDTH (mm)'), disp(W*10^3);

end %end if for coupled line





