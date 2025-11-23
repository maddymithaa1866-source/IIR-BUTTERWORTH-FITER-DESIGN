# EXP 3 : IIR-BUTTERWORTH-FITER-DESIGN

## AIM: 

 To design an IIR Butterworth filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF):
```
clc ; 
close ; 
wp=input('Enter the pass band frequency (Radians )= ' ); 
ws=input('Enter the stop band frequency (Radians )= ' ); 
alphap=input( ' Enter the pass band attenuation (dB)=' ); 
alphas=input( ' Enter the stop band attenuation(dB)=' ); 
T=input('Enter the Value of sampling Time='); 
 
//Pre warping- Bilinear Transformation 
omegap=(2/T)*tan(wp/2); 
disp(omegap,'omegap='); 
omegas=(2/T)*tan(ws/2); 
disp(omegas,'omegas='); 
 
//Order of the filter 
 
N=log10(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1))/(2*log10(omegas/omegap)); 
disp(N,'N='); 
13 
 
N=ceil(N); 
 
disp(N,'Round off value of N='); 
 //Cut off frequency 
omegac=omegap/(((10^(0.1*alphap)) -1)^(1/(2* N))); 
disp(omegac,'omegac='); 
 
disp('Normalised Analog LPF Transfer function H(S)='); 
hs_Normalised = analpf(N,'butt',[0,0],1); 
disp(hs_Normalised); 
 
disp('Analog LPF Transfer function H(S)='); 
hs= analpf(N,'butt',[0,0],omegac); 
disp(hs); 
 
 
z=poly(0,'z');//Defining variable z 
 
Hz=horner(hs,(2/ T)*((z -1)/(z+1)))// Bilinear Transformation 
disp('Digital LPF Transfer function H(Z)='); 
disp(Hz); 
 
 
HW=frmag(Hz,512); // Frequency response 
w=0:%pi/511:%pi ; 
plot(w/%pi,abs(HW)); 
 
xlabel(' Normalized Digital Frequency w'); 
ylabel('Magnitude ');
title(' Frequency Response of Butterworth IIR LPF'); 

```
## PROGRAM (HPF): 
```
clc;
close;
wp = input('Enter the pass band frequency (Radians )= ');
ws = input('Enter the stop band frequency (Radians )= ');
alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');
T = input('Enter the Value of sampling Time= ');

// Pre-warping (Bilinear Transformation)
omegap = (2/T) * tan(wp/2);
disp(omegap, 'omegap=');
omegas = (2/T) * tan(ws/2);
disp(omegas, 'omegas=');

// Order of the filter
N = log10(((10^(0.1*alphas))-1) / ((10^(0.1*alphap))-1)) / (2*log10(omegas/omegap));
disp(N,'N=');
N = ceil(N);
disp(N,'Round off value of N=');

// Cut off frequency
omegac = omegap / (((10^(0.1*alphap))-1)^(1/(2*N)));
disp(omegac,'omegac=');

// Normalised Analog LPF Transfer function
disp('Normalised Analog LPF Transfer function H(S)=');
hs_Normalised = analpf(N,'butt',[0,0],1);
disp(hs_Normalised);

// Analog LPF Transfer function
disp('Analog LPF Transfer function H(S)=');
hs = analpf(N,'butt',[0,0],omegac);
disp(hs);

s = poly(0,'s');
hpf_s = horner(hs, omegac/s);   // substitute s → omegac/s
disp('Analog HPF Transfer function H(S)=');
disp(hpf_s);
// Bilinear Transformation to Digital
z = poly(0,'z'); // Defining variable z
Hz = horner(hpf_s,(2/T)*((z-1)/(z+1))); 
disp('Digital HPF Transfer function H(Z)=');
disp(Hz);

// Frequency Response
HW = frmag(Hz,512);
w = 0:%pi/511:%pi;
plot(w/%pi, abs(HW));

xlabel(' Normalized Digital Frequency w');
ylabel('Magnitude');
title(' Frequency Response of Butterworth IIR HPF');
```



## OUTPUT (LPF) :
<img width="1919" height="1079" alt="image" src="https://github.com/user-attachments/assets/767cf8b6-78ab-4450-8569-254d24d5d04b" />
<img width="1919" height="1079" alt="image" src="https://github.com/user-attachments/assets/653a4a2f-e7a4-49c7-9790-150070400f25" />




## OUTPUT (HPF) : 
<img width="1919" height="1079" alt="image" src="https://github.com/user-attachments/assets/09839817-aca2-4310-9817-909fe2c69839" />
<img width="1919" height="1079" alt="image" src="https://github.com/user-attachments/assets/a296b582-3e99-4515-aba6-8782d4e6d2a4" />


## RESULT:
HENCE THE OUTPUT WAS VERIFIED FOR THE IIR BUTTERWORTH IIR FILTER.
