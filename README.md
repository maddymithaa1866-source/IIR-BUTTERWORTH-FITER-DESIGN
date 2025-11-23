# EXP 3 : IIR-BUTTERWORTH-FITER-DESIGN

## AIM: 

 To design an IIR Butterworth filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF):
'''
clc;
clear;
close;

// ---- Given specifications ----
wp = 0.3 * %pi;     // Passband frequency (radians)
ws = 0.6 * %pi;     // Stopband frequency (radians)
alphap = 3;         // Passband attenuation (dB)
alphas = 40;        // Stopband attenuation (dB)
T = 1;              // Sampling time

// ---- Pre-warping (for Bilinear Transformation) ----
omegap = (2 / T) * tan(wp / 2);
omegas = (2 / T) * tan(ws / 2);

disp(omegap, "Prewarped Passband Frequency (omegap) =");
disp(omegas, "Prewarped Stopband Frequency (omegas) =");

// ---- Filter Order Calculation ----
N = log10(((10^(0.1 * alphas)) - 1) / ((10^(0.1 * alphap)) - 1)) / (2 * log10(omegas / omegap));
N = ceil(N);  // Round off to next integer
disp(N, "Filter Order (N) =");

// ---- Cutoff Frequency ----
omegac = omegap / (((10^(0.1 * alphap)) - 1)^(1 / (2 * N)));
disp(omegac, "Cutoff Frequency (omegac) =");

// ---- Analog Butterworth LPF ----
disp("Analog Butterworth LPF Transfer Function H(s):");
hs = analpf(N, 'butt', [0, 0], omegac);
disp(hs);

// ---- Bilinear Transformation to Digital Filter ----
z = poly(0, 'z');
Hz = horner(hs, (2 / T) * ((1 - z^-1) / (1 + z^-1)));  // Bilinear transform
disp("Digital LPF Transfer Function H(z):");
disp(Hz);

// ---- Frequency Response ----
[Hf, fr] = frmag(Hz, 512);

plot(fr / %pi, abs(Hf));
xlabel('Normalized Digital Frequency (\omega / \pi)');
ylabel('Magnitude');
title('Frequency Response of Butterworth IIR Low Pass Filter');
xgrid();
'''




## PROGRAM (HPF): 




## OUTPUT (LPF) : <img width="1920" height="1200" alt="image" src="https://github.com/user-attachments/assets/c558c617-3952-45df-b555-5c064792147d" />


## OUTPUT (HPF) : ![WhatsApp Image 2025-10-30 at 11 32 13_d80cd448](https://github.com/user-attachments/assets/37a06fa3-32d3-4fd8-bdae-b30fcc3e70af)



## RESULT: HENCE THE OUTPUT WAS VERIFIED FOR THE IIR BUTTERWORTH IIR FILTER.
