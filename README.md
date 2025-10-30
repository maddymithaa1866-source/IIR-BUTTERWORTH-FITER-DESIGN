# EXP 3 : IIR-BUTTERWORTH-FITER-DESIGN

## AIM: 

 To design an IIR Butterworth filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF): clc;
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




## PROGRAM (HPF): // ------------------------------
//   High-Pass Filter in Scilab
// ------------------------------
clc;
clear;

// Filter Design Parameters
fs    = 1000;     // Sampling Frequency
fc    = 200;      // Cutoff Frequency
order = 4;        // Filter Order

// Normalized Cutoff Frequency (fc/fs)
Wn = fc / fs;

// Design High-Pass Butterworth Filter
hz = iir(order, 'hp', 'butt', Wn, [0 0]);
b  = coeff(hz, 'num');   // Numerator coefficients
a  = coeff(hz, 'den');   // Denominator coefficients

// Display Filter Details
disp("----- FILTER DESIGN -----");
disp("Sampling Frequency (fs): " + string(fs));
disp("Cutoff Frequency (fc): " + string(fc));
disp("Filter Order: " + string(order));

disp("Numerator Coefficients (b):");
disp(b);
disp("Denominator Coefficients (a):");
disp(a);

// ------------------------------
//   Test Signal Input
// ------------------------------
t  = 0:1/fs:0.2;
f1 = 50;     // Low frequency (should be removed)
f2 = 350;    // High frequency (should pass)
x  = 0.5*sin(2*%pi*f1*t) + sin(2*%pi*f2*t);

// Filter the Signal
y = filter(b, a, x);

// ------------------------------
//   Time Domain Plots
// ------------------------------
scf(0);
subplot(2,1,1);
plot(t, x);
xtitle('Input Signal', 'Time (s)', 'Amplitude');

subplot(2,1,2);
plot(t, y);
xtitle('Filtered Output (High-Pass)', 'Time (s)', 'Amplitude');

// ------------------------------
//   Frequency Response
// ------------------------------
[hm, fr] = frmag(hz, 512);   // Magnitude response

scf(1);
plot(fr * fs, hm);
xtitle('Magnitude Response', 'Frequency (Hz)', 'Magnitude');




## OUTPUT (LPF) : <img width="1920" height="1200" alt="image" src="https://github.com/user-attachments/assets/c558c617-3952-45df-b555-5c064792147d" />


## OUTPUT (HPF) : ![WhatsApp Image 2025-10-30 at 11 29 12_be715f6d](https://github.com/user-attachments/assets/fb45d99b-1f8f-4d2c-adfb-30af7fd18718)


## RESULT: HENCE THE OUTPUT WAS VERIFIED FOR THE IIR BUTTERWORTH IIR FILTER.
