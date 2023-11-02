format long g
close all
clear all

%% Data Import Functions

[FBG2, Piezo2, freq, Crest] = sourcedata1;          % Calls function that imports data set 1 with specified delimiters and removes Nans

[FBG1, Piezo1, FBGFFT, PiezoFFT, Cresttotroph] = sourcedata2;     % Calls function to import data set 2 with delimiters to organise data

%% New Data

% Random Set
a = min(FBG1);          % Creating boundaries for randomised wavelength set
b = max(FBG1);

RandomSet = (b-a)*rand(500, 1) + a;         % Producing random wavelengths within boundary conditions

% Sinusoidal Pattern
Fs = 500;                   % samples per second
dt = 1/Fs;                  % seconds per sample
StopTime = 1;               % seconds
t = (0:dt:StopTime-dt)';    % seconds

Fc = 60;                    % Choose frequency of sine wave in hertz
SinData = abs(RandomSet .* sin(pi*Fc*t));           % Overlaying the original random data set with a Sine wave of chosen frequency

% Sensitivity
NewSens = [];

for i = 20:20:120               % Loop of fast Fourier transform (FFT) at different sine wave frequencies to produce sensitivity plot
SensFC = i;
SensData = abs(RandomSet .* sin(pi*SensFC*t));
d = SensData;
y = d - mean(d);               
L = length(d);                                
NFFT = 2^nextpow2(L);  
Y = fft(y,NFFT)/L;          
SensFFTmodulus = 2*abs(Y(1:NFFT/2+1));
NewSens(end+1) = [max(SensFFTmodulus)];
end

SensFreq = [20, 40, 60, 80, 100, 120];      % Creating vector for frequencies that the sensitivity is plotted against

%% FFT Computation

Fs = 500;  T = 1/Fs;  clf;      % Sampling frequency 500hz
a=19000;                       % Set the starting data (in this case, it is from 3000, and the data before 3000 is not used)
x1 = a:1:(a+499);               % Data length is 500. x1 is a one-row, 500-column matrix (here, the values are from 3000 to 3499)
Data = [FBG1(x1,1), RandomSet, SinData, Piezo1(x1,1)];  % Creates matrix of all data to be put through FFT 
y = Data - mean(Data);               % Getting rid of the average to filter the 0Hz data
L = length(Data);                 % Length of c
t = (0:L-1)*T;               % Set time interval 
NFFT = 2^nextpow2(L);   % The closest power of 2 larger than L
Y = fft(y,NFFT)/L;          %FFT
f = Fs/2*linspace(0,1,NFFT/2+1);   % Frequency from 0 to 250Hz with NFFT/2+1 intervals 
OrigFFTmod = 2*abs(Y(1:NFFT/2+1, 1));  % modulus of Y, times 2 to get both the positive and negative modulus, only use Y from 1 to NFFT/2+1 as Y is symmetrical excluding the first point. 
PiezoFFTmod = 2*abs(Y(1:NFFT/2+1, 4)); % 
RandFFTmod = 2*abs(Y(1:NFFT/2+1, 2));
SinFFTmod = 2*abs(Y(1:NFFT/2+1, 3));

%% Original Plots

% FBG Resonant frequency plot
figure;
subplot(2,2,1)
plot(f, OrigFFTmod);  
xlabel('Frequency (Hz)'); 
ylabel('Resonant Wavelength Shift (pm)'); % Resonant Wavelength Shift (pm) or Acceleration(g)
title('FBG Accelerometer');
grid minor;

%Piezo frequency plot
subplot(2,2,2)
plot(f, PiezoFFTmod);  
xlabel('Frequency (Hz)'); 
ylabel('Resonant Wavelength Shift (pm)');
title('Piezo Accelerometer');
grid minor;

% Resonance of FBG
subplot(2,2,3)
plot(freq, Crest, 'b', freq, Crest, 'bo');
xlabel('Frequency (Hz)'); 
ylabel('Resonant Wavelength Shift (pm)');
legend(['Resonant Frequency of ' [num2str(freq(9))] ' Hz']);
title('Resonance Profile of FBG');
grid minor;

% Sensitivity plot
subplot(2,2,4)
[P, S] = polyfit(PiezoFFT, Cresttotroph, 1);         % Linear fitting to give sensitivity value
fit = polyval(P, PiezoFFT);
errOrig = 1 - (S.normr/norm(Cresttotroph - mean(Cresttotroph)))^2;
plot(PiezoFFT, fit, 'b-', PiezoFFT, Cresttotroph, 'rx');
title('Original Sensitivity Plot');
legend(['Sensitivity = ' [num2str(P(1))] ' pm/g'], ['R^2 = ' [num2str(errOrig)]]);     % Showing sensitivity slop value and R^2 on plot
xlabel('Acceleration from Piezo (g)');
ylabel('Wavelength Change (pm)');
grid minor;

%% New Plots

% Created data plot
figure;
subplot(2,2,1)
plot(f, RandFFTmod);  
xlabel('Frequency (Hz)'); 
ylabel('Resonant Wavelength Shift (pm)'); 
title('Randomised Wavelength');
grid minor;

subplot(2,2,2)
plot(f, SinFFTmod);  
xlabel('Frequency (Hz)'); 
ylabel('Resonant Wavelength Shift (pm)');
title('Random with Sine Wave');
grid minor;

subplot(2,2,3)
[Q, H] = polyfit(SensFreq, NewSens, 1);          %Linear fit
fit = polyval(Q, SensFreq);
errNew = 1 - (H.normr/norm(NewSens - mean(NewSens)))^2;
plot(SensFreq, fit, 'b-', SensFreq, NewSens, 'rx');
legend(['Sensitivity = ' [num2str(Q(1))] ' pm/hz'], ['R^2 = ' [num2str(errNew)]]);
title('New Random Sensitivity');
xlabel('Driving Frequency (Hz)');
ylabel('Wavelength Change (pm)');
grid minor;

%% Functions

function [FBG2, Piezo2, freq, Crest] = sourcedata1

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [8, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Piezo1", "FBG1", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "freq", "Var10", "Var11", "Var12", "Crest"];
opts.SelectedVariableNames = ["Piezo1", "FBG1", "freq", "Crest"];
opts.VariableTypes = ["double", "double", "string", "string", "string", "string", "string", "string", "double", "string", "string", "string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var10", "Var11", "Var12"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var10", "Var11", "Var12"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable("C:\Users\lachl\Desktop\COMP7860\Data Sets\Frequency response at different frequencies.csv", opts);

% Convert to output type
Piezo2 = tbl.Piezo1;
FBG2 = tbl.FBG1;
freq = tbl.freq;
freq = freq(~isnan(freq));
Crest = tbl.Crest;
Crest = Crest(~isnan(Crest));

% Clear temporary variables
clear opts tbl

end

function [FBG1, Piezo1, FBGFFT, PiezoFFT, Cresttotroph] = sourcedata2

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Piezo1", "FBG1", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "PiezoFFT", "FBGFFT", "Cresttotroph"];
opts.SelectedVariableNames = ["Piezo1", "FBG1", "PiezoFFT", "FBGFFT", "Cresttotroph"];
opts.VariableTypes = ["double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable("C:\Users\lachl\Desktop\COMP7860\Data Sets\Sensitivity at 5Hz.csv", opts);

% Convert to output type
Piezo1 = tbl.Piezo1;
FBG1 = tbl.FBG1;
PiezoFFT = tbl.PiezoFFT;
PiezoFFT = PiezoFFT(~isnan(PiezoFFT));
FBGFFT = tbl.FBGFFT;
FBGFFT = FBGFFT(~isnan(FBGFFT));
Cresttotroph = tbl.Cresttotroph;
Cresttotroph = Cresttotroph(~isnan(Cresttotroph));

% Clear temporary variables
clear opts tbl

end
