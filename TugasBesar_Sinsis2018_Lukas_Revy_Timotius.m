%%1)   Load the 'SoundSample.wav' file in MATLAB
% [sample_data, sample_rate] = audioread('SoundSample.m4a');
[sample_data, sample_rate] = audioread('SoundSampleWoman.aac');

% a.  Plot the signal in time and the amplitude of its frequency components
% using the FFT.
sample_period = 1/sample_rate; %T=1/f = 2.0833e-051
t = (0:sample_period:(length(sample_data)-1)/sample_rate);
figure(2);
subplot(2,2,1)
plot(t,sample_data)
title('Time Domain Representation - Unfiltered Sound')
xlabel('Time (seconds)')
ylabel('Amplitude')
xlim([0 t(end)])


m = length(sample_data); %Panjang sinyal dari sample_data = 165887
%POW
% Mencari panjang sample asli.
% n adalah panjang vektor hasil perhitungan FFT. 
% idealnya adalah bilangan dari angka 2 yang dipangkatkan, 
% dan lebih panjang dari sinyal informasi (m). Sehingga digunakanlah fungsi nextpow2.

n = pow2(nextpow2(m)); %nextpow2(m) = 17 >> pow2(17) >> 2^17 = 131072
% Mengubah panjang sinyal, jadi jumlah sample = pow2
% Membuat perhitungan FFT menjadi lebih cepat
% Terutama untuk ukuran sample data dengan faktor prima yang lebih besar

% Untuk mendapatkan normalisasi magnitud sinyal, 
% nilai hasil FFT harus dibagi dengan L (panjang data)
y = fft(sample_data, n);
f = (0:n-1)*(sample_rate/n); %Membuat sumbu x dalam domain frekuensi (hasil ft)
amplitude = abs(y)/n; %Meng-absolute kan nilai dari y (fft dari sample_data)
subplot(2,2,2)
%semilogy(f(1:floor(n/2)),amplitude(1:floor(n/2)))
plot(f(1:floor(n/2)),amplitude(1:floor(n/2)))
title('Frequency Domain Representation - Unfiltered Sound')
xlabel('Frequency')
% ylabel('dB')
ylabel('Amplitude')
% b.  Listen to the audio file.
% sound(sample_data, sample_rate)

% max amplitude sinyal lama/max amplitude sinyal baru * array dari sinyal baru
% dalam domain waktu

%2)Filter the audio sample data to remove noise from the signal.
order = 2; %Order dalam Butterworth Filter, Semakin tinggi ordernya maka semakin bagus tapi semakin kompleks
% order, cutoff frequency (Wn)= 1200 and return as vector or scalar, lowpass / highpass
%[b,a] = butter(order,1200/(sample_rate/2),'low'); %Filter Butterworth Untuk menghitung filternya menggunakan lowpass

%Wn=[4000 8000]/(sample_rate/2);
Wn=[2.5e6 29e6]/450e6;
[b,a] = butter(order,Wn,'bandpass'); %Filter Butterworth Untuk menghitung filternya menggunakan lowpass

% beginFreq = 400 / (sample_rate/2);
% endFreq = 20000 / (sample_rate/2);
% [b,a] = butter(order, [beginFreq, endFreq], 'bandpass');

filtered_sound = filter(b,a,sample_data); %Untuk memasukkan hasil filter ke sinyal
sound(filtered_sound, sample_rate)  
t1 = (0:sample_period:(length(filtered_sound)-1)/sample_rate);
subplot(2,2,3)
plot(t1,filtered_sound)
title('Timedomain Representation - Filtered Sound')
xlabel('Time (seconds)')
ylabel('Amplitude')

xlim([0 t1(end)])
m1 = length(sample_data); % Original sample length.
n1 = pow2(nextpow2(m1)); % Transforming the length so that the number of 
% samples is a power of 2. This can make the transform computation 
% significantly faster,particularly for sample sizes with large prime 
% factors.
y1 = fft(filtered_sound, n1);
f = (0:n1-1)*(sample_rate/n1);
amplitude = abs(y1)/n1;
subplot(2,2,4)
%semilogy(f(1:floor(n1/2)),amplitude(1:floor(n1/2)))
plot(f(1:floor(n1/2)),amplitude(1:floor(n1/2)))
title('Frequency Domain Representation - Filtered Sound')
xlabel('Frequency')
%ylabel('dB')
ylabel('Amplitude')

L=1024;
f_d=(L*sample_period);
figure(3);
numberF=11;

f_size=round(f_d*sample_rate); %1024
n_f=floor(m/f_size); %161
temp=0;
q=0;
for i = 1:n_f
    frames(i,:) = sample_data(temp+1:temp+f_size);
    temp=temp+f_size;
    tList(i,:)=q:sample_period:(L*sample_period)-sample_period+q; %Sumbu Horizontal
    q=q+0.025;
end
subplot(3,1,1);
plot(tList(numberF,:),frames(numberF,:)); xlabel('Time (seconds)'); ylabel('Amplitude');
title('Time Domain Representation - Unfiltered Sound (USING FRAMING)')
hold on;

f_size=round(f_d*sample_rate);
n_f=floor(m/f_size);
temp=0;
q=0;
for i = 1:n_f
    frames(i,:) = filtered_sound(temp+1:temp+f_size);
    temp=temp+f_size;
    tList(i,:)=q:sample_period:(L*sample_period)-sample_period+q;
    q=q+0.025;
end
subplot(3,1,2);
plot(tList(numberF,:),frames(numberF,:)); xlabel('Time (seconds)'); ylabel('Amplitude');
title('Time Domain Representation - Filtered Sound (FRAME)')
hold on;

for i = 1:n_f
    aa=reshape(frames(i,:),1024,1);
    %%FFT
    % L=length(y);

    y1_hamming=aa.*hamming(L);
    NFFT = 2^nextpow2(L); 
    Y = fft(y1_hamming,NFFT)/L; 
    f = sample_rate/2*linspace(0,1,NFFT/2); 
    y_hammingList(i,:)=y1_hamming;
    NFFTList(i,:)=NFFT;
    fList(i,:)=f;
end

% Plot single-sided amplitude spectrum.
subplot(3,1,3);
plot(tList(numberF,:),y_hammingList(numberF,:)); xlabel('Time (seconds)'); ylabel('Amplitude');
title('Time Domain Representation - Filtered Sound (USING HAMMING)')

figure(4);
% y2 = fft(frames(numberF,:));
% subplot(3,1,1);
% plot(tList(numberF,:),real(y2),tList(numberF,:),imag(y2)); xlabel('Time (seconds)'); ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(tList(numberF,:),y2); xlabel('Time (seconds)'); ylabel('Amplitude');

% x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
% y = x + 2*randn(size(t));
y1 = fft(filtered_sound, n1);

% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(numberF,NFFT)/L;
% f = sample_rate/2*linspace(0,1,NFFT/2+1);
% % Plot single-sided amplitude spectrum.
% % subplot(3,1,3);
% plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')

Y = fft(y_hammingList(numberF,:),NFFTList(numberF,:))/L;
Yf = 2*abs(Y(1:NFFTList(numberF,:)/2));
fY = fList(numberF,:);
plot(fY, Yf);
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('Amplitude')