% Generate the TX waveform in frequency and time domain
clear; clc;
nFFT = 1024;
constellation = [1+1j 1-1j -1+1j -1-1j];
tx_fd = zeros(nFFT, 1);
for sc = -nFFT/2 : nFFT/2-1
    tx_fd(sc + nFFT/2 + 1) = constellation(randi(4));
end
tx_fd = fftshift(tx_fd);
tx_td = ifft(tx_fd);

% Create a "packet" with 10 such sumbols. Ignore the cyclic prefix for now
nSym = 40;
tx_td_app = zeros(nFFT*nSym,1);
for iSym = 1:nSym
    tx_td_app((iSym-1)*nFFT + 1 : iSym*nFFT) = tx_td;
end
tx_td = tx_td_app;
clear tx_td_app;

%%%%%%%%%% Start Channel Emulator %%%%%%%%%%

% A: "Do nothing" channel. Do not comment this out.
rx_td = tx_td;

% B: Timing Offset. Comment this out for no timing offset.
k = 3;
n = size(rx_td, 1);
rx_td = [rx_td(k+1:n); rx_td(1:k)];
clear n k;

% C: Residual Frequency Offset. Comment this out for no frequency offset.
    % We assume a sampling frequency of 1 GSps complex, for a 1 GHz bandwidth.
    % The sample duration is therefore 1ns.
fs = 20e3; % Frequency offset in Hertz
n = 1e9 / fs; % In 'n' samples, the channel will rotate a full circle
p = 2*pi/n; % The channel will rotate by 'p' radian every sample
clear n;
n = size(rx_td, 1); % Now, 'n' is the number of samples in the waveform
phasor = exp(1j*linspace(0, p*(n-1), n));
phasor = phasor';
rx_td = rx_td .* phasor;
clear fs n phasor;

% D: Add Noise
target_snr = 18; % in dB
rx_td = awgn(rx_td, target_snr, 'measured');

%%%%%%%%%%  End Channel Emulator %%%%%%%%%%


% Calculate and plot the 2D Channel
h_fd = zeros(nFFT, nSym);
rx_fd = zeros(nFFT, nSym);
for iSym = 1:nSym
    sym_td = rx_td((iSym-1)*nFFT + 1 : iSym*nFFT);
    sym_fd = (fft(sym_td));
    rx_fd(:, iSym) = sym_fd;
    h_fd(:, iSym) = sym_fd ./ tx_fd;
end
figure(1); clf;

subplot(3,2,1);
s = surf(abs(h_fd));
%caxis([0 2]);
s.EdgeColor = 'none'; view(2); xlim([1 nSym]); ylim([1 nFFT]);
colorbar;
xlabel('Symbol Index');
ylabel('Subcarrier Index');
title('Abs(H)');

subplot(3,2,2);
s = surf(angle(h_fd));
%caxis([-pi pi]);
s.EdgeColor = 'none'; view(2); xlim([1 nSym]); ylim([1 nFFT]); colorbar;
xlabel('Symbol Index');
ylabel('Subcarrier Index');
title('Angle(H)');

subplot(3,2,3);
s = surf(real(h_fd));
%caxis([-2 2]);
s.EdgeColor = 'none'; view(2); xlim([1 nSym]); ylim([1 nFFT]); colorbar;
xlabel('Symbol Index');
ylabel('Subcarrier Index');
title('Real(H)');

subplot(3,2,4);
s = surf(imag(h_fd));
%caxis([-2 2]);
s.EdgeColor = 'none'; view(2); xlim([1 nSym]); ylim([1 nFFT]); colorbar;
xlabel('Symbol Index');
ylabel('Subcarrier Index');
title('Imag(H)');

% Plot the unequalized frequency-domain samples
subplot(3,2,5);
for iSym = 1:nSym
    scatter(real(rx_fd(:, iSym)), imag(rx_fd(:, iSym)), 'b.'); hold on;
end
xlim([-2 2]); ylim([-2 2]);
title('Unequalized');

% Now, plot the received constellation. For this, we will assume that the
% pilots are located every 'k' subcarriers. The remaining subcarriers
% carry the unknown user data that we're trying to demodulate.
k = 4;
subplot(3,2,6);
eq_sym = zeros(nFFT, nSym);
for iSym = 1:nSym
    for sc = 1:nFFT
        if mod(sc-1, k) == 0
            % Pilot
            h = h_fd(sc, iSym);
        else
            % Data
            eq_sym(sc, iSym) = rx_fd(sc, iSym) / h;
        end
    end
    scatter(real(eq_sym(:, iSym)), imag(eq_sym(:, iSym)), '.b'); hold on;
end
xlim([-2 2]); ylim([-2 2]);
title('Equalized');