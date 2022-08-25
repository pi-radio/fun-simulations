clear;
nFFT = 1024;
constellation = [1+1j 1-1j -1+1j -1-1j];
tx_fd = zeros(nFFT, 1);
scMin = -nFFT/2 + 100;
scMax = nFFT/2-1 - 100;
for sc = scMin:scMax
    %if (sc ~= 0)
        tx_fd(sc+nFFT/2 + 1) = constellation(randi(4));
    %end
end
tx_fd = fftshift(tx_fd);
tx_td = ifft(tx_fd);

% Add noise
target_snr = 10; % dB
tx_td = awgn(tx_td, target_snr, 'measured');

tx_td = tx_td / max(abs(tx_td)); % Move into the range of -1 to +1 (approximately)

numEffectiveBits = 5;
tx_td = tx_td * 2^(numEffectiveBits-1);
tx_td = complex(double(int16(real(tx_td))) , double((int16(imag(tx_td)))));

rx_td = tx_td;
rx_td = [rx_td(100:nFFT); rx_td(1:99)];
rx_fd = fft(rx_td);

corr_fd = tx_fd .* conj(rx_fd);
corr_td = ifft(corr_fd);
figure(1); clf;
plot(mag2db(abs(corr_td)));
%ylim([-50 10]);
grid on;

noise = ((sum(abs(corr_td)) - abs(corr_td(100))))/nFFT;
noise_dB = mag2db(noise);
SNR = max(mag2db(abs(corr_td))) - noise_dB

figure(2); clf;
plot(real(tx_td), 'r'); hold on;
plot(imag(tx_td), 'g'); hold on;
%plot(abs(tx_td), 'b'); hold off;