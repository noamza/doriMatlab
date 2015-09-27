function Clean = WienerFilter(Signal,Noisy,Sigma)
%   Wiener Filter for denoising signals.
%
%   INPUT:
%   Signal - Original signal
%   Noisy - Noisy signal
%   Sigma - Noise level
%
%   OUTPUT:
%   Clean_Signal - Denoised signal

    n = length(Signal);

    Signal_FFT = fft(Signal);
    Noisy_FFT = fft(Noisy);
    
    P_Sig = abs(Signal_FFT).^2; % Spectral power

    Wie_Fltr = P_Sig./(P_Sig + n*Sigma^2);

    Clean = real(ifft(Noisy_FFT.*Wie_Fltr));
    
end