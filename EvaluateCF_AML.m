function CF = EvaluateCF_AML(p,y,NoB)
%p(1) = A; p(2) = B; p(3) = C; p(4) = theta; p(5) = sigma; p(6) =
%INLFFT(1); p(7) = real(INLFFT(2)); p(8) = imag(INLFFT(2)); ...

p = p(:);
NUM_OF_FOURIER_COEFFS = length(p)-5;
INLFFT = zeros(1,2^NoB-1);
INLFFT(1) = p(6);
for k = 1:(NUM_OF_FOURIER_COEFFS-1)/2
    INLFFT(k+1) = p(6+2*k-1) + 1i*p(6+2*k);
end

for k = 1:NUM_OF_FOURIER_COEFFS
    INLFFT(end-k+1) = conj(INLFFT(k+1));
end

INL = ifft(INLFFT);

if (max(imag(INL)) > 1e-14)
    warning('Calculated INL is incorrect (the imaginary part is larger than quantization noise)');
end

INL = real(INL);    %Discarding imagiary parts owing to quantization noise

[~,CF] = EvaluateCF(y,p(1:5),NoB,INL);

end