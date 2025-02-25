// simpleSignalProcessing.js

// Compute FFT using Discrete Fourier Transform (DFT)
function computeFFT(signal, sampleRate) {
    let N = signal.length;
    let real = new Array(N).fill(0);
    let imag = new Array(N).fill(0);
    let frequencies = Array.from({ length: N / 2 }, (_, i) => (i * sampleRate) / N);

    for (let k = 0; k < N; k++) {
        for (let n = 0; n < N; n++) {
            let angle = (-2 * Math.PI * k * n) / N;
            real[k] += signal[n] * Math.cos(angle);
            imag[k] += signal[n] * Math.sin(angle);
        }
    }

    let magnitudes = real.map((re, i) => Math.sqrt(re * re + imag[i] * imag[i]));

    return { frequencies, magnitudes, real, imag };
}

// Compute Inverse FFT (IFT)
function inverseFFT(real, imag) {
    let N = real.length;
    let signal = new Array(N).fill(0);

    for (let n = 0; n < N; n++) {
        for (let k = 0; k < N; k++) {
            let angle = (2 * Math.PI * k * n) / N;
            signal[n] += real[k] * Math.cos(angle) - imag[k] * Math.sin(angle);
        }
        signal[n] /= N; // Normalize
    }

    return signal;
}

// Compute Hilbert Transform using FFT and IFT
function hilbertTransform(signal, sampleRate) {
    let [, , real, imag] = computeFFT(signal, sampleRate);
    let N = real.length;

    // Zero out negative frequencies and double positive ones (analytic signal)
    let hilbertImag = new Array(N).fill(0);
    for (let k = 1; k < N / 2; k++) {
        hilbertImag[k] = imag[k] * 2;
    }

    return inverseFFT(new Array(N).fill(0), hilbertImag); // Return only the imaginary part
}

function analyticSignal(signal, sampleRate) {
    let [, , real, imag] = computeFFT(signal, sampleRate);
    let N = real.length;

    let halfN = Math.floor(N / 2); // Integer division
    let isOdd = N % 2 !== 0;

    // Apply Hilbert Transform in frequency domain
    for (let k = 1; k < halfN; k++) {
        real[k] *= 2;
        imag[k] *= 2;
    }

    if (!isOdd) {
        real[halfN] *= 2;  // Double the Nyquist frequency for even N
        imag[halfN] *= 2;
    }

    for (let k = halfN + 1; k < N; k++) {
        real[k] = 0;
        imag[k] = 0;
    }

    let hilbertImag = inverseFFT(new Array(N).fill(0), imag); // Compute Hilbert transform
    return signal.map((val, i) => [val, hilbertImag[i]]); // (Real part, Imaginary part)
}

function pearsonCorrelation(x, y) {
	let n = x.length;
	let sumX = x.reduce((a, b) => a + b, 0);
	let sumY = y.reduce((a, b) => a + b, 0);
	let sumXY = x.map((xi, i) => xi * y[i]).reduce((a, b) => a + b, 0);
	let sumX2 = x.map(xi => xi * xi).reduce((a, b) => a + b, 0);
	let sumY2 = y.map(yi => yi * yi).reduce((a, b) => a + b, 0);

	let numerator = (n * sumXY) - (sumX * sumY);
	let denominator = Math.sqrt((n * sumX2 - sumX ** 2) * (n * sumY2 - sumY ** 2));
	return (denominator === 0) ? 0 : (numerator / denominator);
}

// Export functions for use in another script
export { computeFFT, inverseFFT, hilbertTransform, analyticSignal, pearsonCorrelation};
