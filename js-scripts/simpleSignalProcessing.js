// simpleSignalProcessing.js

function computeDFT(signal, sampleRate) {
	let N = signal.length;
	let real = new Array(N).fill(0);
	let imag = new Array(N).fill(0);
	
	// Handle odd N for frequency array indexing
	let halfN = Math.floor(N / 2);
	let frequencies = Array.from({ length: halfN + 1 }, (_, i) => (i * sampleRate) / N);

	for (let k = 0; k < N; k++) {
		for (let n = 0; n < N; n++) {
			let angle = (-2 * Math.PI * k * n) / N;
			real[k] += signal[n] * Math.cos(angle);
			imag[k] += signal[n] * Math.sin(angle);
		}
	}

	let psd = real.map((re, i) => 2/N*Math.sqrt(re * re + imag[i] * imag[i]));
	let phase = real.map((re, i) => Math.atan2(imag[i], re));

	return { frequencies, psd, phase, real, imag };
}

// Compute Inverse DFT (IDFT)
function inverseDFT(real, imag) {
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

function analyticSignal(signal, sampleRate) {
	let dftResults = computeDFT(signal, sampleRate);
	let N = signal.length;

	let halfN = Math.floor(N / 2); // Integer division
	let isOdd = N % 2 !== 0;

	// Apply Hilbert Transform in frequency domain
	for (let k = 1; k < halfN; k++) {
		dftResults.real[k] *= 2;
		dftResults.imag[k] *= 2;
	}

	if (!isOdd) {
		dftResults.real[halfN] *= 2;  // Double the Nyquist frequency for even N
		dftResults.imag[halfN] *= 2;
	}

	for (let k = halfN + 1; k < N; k++) {
		dftResults.real[k] = 0;
		dftResults.imag[k] = 0;
	}

	let hilbertImag = inverseFFT(new Array(N).fill(0), dftResults.imag); // Compute Hilbert transform
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
	let correlation = (denominator === 0) ? 0 : (numerator / denominator);
	
	// Compute Linear Regression inside the same function
    let slope = numerator / (n * sumX2 - sumX ** 2);
    let intercept = (sumY - slope * sumX) / n;

    return { correlation, slope, intercept };
}

function generateGaussianNoise(noiseLvl, length){
	const noise = [];
	for (let i = 0; i < length; i++) {
		const u1 = Math.random();
		const u2 = Math.random();
		const z0 = noiseLvl*Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2); // Generate standard normal variable
		noise.push(z0); // Scale and shift to desired mean and standard deviation
	}
	return noise
}

// Export functions for use in another script
export { computeDFT, inverseDFT, analyticSignal, pearsonCorrelation, generateGaussianNoise};
