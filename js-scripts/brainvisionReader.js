// brainVisionReader.js
export async function readBrainvisionEEGData(url) {
    const sampleSize = 4; // 4 bytes for a 32-bit float
    const numberOfChannels = 204; // Equivalent to hdr.NumberOfChannels in MATLAB
    const calib = 1; // Calibration factor
    let data = [];

    try {
        // Fetch the file from the URL
        const response = await fetch(url);
        if (!response.ok) {
            throw new Error(`HTTP error! Status: ${response.status}`);
        }

        // Get the file as an ArrayBuffer
        const arrayBuffer = await response.arrayBuffer();
        const nSamples = arrayBuffer.byteLength / sampleSize/numberOfChannels;

        // Create a DataView for reading binary data
        const dataView = new DataView(arrayBuffer);

        // Iterate through the data to extract samples
        for (let i = 0; i < nSamples; i++) {
            const sample = [];
            for (let j = 0; j < numberOfChannels; j++) {
                const offset = (i * numberOfChannels + j) * sampleSize;
                const value = dataView.getFloat32(offset, true); // true for little-endian
                sample.push(value * calib); // Apply calibration
            }
            data.push(sample);
        }

        console.log("EEG Data:", data); // The resulting 2D array
		return data
    } catch (error) {
        console.error("Error reading EEG file:", error);
		return null; // Return null in case of an error
    }
}

export async function readBrainvisionEEGHeader(url) {
}