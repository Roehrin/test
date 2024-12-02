// brainVisionReader.js
export async function readBrainvisionEEGData(url) {
    const sampleSize = 4; // 4 bytes for a 32-bit float
    const numberOfChannels = 204; // Equivalent to hdr.NumberOfChannels in MATLAB
    const calib = 1; // Calibration factor

    try {
        // Fetch the file from the URL
        const response = await fetch(url);
        if (!response.ok) {
            throw new Error(`HTTP error! Status: ${response.status}`);
        }

        // Get the file as an ArrayBuffer
        const arrayBuffer = await response.arrayBuffer();
        const nSamples = arrayBuffer.byteLength / sampleSize/numberOfChannels;
		
		// create time vector
		let timeVec = Array.from({length: nSamples}, (_, i) => (i + 1)/Fs);
        // Create a DataView for reading binary data
        const dataView = new DataView(arrayBuffer);
		
		// write the data so they can be directly used by plot.ly
		// init data by pre-allocating ys as vector of zeros
		let data = [];
		// important because data are multiplex
		for (let i = 0; i < numberOfChannels; i++) {
			let dummyLine = {x:timeVec, y:Array(nSamples).fill(0), mode: "lines", line: {color: 'rgb(0,0,0)'}}; 
			data.push(dummyLine);
		}
		for (let sampId = 0; sampId < nSamples; sampId++) {
			for (let elId = 0; elId < numberOfChannels; elId++) {
				const offset = (sampId * numberOfChannels + elId) * sampleSize;
                const value = dataView.getFloat32(offset, true); // true for little-endian
				data[elId].y[sampId] = value;
			}
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