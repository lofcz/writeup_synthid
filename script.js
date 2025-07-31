import * as rs from "./reedsolomon.js";
import { dctApply, idctApply } from "./dct.js";
import { dwtApply, idwtApply } from "./dwt.js";
import { cyrb128, sfc32 } from "./rand.js";

document.addEventListener('DOMContentLoaded', () => {
    // #region UI Elements and Constants
    const canvas = document.getElementById('canvas');
    const ctx = canvas.getContext('2d', { willReadFrequently: true });
    const DELIMITER = "####";

    const messageInput = document.getElementById('message-input');
    const dropZoneEncode = document.getElementById('drop-zone-encode');
    const fileInputEncode = document.getElementById('file-input-encode');
    const previewEncode = document.getElementById('preview-encode');
    const encodeButton = document.getElementById('encode-button');
    const resultEncode = document.getElementById('result-encode');
    let encodeFile = null;

    const dropZoneDecode = document.getElementById('drop-zone-decode');
    const fileInputDecode = document.getElementById('file-input-decode');
    const previewDecode = document.getElementById('preview-decode');
    const resultDecode = document.getElementById('result-decode');

    const exampleImagesContainer = document.querySelector('#example-images .flex');

    const picaInstance = pica();
    const exampleImageURLs = ['terry.jpg', 'lion.jpg', 'lena.jpg'];
    const exampleImageBlobs = {};
    const exampleMessages = {
        'terry.jpg': 'Hello Terry!',
        'lion.jpg': 'Hello lion!',
        'lena.jpg': 'Hello Lena!'
    };

    const gain = 5.0;
    const seed = 'secret';
    const mainSize = 256;
    const N = 4; // DCT block size
    // #endregion

    // #region Redundancy Config
    const RS_CONFIG = {
        messageLength: 32, errorCorrectionLength: 16,
        dataLength: 16,
        encoder: new rs.Encoder(rs.GF),
        decoder: new rs.Decoder(rs.GF)
    };
    const CHANNELS = 2; // Embed the message twice for redundancy
    const TOTAL_BLOCKS = Math.pow(mainSize / 2 / N, 2); // 1024
    const BLOCKS_PER_CHANNEL = TOTAL_BLOCKS / CHANNELS; // 512
    const BITS_PER_CHANNEL = BLOCKS_PER_CHANNEL; // 512
    const BYTES_PER_CHANNEL = BITS_PER_CHANNEL / 8; // 64
    const NUM_CHUNKS = BYTES_PER_CHANNEL / RS_CONFIG.messageLength; // 2
    RS_CONFIG.contentLength = NUM_CHUNKS * RS_CONFIG.dataLength; // 32
    messageInput.maxLength = RS_CONFIG.contentLength;
    // #endregion

    // #region Image Pre-processing
    async function processExampleImages() {
        for (const url of exampleImageURLs) {
            const img = new Image();
            img.crossOrigin = 'anonymous'; 
            
            const response = await fetch(url);
            const blob = await response.blob();
            img.src = URL.createObjectURL(blob);

            await new Promise(resolve => { img.onload = resolve; });

            const offScreenCanvas512 = document.createElement('canvas');
            const aspectRatio = img.width / img.height;
            offScreenCanvas512.width = 512;
            offScreenCanvas512.height = 512 / aspectRatio;
            
            await picaInstance.resize(img, offScreenCanvas512);
            const blob512 = await picaInstance.toBlob(offScreenCanvas512, 'image/jpeg', 0.9);
            exampleImageBlobs[url] = { main: blob512 };

            const offScreenCanvas64 = document.createElement('canvas');
            offScreenCanvas64.width = 64;
            offScreenCanvas64.height = 64 / aspectRatio;
            await picaInstance.resize(img, offScreenCanvas64);
            const thumbBlob = await picaInstance.toBlob(offScreenCanvas64, 'image/jpeg', 0.8);
            
            const thumbUrl = URL.createObjectURL(thumbBlob);
            const thumbImg = document.createElement('img');
            thumbImg.src = thumbUrl;
            thumbImg.className = 'w-16 h-16 object-cover rounded-md cursor-pointer hover:opacity-75 transition';
            thumbImg.title = `Use ${url}`;
            thumbImg.addEventListener('click', () => {
                const mainBlobUrl = URL.createObjectURL(exampleImageBlobs[url].main);
                previewEncode.src = mainBlobUrl;
                previewEncode.classList.remove('hidden');
                encodeFile = new File([exampleImageBlobs[url].main], url, { type: 'image/jpeg' });
                if (messageInput.value.trim() === '') {
                    messageInput.value = exampleMessages[url];
                }
            });
            exampleImagesContainer.appendChild(thumbImg);
        }
    }
    // #endregion
    
    // #region Core Stego Logic
    let cbArray, crArray;

    function getPo2Size(width, height) {
        let size = Math.pow(2, Math.round(Math.log(Math.min(width, height)) / Math.log(2)));
        if (size < mainSize) size = mainSize;
        return [size, Math.round(width * size / height), size];
    }
    
    function convertYCbCr(imageData, dataArray) {
        const size = dataArray.length;
        cbArray = new Float32Array(size);
        crArray = new Float32Array(size);
        for (let i = 0; i < size; i++) {
            const R = imageData.data[i*4], G = imageData.data[i*4+1], B = imageData.data[i*4+2];
            dataArray[i] = 0.299 * R + 0.587 * G + 0.114 * B;
            cbArray[i] = 128 - 0.168736 * R - 0.331264 * G + 0.5 * B;
            crArray[i] = 128 + 0.5 * R - 0.418688 * G - 0.081312 * B;
        }
    }

    function convertRGB(imageData, dataArray) {
        for (let i = 0; i < dataArray.length; i++) {
            const Y = dataArray[i];
            const Cb = cbArray[i];
            const Cr = crArray[i];
            imageData.data[i*4+0] = Math.max(0, Math.min(255, Math.round(Y + 1.402 * (Cr - 128))));
            imageData.data[i*4+1] = Math.max(0, Math.min(255, Math.round(Y - 0.344136 * (Cb - 128) - 0.714136 * (Cr - 128))));
            imageData.data[i*4+2] = Math.max(0, Math.min(255, Math.round(Y + 1.772 * (Cb - 128))));
        }
    }

    function recursiveDwt(dataArray, currentSize, desiredSize, binaryString, isEncode) {
        if (currentSize !== desiredSize) {
            const [ LL, LH, HL, HH ] = dwtApply(dataArray, currentSize);
            const result = recursiveDwt(LL, currentSize / 2, desiredSize, binaryString, isEncode);
            if (isEncode) idwtApply(dataArray, currentSize, LL, LH, HL, HH);
            return result;
        } else {
            const [ LL, LH, HL, HH ] = dwtApply(dataArray, currentSize);
            const result = dctEmbed(HH, currentSize / 2, binaryString, isEncode);
            if (isEncode) idwtApply(dataArray, currentSize, LL, LH, HL, HH);
            return result;
        }
    }
    
    function pearsonCorrelation(X, Y) {
        const len = X.length;
        let x_bar = 0, y_bar = 0;
        for (let i = 0; i < len; i++) { x_bar += X[i]; y_bar += Y[i]; }
        x_bar /= len; y_bar /= len;
        let sumDiffProd = 0, sumSquareDiffX = 0, sumSquareDiffY = 0;
        for (let i = 0; i < len; i++) {
            sumDiffProd += (X[i] - x_bar) * (Y[i] - y_bar);
            sumSquareDiffX += Math.pow(X[i] - x_bar, 2);
            sumSquareDiffY += Math.pow(Y[i] - y_bar, 2);
        }
        const denominator = Math.sqrt(sumSquareDiffX * sumSquareDiffY);
        if (denominator === 0) return 0;
        return sumDiffProd / denominator;
    }
    
        function dctEmbed(dataArray, size, binaryString, isEncode) {
        const prngs = [];
        for (let i = 0; i < CHANNELS; i++) {
            const channelSeed = seed + `_ch${i}`;
            prngs.push({
                p0: sfc32(...cyrb128(channelSeed + "_0")),
                p1: sfc32(...cyrb128(channelSeed + "_1")),
            });
        }

        // --- Phase 1: Data Shuffling ---
        const blockCoords = [];
        for (let y = 0; y < size; y += N) {
            for (let x = 0; x < size; x += N) {
                blockCoords.push({ x, y });
            }
        }

        const shuffleSeed = cyrb128(seed + "_shuffle");
        const shufflePrng = sfc32(...shuffleSeed);
        for (let i = blockCoords.length - 1; i > 0; i--) {
            const j = Math.floor(shufflePrng() * (i + 1));
            [blockCoords[i], blockCoords[j]] = [blockCoords[j], blockCoords[i]];
        }
        // --- End of Shuffling ---
    
        const decodedChannels = Array.from({ length: CHANNELS }, () => []);
    
        for (let blockIndex = 0; blockIndex < blockCoords.length; blockIndex++) {
            const { x, y } = blockCoords[blockIndex];
            const channelIndex = blockIndex % CHANNELS;

            const dctSquare = new Float32Array(N * N);
            for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) dctSquare[i*N+j] = dataArray[(y+i)*size + (x+j)];
            dctApply(dctSquare);

            const dataVector = new Float32Array(10);
            let k = 0;
            for (let i = 2; i < N; i++) for (let j = 0; j < i + 1; j++) dataVector[k++] = dctSquare[j * N + (i - j)];
            for (let j = 0; j < N - 1; j++) dataVector[k++] = dctSquare[(1 + j) * N + (N - 1 - j)];
            
            const p0 = prngs[channelIndex].p0;
            const p1 = prngs[channelIndex].p1;

            if (isEncode) {
                // --- Adaptive Gain ---
                let energy = 0;
                for (let i = 1; i < 16; i++) energy += Math.abs(dctSquare[i]);
                const blockGain = gain * (1 + Math.log1p(energy) * 0.25); // Refined logarithmic gain
                // ---

                const bitIndex = Math.floor(blockIndex / CHANNELS);
                const bit = binaryString[bitIndex];
                if (bit === undefined) { continue; }
                
                const noise0 = new Float32Array(10).map(() => (p0() * 2.0 - 1.0) * blockGain);
                const noise1 = new Float32Array(10).map(() => (p1() * 2.0 - 1.0) * blockGain);
                const noiseVector = (bit === '0' ? noise0 : noise1);

                k = 0;
                for (let i = 2; i < N; i++) for (let j = 0; j < i + 1; j++) dctSquare[j*N+(i-j)] += noiseVector[k++];
                for (let j = 0; j < N - 1; j++) dctSquare[(1+j)*N+(N-1-j)] += noiseVector[k++];
            } else {
                // --- Adaptive Gain ---
                let energy = 0;
                for (let i = 1; i < 16; i++) energy += Math.abs(dctSquare[i]);
                const blockGain = gain * (1 + Math.log1p(energy) * 0.25); // Refined logarithmic gain
                // ---
                const noise0 = new Float32Array(10).map(() => (p0() * 2.0 - 1.0) * blockGain);
                const noise1 = new Float32Array(10).map(() => (p1() * 2.0 - 1.0) * blockGain);
                const c0 = pearsonCorrelation(dataVector, noise0);
                const c1 = pearsonCorrelation(dataVector, noise1);
                const confidence = Math.abs(c1 - c0);
                decodedChannels[channelIndex].push({ bit: c0 < c1 ? 1 : 0, confidence });
            }

            if (isEncode) {
                idctApply(dctSquare);
                for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) dataArray[(y+i)*size + (x+j)] = dctSquare[i*N+j];
            }
        }
        return isEncode ? null : decodedChannels;
    }
    
    function encode(image, message) {
        console.log('--- ENCODING ---');
        const [ baseSize, newWidth, newHeight ] = getPo2Size(image.width, image.height);
        canvas.width = newWidth; canvas.height = newHeight;
        ctx.drawImage(image, 0, 0, newWidth, newHeight);
        
        const fullMessage = message + DELIMITER;
        console.log(`Message: "${fullMessage}" (length: ${fullMessage.length}) capacity: ${RS_CONFIG.contentLength}`);
        
        const contentArray = new Uint8Array(RS_CONFIG.contentLength);
        for (let i = 0; i < fullMessage.length; i++) contentArray[i] = fullMessage.charCodeAt(i);
        
        const totalBytes = NUM_CHUNKS * RS_CONFIG.messageLength;
        const fullMessageEncoded = new Uint8Array(totalBytes);

        for (let i = 0; i < NUM_CHUNKS; i++) {
            const dataChunk = contentArray.subarray(i * RS_CONFIG.dataLength, (i + 1) * RS_CONFIG.dataLength);
            const buffer = new Uint8Array(RS_CONFIG.messageLength);
            buffer.set(dataChunk);
            RS_CONFIG.encoder.encode(buffer, RS_CONFIG.errorCorrectionLength);
            fullMessageEncoded.set(buffer, i * RS_CONFIG.messageLength);
        }
        const binaryString = Array.from(fullMessageEncoded).map(n => n.toString(2).padStart(8, '0')).join('');
        console.log(`Binary string length: ${binaryString.length} (capacity per channel: ${BITS_PER_CHANNEL})`);

        const imageData = ctx.getImageData(0, 0, baseSize, baseSize);
        const dataArray = new Float32Array(baseSize * baseSize);
        convertYCbCr(imageData, dataArray);
        recursiveDwt(dataArray, baseSize, mainSize, binaryString, true);
        convertRGB(imageData, dataArray);
        ctx.putImageData(imageData, 0, 0);
        console.log('Encoding complete.');
        return canvas.toDataURL('image/jpeg', 0.95);
    }

    function decode(image) {
        console.log('--- DECODING ---');
        const [ baseSize, newWidth, newHeight ] = getPo2Size(image.width, image.height);
        canvas.width = newWidth; canvas.height = newHeight;
        ctx.drawImage(image, 0, 0, newWidth, newHeight);
        
        const imageData = ctx.getImageData(0, 0, baseSize, baseSize);
        const dataArray = new Float32Array(baseSize * baseSize);
        convertYCbCr(imageData, dataArray);
        const decodedChannels = recursiveDwt(dataArray, baseSize, mainSize, null, false);
        console.log(`Extracted ${decodedChannels.length} channels of length ${decodedChannels[0].length}`);

        // Phase 2: Confidence-Based Voting
        const votedBitstream = [];
        for (let i = 0; i < BITS_PER_CHANNEL; i++) {
            const bit0 = decodedChannels[0][i];
            const bit1 = decodedChannels[1][i];

            if (bit0.bit === bit1.bit) {
                votedBitstream.push(bit0.bit); // They agree, simple case
            } else {
                // They disagree, vote based on confidence
                votedBitstream.push(bit0.confidence > bit1.confidence ? bit0.bit : bit1.bit);
            }
        }
        const binaryString = votedBitstream.join('');

        const intArray = new Uint8Array(binaryString.length / 8);
        for (let i = 0; i < intArray.length; i++) intArray[i] = parseInt(binaryString.substr(i*8, 8), 2);
        
        let res = "";
        for (let i = 0; i < NUM_CHUNKS; i++) {
            const chunkOffset = i * RS_CONFIG.messageLength;
            const chunk = intArray.slice(chunkOffset, chunkOffset + RS_CONFIG.messageLength);
            if (chunk.length < RS_CONFIG.messageLength) continue;
            
            try {
                RS_CONFIG.decoder.decode(chunk, RS_CONFIG.errorCorrectionLength);
                console.log(`Chunk ${i} decoded successfully.`);
                for (let j = 0; j < RS_CONFIG.dataLength; j++) {
                     if (chunk[j] === 0) break;
                    res += String.fromCharCode(chunk[j]);
                }
            } catch (e) {
                console.error(`Could not recover chunk ${i}. Message recovery failed.`, e.message);
                return "";
            }
        }
        console.log('Decoded raw text:', res);

        const delimiterIndex = res.indexOf(DELIMITER);
        if (delimiterIndex === -1) {
            console.log('Delimiter not found.');
            return "";
        }
        
        const finalMessage = res.substring(0, delimiterIndex);
        console.log('Final message:', finalMessage);
        return finalMessage;
    }
    // #endregion

    // #region UI Event Handlers
    function handleFileSelect(file, previewElement) {
        if (!file || !file.type.startsWith('image/')) { alert('Please select a valid image file.'); return; }
        const reader = new FileReader();
        reader.onload = (e) => { previewElement.src = e.target.result; previewElement.classList.remove('hidden'); };
        reader.readAsDataURL(file);
    }
    
    dropZoneEncode.addEventListener('click', () => fileInputEncode.click());
    fileInputEncode.addEventListener('change', (e) => { encodeFile = e.target.files[0]; handleFileSelect(encodeFile, previewEncode); });
    dropZoneEncode.addEventListener('dragover', (e) => { e.preventDefault(); dropZoneEncode.classList.add('dragover'); });
    dropZoneEncode.addEventListener('dragleave', (e) => { e.preventDefault(); dropZoneEncode.classList.remove('dragover'); });
    dropZoneEncode.addEventListener('drop', (e) => { e.preventDefault(); dropZoneEncode.classList.remove('dragover'); encodeFile = e.dataTransfer.files[0]; handleFileSelect(encodeFile, previewEncode); });

    encodeButton.addEventListener('click', () => {
        if (!encodeFile || !messageInput.value) { alert("Please select a file and enter a message."); return; }
        const img = new Image();
        img.onload = () => {
            const dataUrl = encode(img, messageInput.value);
            if (dataUrl) {
                const imgElement = document.createElement('img');
                imgElement.src = dataUrl;
                imgElement.className = 'rounded-md mt-4 max-w-full h-auto';
                const p = document.createElement('p');
                p.className = 'text-sm text-gray-400 mb-2';
                p.textContent = 'Encoded image (drag to decoder):';
                resultEncode.innerHTML = ''; resultEncode.appendChild(p); resultEncode.appendChild(imgElement);
            }
        };
        img.src = URL.createObjectURL(encodeFile);
    });

    function handleDecodeFile(file) {
        if (!file) return;
        const reader = new FileReader();
        reader.onload = function(e) {
            previewDecode.src = e.target.result; previewDecode.classList.remove('hidden');
            const img = new Image();
            img.onload = () => {
                const decodedMessage = decode(img);
                if (decodedMessage && decodedMessage.length > 0) {
                    resultDecode.innerHTML = `<p class="text-green-400">Success!</p><p class="font-mono p-2 bg-gray-600 rounded mt-2 break-all">${decodedMessage}</p>`;
                } else {
                    resultDecode.innerHTML = `<p class="text-red-400">No hidden message found.</p>`;
                }
            };
            img.src = e.target.result;
        };
        reader.readAsDataURL(file);
    }
    
    dropZoneDecode.addEventListener('click', () => fileInputDecode.click());
    fileInputDecode.addEventListener('change', (e) => handleDecodeFile(e.target.files[0]));
    dropZoneDecode.addEventListener('dragover', (e) => { e.preventDefault(); dropZoneEncode.classList.add('dragover'); });
    dropZoneDecode.addEventListener('dragleave', (e) => { e.preventDefault(); dropZoneEncode.classList.remove('dragover'); });
    dropZoneDecode.addEventListener('drop', (e) => { e.preventDefault(); dropZoneEncode.classList.remove('dragover'); handleDecodeFile(e.dataTransfer.files[0]); });
    
    processExampleImages().catch(console.error);
    // #endregion
});
