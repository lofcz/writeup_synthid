var GenericGF = function () { this.init.apply(this, arguments) };
GenericGF.prototype = {
	init : function (primitive, size, b) {
		this.primitive = primitive;
		this.size = size;
		this.generatorBase = b;

		this.expTable = new Int32Array(size);
		this.logTable = new Int32Array(size);

		var x = 1;
		for (var i = 0; i < size; i++) {
			this.expTable[i] = x;
			x *= 2; // we're assuming the generator alpha is 2
			if (x >= size) {
				x ^= primitive;
				x &= size-1;
			}
		}
		for (var i = 0; i < size-1; i++) {
			this.logTable[this.expTable[i]] = i;
		}
		// logTable[0] == 0 but this should never be used

		this.zero = new GenericGFPoly(this, new Int32Array([0]));
		this.one = new GenericGFPoly(this, new Int32Array([1]));
	},

	buildMonomial : function (degree, coefficient) {
		if (degree < 0) {
			throw new Error("IllegalArgumentException()");
		}
		if (coefficient === 0) {
			return this.zero;
		}
		var coefficients = new Int32Array(degree + 1);
		coefficients[0] = coefficient;
		return new GenericGFPoly(this, coefficients);
	},
    exp : function (a) { return this.expTable[a]; },
    log : function (a) { if (a === 0) { throw new Error("IllegalArgumentException()"); } return this.logTable[a]; },
    inverse: function (a) { if (a === 0) { throw new Error("ArithmeticException()"); } return this.expTable[this.size - this.logTable[a] - 1]; },
    multiply: function (a, b) { if (a === 0 || b === 0) { return 0; } return this.expTable[(this.logTable[a] + this.logTable[b]) % (this.size - 1)]; },
};
GenericGF.addOrSubtract = function (a, b) { return a ^ b };

var GenericGFPoly = function () { this.init.apply(this, arguments) };
GenericGFPoly.prototype = {
	init : function (field, coefficients) {
		if (coefficients.length === 0) {
			throw new Error("IllegalArgumentException()");
		}
		this.field = field;
		var coefficientsLength = coefficients.length;
		if (coefficientsLength > 1 && coefficients[0] === 0) {
			// Leading term must be non-zero for anything except the constant polynomial "0"
			var firstNonZero = 1;
			while (firstNonZero < coefficientsLength && coefficients[firstNonZero] === 0) {
				firstNonZero++;
			}
			if (firstNonZero == coefficientsLength) {
				this.coefficients = new Int32Array([0]);
			} else {
				this.coefficients = coefficients.subarray(firstNonZero, coefficientsLength);
			}
		} else {
			this.coefficients = coefficients;
		}
		this.degree = this.coefficients.length - 1;
	},

	getDegree : function () { return this.degree; },
	isZero : function () { return this.coefficients[0] === 0; },
	getCoefficient : function (degree) { return this.coefficients[this.coefficients.length - 1 - degree]; },
	evaluateAt : function (a) {
		if (a === 0) { return this.getCoefficient(0); }
		var coefficients = this.coefficients;
		var size = coefficients.length;
		var result;
		if (a == 1) {
			result = 0;
			for (var i = 0, len = coefficients.length; i < len; i++) {
				result = GenericGF.addOrSubtract(result, coefficients[i]);
			}
			return result;
		}
		result = coefficients[0];
		for (var i = 1; i < size; i++) {
			result = GenericGF.addOrSubtract(this.field.multiply(a, result), coefficients[i]);
		}
		return result;
	},
	addOrSubtract : function (other) {
		if (this.isZero()) { return other; }
		if (other.isZero()) { return this; }
		var smallerCoefficients = this.coefficients;
		var largerCoefficients = other.coefficients;
		if (smallerCoefficients.length > largerCoefficients.length) {
			[smallerCoefficients, largerCoefficients] = [largerCoefficients, smallerCoefficients];
		}
		var sumDiff = new Int32Array(largerCoefficients.length);
		var lengthDiff = largerCoefficients.length - smallerCoefficients.length;
        sumDiff.set(largerCoefficients.subarray(0, lengthDiff));
		for (var i = lengthDiff; i < largerCoefficients.length; i++) {
			sumDiff[i] = GenericGF.addOrSubtract(smallerCoefficients[i - lengthDiff], largerCoefficients[i]);
		}
		return new GenericGFPoly(this.field, sumDiff);
	},
	multiply : function (other) {
        if(typeof other === 'number') { return this.multiplyScalar(other); }
		if (this.isZero() || other.isZero()) { return this.field.zero; }
		var aCoefficients = this.coefficients;
		var aLength = aCoefficients.length;
		var bCoefficients = other.coefficients;
		var bLength = bCoefficients.length;
		var product = new Int32Array(aLength + bLength - 1);
		for (var i = 0; i < aLength; i++) {
			var aCoeff = aCoefficients[i];
			for (var j = 0; j < bLength; j++) {
				product[i + j] = GenericGF.addOrSubtract(product[i + j], this.field.multiply(aCoeff, bCoefficients[j]));
			}
		}
		return new GenericGFPoly(this.field, product);
	},
	multiplyScalar : function (scalar) {
		if (scalar === 0) { return this.field.zero; }
		if (scalar == 1) { return this; }
		var size = this.coefficients.length;
		var product = new Int32Array(size);
		for (var i = 0; i < size; i++) {
			product[i] = this.field.multiply(this.coefficients[i], scalar);
		}
		return new GenericGFPoly(this.field, product);
	},
	multiplyByMonomial : function (degree, coefficient) {
		if (degree < 0) { throw new Error('IllegalArgumentException()'); }
		if (coefficient === 0) { return this.field.zero; }
		var size = this.coefficients.length;
		var product = new Int32Array(size + degree);
		for (var i = 0; i < size; i++) {
			product[i] = this.field.multiply(this.coefficients[i], coefficient);
		}
		return new GenericGFPoly(this.field, product);
	},
	divide : function (other) {
		if (other.isZero()) { throw new Error('IllegalArgumentException("Divide by 0")'); }
		var quotient = this.field.zero;
		var remainder = this;
		var denominatorLeadingTerm = other.getCoefficient(other.getDegree());
		var inverseDenominatorLeadingTerm = this.field.inverse(denominatorLeadingTerm);
		while (remainder.getDegree() >= other.getDegree() && !remainder.isZero()) {
			var degreeDifference = remainder.getDegree() - other.getDegree();
			var scale = this.field.multiply(remainder.getCoefficient(remainder.getDegree()), inverseDenominatorLeadingTerm);
			var term = other.multiplyByMonomial(degreeDifference, scale);
			var iterationQuotient = this.field.buildMonomial(degreeDifference, scale);
			quotient = quotient.addOrSubtract(iterationQuotient);
			remainder = remainder.addOrSubtract(term);
		}
		return [ quotient, remainder ];
	}
};

var ReedSolomonEncoder = function () { this.init.apply(this, arguments) };
ReedSolomonEncoder.prototype = {
	init : function (field) { this.field = field; this.cachedGenerators = [new GenericGFPoly(field, new Int32Array([1]))]; },
	buildGenerator : function (degree) {
		if (degree >= this.cachedGenerators.length) {
			var lastGenerator = this.cachedGenerators[this.cachedGenerators.length - 1];
			for (var d = this.cachedGenerators.length; d <= degree; d++) {
				var nextGenerator = lastGenerator.multiply(new GenericGFPoly(this.field, new Int32Array([ 1, this.field.exp(d - 1 + this.field.generatorBase) ]) ));
				this.cachedGenerators.push(nextGenerator);
				lastGenerator = nextGenerator;
			}
		}
		return this.cachedGenerators[degree];
	},
	encode : function (toEncode, ecBytes) {
		var dataBytes = toEncode.length - ecBytes;
		var generator = this.buildGenerator(ecBytes);
		var infoCoefficients = new Int32Array(dataBytes);
		infoCoefficients.set(toEncode.subarray(0, dataBytes));
		var info = new GenericGFPoly(this.field, infoCoefficients);
		info = info.multiplyByMonomial(ecBytes, 1);
		var remainder = info.divide(generator)[1];
		var coefficients = remainder.coefficients;
		var numZeroCoefficients = ecBytes - coefficients.length;
        toEncode.set(new Int32Array(numZeroCoefficients), dataBytes);
		toEncode.set(coefficients, dataBytes + numZeroCoefficients);
	}
};

var ReedSolomonDecoder = function (field) { this.field = field; };
ReedSolomonDecoder.prototype = {
	decode : function (received, twoS) {
		var poly = new GenericGFPoly(this.field, received);
		var syndromeCoefficients = new Int32Array(twoS);
		var noError = true;
		for (var i = 0; i < twoS; i++) {
			var eval_ = poly.evaluateAt(this.field.exp(i + this.field.generatorBase));
			syndromeCoefficients[syndromeCoefficients.length - 1 - i] = eval_;
			if (eval_ !== 0) { noError = false; }
		}
		if (noError) { return; }
		var syndrome = new GenericGFPoly(this.field, syndromeCoefficients);
		var sigmaOmega = this.runEuclideanAlgorithm(this.field.buildMonomial(twoS, 1), syndrome, twoS);
		var sigma = sigmaOmega[0];
		var omega = sigmaOmega[1];
		var errorLocations = this.findErrorLocations(sigma);
		var errorMagnitudes = this.findErrorMagnitudes(omega, errorLocations);
		for (var i = 0; i < errorLocations.length; i++) {
			var position = received.length - 1 - this.field.log(errorLocations[i]);
			if (position < 0) { throw new Error('ReedSolomonException("Bad error location")'); }
			received[position] = GenericGF.addOrSubtract(received[position], errorMagnitudes[i]);
		}
	},
	runEuclideanAlgorithm : function (a, b, R) {
		if (a.getDegree() < b.getDegree()) { [a,b] = [b,a]; }
		var rLast = a; var r = b; var tLast = this.field.zero; var t = this.field.one;
		while (r.getDegree() >= R / 2) {
			var rLastLast = rLast; var tLastLast = tLast; rLast = r; tLast = t;
			if (rLast.isZero()) { throw new Error('ReedSolomonException("r_{i-1} was zero")'); }
			r = rLastLast;
			var q = this.field.zero;
			var denominatorLeadingTerm = rLast.getCoefficient(rLast.getDegree());
			var dltInverse = this.field.inverse(denominatorLeadingTerm);
			while (r.getDegree() >= rLast.getDegree() && !r.isZero()) {
				var degreeDiff = r.getDegree() - rLast.getDegree();
				var scale = this.field.multiply(r.getCoefficient(r.getDegree()), dltInverse);
				q = q.addOrSubtract(this.field.buildMonomial(degreeDiff, scale));
				r = r.addOrSubtract(rLast.multiplyByMonomial(degreeDiff, scale));
			}
			t = q.multiply(tLast).addOrSubtract(tLastLast);
		}
		var sigmaTildeAtZero = t.getCoefficient(0);
		if (sigmaTildeAtZero === 0) { throw new Error('ReedSolomonException("sigmaTilde(0) was zero")'); }
		var inverse = this.field.inverse(sigmaTildeAtZero);
		var sigma = t.multiply(inverse);
		var omega = r.multiply(inverse);
		return [ sigma, omega ];
	},
	findErrorLocations : function (errorLocator) {
		var numErrors = errorLocator.getDegree();
		if (numErrors == 1) { return new Int32Array([  errorLocator.getCoefficient(1)  ]); }
		var result = new Int32Array(numErrors);
		var e = 0;
		for (var i = 1; i < this.field.size && e < numErrors; i++) {
			if (errorLocator.evaluateAt(i) === 0) {
				result[e] = this.field.inverse(i);
				e++;
			}
		}
		if (e != numErrors) { throw new Error('ReedSolomonException("Error locator degree does not match number of roots")'); }
		return result;
	},
	findErrorMagnitudes : function (errorEvaluator, errorLocations) {
		var s = errorLocations.length;
		var result = new Int32Array(s);
		for (var i = 0; i < s; i++) {
			var xiInverse = this.field.inverse(errorLocations[i]);
			var denominator = 1;
			for (var j = 0; j < s; j++) { if (i != j) { denominator = this.field.multiply(denominator, GenericGF.addOrSubtract(1, this.field.multiply(errorLocations[j], xiInverse))); } }
			result[i] = this.field.multiply(errorEvaluator.evaluateAt(xiInverse), this.field.inverse(denominator));
			if (this.field.generatorBase !== 0) { result[i] = this.field.multiply(result[i], xiInverse); }
		}
		return result;
	}
};

export const GF = new GenericGF(0x011D, 256, 0); // x^8 + x^4 + x^3 + x^2 + 1
export const Encoder = ReedSolomonEncoder;
export const Decoder = ReedSolomonDecoder;
