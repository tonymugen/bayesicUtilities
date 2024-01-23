/*
 * Copyright (c) 2020 Anthony J. Greenberg
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Numerical utilities implementation
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2020 Anthony J. Greenberg
 * \version 1.0
 *
 * Class implementation for a set of numerical utilities.
 * Implemented as a class because this seems to be the only way for these methods to be included using Rcpp with no compilation errors.
 *
 */

#include <cstdint>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <limits>
#include <cassert>

#include "utilities.hpp"


void BayesicSpace::swapXOR(size_t &first, size_t &second) noexcept {
	if (&first != &second) { // no move needed if this is actually the same variable
		first  ^= second;
		second ^= first;
		first  ^= second;
	}
}

double BayesicSpace::logit(const double &probability) noexcept {
	return log(probability) - log(1.0 - probability);
}

double BayesicSpace::logistic(const double &input) noexcept {
	// 35.0 is the magic number because logistic(-35) ~ EPS
	constexpr double cutOff{35.0};
	constexpr double smallishCut{7.0};
	// the other cut-offs have been empirically determined
	if (input <= -cutOff) {
		return 0.0;
	}
	if (input >= cutOff) {
		return 1.0;
	}
	if (input <= -smallishCut) { // approximation for smallish x
		return exp(input);
	}
	if (input >= smallishCut/2.0) {  // approximation for largish x
		return 1.0 - exp(-input);
	}
	return 1.0 / ( 1.0 + exp(-input) );
}

double BayesicSpace::lnGamma(double input) noexcept {
	if (input <= 0.0) {
		return nan("");
	}
	// alternative method for small input values
	constexpr double inputCutOff{0.02};
	if (input < inputCutOff) {
		constexpr std::array<double, 10> smallInputCoeff { // from the GSL lngamma_sgn_0 function
			-0.07721566490153286061,
			-0.01094400467202744461,
			 0.09252092391911371098,
			-0.01827191316559981266,
			 0.01800493109685479790,
			-0.00685088537872380685,
			 0.00399823955756846603,
			-0.00189430621687107802,
			 0.00097473237804513221,
			-0.00048434392722255893
		};
		double lancsozG = std::accumulate(
			std::next( smallInputCoeff.crbegin() ),
			smallInputCoeff.crend(),
			smallInputCoeff.back(),
			[input](const double &cumLG, const double &coeff){return coeff + input * cumLG;}
		);
		lancsozG *= input;
		return log( (lancsozG + 1.0 / (1.0 + input) + 0.5 * input) / input );
	}

	// main Lanczos method
	constexpr double logRootTwoPi{0.91893853320467267};       // 0.5*ln(2.0*pi)
	constexpr double eulerE{2.71828182845904523536028747135};
	constexpr double lanczosOffset{7.5};                      // from the GSL lngamma_lanczos function
	constexpr double term2offset{7.0};
	constexpr std::array<double, 9> lanczos7coeff {
		0.99999999999980993227684700473478,
		676.520368121885098567009190444019,
		-1259.13921672240287047156078755283,
		771.3234287776530788486528258894,
		-176.61502916214059906584551354,
		12.507343278686904814458936853,
		-0.13857109526572011689554707,
		9.984369578019570859563e-6,
		1.50563273514931155834e-7
	};
	input -= 1.0;
	double lanczosAg{lanczos7coeff[0]};
	double sumIdx{1.0};
	std::for_each(
		std::next( lanczos7coeff.cbegin() ),
		lanczos7coeff.cend(),
		[&lanczosAg, &sumIdx, input](double cSubk){
			lanczosAg += cSubk / (input + sumIdx);
			sumIdx += 1.0;}
	);

	return (input + 0.5) * log( (input + lanczosOffset) / eulerE ) + (logRootTwoPi + log(lanczosAg) - term2offset);
}

double BayesicSpace::digamma(const double &input) noexcept {
	constexpr int32_t nMax{100};
	constexpr int32_t nCrit = std::min(-std::numeric_limits<double>::min_exponent, std::numeric_limits<double>::max_exponent);
	constexpr double r1m4{0.5 * std::numeric_limits<double>::epsilon()};
	constexpr double r1m5{0.301029995663981195213738894724};            // log_10(2)
	constexpr double wdtol = std::max(r1m4, 0.5e-18);
	constexpr double elim  = 2.302 * (static_cast<double>(nCrit) * r1m5 - 3.0);  // = 700.6174...
	constexpr double lnRcrit{18.06};
	constexpr double lnFcrit{3.0};
	constexpr double lrg = 1.0 / ( 2.0 * std::numeric_limits<double>::epsilon() );
	constexpr double rln = std::min(r1m5 * static_cast<double>(std::numeric_limits<double>::digits), lnRcrit);
	constexpr double fln = 3.50 + 0.40 * std::max(rln, lnFcrit);

	constexpr std::array<double, 22> bvalues {
		1.00000000000000000e+00, -5.00000000000000000e-01,
		1.66666666666666667e-01, -3.33333333333333333e-02,
		2.38095238095238095e-02, -3.33333333333333333e-02,
		7.57575757575757576e-02, -2.53113553113553114e-01,
		1.16666666666666667e+00, -7.09215686274509804e+00,
		5.49711779448621554e+01, -5.29124242424242424e+02,
		6.19212318840579710e+03, -8.65802531135531136e+04,
		1.42551716666666667e+06, -2.72982310678160920e+07,
		6.01580873900642368e+08, -1.51163157670921569e+10,
		4.29614643061166667e+11, -1.37116552050883328e+13,
		4.88332318973593167e+14, -1.92965793419400681e+16
	};

	// very large input
	const double lnInput = log(input);
	if (std::isnan(input) || input <= 0.0 || lnInput < -elim) {
		return nan("");
	}
	if(input * lnInput > lrg) {
		return lnInput;
	}
	if (input < wdtol) {
		return -1.0 / input;
	}

	// regular calculations
	const double minInput{ceil(fln)};
	double inputCopy{input};          // make a copy of the input
	double lnInputCopy{lnInput};
	double incInput{0.0};
	if (input < minInput) {
		incInput    = minInput - floor(input);
		inputCopy   = input + incInput;
		lnInputCopy = log(inputCopy);
	}

	double seriesElement = 2.0 * lnInputCopy;
	if (seriesElement <= elim) { // for input not large
		double firstSeriesElement{0.5 / inputCopy};
		double testValue  = wdtol * firstSeriesElement;
		double rxsq       = 1.0 / (inputCopy * inputCopy);
		double seriesProd = 0.5 * rxsq;
		double seriesSum  = seriesProd * bvalues[2];
		if (fabs(seriesSum) >= testValue) {
			seriesElement = 2.0;
			for(size_t bernoulliIndex = 3; bernoulliIndex < bvalues.size(); ++bernoulliIndex) {
				seriesProd *= ( (seriesElement + 1.0) / (seriesElement + 1.0) ) * ( seriesElement / (seriesElement + 2.0) ) * rxsq;
				double tmp  = seriesProd * bvalues.at(bernoulliIndex);
				if (fabs(tmp) < testValue) {
					break;
				}
				seriesSum     += tmp;
				seriesElement += 2.0;
			}
		}
		seriesSum += firstSeriesElement;
		if (incInput > 0.0) {
			// backward recursion from inputCopy to input
			const auto integerInput = static_cast<int32_t>(incInput);
			assert( (integerInput <= nMax) && "ERROR: Increment nx too large in digamma()" );
			for(int32_t idx = 1; idx <= integerInput; ++idx) {
				seriesSum += 1.0 / ( input + static_cast<double>(integerInput - idx) ); // avoid disastrous cancellation, according to the comment in the R code
			}
		}
		return lnInputCopy - seriesSum;
	}
	double seriesSum = -input;
	double den       = input;
	for(uint32_t i = 0; i < static_cast<uint32_t>(fln) + 1; ++i) { // fln is derived from positive constants
		den       += 1.0;
		seriesSum += 1.0 / den;
	}
	return -seriesSum;
}

double BayesicSpace::dotProd(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end) noexcept {
	const double dotProd = std::inner_product(begin, end, begin, 0.0);
	return dotProd;
}

double BayesicSpace::dotProd(std::vector<double>::const_iterator firstBegin,
								std::vector<double>::const_iterator firstEnd,
								std::vector<double>::const_iterator secondBegin,
								const std::vector<double> &second) {
	const auto firstSize  = std::distance(firstBegin, firstEnd);
	const auto secondSize = std::distance(secondBegin, second.cend());
	if (secondSize < firstSize) {
		throw std::string("ERROR: second vector must have enough elements to complete inner product in ") +
			std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	const double dotProd = std::inner_product(firstBegin, firstEnd, secondBegin, 0.0); // NOLINT

	return dotProd;
}

double BayesicSpace::stableMean(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end) noexcept {
	double index{1.0};

	const double mean = std::accumulate(
			begin,
			end,
			0.0, 
			[&index](double currentSum, double currentValue){
				currentSum += (currentValue - currentSum) / index;
				index      += 1.0;
				return currentSum;
			});
	return mean;
}

BayesicSpace::ValueWithWeight BayesicSpace::updateWeightedMean(const BayesicSpace::ValueWithWeight &nextDataPoint,
					const BayesicSpace::ValueWithWeight &currentMean) noexcept {
	BayesicSpace::ValueWithWeight result{currentMean};
	if ( nextDataPoint.weight > std::numeric_limits<double>::epsilon() ) {
		const double tmp = currentMean.value * currentMean.weight;
		result.weight += nextDataPoint.weight;
		result.value = (tmp + nextDataPoint.weight * nextDataPoint.value) / result.weight;
	}
	return result;
}
