/*
 * Copyright (c) 2023 Anthony J. Greenberg
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

/// Unit tests
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023 Anthony J. Greenberg
 * \version 0.1
 *
 * Unit tests using Catch2.
 *
 */

#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <ratio>
#include <vector>
#include <array>
#include <limits>

#include <iostream>

#include "random.hpp"
#include "index.hpp"
#include "utilities.hpp"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"

// Number of times tests of random events will be run
static constexpr uint16_t N_RAN_ITERATIONS{10};
// precision for float comparisons
static constexpr float FPREC{1e-4F};
// precision for double comparisons
static constexpr double DPREC{1e-4};
// number of distribution deviates to generate
static constexpr size_t N_DEVIATES{1000};
// 10X the double epsilon
static constexpr double MUL_DOUBLE_EPS{10.0 * std::numeric_limits<double>::epsilon()};

TEST_CASE("Random number generator works properly", "[prng]") {
	constexpr uint64_t seed{2153025619};
	BayesicSpace::RanDraw unseededPRNG;
	BayesicSpace::RanDraw seededPRNG(seed);
	SECTION("Pseudo-random integer generation") {
		constexpr uint64_t maxInt{13};
		constexpr uint64_t minInt{3};
		constexpr uint64_t correctPranInt{13307250832764291096UL};
		const uint64_t pranInt{seededPRNG.ranInt()};
		REQUIRE(pranInt == correctPranInt);
		REQUIRE(unseededPRNG.sampleInt(maxInt) < maxInt);
		REQUIRE(unseededPRNG.sampleInt(0UL) == 0UL);
		const uint64_t mmInt{unseededPRNG.sampleInt(minInt, maxInt)};
		REQUIRE(mmInt >= minInt);
		REQUIRE(mmInt <= maxInt);
		REQUIRE(unseededPRNG.sampleInt(maxInt, maxInt) == maxInt);
		for (uint16_t iIter = 0; iIter < N_RAN_ITERATIONS; ++iIter) {
			const std::vector<size_t> fyUp{unseededPRNG.fyIndexesUp(maxInt)};
			REQUIRE(fyUp.size() == maxInt - 1);
			const auto maxEl = std::max_element( fyUp.cbegin(), fyUp.cend() );
			REQUIRE(*maxEl < maxInt);
		}
		for (uint16_t iIter = 0; iIter < N_RAN_ITERATIONS; ++iIter) {
			std::vector<size_t> fyDown{unseededPRNG.fyIndexesDown(maxInt)};
			REQUIRE(fyDown.size() == maxInt - 1);
			const auto maxEl = std::max_element( fyDown.cbegin(), fyDown.cend() );
			REQUIRE(*maxEl < maxInt);
		}
	}
	SECTION("Distribution deviates") {
		std::array<double, N_DEVIATES> deviateArray{0.0};

		// uniforms
		constexpr double correctUnif{0.721388};
		const double unifVal = seededPRNG.runif();
		REQUIRE(fabs(correctUnif - unifVal) <= DPREC);
		std::for_each(deviateArray.begin(), deviateArray.end(), [&unseededPRNG](double &val){val = unseededPRNG.runif();});
		const auto mmRunif = std::minmax_element( deviateArray.cbegin(), deviateArray.cend() );
		REQUIRE((*mmRunif.second) <= 1.0);
		REQUIRE((*mmRunif.first)  >= 0.0);
		constexpr double correctUnifno{0.817305};
		const double unifValNO = seededPRNG.runifno();
		REQUIRE(fabs(correctUnifno - unifValNO) <= DPREC);
		std::for_each(deviateArray.begin(), deviateArray.end(), [&unseededPRNG](double &val){val = unseededPRNG.runifno();});
		const auto mmRunifNO = std::minmax_element( deviateArray.cbegin(), deviateArray.cend() );
		REQUIRE((*mmRunifNO.second) <  1.0);
		REQUIRE((*mmRunifNO.first)  >= 0.0);
		constexpr double correctUnifnz{0.807286};
		const double unifValNZ = seededPRNG.runifnz();
		REQUIRE(fabs(correctUnifnz - unifValNZ) <= DPREC);
		std::for_each(deviateArray.begin(), deviateArray.end(), [&unseededPRNG](double &val){val = unseededPRNG.runifnz();});
		const auto mmRunifNZ = std::minmax_element( deviateArray.cbegin(), deviateArray.cend() );
		REQUIRE((*mmRunifNZ.second) <= 1.0);
		REQUIRE((*mmRunifNZ.first)  >  0.0);
		constexpr double correctUnifop{0.737129};
		const double unifValOP = seededPRNG.runifop();
		REQUIRE(fabs(correctUnifop - unifValOP) <= DPREC);
		std::for_each(deviateArray.begin(), deviateArray.end(), [&unseededPRNG](double &val){val = unseededPRNG.runifop();});
		const auto mmRunifOP = std::minmax_element( deviateArray.cbegin(), deviateArray.cend() );
		REQUIRE((*mmRunifOP.second) < 1.0);
		REQUIRE((*mmRunifOP.first)  > 0.0);

		// normals
		constexpr double correctNorm{-0.166741};
		REQUIRE(fabs( correctNorm - seededPRNG.rnorm() ) <= DPREC);
		constexpr double sigma{2.9};
		constexpr double correctNormSig{0.120846};
		REQUIRE(fabs( correctNormSig - seededPRNG.rnorm(sigma) ) <= DPREC);
		constexpr double mean{101.245};
		constexpr double correctNormMuSig{96.9874};
		REQUIRE(fabs( correctNormMuSig - seededPRNG.rnorm(mean, sigma) ) <= DPREC);

		// gammas
		constexpr double negAlpha{-0.1};
		REQUIRE( std::isnan( unseededPRNG.rgamma(negAlpha) ) );
		REQUIRE( std::isnan( unseededPRNG.rgamma(0.0) ) );
		constexpr double smallAlpha{0.77};
		constexpr double correctRgammaSm{1.23799};
		REQUIRE(fabs( correctRgammaSm - seededPRNG.rgamma(smallAlpha) ) <= DPREC);
		constexpr double bigAlpha{11.3};
		constexpr double correctRgammaBg{11.0797};
		REQUIRE(fabs( correctRgammaBg - seededPRNG.rgamma(bigAlpha) ) <= DPREC);
		std::for_each(deviateArray.begin(), deviateArray.end(), [&unseededPRNG, smallAlpha](double &val){val = unseededPRNG.rgamma(smallAlpha);});
		const auto *const minGammSm = std::min_element( deviateArray.cbegin(), deviateArray.cend() );
		REQUIRE( (*minGammSm) > 0.0 );
		std::for_each(deviateArray.begin(), deviateArray.end(), [&unseededPRNG, bigAlpha](double &val){val = unseededPRNG.rgamma(bigAlpha);});
		const auto *const minGammBg = std::min_element( deviateArray.cbegin(), deviateArray.cend() );
		REQUIRE( (*minGammBg) > 0.0 );
		REQUIRE( std::isnan( unseededPRNG.rgamma(smallAlpha, 0.0) ) );
		REQUIRE( std::isnan( unseededPRNG.rgamma(bigAlpha, 0.0) ) );
		constexpr double smallBeta{0.5};
		constexpr double correctRgammaSmBeta{0.541837};
		REQUIRE(fabs( correctRgammaSmBeta - seededPRNG.rgamma(smallAlpha, smallBeta) ) <= DPREC);
		constexpr double bigBeta{2.0};
		constexpr double correctRgammaBgBeta{0.375599};
		REQUIRE(fabs( correctRgammaBgBeta - seededPRNG.rgamma(smallAlpha, bigBeta) ) <= DPREC);
		std::for_each(deviateArray.begin(), deviateArray.end(),
			[&unseededPRNG, smallAlpha, smallBeta](double &val){val = unseededPRNG.rgamma(smallAlpha, smallBeta);});
		const auto *const minGammSmSm = std::min_element( deviateArray.cbegin(), deviateArray.cend() );
		REQUIRE( (*minGammSmSm) > 0.0 );
		std::for_each(deviateArray.begin(), deviateArray.end(),
			[&unseededPRNG, smallAlpha, bigBeta](double &val){val = unseededPRNG.rgamma(smallAlpha, bigBeta);});
		const auto *const minGammSmBg = std::min_element( deviateArray.cbegin(), deviateArray.cend() );
		REQUIRE( (*minGammSmBg) > 0.0 );

		// Dirichlet
		constexpr size_t nClasses{11};
		constexpr double eachConc{0.45};
		constexpr double tinyConc{1e-5};

		std::vector<double> concentrations(nClasses, eachConc);
		const std::vector<double> dirProbs{unseededPRNG.rdirichlet(concentrations)};
		REQUIRE(std::all_of(dirProbs.cbegin(), dirProbs.cend(), [](double prob){return prob <= 1.0;}));
		REQUIRE(std::all_of(dirProbs.cbegin(), dirProbs.cend(), [](double prob){return prob >= 0.0;}));
		const auto pSum = std::accumulate(dirProbs.cbegin(), dirProbs.cend(), 0.0);
		REQUIRE(fabs(1.0 - pSum) < DPREC);

		for (size_t iConc = 0; iConc < concentrations.size(); ++iConc) { // if any alpha is 0.0, should return all NaNs
			std::vector<double> locConc = concentrations;
			locConc.at(iConc) = 0.0;
			const std::vector<double> locProb{unseededPRNG.rdirichlet(locConc)};
			REQUIRE(std::all_of(locProb.cbegin(), locProb.cend(), [](double prob){return std::isnan(prob);}));
		}

		// testing the underflow situation
		std::fill(concentrations.begin(), concentrations.end(), tinyConc);
		std::vector<double> dirProbsTiny{seededPRNG.rdirichlet(concentrations)};
		REQUIRE(std::all_of(dirProbsTiny.cbegin(), dirProbsTiny.cend(), [](double prob){return prob <= 1.0;}));
		REQUIRE(std::all_of(dirProbsTiny.cbegin(), dirProbsTiny.cend(), [](double prob){return prob >= 0.0;}));
		const auto tPsum = std::accumulate(dirProbsTiny.cbegin(), dirProbsTiny.cend(), 0.0);
		REQUIRE(fabs(1.0 - tPsum) < DPREC);

		// Chi-squared
		constexpr double degFreedom{3.0};
		constexpr double correctChiSq{1.89286};
		REQUIRE(fabs( correctChiSq - seededPRNG.rchisq(degFreedom) ) < DPREC);
		REQUIRE(std::isnan( unseededPRNG.rchisq(0.0) ));

		// Vitter
		constexpr double nLeft{100.0};
		constexpr uint64_t nLeftUI{100UL};
		constexpr double nToPic{5.0};
		constexpr double nToPicA{20.0}; // forces VitterA
		constexpr uint64_t correctVitter{46};
		REQUIRE( correctVitter == seededPRNG.vitter(nToPic, nLeft) );
		constexpr uint64_t correctVitterA{1};
		REQUIRE( correctVitterA == seededPRNG.vitter(nToPicA, nLeft) );

		std::array<uint64_t, N_DEVIATES> vitterArray{0UL};
		std::for_each(vitterArray.begin(), vitterArray.end(), 
				[&unseededPRNG, nToPic, nLeft](uint64_t &sample){sample = unseededPRNG.vitter(nToPic, nLeft);});
		REQUIRE(std::all_of(vitterArray.cbegin(), vitterArray.cend(), [nLeftUI](uint64_t sample){return sample < nLeftUI;}));
		std::for_each(vitterArray.begin(), vitterArray.end(), 
				[&unseededPRNG, nToPicA, nLeft](uint64_t &sample){sample = unseededPRNG.vitter(nToPicA, nLeft);});
		REQUIRE(std::all_of(vitterArray.cbegin(), vitterArray.cend(), [nLeftUI](uint64_t sample){return sample < nLeftUI;}));
	}
}

TEST_CASE("Utilities work properly", "[util]") {
	// swap
	constexpr size_t testInt1{13};
	constexpr size_t testInt2{47};
	size_t swapTest1{testInt1};
	size_t swapTest2{testInt2};
	BayesicSpace::swapXOR(swapTest1, swapTest2);
	REQUIRE(swapTest1 == testInt2);
	REQUIRE(swapTest2 == testInt1);
	BayesicSpace::swapXOR(swapTest1, swapTest1);
	REQUIRE(swapTest1 == testInt2);

	// logit
	constexpr double negLogit{-0.1};
	constexpr double gOneLogit{1.3};
	constexpr double smallLogit{0.25};
	constexpr double bigLogit{0.75};
	constexpr double testLogit{0.815};
	constexpr double correctTestLogit{1.48283};
	REQUIRE(fabs( BayesicSpace::logit(0.5) ) <= DPREC);
	REQUIRE(fabs( BayesicSpace::logit(smallLogit) + BayesicSpace::logit(bigLogit) ) <= DPREC);
	REQUIRE(fabs(BayesicSpace::logit(testLogit) - correctTestLogit) <= DPREC);
	REQUIRE( std::isnan( BayesicSpace::logit(negLogit) ) );
	REQUIRE( std::isnan( BayesicSpace::logit(gOneLogit) ) );

	// logistic
	constexpr double hugeLogistic{40.0};
	constexpr double bigLogistic{10.0};
	constexpr double niceLogistic{2.3};
	REQUIRE(fabs(BayesicSpace::logistic(0.0) - 0.5) <= DPREC);
	REQUIRE(BayesicSpace::logistic(hugeLogistic) == 1.0);
	REQUIRE(BayesicSpace::logistic(-hugeLogistic) == 0.0);
	REQUIRE(fabs(BayesicSpace::logistic(bigLogistic) + BayesicSpace::logistic(-bigLogistic) - 1.0) <= DPREC);
	REQUIRE(fabs(BayesicSpace::logistic(niceLogistic) + BayesicSpace::logistic(-niceLogistic) - 1.0) <= DPREC);
	REQUIRE(fabs(BayesicSpace::logistic( BayesicSpace::logit(testLogit) ) - testLogit) < DPREC);

	// lnGamma
	REQUIRE( std::isnan( BayesicSpace::lnGamma(-1.0) ) );
	REQUIRE(fabs(BayesicSpace::lnGamma(1.0)) <= MUL_DOUBLE_EPS);
	constexpr double testLG1{0.013};
	constexpr double testLG2{7.35};
	constexpr double correctLG1{4.33544};
	constexpr double correctLG2{7.24397};
	REQUIRE(fabs(BayesicSpace::lnGamma(testLG1) - correctLG1) <= DPREC);
	REQUIRE(fabs(BayesicSpace::lnGamma(testLG2) - correctLG2) <= DPREC);
	constexpr uint16_t nLGtests{10};
	double currZ{testLG1};
	for (uint16_t iLG = 0; iLG < nLGtests; ++iLG) {
		// testing the log-Gamma recursion
		REQUIRE(fabs( BayesicSpace::lnGamma(currZ) - BayesicSpace::lnGamma(currZ + 1.0) + log(currZ) ) <= DPREC);
		currZ += 0.5;
	}

	// digamma
	REQUIRE( std::isnan( BayesicSpace::digamma(-1.0) ) );
	constexpr double largeInput{3.1e21};
	constexpr double tinyInput{1.0e-17};
	constexpr double rootVal{1.461632144968362};
	constexpr double testDGalue1{0.13};
	constexpr double testDGalue2{11.34};
	constexpr double correctDGalue1{-8.07388};
	constexpr double correctDGalue2{2.3836};
	constexpr double eulerConstant{0.5772156649};
	REQUIRE(fabs( BayesicSpace::digamma(largeInput) - log(largeInput) ) <= DPREC);
	REQUIRE(fabs(BayesicSpace::digamma(tinyInput) + 1.0 / tinyInput) <= DPREC);
	REQUIRE(fabs( BayesicSpace::digamma(rootVal) ) <= MUL_DOUBLE_EPS);
	REQUIRE(fabs(BayesicSpace::digamma(testDGalue1) - correctDGalue1) <= DPREC);
	REQUIRE(fabs(BayesicSpace::digamma(testDGalue2) - correctDGalue2) <= DPREC);
	REQUIRE(fabs(BayesicSpace::digamma(1.0) + eulerConstant) <= DPREC);
	constexpr uint16_t nDGtests{10};
	double currX{testLG1};
	for (uint16_t iDG = 0; iDG < nDGtests; ++iDG) {
		// testing the digamma recursion
		REQUIRE(fabs( BayesicSpace::digamma(currX + 1.0) - BayesicSpace::digamma(currX)  - 1.0 / currX ) <= DPREC);
		currX += 0.5;
	}

	// dot products
	constexpr std::array<double, 11> testArray1 {
		0.9599027, 0.5084086, 1.1381224, -0.5632494, -0.5972133, 1.2555985,
		-0.2760020, -1.1939086, -0.7454797, 1.7004699, -1.4847591
	};
	constexpr std::array<double, 11> testArray2 {
		1.83941715332874, -1.01071106140829, -1.37696213968154, 0.227406257382323,
		2.08788585893456, 1.56032485376061, -1.01659801265071, 0.329786842248178,
		-0.851614769872224, -0.976615503871935, 0.46038209714626
	};
	constexpr double correctDP1{11.8791};
	constexpr double correctDP2{-1.553755};
	const std::vector<double> testVector1( testArray1.cbegin(), testArray1.cend() );
	const std::vector<double> testVector2( testArray2.cbegin(), testArray2.cend() );
	REQUIRE(fabs(BayesicSpace::dotProd( testVector1.cbegin(), testVector1.cend() ) - correctDP1) <= DPREC);
	REQUIRE(fabs(BayesicSpace::dotProd(testVector1.cbegin(), testVector1.cend(), testVector2.cbegin(), testVector2) - correctDP2) <= DPREC);
	const std::vector<double> testVectorShort( testArray2.cbegin(), testArray2.cend() - 1 );
	REQUIRE_THROWS_WITH(BayesicSpace::dotProd(testVector1.cbegin(), testVector1.cend(), testVectorShort.cbegin(), testVectorShort),
		Catch::Matchers::StartsWith("ERROR: second vector must have enough elements to complete inner product in ") );

	// means
	constexpr std::array<double, 10> badSequence{
		2.0e16, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, -2.0e16
	};
	std::vector<double> badVec( badSequence.cbegin(), badSequence.cend() );
	std::cout << BayesicSpace::stupidMean( badVec.cbegin(), badVec.cend() ) << "; " 
		<< BayesicSpace::stupidMean( badVec.cbegin() + 1, badVec.cend() - 1 ) << "; " 
		<< BayesicSpace::stableMean( badVec.cbegin(), badVec.cend() ) << "; "
		<< BayesicSpace::stableMean( badVec.cbegin() + 1, badVec.cend() - 1 ) << "\n";
}
