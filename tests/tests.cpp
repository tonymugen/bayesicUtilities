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
#include <ratio>
#include <vector>
#include <array>

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
// precision for float comparisons
static constexpr double DPREC{1e-4};
// number of distribution deviates to generate
static constexpr size_t N_DEVIATES{1000};

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
		std::vector<double> concentrations{nClasses, eachConc};
	}
}
