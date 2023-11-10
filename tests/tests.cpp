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
#include <algorithm>
#include <vector>

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
}
