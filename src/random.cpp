/*
 * Copyright (c) 2022 Anthony J. Greenberg
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

/// Random number generation
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 -- 2022 Anthony J. Greenberg
 * \version 1.0
 *
 * Class implementation for facilities that generate random draws from various distributions.
 *
 */

#include <array>
#include <vector>
#include <cmath>
#include <numeric>
#include <cassert>
#include <random>

#include "random.hpp"

using namespace BayesicSpace;

// RanDraw static members
constexpr uint64_t RanDraw::seedIncrement_{0x9e3779b97f4a7c15};
constexpr std::array<uint64_t, 2> RanDraw::initMultipliers_{0xbf58476d1ce4e5b9, 0x94d049bb133111eb};
constexpr std::array<uint64_t, 2> RanDraw::initShifts_{30, 27};
constexpr std::array<uint64_t, 3> RanDraw::riShifts_{17, 45, 19};

constexpr uint64_t RanDraw::llWordLen_{64};
constexpr std::array<double, 4> RanDraw::unifDivisors_{9007199254740991.0, 4503599627370495.5, 9007199254740992.0, 4503599627370496.0};
constexpr std::array<uint64_t, 2> RanDraw::unifShifts_{11, 12};
constexpr double RanDraw::paramR_{3.44428647676};
constexpr uint64_t RanDraw::signMask_{0x80};
constexpr std::array<double, 128> RanDraw::ytab_ = {
	1.0, 0.963598623011, 0.936280813353, 0.913041104253,
	0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
	0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
	0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
	0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
	0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
	0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
	0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
	0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
	0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
	0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
	0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
	0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
	0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
	0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
	0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
	0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
	0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
	0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
	0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
	0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
	0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
	0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
	0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
	0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
	0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
	0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
	0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
	0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
	0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
	0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
	0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};
constexpr std::array<uint64_t, 128> RanDraw::ktab_ = {
	0, 12590644, 14272653, 14988939,
	15384584, 15635009, 15807561, 15933577,
	16029594, 16105155, 16166147, 16216399,
	16258508, 16294295, 16325078, 16351831,
	16375291, 16396026, 16414479, 16431002,
	16445880, 16459343, 16471578, 16482744,
	16492970, 16502368, 16511031, 16519039,
	16526459, 16533352, 16539769, 16545755,
	16551348, 16556584, 16561493, 16566101,
	16570433, 16574511, 16578353, 16581977,
	16585398, 16588629, 16591685, 16594575,
	16597311, 16599901, 16602354, 16604679,
	16606881, 16608968, 16610945, 16612818,
	16614592, 16616272, 16617861, 16619363,
	16620782, 16622121, 16623383, 16624570,
	16625685, 16626730, 16627708, 16628619,
	16629465, 16630248, 16630969, 16631628,
	16632228, 16632768, 16633248, 16633671,
	16634034, 16634340, 16634586, 16634774,
	16634903, 16634972, 16634980, 16634926,
	16634810, 16634628, 16634381, 16634066,
	16633680, 16633222, 16632688, 16632075,
	16631380, 16630598, 16629726, 16628757,
	16627686, 16626507, 16625212, 16623794,
	16622243, 16620548, 16618698, 16616679,
	16614476, 16612071, 16609444, 16606571,
	16603425, 16599973, 16596178, 16591995,
	16587369, 16582237, 16576520, 16570120,
	16562917, 16554758, 16545450, 16534739,
	16522287, 16507638, 16490152, 16468907,
	16442518, 16408804, 16364095, 16301683,
	16207738, 16047994, 15704248, 15472926
};
constexpr std::array<double, 128> RanDraw::wtab_ = {
	1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
	3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
	3.8950989572e-08,  4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
	4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
	5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
	5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
	5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
	6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
	6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
	6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
	7.2824062723e-08,  7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
	7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
	7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
	8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
	8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
	8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
	9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
	9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
	9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
	1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
	1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
	1.09477300508e-07, 1.1042504257e-07,  1.11384564771e-07, 1.12356564007e-07,
	1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
	1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
	1.21808209468e-07, 1.2295639141e-07,  1.24129212952e-07, 1.25328445797e-07,
	1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
	1.31797105598e-07, 1.3320433736e-07,  1.34657379914e-07, 1.36160594606e-07,
	1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
	1.44635331499e-07, 1.4657889173e-07,  1.48632138436e-07, 1.50811780719e-07,
	1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
	1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
	1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};

RanDraw::RanDraw() noexcept : state_{0, 0, 0, 0} {
	std::random_device cppRandom;
	uint64_t seed = cppRandom();
	for (auto &eachS : state_){
		seed += seedIncrement_;
		uint64_t currentVal{seed};
		currentVal = ( currentVal ^ (currentVal >> initShifts_[0]) ) * initMultipliers_[0];
		currentVal = ( currentVal ^ (currentVal >> initShifts_[1]) ) * initMultipliers_[1];
		eachS      = currentVal;
	}
}

RanDraw::RanDraw(const uint64_t &seed) noexcept : state_{0, 0, 0, 0} {
	uint64_t locSeed{seed};
	for (auto &eachS : state_){
		locSeed += seedIncrement_;
		uint64_t lsCopy{locSeed};
		lsCopy = (lsCopy ^ (lsCopy >> initShifts_[0])) * initMultipliers_[0];
		lsCopy = (lsCopy ^ (lsCopy >> initShifts_[1])) * initMultipliers_[1];
		eachS  = lsCopy;
	}
}

RanDraw::RanDraw(RanDraw &&old) noexcept : state_{old.state_} {
	*this = std::move(old);
}

RanDraw & RanDraw::operator= (RanDraw &&old) noexcept {
	if (this != &old){
		state_ = old.state_;
	}
	return *this;
}

uint64_t RanDraw::ranInt() noexcept {
	const uint64_t stateSum  = state_[0] + state_[3];
	const uint64_t sumRotSum = ( (stateSum << 23) | (stateSum >> 41) ) + state_[0]; // this is rotl and sum; compilers will generate the right instruction

	state_[2] ^= state_[0];
	state_[3] ^= state_[1];
	state_[1] ^= state_[2];
	state_[0] ^= state_[3];

	state_[2] ^= state_[1] << riShifts_[0];

	state_[3]  = (state_[3] << riShifts_[1]) | (state_[3] >> riShifts_[2]); // rotl again

	return sumRotSum;
}

uint64_t RanDraw::sampleInt(const uint64_t &max) noexcept {
	const auto max128  = static_cast<__uint128_t>(max);
	uint64_t bigRanInt = this->ranInt();
	__uint128_t prod   = static_cast<__uint128_t>(bigRanInt) * max128;
	auto prod64        = static_cast<uint64_t>(prod);
	if (prod64 < max){
		const uint64_t test = static_cast<uint64_t>(-max) % max;
		while (prod64 < test){
			bigRanInt = this->ranInt();
			prod      = static_cast<__uint128_t>(bigRanInt) * max128;
			prod64    = static_cast<uint64_t>(prod);
		}
	}
	return static_cast<uint64_t>(prod >> llWordLen_);
}

uint64_t RanDraw::sampleInt(const uint64_t &min, const uint64_t &max) noexcept {
	assert( (min < max) && "ERROR: Lower bound not smaller than upper bound in RanDraw::sampleInt()" );

	return min + this->ranInt() % (max - min);
}

std::vector<size_t> RanDraw::fyIndexesDown(const size_t &Nidx){
	std::vector<size_t> perInd(Nidx, 0);
	for (size_t i = Nidx - 1; i > 0; --i){
		perInd[i] = this->sampleInt(i + 1); // sampleInt(max) samples j < max
	}
	return perInd;
}

std::vector<size_t> RanDraw::fyIndexesUp(const size_t &Nidx){
	std::vector<size_t> perInd;
	perInd.reserve(Nidx - 1);
	for (size_t i = 0; i < Nidx - 1; ++i){
		perInd.push_back( this->sampleInt(i, Nidx) );
	}
	return perInd;
}

double RanDraw::rnorm() noexcept {
	double absValue{0.0};
	double sign{0.0};
	uint64_t iIdx{0};
	uint64_t jIdx{0};
	constexpr uint64_t mask255{255};
	constexpr uint64_t maxJind{16777216};
	constexpr uint64_t mask127{0x7f};

	while (true){
		iIdx  = this->ranInt() & mask255;        // choose the step
		jIdx  = this->ranInt() % maxJind;        // sample from 2^24

		sign  = ( (iIdx & signMask_) > 0 ) ? 1.0 : -1.0;
		iIdx &= mask127;                               // convert iIdx to [0, 127], an index into the ziggurat slices

		absValue = static_cast<double>(jIdx) * wtab_.at(iIdx);
		if ( jIdx < ktab_.at(iIdx) ){ // ktab_ is used to test for acceptance without floating point operations
			break;
		}
		double ySample{0.0};
		if ( iIdx < (ytab_.size() - 1) ){
			ySample = ytab_.at(iIdx + 1) + ( ytab_.at(iIdx) - ytab_.at(iIdx + 1) ) * ( this->runifno() );
		} else {
			double unif1 = this->runifnz(); // (0,1]
			double unif2 = this->runifno(); // [0,1)
			absValue = paramR_ - log(unif1) / paramR_;
			ySample = exp( -paramR_ * (absValue - 0.5 * paramR_) ) * unif2;
		}
		if ( ySample < exp(-0.5 * absValue * absValue) ){
			break;
		}
	}

	return sign * absValue;
}

double RanDraw::rgamma(const double &alpha) noexcept {
	if (alpha <= 0.0){
		return nan("");
	}
	if (alpha < 1.0){
		return this->rgamma(alpha + 1.0) * pow(this->runifop(), 1.0 / alpha);
	}
	constexpr double oneThird{0.3333333333};
	constexpr double smallMultiplier{0.0331};
	double stdNorm{0.0};
	double scaledNorm{0.0};
	double unifSample{0.0};
	double shiftedAlpha{alpha - oneThird};
	const double normSD{oneThird / sqrt(shiftedAlpha)};

	while (true){
		do {
			stdNorm    = this->rnorm();
			scaledNorm = 1.0 + normSD * stdNorm;
		} while (scaledNorm <= 0.0);

		scaledNorm = scaledNorm * scaledNorm * scaledNorm;
		unifSample = this->runifop();
		const double stdNormSq = {stdNorm * stdNorm};

		if (unifSample < 1.0 - smallMultiplier * stdNormSq * stdNormSq){
			break;
		}
		if ( log(unifSample) < 0.5 * stdNormSq + shiftedAlpha * ( 1.0 - scaledNorm + log(scaledNorm) ) ){
			break;
		}
	}

	return shiftedAlpha * scaledNorm;
}

void RanDraw::rdirichlet(const std::vector<double> &alpha, std::vector<double> &probabilities) noexcept {
	assert( ( alpha.size() == p.size() ) && "ERROR: length of alpha vector not the same as length of p vector in RanDraw::rdirichlet()" );

	double sum = 0.0;
	for (size_t k = 0; k < alpha.size(); ++k){
		probabilities[k] = this->rgamma(alpha[k]);
		sum += probabilities[k];
	}
	for (auto &eachP : probabilities){
		eachP = eachP / sum;
	}
}

uint64_t RanDraw::vitterA(const double &nToPick, const double &Nremain) noexcept {
	uint64_t sample{0};
	auto top{static_cast<double>(Nremain - nToPick)};
	double quot{top / Nremain};

	// some trivial conditions first
	if ( (nToPick == 0) || (nToPick > Nremain) ){
		return sample;
	}
	if (nToPick == 1){
		return static_cast<uint64_t>( floor( Nremain * this->runifop() ) );
	}
	const double unifSample{this->runifop()};
	double Nloc{Nremain};
	while (quot > unifSample){
		++sample;
		--top;
		--Nloc;
		quot = quot * top / Nloc;
	}

	return sample;
}

uint64_t RanDraw::vitter(const double &nToPick, const double &Nremain) noexcept {
	// Notation is as close as possible to Vitter (1987), Appendix A2
	uint64_t sample{0};
	const double alphaInv{13.0};
	if (nToPick >= Nremain / alphaInv){ // if the threshold is not satisfied, use Vitter's A algorithm
		return this->vitterA(nToPick, Nremain);
	}
	if (nToPick == 1){ // trivial case
		return static_cast<uint64_t>( floor( Nremain * this->runifop() ) );
	}

	// if we pass all thresholds, we use Vitter's rejection scheme
	const double nInv{1.0 / nToPick};
	const double nMin1inv{1.0 / (nToPick - 1.0)};
	const double qu1db{1.0 + Nremain - nToPick};
	const auto qu1{static_cast<uint64_t>(qu1db)};

	// outer loop in Vitter (1987) A2
	while (true){
		double roofSample{0.0};                  // covering distribution sample for rejection
		double Vprime{0.0};
		do {  // step D2; generate U and X
			Vprime     = pow(this->runifop(), nInv); // same as exp(log(u)*nInv)
			roofSample = Nremain * (1.0 - Vprime);
			sample     = static_cast<uint64_t>( floor(roofSample) );
		} while (sample >= qu1);

		const auto dSample{static_cast<double>(sample)};
		const double y1sample{pow(this->runifop() * Nremain / qu1db, nMin1inv)};
		Vprime = y1sample * (1.0 - roofSample / Nremain) * ( qu1db / (qu1db - dSample) );
		if (Vprime < 1.0){ // Step D3: accept test 2.8 (Vitter 1987)
			break;
		}
		// moving on to Step D4
		double y2sample{1.0};
		double top{Nremain - 1.0};
		double bottom{0.0};
		uint64_t limit{0};
		if (dSample < nToPick - 1.0){
			bottom = Nremain - nToPick;
			limit  = static_cast<uint64_t>(Nremain - dSample);
		} else {
			bottom = Nremain - dSample - 1.0;
			limit  = qu1;
		}

		// calculate f(|_X_|)
		auto loopVar = static_cast<uint64_t>(Nremain - 1.0);
		while (loopVar >= limit){
			y2sample *= top / bottom;
			top      -= 1.0;
			bottom   -= 1.0;
			--loopVar;
		}
		if ( Nremain / (Nremain - roofSample) >= y1sample * pow(y2sample, nMin1inv) ){ // Accept D4 condition
			break;
		}
		// reject everything, go back to the start
	}
	return sample;
}



