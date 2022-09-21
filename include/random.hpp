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
 * Class definition and interface documentation for facilities that generate random draws from various distributions.
 *
 */

#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <cmath>

namespace BayesicSpace {

	class RanDraw;


	/** \brief Random number generating class
	 *
	 * Generates pseudo-random deviates from a number of distributions.
	 * Uses an implementation of the 64-bit MT19937 ("Mersenne Twister")  \cite matsumoto98a pseudo-random number generator (PRNG) for random integers.
	 * This implementation of MT is ~35% faster than in `std::random` and ~250-fold faster than hardware RDRAND.
	 *
	 */
	class RanDraw {
	public:
		/** \brief Default constructor
		 *
		 * Seeded internally with a random number.
		 */
		RanDraw() : RanDraw( randomSeed_() ) {};
		/** \brief Constructor with seed
		 *
		 * Sets the provided seed.
		 *
		 * \param[in] seed seed value
		 */
		RanDraw(const uint64_t &seed);
		/** \brief Destructor */
		~RanDraw(){ };
		/** \brief Copy constructor (deleted)
		 *
		 * \param[in] old object to be copied
		 */
		RanDraw(const RanDraw &old) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] old object to be moved
		 */
		RanDraw(RanDraw &&old);

		/** \brief Copy assignment (deleted)
		 *
		 * \param[in] old object to be copied
		 */
		RanDraw & operator= (const RanDraw &old) = delete;
		/** \brief Move assignment
		 *
		 * \param[in] old object to be moved
		 */
		RanDraw & operator= (RanDraw &&old);
		/** \brief Generate random integer
		 *
		 * \return An unsigned random 64-bit integer
		 */
		uint64_t ranInt() noexcept;
		/** \brief Sample and integer from the \f$ [0, n) \f$ interval
		 *
		 * \param[in] max the maximal value \f$n\f$ (does not appear in the sample)
		 * \return sampled value
		 */
		uint64_t sampleInt(const uint64_t &max) noexcept { return this->ranInt() % max; };
		/** \brief Sample and integer from the \f$ [m, n) \f$ interval
		 *
		 *
		 * \param[in] min the minimal value \f$m\f$ (can appear in the sample)
		 * \param[in] max the maximal value \f$n\f$ (does not appear in the sample)
		 * \return sampled value
		 */
		uint64_t sampleInt(const uint64_t &min, const uint64_t &max) noexcept;
		/** \brief Draw non-negative intergers in random order
		 *
		 * Uses the Fisher-Yates-Durstenfeld algorithm to produce a random shuffle of integers in \f$ [0, N) \f$.
		 *
		 * \param[in] Nmax the upper bound of the integer sequence
		 *
		 * \return vector of \f$ N \f$ shuffled integers
		 */
		std::vector<uint64_t> shuffleUint(const uint64_t &N);

		/** \brief Generate a uniform deviate
		 *
		 * \return A double-precision value from the \f$ U[0,1]\f$ distribution
		 */
		double runif() noexcept {return 5.42101086242752217E-20 * static_cast<double>( this->ranInt() ); };
		/** \brief Generate a non-zero uniform deviate
		 *
		 * \return A double-precision value from the \f$ U(0,1]\f$ distribution
		 */
		double runifnz() noexcept;
		/** \brief Generate a non-one uniform deviate
		 *
		 * \return A double-precision value from the \f$ U[0,1)\f$ distribution
		 */
		double runifno() noexcept;
		/** \brief Generate an open-interval uniform deviate
		 *
		 * \return A double-precision value from the \f$ U(0,1)\f$ distribution
		 */
		double runifop() noexcept;
		/** \brief A standard Gaussian deviate
		 *
		 * Generates a Gaussian random value with mean \f$ \mu = 0.0 \f$ and standard deviation \f$ \sigma = 1.0 \f$. Implemented using a version of the Marsaglia and Tsang (2000) ziggurat algorithm, modified according to suggestions in the GSL implementation of the function.
		 *
		 * \return a sample from the standard Gaussian distribution
		 */
		double rnorm() noexcept;
		/** \brief A zero-mean Gaussian deviate
		 *
		 * Generates a Gaussian random value with mean \f$ \mu = 0.0 \f$ and standard deviation \f$ \sigma \f$. Implemented using a version of the Marsaglia and Tsang (2000) ziggurat algorithm, modified according to suggestions in the GSL implementation of the function.
		 *
		 * \param[in] sigma standard deviation
		 * \return a sample from the zero-mean Gaussian distribution
		 */
		double rnorm(const double &sigma) noexcept { return this->rnorm() * sigma; };
		/** \brief A Gaussian deviate
		 *
		 * Generates a Gaussian random value with mean \f$ \mu \f$ and standard deviation \f$ \sigma \f$. Implemented using a version of the Marsaglia and Tsang (2000) ziggurat algorithm, modified according to suggestions in the GSL implementation of the function.
		 *
		 * \param[in] mu standard deviation
		 * \param[in] sigma standard deviation
		 * \return a sample from the Gaussian distribution
		 */
		double rnorm(const double &mu, const double &sigma) noexcept { return mu + this->rnorm() * sigma; };
		/** \brief A standard Gamma deviate
		 *
		 * Generates a Gamma random variable with shape \f$ \alpha > 0 \f$ and standard scale \f$ \beta = 1.0 \f$. Implements the Marsaglia and Tsang (2000) method.
		 *
		 * \param[in] alpha shape parameter \f$ \alpha \f$
		 * \return a sample from the standard Gamma distribution
		 */
		double rgamma(const double &alpha) noexcept;
		/** \brief A general Gamma deviate
		 *
		 * Generates a Gamma random variable with shape \f$ \alpha > 0 \f$ and scale \f$ \beta > 0 \f$.
		 *
		 * \param[in] alpha shape parameter \f$ \alpha \f$
		 * \param[in] beta scale parameter \f$ \beta \f$
		 * \return a sample from the general Gamma distribution
		 */
		double rgamma(const double &alpha, const double &beta) noexcept { return beta > 0.0 ? ( this->rgamma(alpha) ) / beta : nan(""); };
		/** \brief A Dirichlet deviate
		 *
		 * Generates a vector of probabilities, given a vector of concentration parameters \f$ \alpha_K > 0 \f$.
		 *
		 * \param[in] alpha vector of concentration parameters
		 * \param[out] p vector of probabilities, must be the same length as \f$ \alpha \f$.
		 */
		void rdirichlet(const std::vector<double> &alpha, std::vector<double> &p) noexcept;
		/** \brief A chi-square deviate
		 *
		 * Generates a \f$ \chi^2 \f$ random variable with degrees of freedom \f$ \nu > 0.0 \f$.
		 *
		 * \param[in] nu degrees of freedom
		 * \return a sample from the \f$ \chi^2 \f$ distribution
		 */
		double rchisq(const double &nu) noexcept { return 2.0 * this->rgamma(nu / 2.0); };
		/** \brief Sample from Vitter's distribution, method A
		 *
		 * Given the number of remaining records in a file \f$N\f$ and the number of records \f$n\f$ remaining to be selected, sample the number of records to skip over. This function implements Vitter's \cite vitter84a \cite vitter87a method A. It is useful for online one-pass sampling of records from a file.
		 * While the inputs are integer, we pass them in as _double_ because that is more efficient for calculations.
		 *
		 * \param[in] n number of records remaining to be picked
		 * \param[in] N number of remaining records in the file
		 *
		 * \return the number of records to skip
		 */
		uint64_t vitterA(const double &n, const double &N) noexcept;
		/** \brief Sample from Vitter's distribution, method D
		 *
		 * Given the number of remaining records in a file \f$N\f$ and the number of records \f$n\f$ remaining to be selected, sample the number of records to skip over. This function implements Vitter's \cite vitter84a \cite vitter87a method D. It is useful for online one-pass sampling of records from a file.
		 * While the inputs are integer, we pass them in as _double_ because that is more efficient for calculations.
		 *
		 * \param[in] n number of records remaining to be picked
		 * \param[in] N number of remaining records in the file
		 *
		 * \return the number of records to skip
		 */
		uint64_t vitter(const double &n, const double &N) noexcept;
	private:
		// Mersenne twister constants
		/** \brief MT degree of recurrence */
		static const uint16_t n_;
		/** \brief MT middle word */
		static const uint16_t m_;
		/** \brief MT most significant 33 bits */
		static const uint64_t um_;
		/** \brief MT least significant 31 bits */
		static const uint64_t lm_;
		/** \brief MT tempering bitmask */
		static const uint64_t b_;
		/** \brief MT tempering bitmask */
		static const uint64_t c_;
		/** \brief MT tempering bitmask */
		static const uint64_t d_;
		/** \brief MT tempering shift */
		static const uint32_t l_;
		/** \brief MT tempering shift */
		static const uint32_t s_;
		/** \brief MT tempering shift */
		static const uint32_t t_;
		/** \brief MT tempering shift */
		static const uint32_t u_;
		/** \brief MT array of alternative values for the twist */
		static const std::array<uint64_t, 2> alt_;
		/** \brief MT generator state array */
		std::array<uint64_t, 312> mt_;
		/** \brief MT state of the array index */
		size_t mti_;
		/** \brief MT current state */
		uint64_t x_;
		/** \brief Right-most step position
		 *
		 * Used for the Gaussian ziggurat algorithm.
		 */
		static const double paramR_;

		/** \brief Tabulated values for the height of the ziggurat levels */
		static const std::array<double, 128> ytab_;
		/** \brief Tabulated \f$2^{24} \times x_i/x_{i+1}\f$ values
		 *
		 * Used for the Gaussian ziggurat algorithm.
		 */
		static const std::array<uint64_t, 128> ktab_;

		/** \brief Tabulated \f$2^{-24} \times x_i \f$ values
		 *
		 * Used for the Gaussian ziggurat algorithm.
		 */
		static const std::array<double, 128> wtab_;
		/** \brief Random seed for internal PRNG seeding */
		uint64_t randomSeed_() const ;
	};

}


