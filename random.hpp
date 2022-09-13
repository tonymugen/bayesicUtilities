/*
 * Copyright (c) 2018 Anthony J. Greenberg
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
 * \copyright Copyright (c) 2017 Anthony J. Greenberg
 * \version 1.0
 *
 * Class definition and interface documentation for facilities that generate random draws from various distributions.
 *
 */

#ifndef random_hpp
#define random_hpp

#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <memory>

namespace BayesicSpace {
	class RanDraw;
	class Generate;
	class GenerateHR;
	class GenerateMT;

	/** \brief Abstract base random number class
	 *
	 * Provides the interface for random or pseudorandom (depending on derived class) generation. For internal use by the `RanDraw` interface class.
	 */
	class Generate {
	public:
		/** \brief Destructor */
		virtual ~Generate(){};
		/** \brief Generate a (pseudo-)random 64-bit unsigned integer
		 *
		 * \return random or pseudo-random 64-bit unsigned integer
		 */
		virtual uint64_t ranInt() noexcept = 0;
		/** \brief Generate a (pseudo-)random 64-bit unsigned integer
		 *
		 * \param[out] ranInt (pseudo-)ranom number
		 * \return (P)RNG status report
		 */
		virtual int32_t ranInt(uint64_t &ranInt) noexcept = 0;
	protected:
		/** \brief Protected default constructor */
		Generate(){};
		/** \brief Protected copy constructor (deleted)
		 *
		 * \param[in] old object to copy
		 */
		Generate(const Generate &old) = delete;
		/** \brief Protected move constructor
		 *
		 * \param[in] old object to move
		 */
		Generate(Generate &&old) = default;
		/** \brief Protected copy assignment operator (deleted)
		 *
		 * \param[in] old object to copy
		 */
		Generate & operator= (const Generate &old) = delete;
		/** \brief Protected move assignment
		 *
		 * \param[in] old object to move
		 */
		Generate & operator= (Generate &&old) = default;
	};

	/** \brief Hardware random number generating class
	 *
	 * Generates random deviates from a number of distributions, using hardware random numbers (_RDRAND_ processor instruction). Health of the RDRAND generator is tested every time a new number is required. Throws a `string` object "RDRAND_failed" if the test fails.
	 * The implementation of random 64-bit integer generation follows [Intel's suggestions](https://software.intel.com/en-us/articles/intel-digital-random-number-generator-drng-software-implementation-guide ).
	 */
	class GenerateHR final : public Generate {
	public:
		/** \brief Default constructor */
		GenerateHR(){};
		/** \brief Destructor */
		~GenerateHR(){};
		/** \brief Copy constructor (deleted)
		 *
		 * \param[in] old object to copy
		 */
		GenerateHR(const GenerateHR &old) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] old object to move
		 */
		GenerateHR(GenerateHR &&old) = default;
		/** \brief Copy assignment operator (deleted)
		 *
		 * \param[in] old object to copy
		 */
		GenerateHR & operator= (const GenerateHR &old) = delete;
		/** \brief Move assignment
		 *
		 * \param[in] old object to move
		 */
		GenerateHR & operator= (GenerateHR &&old) = default;

		/** \brief Generate a random 64-bit unsigned integer
		 *
		 * Health of the RNG not checked to increase efficiency.
		 * Extensive testing has not revealed any instances of RNG failure in practice.
		 *
		 * \return digital random 64-bit unsigned integer
		 */
		uint64_t ranInt() noexcept override;
		/** \brief Generate a random 64-bit unsigned integer with health check
		 *
		 * RNG generation success returns 1, failure returns 0.
		 * In case of failure, the random number value cannot be trusted (should be 0).
		 *
		 * \param[out] ranInt random number
		 * \return RNG status report
		 */
		int32_t ranInt(uint64_t &ranInt) noexcept override;
	protected:
		// no protected members
	};

	/** \brief Pseudo-random number generator
	 *
	 * An implementation of the 64-bit MT19937 ("Mersenne Twister")  \cite matsumoto98a pseudo-random number generator (PRNG). The constructor automatically seeds the PRNG. The implementation was guided by the reference code [posted by the authors](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html ).
	 */
	class GenerateMT final : public Generate {
	public:
		/** \brief Default constructor
		 *
		 * Seeds the PRNG with a call to the _RDTSC_ instruction.
		 */
		GenerateMT();
		/** \brief Constructor with seed
		 *
		 * \param[in] seed seed value
		 */
		GenerateMT(const uint64_t &seed);
		/** \brief Protected destructor */
		~GenerateMT(){};
		/** \brief Copy constructor (deleted)
		 *
		 * \param[in] old object to copy
		 */
		GenerateMT(const GenerateMT &old) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] old object to move
		 */
		GenerateMT(GenerateMT &&old);
		/** \brief Copy assignment operator (deleted)
		 *
		 * \param[in] old object to copy
		 */
		GenerateMT & operator= (const GenerateMT &old) = delete;
		/** \brief Move assignment
		 *
		 * \param[in] old object to move
		 */
		GenerateMT & operator= (GenerateMT &&old);

		/** \brief Generate a pseudo-random 64-bit unsigned integer
		 *
		 *
		 * \return pseudo-random 64-bit unsigned integer
		 */
		uint64_t ranInt() noexcept override;
		/** \brief Generate a pseudo-random 64-bit unsigned integer with health check
		 *
		 * Health check always returns 1 at moment and can be ignored.
		 *
		 * \param[out] ranInt pseudo-random number
		 * \return PRNG status report
		 */
		int32_t ranInt(uint64_t &ranInt) noexcept override;
	protected:
		/** \brief Degree of recurrence */
		static const uint16_t n_;
		/** \brief Middle word */
		static const uint16_t m_;
		/** \brief Most significant 33 bits */
		static const uint64_t um_;
		/** \brief Least significant 31 bits */
		static const uint64_t lm_;
		/** \brief Tempering bitmask */
		static const uint64_t b_;
		/** \brief Tempering bitmask */
		static const uint64_t c_;
		/** \brief Tempering bitmask */
		static const uint64_t d_;
		/** \brief Tempering shift */
		static const uint32_t l_;
		/** \brief Tempering shift */
		static const uint32_t s_;
		/** \brief Tempering shift */
		static const uint32_t t_;
		/** \brief Tempering shift */
		static const uint32_t u_;
		/** \brief Array of alternative values for the twist */
		static const std::array<uint64_t, 2> alt_;
		/** \brief Generator state array */
		std::array<uint64_t, 312> mt_;
		/** \brief State of the array index */
		size_t mti_;
		/** \brief Current state */
		uint64_t x_;
	};

	/** \brief Random number generating class
	 *
	 * Generates (pseudo-)random deviates from a number of distributions. If hardware random numbers are supported, uses them. Otherwise, falls back to 64-bit MT19937 ("Mersenne Twister").
	 *
	 */
	class RanDraw {
	public:
		/** \brief Default constructor
		 *
		 * Checks if the processor provides hardware random number support. Seeds the Mersenne Twister if not.
		 * Throws "CPU_unsupported" string object if the CPU is not AMD or Intel.
		 */
		RanDraw();
		/** \brief Constructor with seed
		 *
		 * Sets the provided seed and uses the Mersenne Twister. Is CPU agnostic.
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
		/** \brief Query RNG kind
		 *
		 * Find out the kind of (P)RNG in use.
		 *
		 * \return String reflecting the RNG type
		 */
		std::string type() const {return kind_ == 'h' ? "hardware" : "mersenne_twister"; };
		/** \brief Generate random integer
		 *
		 * \return An unsigned random 64-bit integer
		 */
		uint64_t ranInt() const noexcept { return rng_->ranInt(); };
		/** \brief Sample and integer from the \f$ [0, n) \f$ interval
		 *
		 * \param[in] max the maximal value \f$n\f$ (does not appear in the sample)
		 * \return sampled value
		 */
		uint64_t sampleInt(const uint64_t &max) const noexcept { return rng_->ranInt() % max; };
		/** \brief Sample and integer from the \f$ [m, n) \f$ interval
		 *
		 * Throws `string` "Lower bound not smaller than upper bound" if \f$ m \ge n \f$.
		 *
		 * \param[in] min the minimal value \f$m\f$ (can appear in the sample)
		 * \param[in] max the maximal value \f$n\f$ (does not appear in the sample)
		 * \return sampled value
		 */
		uint64_t sampleInt(const uint64_t &min, const uint64_t &max) const;
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
		double runif() const noexcept {return 5.42101086242752217E-20 * static_cast<double>( rng_->ranInt() ); };
		/** \brief Generate a non-zero uniform deviate
		 *
		 * \return A double-precision value from the \f$ U(0,1]\f$ distribution
		 */
		double runifnz() const noexcept;
		/** \brief Generate a non-one uniform deviate
		 *
		 * \return A double-precision value from the \f$ U[0,1)\f$ distribution
		 */
		double runifno() const noexcept;
		/** \brief Generate an open-interval uniform deviate
		 *
		 * \return A double-precision value from the \f$ U(0,1)\f$ distribution
		 */
		double runifop() const noexcept;
		/** \brief A standard Gaussian deviate
		 *
		 * Generates a Gaussian random value with mean \f$ \mu = 0.0 \f$ and standard deviation \f$ \sigma = 1.0 \f$. Implemented using a version of the Marsaglia and Tsang (2000) ziggurat algorithm, modified according to suggestions in the GSL implementation of the function.
		 *
		 * \return a sample from the standard Gaussian distribution
		 */
		double rnorm() const noexcept;
		/** \brief A zero-mean Gaussian deviate
		 *
		 * Generates a Gaussian random value with mean \f$ \mu = 0.0 \f$ and standard deviation \f$ \sigma \f$. Implemented using a version of the Marsaglia and Tsang (2000) ziggurat algorithm, modified according to suggestions in the GSL implementation of the function.
		 *
		 * \param[in] sigma standard deviation
		 * \return a sample from the zero-mean Gaussian distribution
		 */
		double rnorm(const double &sigma) const noexcept { return ( this->rnorm() ) * sigma; };
		/** \brief A Gaussian deviate
		 *
		 * Generates a Gaussian random value with mean \f$ \mu \f$ and standard deviation \f$ \sigma \f$. Implemented using a version of the Marsaglia and Tsang (2000) ziggurat algorithm, modified according to suggestions in the GSL implementation of the function.
		 *
		 * \param[in] mu standard deviation
		 * \param[in] sigma standard deviation
		 * \return a sample from the Gaussian distribution
		 */
		double rnorm(const double &mu, const double &sigma) const noexcept { return mu + ( this->rnorm() ) * sigma; };
		/** \brief A standard Gamma deviate
		 *
		 * Generates a Gamma random variable with shape \f$ \alpha > 0 \f$ and standard scale \f$ \beta = 1.0 \f$. Implements the Marsaglia and Tsang (2000) method.
		 *
		 * \param[in] alpha shape parameter \f$ \alpha \f$
		 * \return a sample from the standard Gamma distribution
		 */
		double rgamma(const double &alpha) const noexcept;
		/** \brief A general Gamma deviate
		 *
		 * Generates a Gamma random variable with shape \f$ \alpha > 0 \f$ and scale \f$ \beta > 0 \f$.
		 *
		 * \param[in] alpha shape parameter \f$ \alpha \f$
		 * \param[in] beta scale parameter \f$ \beta \f$
		 * \return a sample from the general Gamma distribution
		 */
		double rgamma(const double &alpha, const double &beta) const noexcept { return beta > 0.0 ? ( this->rgamma(alpha) ) / beta : nan(""); };
		/** \brief A Dirichlet deviate
		 *
		 * Generates a vector of probabilities, given a vector of concentration parameters \f$ \alpha_K > 0 \f$.
		 *
		 * \param[in] alpha vector of concentration parameters
		 * \param[out] p vector of probabilities, must be the same length as \f$ \alpha \f$.
		 */
		void rdirichlet(const std::vector<double> &alpha, std::vector<double> &p) const;
		/** \brief A chi-square deviate
		 *
		 * Generates a \f$ \chi^2 \f$ random variable with degrees of freedom \f$ \nu > 0.0 \f$.
		 *
		 * \param[in] nu degrees of freedom
		 * \return a sample from the \f$ \chi^2 \f$ distribution
		 */
		double rchisq(const double &nu) const noexcept { return 2.0*this->rgamma(nu/2.0); };
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
		uint64_t vitterA(const double &n, const double &N) const noexcept;
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
		uint64_t vitter(const double &n, const double &N) const noexcept;
	private:
		/** \brief (Pseudo-)random number generator */
		std::unique_ptr<Generate> rng_;
		/** \brief Which generator
		 *
		 * Set to "h" if the hardware generator is supported, "m" is had to use the MT.
		 */
		char kind_;
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
	};

}

#endif /* random_hpp */

