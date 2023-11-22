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

/// Numerical utilities
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2020 -- 2022 Anthony J. Greenberg
 * \version 1.0
 *
 * Function definitions for a set of numerical utilities.
 *
 */

#pragma once

#include <cmath>
#include <vector>
#include <array>

namespace BayesicSpace {
	/** \brief Value and its weight
	 *
	 * Set for calculating weighted means.
	 */
	struct ValueWithWeight {
		double value;
		double weight;
	};
	/** \brief Swap two `size_t` values
	 *
	 * Uses the three XORs trick to swap two integers. Safe if the variables happen to refer to the same address.
	 *
	 * \param[in,out] first first integer
	 * \param[in,out] second second integer
	 */
	void swapXOR(size_t &first, size_t &second) noexcept;
	/** \brief Logit function
	 *
	 * \param[in] probability probability in the (0, 1) interval
	 * \return logit transformation
	 */
	[[nodiscard]] double logit(const double &probability) noexcept;
	/** \brief Logistic function
	 *
	 * There is a guard against under- and overflow: the function returns 0.0 for \f$ x \le -35.0\f$ and 1.0 for \f$x \ge 35.0\f$.
	 *
	 * \param[in] input value to be projected to the (0, 1) interval
	 * \return logistic transformation
	 */
	[[nodiscard]] double logistic(const double &input) noexcept;
	/** \brief Logarithm of the Gamma function
	 *
	 * The log of the \f$ \Gamma(x) \f$ function for real \f$ x > 0 \f$ (non-positive input values produce `NaN`).
	 * Implementing the Lanczos algorithm following GSL.
	 *
	 * \param[in] input input, passed by value
	 * \return \f$ \log \Gamma(input) \f$
	 *
	 */
	[[nodiscard]] double lnGamma(double input) noexcept;
	/** \brief Digamma function
	 *
	 * Defined only for \f$ x > 0 \f$, will return `NaN` otherwise. Adopted from the `dpsifn` function in R.
	 *
	 * \param[in] input function argument (must be positive)
	 * \return value of the digamma function
	 */
	[[nodiscard]] double digamma(const double &input) noexcept;
	/** \brief Vector self-dot-product
	 *
	 * \param[in] begin first element iterator
	 * \param[in] end one past the last element iterator
	 * \return dot-product value
	 */
	[[nodiscard]] double dotProd(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end) noexcept;
	/** \brief Dot-product of two `double` vectors
	 *
	 * Throws if the second vector range goes past its end (the second vector reference is needed for the check).
	 *
	 * \param[in] firstBegin first vector start iterator
	 * \param[in] firstEnd first vector end iterator 
	 * \param[in] secondBegin second vector start iterator
	 * \param[in] second second vector
	 * \return dot-product value
	 */
	[[nodiscard]] double dotProd(std::vector<double>::const_iterator firstBegin,
								std::vector<double>::const_iterator firstEnd,
								std::vector<double>::const_iterator secondBegin,
								const std::vector<double> &second);
	/** \brief Mean of a vector of `double`
	 *
	 * Uses the numerically stable recursive algorithm.
	 *
	 * \param[in] begin iterator to first element
	 * \param[in] end iterator to one past the last element
	 * \return mean value
	 */
	[[nodiscard]] double stableMean(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end) noexcept;
	[[nodiscard]] double stupidMean(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end) noexcept;
	/** \brief Weighted mean update
	 *
	 * Takes the current weighted mean and updates using the new data point and weight. The formula is
	 *
	 * \f$
	 *     \bar{\mu}_n = \frac{\bar{\mu}_{n-1}\sum_{i=1}^{n-1}w_i + w_n x_n}{\sum_{i=1}^{n-1}w_i + w_n}
	 * \f$
	 *
	 * \param[in] nextDataPoint new point \f$ x_n \f$ and \f$ w_n \f$
	 * \param[in] currentMean mean and sum of weights up to the next point
	 * \return new mean and sum of weights
	 */
	 [[nodiscard]] ValueWithWeight updateWeightedMean(const ValueWithWeight &nextDataPoint, const ValueWithWeight &currentMean) noexcept;
}

