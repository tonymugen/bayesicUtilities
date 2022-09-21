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
 * Class definition for a set of numerical utilities.
 * Implemented as a class because this seems to be the only way for these methods to be included using Rcpp with no compilation errors.
 *
 */

#pragma once

#include <math.h>
#include <vector>

namespace BayesicSpace {
	/** \brief Numerical utilities collection
	 *
	 * Implements numerical functions for use throughout the project.
	 *
	 */
	class NumerUtil {
	public:
	/** \brief Swap two `size_t` values
	 *
	 * Uses the three XORs trick to swap two integers. Safe if the variables happen to refer to the same address.
	 *
	 * \param[in,out] i first integer
	 * \param[in,out] j second integer
	 */
	void swapXOR(size_t &i, size_t &j) const noexcept;
	/** \brief Logit function
	 *
	 * \param[in] p probability in the (0, 1) interval
	 * \return logit transformation
	 */
	double logit(const double &p) const noexcept { return log(p) - log(1.0 - p); }
	/** \brief Logistic function
	 *
	 * There is a guard against under- and overflow: the function returns 0.0 for \f$ x \le -35.0\f$ and 1.0 for \f$x \ge 35.0\f$.
	 *
	 * \param[in] x value to be projected to the (0, 1) interval
	 * \return logistic transformation
	 */
	double logistic(const double &x) const noexcept;
	/** \brief Logarithm of the Gamma function
	 *
	 * The log of the \f$ \Gamma(x) \f$ function. Implementing the Lanczos algorithm following Numerical Recipes in C++.
	 *
	 * \param[in] x value
	 * \return \f$ \log \Gamma(x) \f$
	 *
	 */
	double lnGamma(const double &x) const noexcept;
	/** \brief Digamma function
	 *
	 * Defined only for \f$ x > 0 \f$, will return _NaN_ otherwise. Adopted from the `dpsifn` function in R.
	 *
	 * \param[in] x function argument (must be positive)
	 * \return value of the digamma function
	 */
	double digamma(const double &x) const noexcept;
	/** \brief Vector self-dot-product
	 *
	 * \param[in] v vector
	 * \return dot-product value
	 */
	double dotProd(const std::vector<double> &v) const noexcept;
	/** \brief Dot-product of two vectors
	 *
	 * \param[in] v1 vector 1
	 * \param[in] v2 vector 2
	 * \return dot-product value
	 */
	double dotProd(const std::vector<double> &v1, const std::vector<double> &v2) const noexcept;
	/** \brief Weighted mean update
	 *
	 * Takes the current weighted mean and updates using the new data point and weight. The formula is
	 *
	 * \f$
	 *     \bar{\mu}_n = \dfrac{\bar{\mu}_{n-1}\sum_{i=1}^{n-1}w_i + w_n x_n}{\sum_{i=1}^{n-1}w_i + w_n}
	 * \f$
	 *
	 * \param[in] xn new point \f$ x_n \f$
	 * \param[in] wn weight \f$ w_n \f$
	 * \param[out] mu new mean
	 * \param[out] w new weight
	 *
	 */
	 void updateWeightedMean(const double &xn, const double &wn, double &mu, double &w) const noexcept;
	/** \brief Mean of an array
	 *
	 * Uses the numerically stable recursive algorithm.
	 *
	 * \param[in] arr c-style array of values
	 * \param[in] len array length
	 * \return mean value
	 */
	double mean(const double arr[], const size_t &len) const noexcept;

	private:
		/** \brief Gamma function magical coefficients */
		static const double gCoeff_[14];
		/** \brief Bernoulli numbers */
		static const double bvalues_[22];
	};
}

