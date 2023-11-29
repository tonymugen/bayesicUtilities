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

/// Connect lines with groups
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2017 -- 2022 Anthony J. Greenberg
 * \version 1.1
 *
 * Definitions and interface documentation for a class that relates individuals to groups, similar to an factor in R.
 *
 */

#pragma once

#include <vector>
#include <string>

namespace BayesicSpace {
	/** \brief Group index
	 *
	 * For each group, contains indexes of the lines that belong to it. Can also identify the group a given element belongs to.
	 * Group numbers need not be consecutive. Although group IDs are assumed to be base-0, everything should work even if they are not.
	 */
	class Index {
	public:
		/** \brief Default constructor */
		Index() = default;
		/** \brief Group constructor
		 *
		 * Sets up empty groups.
		 *
		 * \param[in] Ngroups number of groups to set up
		 */
		Index(const size_t &nGroups);
		/** \brief Vector iterator constructor
		 *
		 * The input vector has an element for each line, and the value of that element is the base-0 group ID (i.e., if line _n_ is in the first group, then `vec[n] == 0`).
		 *
		 * \param[in] groupVecBegin start iterator of the group ID vector
		 * \param[in] groupVecEnd end iterator of the group ID vector
		 */
		Index(std::vector<size_t>::const_iterator groupVecBegin, std::vector<size_t>::const_iterator groupVecEnd);
		/** \brief File read constructor
		 *
		 * The input file has an entry for each line (separated by white space), and the value of that entry is the base-0 group ID.
		 * If the file cannot be opened, throws "Cannot open file file_name". If a negative group value is detected, thorws "Negative group ID".
		 *
		 * \param[in] inFileName input file name
		 */
		Index(const std::string &inFileName);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy `Index` to be copied
		 */
		Index(const Index &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to be copied
		 * \return an `Index` object
		 */
		Index &operator=(const Index &toCopy) = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove `Index` object to be moved
		 */
		Index(Index &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to be moved
		 * \return an `Index` object
		 */
		Index &operator=(Index &&toMove) noexcept = default;
		/** \brief Destructor */
		~Index() = default;

		/** \brief Vector subscript operator
		 *
		 * Returns the index of group _i_.
		 *
		 * \param[in] index group index
		 * \return index of line IDs
		 */
		[[nodiscard]] const std::vector<size_t> & operator[] (const size_t &index) const { return index_[index]; };

		/** \brief Group size
		 *
		 * \param[in] index group index
		 * \return size of the _i_th group
		 */
		[[nodiscard]] size_t groupSize(const size_t &index) const {return index_[index].size(); };

		/** \brief Total sample size
		 *
		 * \return total sample size
		 */
		[[nodiscard]] size_t size() const noexcept {return groupVal_.size(); };

		/** \brief Number of groups
		 *
		 * \return number of groups
		 */
		[[nodiscard]] size_t groupNumber() const noexcept {return index_.size(); };

		/** \brief Number of non-empty groups
		 *
		 * \return number of non-empty groups
		 */
		[[nodiscard]] size_t neGroupNumber() const noexcept;

		/** \brief Group ID
		 *
		 * Returns the group ID for a given individual.
		 *
		 * \param[in] ind index of an individual
		 *
		 * \return group ID
		 */
		[[nodiscard]] size_t groupID(const size_t &ind) const {return groupVal_[ind]; };

		/** \brief Update the index
		 *
		 * Updates the groups with a new index. If a group is not present in the new vector, it is left empty but still exists.
		 * The new vector must not contain any new groups. This is checked only in debug mode via `assert()`.
		 *
		 * \param[in] newGrpVecBegin start iterator of the group ID vector
		 * \param[in] newGrpVecEnd end iterator of the group ID vector
		 */
		void update(std::vector<size_t>::const_iterator newGrpVecBegin, std::vector<size_t>::const_iterator newGrpVecEnd);

	private:
		/** \brief Vector of index vectors
		 *
		 * The outside vector is the same length as the number of groups. Each inside vector has the line indexes.
		 */
		std::vector< std::vector<size_t> > index_;
		/** \brief Vector of group IDs
		 *
		 * Each element of the vector stores the corresponding (base-0) group ID.
		 */
		std::vector<size_t> groupVal_;

	};

}

