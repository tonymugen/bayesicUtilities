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

/// Connect lines with populations
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2018 Anthony J. Greenberg
 * \version 1.1
 *
 * Implementation of a class that relates individuals to groups, similar to an factor in R.
 *
 */

#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>

#include "index.hpp"

using namespace BayesicSpace;

/*
 * Index methods
 */

Index::Index(const size_t &Ngroups) {
	index_.resize(Ngroups);
}

Index::Index(std::vector<size_t>::const_iterator groupVecBegin, std::vector<size_t>::const_iterator groupVecEnd) {
	groupVal_.reserve( static_cast<size_t>( std::distance(groupVecBegin, groupVecEnd) ) );
	size_t elInd{0};
	std::for_each(
		groupVecBegin,
		groupVecEnd,
		[this, &elInd](size_t currentGroupID) {
			groupVal_.push_back(currentGroupID);
			if ( currentGroupID >= index_.size() ) {
				index_.resize(currentGroupID + 1); // maybe this groupID is a few steps ahead of the current
				index_[currentGroupID].push_back(elInd);
			} else {
				index_[currentGroupID].push_back(elInd);
			}
			++elInd;
		}
	);
	// now remove all empty elements from index_
	index_.erase(
		std::remove_if(
			index_.begin(), index_.end(), [](const std::vector<size_t> &idxVector){return idxVector.empty();}
		),
		index_.end()
	);
	index_.shrink_to_fit();
}

Index::Index(const std::string &inFileName) {
	std::ifstream idxFl(inFileName);

	if (!idxFl) {
		throw std::string("ERROR: Cannot open file ") + inFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	int tmpIn{0};
	size_t iEl{0};
	while (idxFl >> tmpIn){
		if (tmpIn < 0) {
			idxFl.close();
			throw std::string("ERROR: Negative group ID in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
		}
		const auto tmpUI = static_cast<size_t>(tmpIn);
		groupVal_.push_back(tmpUI);
		if ( tmpUI >= index_.size() ) {
			index_.resize(tmpUI + 1);                       // maybe this goupID is a few steps ahead of the current
			index_[tmpUI].push_back(iEl);
			++iEl;
		}
		else {
			index_[tmpUI].push_back(iEl);
			++iEl;
		}
	}
	idxFl.close();
	// now remove all empty elements from index_
	index_.erase(
		std::remove_if(
			index_.begin(), index_.end(), [](const std::vector<size_t> &idxVector){return idxVector.empty();}
		),
		index_.end()
	);
}

size_t Index::neGroupNumber() const noexcept {
	return static_cast<size_t>( std::count_if(index_.cbegin(), index_.cend(), [](const std::vector<size_t> &group){return !group.empty();}) );
}

void Index::update(std::vector<size_t>::const_iterator newGrpVecBegin, std::vector<size_t>::const_iterator newGrpVecEnd) {
	std::for_each(index_.begin(), index_.end(), [](std::vector<size_t> &groupVec){groupVec.clear();});

	size_t elInd{0};
	std::for_each(
		newGrpVecBegin,
		newGrpVecEnd,
		[this, &elInd](size_t groupID) {
			assert( ( groupID < index_.size() ) && "ERROR: updating vector has an extra group ID in Index::update()" );
			index_[groupID].push_back(elInd);
			++elInd;
		}
	);
	groupVal_.resize( static_cast<size_t>( std::distance(newGrpVecBegin, newGrpVecEnd) ) );
	std::copy( newGrpVecBegin, newGrpVecEnd, groupVal_.begin() );
}
