/**
 ** This file is based on the Consensus Assessment of Sequence And VAriation
 ** (CASAVA) software package.
 **
 ** \file ./CASAVA_v1.8.2/src/c++/lib/demultiplex/MaskQvalsByEamss.cpp
 **
 ** \brief Masks quality values of a read, using the EAMSS.
 **
 ** \author Come Raczy
 **
 ** \rewrite Eric Biggers
 **
 ** \changed Use vectors instead fo strings; also got rid of the boost header.
 **/

#include <string>
#include <vector>
#include <string.h>
#include <assert.h>

using std::string;
using std::vector;
using std::min;
using std::max;
#include "MaskQvalsByEamss.h"

namespace casava
{
namespace demultiplex
{

const char *get_quality_mask_mode_name(enum QualityMaskMode mode) 
{
	switch (mode) {
	case MODE_QUALITY:
		return "quality";
	case MODE_LOWERCASE:
		return "lowercase";
	default:
		return "None";
	}
}

const int MaskQvalsByEamss::lowScore = 1;
const int MaskQvalsByEamss::mediumScore = 0;
const int MaskQvalsByEamss::highScore = -2;
//const char MaskQvalsByEamss::mediumThreshold = 'O';
const char MaskQvalsByEamss::mediumThreshold = 15;
//const char MaskQvalsByEamss::highThreshold = '^';
const char MaskQvalsByEamss::highThreshold = 30;
const int MaskQvalsByEamss::minScore = 1;
const std::vector<std::string> MaskQvalsByEamss::motifList = std::vector<std::string>(0);

/**
 ** \brief Mask quality using the EAMSS algorithm
 **/
int MaskQvalsByEamss::operator()(vector<char> &qValues,
				 vector<char> &baseCalls,
				 enum QualityMaskMode mode) const
{
	std::pair<int, int> tmp = eamss(qValues);
	const int score = tmp.first;
	int position = tmp.second;

	if (score < minScore)
		return 0;

	// look for troublemaker motifs
	int extendedPosition = position;
	std::vector<std::string>::const_iterator motif;
	for (motif = motifList.begin(); motifList.end() != motif; ++motif) {
		const int troublemakerPosition = findStr(&baseCalls[0], *motif,
							 max(0, position - 15),
							 min(baseCalls.size(), (size_t)position));
		if (troublemakerPosition != -1)
			extendedPosition = min(extendedPosition, troublemakerPosition);
	}
	position = min(position, extendedPosition);

	// look for extended run of polyG, with at least 90% G bases
	unsigned int numG = 0;
	unsigned int numBases = 0;
	bool maskPolyG = false;
	int maskStart = position;
	for (int curPos = position; 0 <= curPos; --curPos) {
		if ('G' == baseCalls[curPos])
			++numG;

		if (++numBases < 10)
			continue;
		const float gFrac = static_cast<float>(numG)/static_cast<float>(numBases);
		if (gFrac >= 0.9f && 'G' == baseCalls[curPos]) {
			// only start masking at G
			maskPolyG = true;
			maskStart = curPos;
		}
		else if (gFrac < 0.9f)
			break;
	}
	if (maskPolyG)
		position = maskStart;
	
	switch (mode) {
	case MODE_QUALITY:
		for (size_t idx = position; idx < baseCalls.size(); ++idx)
			qValues[idx] = 2; // will become 'B' or '#'
		break;
	case MODE_LOWERCASE:
		for (size_t idx = position; idx < baseCalls.size(); ++idx)
			baseCalls[idx] += ('a' - 'A'); // to lowercase
		break;
	default:
		assert(0);
	}

	return baseCalls.size() - position;
	
}

/**
 ** \brief Compute the position with the highest EAMSS score
 **
 ** \return the pair (bestScore, bestPosition) or (highScore-1, -1) if qValues
 ** is empty
 **/
std::pair<int, int> MaskQvalsByEamss::eamss(const vector<char>& qValues) const
{
	int curScore = 0;
	// initialize the bestscore to something lower than the first value of curScore
	int bestScore = min(min(highScore, mediumScore), lowScore) - 1;
	int bestPosition = -1;
	for (int idx = qValues.size() - 1; 0 <= idx; --idx) {
		if (qValues[idx] >= highThreshold)
			curScore += highScore;
		else if (qValues[idx] >= mediumThreshold)
			curScore += mediumScore;
		else
			curScore += lowScore;
		if (curScore >= bestScore) {
			bestScore = curScore;
			bestPosition = idx;
		}
	}
	return std::pair<int, int>(bestScore, bestPosition);
}

int MaskQvalsByEamss::findStr(const char targetString[],
			      const std::string &queryString,
			      int start, int stop) const
{
	void *p = memmem(targetString + start, stop - start,
			 queryString.c_str(), queryString.size());
	if (p)
		return (const char*)p - targetString;
	else
		return -1;
}

} // namespace demultiplex
} // namespace casava
