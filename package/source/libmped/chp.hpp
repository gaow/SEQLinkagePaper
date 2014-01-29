// $File: chp.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#ifndef _CHP_HPP_
#define _CHP_HPP_

#include <string>
#include "Core.hpp"
namespace SEQLinco {
class CHP
{
public:
	CHP(const int wsize, const double padj = 0.01) : __windowSize(wsize), __positionAdjustment(padj) {};
	~CHP() {};
	CHP * clone() const { return new CHP(*this); }
	VecVecString Apply(const std::string & chrom, const VecString & marker_names,
	                   const VecString & marker_positions, const VecVecString & samples)
	{
		// input chrom must be 1 .. 22 and X
		// input samples follows PED convention, e.g., { "1", "1", "0", "0", "1", "21", "21", "21" }
		// first 5 cols are fid, sid, pid, mid and sex (0/1/2 coding), followed by genotypes (1/2 coding, 0 for missing)
		PedigreeData ped;

		ped.LoadVariants(marker_names, marker_positions, chrom, __positionAdjustment);
		ped.LoadSamples(samples);
		MendelianErrorChecker mc;
		mc.Apply(ped.data);
		__mendelianErrorCount = mc.errorCount;
		GeneticHaplotyper gh(chrom);
		gh.Apply(ped.data);
		HaplotypeCoder hc(__windowSize);
		hc.Apply(gh.data);
		__recombCount = hc.recombCount;
		return hc.data;
	}


	int countMendelianErrors() { return __mendelianErrorCount; }
	int countRecombs() { return __recombCount; }

private:
	int __mendelianErrorCount;
	int __recombCount;
	int __windowSize;
	double __positionAdjustment;
};
}
#endif
