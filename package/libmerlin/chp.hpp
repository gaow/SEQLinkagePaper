// $File: chp.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#ifndef _CHP_HPP_
#define _CHP_HPP_

#include <vector>
#include <string>
#include "Core.hpp"
using namespace SEQLinco;

class CHP
{
public:
	CHP(int wsize) : __windowSize(wsize) {};
	~CHP() {};
	std::vector< std::vector<std::string> >
	Apply(std::string chrom, const std::vector<std::string> & marker_names,
	      const std::vector<int> & marker_positions, const std::vector< std::vector<std::string> > & samples)
	{
		// input chrom must be 1 .. 22 and X
		// input samples follows PED convention, e.g., { "1", "1", "0", "0", "1", "21", "21", "21" }
		// first 5 cols are fid, sid, pid, mid and sex (0/1/2 coding), followed by genotypes (1/2 coding, 0 for missing)
		PedigreeData ped;

		ped.LoadVariants(marker_names, marker_positions, chrom);
		ped.LoadSamples(samples);
		MendelianErrorChecker mc;
		mc.Apply(ped.data);
		__mendelianErrorCount = mc.errorCount;
		GeneticHaplotyper gh(chrom);
		gh.Apply(ped.data);
		HaplotypeCoder hc(__windowSize);
		hc.apply(gh.data);
		__recombCount = hc.recombCount;
		return hc.data;
	}


	int countMendelianErrors() { return __mendelianErrorCount; }
	int countRecombs() { return __recombCount; }

private:
	int __mendelianErrorCount;
	int __recombCount;
	int __windowSize;
};
#endif
