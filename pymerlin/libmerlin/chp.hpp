// $File: chp.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#ifndef _CHP_HPP_
#define _CHP_HPP_

#include "Pedigree.h"
#include "MerlinFamily.h"
#include "MerlinHaplotype.h"
#include "MerlinSort.h"

#include <algorithm>
#include <vector>
#include <string>
#include <iterator>
#include <iostream>

namespace SEQLinco {
class PedigreeData
{
public:
	PedigreeData() {};
	~PedigreeData() {};
	Pedigree data;
	void LoadVariants(std::vector<std::string> & names,
		std::vector<int> & positions,
		std::string chrom);

	void LoadSamples(std::vector< std::vector<std::string> > & samples);

private:
	void __AddPerson(std::vector<std::string> & fam_info,
		std::vector<std::string> & genotypes);

};

class MendelianErrorChecker
{
public:
	MendelianErrorChecker() {};
	~MendelianErrorChecker() {};
	int errorCount;
	void Apply(Pedigree & ped);

};

class GeneticHaplotyper
{
public:
	GeneticHaplotyper(std::string chrom) : __chrom(chrom) {}
	~GeneticHaplotyper() {};
	// [family][sample][haplotypes]
	std::vector< std::vector< std::vector<std::string> > > data;
	void Apply(Pedigree & ped);

private:
	std::string __chrom;
};

class HaplotypeCoder
{
public:
	HaplotypeCoder(int size) : __size(size) {}
	~HaplotypeCoder() {};
	// [family][sample][haplotypes]
	std::vector< std::vector< std::vector<std::string> > > data;
	void Apply(std::vector< std::vector< std::vector<std::string> > > & ghdata);

private:
	int __size;
	int __AdjustSize(int n);

};
}
#endif
