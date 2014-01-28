// $File: Core.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#ifndef _CORE_HPP_
#define _CORE_HPP_

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
	MendelianErrorChecker() : errorCount(0) {};
	~MendelianErrorChecker() {};
	int errorCount;
	void Apply(Pedigree & ped);

};

class GeneticHaplotyper
{
public:
	GeneticHaplotyper(std::string chrom) : __chrom(chrom), data(0) {}
	~GeneticHaplotyper() {};
	// [family][sample][haplotypes]
	std::vector< std::vector< std::vector<std::string> > > data;
	// Apply haplotyping. Missing data are imputed as possible
	void Apply(Pedigree & ped);

private:
	std::string __chrom;
};

class HaplotypeCoder
{
public:
	HaplotypeCoder(int size) : __size(size), recombCount(0), data(0) {}
	~HaplotypeCoder() {};
	// [[familyid, sampleid, hap1, hap2] ...]
	std::vector< std::vector<std::string> > data;
	int recombCount;
	void Apply(std::vector< std::vector< std::vector<std::string> > > & ghdata);

private:
	int __size;
	unsigned __AdjustSize(int n);

	std::string __Collapse(std::string & haplotype);

};
}
#endif
