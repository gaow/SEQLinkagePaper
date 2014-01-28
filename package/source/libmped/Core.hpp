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

#include "Exception.hpp"

namespace SEQLinco {
class PedigreeData
{
public:
	PedigreeData() {};
	~PedigreeData() {};
	PedigreeData * clone() const { return new PedigreeData(*this); }
	Pedigree data;
	void LoadVariants(const std::vector<std::string> & names,
		const std::vector<int> & positions,
		const std::string & chrom);

	void LoadSamples(const std::vector< std::vector<std::string> > & samples);

private:
	void __AddPerson(std::vector<std::string> & fam_info,
		std::vector<std::string> & genotypes);

};

class MendelianErrorChecker
{
public:
	MendelianErrorChecker() : errorCount(0) {};
	~MendelianErrorChecker() {};
	MendelianErrorChecker * clone() const { return new MendelianErrorChecker(*this); }
	int errorCount;
	void Apply(Pedigree & ped);

};

class GeneticHaplotyper
{
public:
	GeneticHaplotyper(const std::string & chrom) : data(0), __chrom(chrom) {}
	~GeneticHaplotyper() {};
	GeneticHaplotyper * clone() const { return new GeneticHaplotyper(*this); }
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
	HaplotypeCoder(const int size) : data(0), recombCount(0), __size(size) {}
	~HaplotypeCoder() {};
	HaplotypeCoder * clone() const { return new HaplotypeCoder(*this); }
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
