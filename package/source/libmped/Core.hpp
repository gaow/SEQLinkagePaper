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
typedef std::vector<std::string> VecString;
typedef std::vector<std::vector<std::string> > VecVecString;
typedef std::vector<std::vector<std::vector<std::string> > > VecVecVecString;

class PedigreeData
{
public:
	PedigreeData() {};
	~PedigreeData() {};
	PedigreeData * clone() const { return new PedigreeData(*this); }
	Pedigree data;
	void LoadVariants(const VecString & names,
		const VecString & positions,
		const std::string & chrom,
		double positionAdjustment = 0.01);

	void LoadSamples(const VecVecString & samples);

private:
	void __AddPerson(VecString & fam_info, VecString & genotypes);

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
	VecVecVecString data;
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
	VecVecString data;
	int recombCount;
	void Apply(VecVecVecString & ghdata);

private:
	int __size;
	unsigned __AdjustSize(int n);

	std::string __Collapse(std::string & haplotype);

};
}
#endif
