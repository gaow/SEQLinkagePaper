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

#include "Exception.hpp"

namespace SEQLinco {
typedef std::vector<int> VecInt;
typedef std::vector<std::string> VecString;
typedef std::vector<std::vector<std::string> > VecVecString;
typedef std::vector<std::vector<std::vector<std::string> > > VecVecVecString;

class DataLoader
{
public:
	DataLoader() {};
	~DataLoader() {};
	DataLoader * clone() const { return new DataLoader(*this); }
	void LoadVariants(Pedigree & ped,
		const VecString & names,
		const VecString & positions,
		const std::string & chrom,
		double positionAdjustment = 0.01);

	void LoadSamples(Pedigree & ped,
		const VecVecString & samples);

private:
	void __AddPerson(Pedigree & ped,
		VecString & fam_info,
		VecString & genotypes);

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

	void Print();

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
	void Apply(VecVecVecString & haploVecs);

	void Print();

private:
	int __size;
	unsigned __AdjustSize(unsigned n);

	std::string __Collapse(VecString & haplotype, unsigned start = 0, unsigned end = 0);

};
}
#endif
