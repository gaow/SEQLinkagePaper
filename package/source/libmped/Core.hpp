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


inline bool hasEnding(std::string const & fullString, std::string const & ending)
{
	if (fullString.length() >= ending.length()) {
		return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
	} else {
		return false;
	}
}


inline unsigned adjustSize(unsigned n, unsigned a)
{
	if (a <= 0) return n;
	if (a == 1) return a;
	div_t divresult = div(n, a);
	// reduce size by x such that rem + res * x = a - x
	return (unsigned)(a - ((a - divresult.rem) / (divresult.quot + 1)));
}


inline std::string collapse(VecString & haplotype, unsigned start, unsigned end, unsigned size)
{
	if (end == 0) end = haplotype.size();
	if (start == end) return "?";
	std::string collapsed_haplotype = "";
	unsigned wsize = adjustSize(haplotype.size(), size);
	unsigned counter = 0;
	std::string code = "1";

	for (unsigned i = start; i < end; ++i) {
		counter += 1;
		if (haplotype[i] == "2") code = "2";
		else code = (code == "?" || code == "2") ? code : haplotype[i];
		if (counter == wsize) {
			collapsed_haplotype += code;
			counter = 0;
			code = "1";
		}
	}
	// make it missing data if "?" is in the haplotype
	if (collapsed_haplotype.find("?") != std::string::npos) return "?";
	else return collapsed_haplotype;
}


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
};
}
#endif
