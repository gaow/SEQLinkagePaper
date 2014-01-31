// $File: Core.cpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#include "Core.hpp"
#include <iostream>
#include <cstdlib>
#include <set>
#include <algorithm>

inline bool hasEnding(std::string const & fullString, std::string const & ending)
{
	if (fullString.length() >= ending.length()) {
		return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
	} else {
		return false;
	}
}


void SEQLinco::DataLoader::LoadVariants(Pedigree & ped,
                                        const VecString & names,
                                        const VecString & positions,
                                        const std::string & chrom,
                                        double positionAdjustment)
{
	VecInt vs(names.size());

	for (unsigned i = 0; i < names.size(); ++i) {
		ped.pd.columnHash.Push(ped.GetMarkerID(names[i].c_str()));
		ped.pd.columns.Push(1);
		ped.pd.columnCount++;
		MarkerInfo * info = ped.GetMarkerInfo(i);
		info->chromosome = (chrom == "X" || chrom == "x") ? 999 : atoi(chrom.c_str());
		int position = atoi(positions[i].c_str());
		if (std::find(vs.begin(), vs.end(), position) != vs.end()) position++;
		vs[i] = position;
		info->positionFemale = info->positionMale = info->position = position * positionAdjustment;
		// std::clog << info->freq.dim << std::endl;
	}
}


void SEQLinco::DataLoader::LoadSamples(Pedigree & ped, const VecVecString & samples)
{

	for (unsigned i = 0; i < samples.size(); ++i) {
		VecString fam_info(samples[i].begin(), samples[i].begin() + 5);
		VecString genotypes(samples[i].begin() + 5, samples[i].end());
		__AddPerson(ped, fam_info, genotypes);
	}
	ped.Sort();
	SortFamilies(ped);
}


void SEQLinco::DataLoader::__AddPerson(Pedigree & ped, VecString & fam_info, VecString & genotypes)
{
	// add person info
	bool sex_failure = false;

	ped.AddPerson(fam_info[0].c_str(), fam_info[1].c_str(),
		fam_info[2].c_str(), fam_info[3].c_str(),
		ped.TranslateSexCode(fam_info[4].c_str(), sex_failure));
	// add person genotypes
	for (unsigned i = 0; i < genotypes.size(); ++i) {
		String c1 = genotypes[i].c_str()[0];
		String c2 = genotypes[i].c_str()[1];
		Alleles new_genotype;
		new_genotype[0] = ped.LoadAllele(ped.GetMarkerInfo(i), c1);
		new_genotype[1] = ped.LoadAllele(ped.GetMarkerInfo(i), c2);
		if (new_genotype.isKnown()) ped[ped.count - 1].markers[i] = new_genotype;
	}
}


void SEQLinco::GeneticHaplotyper::Apply(Pedigree & ped)
{
	String chrom = __chrom.c_str();

	if (chrom == "X") ped.chromosomeX = true;
	//
	ped.EstimateFrequencies(0, true);
	// recode alleles so more frequent alleles have lower allele numbers internally
	ped.LumpAlleles(0.0);
	// remove uninformative family or individuals
	// !! Do not trim here, because if a family is uninformative we can report as is
	// ped.Trim(true);
	FamilyAnalysis engine(ped);
	// activate haplotyping options
	engine.bestHaplotype = true;
	engine.zeroRecombination = false;
	engine.SetupGlobals();
	engine.SetupMap(chrom);
	for (int i = 0; i < ped.familyCount; i++) {
		if (engine.SelectFamily(ped.families[i])) {
			engine.Analyse();
			data.push_back(engine.hapOutput);
		}
	}
	engine.CleanupGlobals();
}


void SEQLinco::GeneticHaplotyper::Print()
{
	for (unsigned f = 0; f < data.size(); f++) {
		for (unsigned p = 0; p < data[f].size(); p++) {
			for (unsigned i = 0; i < data[f][p].size(); i++) {
				std::clog << data[f][p][i] << "\t";
			}
			std::clog << std::endl;
		}
		std::clog << std::endl;
	}
	std::clog << std::endl;
}


void SEQLinco::MendelianErrorChecker::Apply(Pedigree & ped)
{
	// check mendelian error for everyone's every marker in input ped object
	for (int i = 0; i < ped.count; i++) {
		// skip founder
		if (ped[i].isFounder()) continue;
		// identify founders for this person
		Person * mom = ped[i].mother;
		Person * dad = ped[i].father;
		for (int m = 0; m < ped.markerCount; m++) {
			//
			// genotype data missing for both founders
			//
			if (!mom->markers[m].isKnown() && !dad->markers[m].isKnown()) continue;
			//
			// if mother/father missing, substitute with uninformative marker
			// otherwise use as is
			//
			int gdad1, gdad2, gmom1, gmom2;
			if (!mom->markers[m].isKnown()) {
				gmom1 = 1; gmom2 = 2;
			} else {
				gmom1 = mom->markers[m][0]; gmom2 = mom->markers[m][1];
			}
			if (!dad->markers[m].isKnown()) {
				gdad1 = 1; gdad2 = 2;
			} else {
				gdad1 = dad->markers[m][0]; gdad2 = dad->markers[m][1];
			}
			//
			// check for mendelian error
			//
			// person missing data
			if (!ped[i].markers[m].isKnown()) {
				if (dad->markers[m].isHomozygous() && mom->markers[m].isHomozygous() && gdad1 == gmom1) {
					ped[i].markers[m][0] = gdad1; ped[i].markers[m][1] = gdad2; continue;
				}
			}
			// no error
			if (((ped[i].markers[m][0] == gdad1 || ped[i].markers[m][0] == gdad2) &&
			     (ped[i].markers[m][1] == gmom1 || ped[i].markers[m][1] == gmom2)) ||
			    ((ped[i].markers[m][1] == gdad1 || ped[i].markers[m][1] == gdad2) &&
			     (ped[i].markers[m][0] == gmom1 || ped[i].markers[m][0] == gmom2)))
				continue;
			// error found, make missing
			else {
				errorCount += 1;
				ped[i].markers[m][0] = ped[i].markers[m][1] = 0;
			}
		}
	}
}


void SEQLinco::HaplotypeCoder::Apply(VecVecVecString & haploVecs)
{
	// each element of haploVecsis a family's data
	// each element of haploVecs[i] is a haplotype with the first 2 items being fid and sid
	if (!haploVecs.size()) return;
	for (unsigned f = 0; f < haploVecs.size(); f++) {
		//
		// record recombination and adjust data format
		//
		// the nearest recombination event index
		unsigned minRecombPos = 999999999;
		for (unsigned p = 0; p < haploVecs[f].size(); p++) {
			for (unsigned i = 2; i < haploVecs[f][p].size(); i++) {
				// recombination event detected
				if (!hasEnding(haploVecs[f][p][i], ":") && !hasEnding(haploVecs[f][p][i], "|")) {
					recombCount += 1;
					minRecombPos = (i < minRecombPos) ? i : minRecombPos;
				}
				// use one of the likely haplotype configuration
				haploVecs[f][p][i] = (haploVecs[f][p][i].substr(0, 1) == "A") ? haploVecs[f][p][i].substr(1, 1) : haploVecs[f][p][i].substr(0, 1);
			}
		}
		//
		// collapse haplotype vectors to data object as a string
		//
		std::set<std::string> pool;
		unsigned dataStart = data.size();
		for (unsigned p = 0; p < haploVecs[f].size(); p++) {
			if (!data.size() || (data.back()[0] != haploVecs[f][p][0] || data.back()[1] != haploVecs[f][p][1])) {
				VecString newperson(haploVecs[f][p].begin(), haploVecs[f][p].begin() + 2);
				data.push_back(newperson);
			}
			std::string haplotype = __Collapse(haploVecs[f][p], 2,
				(minRecombPos < haploVecs[f][p].size()) ? minRecombPos : (unsigned)haploVecs[f][p].size());
			data[data.size() - 1].push_back(haplotype);
			if (haplotype != "?") pool.insert(haplotype);
		}
		//
		// convert haplotype super strings to haplotype patterns
		//
		VecString patterns(pool.begin(), pool.end());
		std::sort(patterns.begin(), patterns.end());
		for (unsigned p = dataStart; p < data.size(); p++) {
			for (unsigned i = 2; i < 4; i++) {
				if (data[p][i] == "?") data[p][i] = "0";
				else data[p][i] = std::to_string(std::find(patterns.begin(), patterns.end(), data[p][i]) - patterns.begin() + 1);
			}
		}
	}
}


void SEQLinco::HaplotypeCoder::Print()
{
	for (unsigned p = 0; p < data.size(); p++) {
		for (unsigned i = 0; i < data[p].size(); i++) {
			std::clog << data[p][i] << "\t";
		}

		std::clog << std::endl;
	}
	std::clog << std::endl;

}


unsigned SEQLinco::HaplotypeCoder::__AdjustSize(unsigned n)
{
	if (__size <= 0) return n;
	if (__size == 1) return __size;
	div_t divresult;
	divresult = div(n, __size);
	// reduce size by x such that rem + res * x = size - x
	return (unsigned)(__size - ((__size - divresult.rem) / (divresult.quot + 1)));
}


std::string SEQLinco::HaplotypeCoder::__Collapse(VecString & haplotype, unsigned start, unsigned end)
{
	if (end == 0) end = haplotype.size();
	if (start == end) return "?";
	std::string collapsed_haplotype = "";
	unsigned wsize = __AdjustSize(haplotype.size());
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


