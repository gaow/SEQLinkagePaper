// $File: Core.cpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#include "Core.hpp"
#include <iostream>
#include <cstdlib>
#include <set>
inline bool hasEnding(std::string const & fullString, std::string const & ending)
{
	if (fullString.length() >= ending.length()) {
		return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
	} else {
		return false;
	}
}


void SEQLinco::PedigreeData::LoadVariants(const VecString & names,
                                          const VecString & positions,
                                          const std::string & chrom,
                                          double positionAdjustment)
{
	for (unsigned i = 0; i < names.size(); ++i) {
		data.pd.columnHash.Push(data.GetMarkerID(names[i].c_str()));
		data.pd.columns.Push(1);
		data.pd.columnCount++;
		MarkerInfo * info = data.GetMarkerInfo(i);
		info->chromosome = (chrom == "X" || chrom == "x") ? 999 : atoi(chrom.c_str());
		info->positionFemale = info->positionMale = info->position = atoi(positions[i].c_str()) * positionAdjustment;
		// FIXME: should not need to clear but somehow this is the case ...
		// std::clog << info->freq.dim << std::endl;
		info->freq.Clear(); info->alleleLabels.Clear(); info->alleleNumbers.Clear();
	}
}


void SEQLinco::PedigreeData::LoadSamples(const VecVecString & samples)
{

	for (unsigned i = 0; i < samples.size(); ++i) {
		VecString fam_info(samples[i].begin(), samples[i].begin() + 5);
		VecString genotypes(samples[i].begin() + 5, samples[i].end());
		__AddPerson(fam_info, genotypes);
	}
	data.Sort();
	SortFamilies(data);
}


void SEQLinco::PedigreeData::__AddPerson(VecString & fam_info, VecString & genotypes)
{
	// add person info
	bool sex_failure = false;

	data.AddPerson(fam_info[0].c_str(), fam_info[1].c_str(),
		fam_info[2].c_str(), fam_info[3].c_str(),
		data.TranslateSexCode(fam_info[4].c_str(), sex_failure));
	// add person genotypes
	for (unsigned i = 0; i < genotypes.size(); ++i) {
		String c1 = genotypes[i].c_str()[0];
		String c2 = genotypes[i].c_str()[1];
		Alleles new_genotype;
		new_genotype[0] = data.LoadAllele(data.GetMarkerInfo(i), c1);
		new_genotype[1] = data.LoadAllele(data.GetMarkerInfo(i), c2);
		if (new_genotype.isKnown()) data[data.count - 1].markers[i] = new_genotype;
	}
}


void SEQLinco::GeneticHaplotyper::Apply(Pedigree & ped)
{
	// activate these analysis options
	FamilyAnalysis::bestHaplotype = true;
	FamilyAnalysis::zeroRecombination = false;
	MerlinHaplotype::outputHorizontal = true;
	String chrom = __chrom.c_str();
	if (chrom == "X") PedigreeGlobals::chromosomeX = true;
	//
	ped.EstimateFrequencies(0, true);
	// recode alleles so more frequent alleles have lower allele numbers internally
	ped.LumpAlleles(0.0);
	// remove uninformative family or individuals
	// !! Do not trim here, because if a family is uninformative we can report as is
	// ped.Trim(true);
	FamilyAnalysis engine(ped);
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


void SEQLinco::HaplotypeCoder::Apply(VecVecVecString & ghdata)
{
	// each element of hpdata is a family's data
	// each element of hpdata[i] is a haplotype with the first item being sample id
	if (!ghdata.size()) return;
	// recode
	VecVecVecString hpdata(ghdata.size());
	for (unsigned f = 0; f < ghdata.size(); f++) {
		hpdata[f].resize(ghdata[f].size());
		// genotype pattern pool for a family
		std::set<std::string> pool;
		for (unsigned p = 0; p < ghdata[f].size(); p++) {
			hpdata[f][p].resize(3);
			// fam and individual ID
			hpdata[f][p][0] = ghdata[f][p][0];
			hpdata[f][p][1] = ghdata[f][p][1];
			// haplotype super string
			std::string haplotype;
			for (unsigned i = 2; i < ghdata[f][p].size(); i++) {
				// break due to recombination event detected
				if (!hasEnding(ghdata[f][p][i], ":") && !hasEnding(ghdata[f][p][i], "|")) {
					recombCount += 1;
					break;
				}
				// use one of the likely haplotype configuration
				haplotype.push_back((ghdata[f][p][i][0] == 'A') ? ghdata[f][p][i][1] : ghdata[f][p][i][0]);
			}
			haplotype = __Collapse(haplotype);
			hpdata[f][p][2] = haplotype;
			if (haplotype != "?") pool.insert(haplotype);
		}
		// convert haplotype super strings to haplotype patterns
		VecString patterns(pool.begin(), pool.end());
		std::sort(patterns.begin(), patterns.end());
		for (unsigned p = 0; p < hpdata[f].size(); p++) {
			if (hpdata[f][p][2] == "?") hpdata[f][p][2] = "0";
			else hpdata[f][p][2] = std::to_string(std::find(patterns.begin(), patterns.end(), hpdata[f][p][2]) - patterns.begin() + 1);
		}
	}
	// format output
	for (unsigned f = 0; f < hpdata.size(); f++) {
		for (unsigned p = 0; p < hpdata[f].size(); p++) {
			if (!data.size() || (data.back()[0] != hpdata[f][p][0] || data.back()[1] != hpdata[f][p][1])) {
				data.push_back(hpdata[f][p]);
			} else {
				data[data.size() - 1].push_back(hpdata[f][p][2]);
			}
		}
	}
}


unsigned SEQLinco::HaplotypeCoder::__AdjustSize(int n)
{
	if (__size <= 0) return n;
	div_t divresult;
	divresult = div(n, __size);
	// reduce size by x such that rem + res * x = self.size - x
	return __size - (int)((__size - divresult.rem) / (double)(divresult.quot + 1));
}


std::string SEQLinco::HaplotypeCoder::__Collapse(std::string & haplotype)
{
	if (__size == 1) {
		if (haplotype.find("?") != std::string::npos) return "?";
		else return haplotype;
	}
	unsigned wsize = __AdjustSize((int)haplotype.size());
	std::string collapsed_haplotype;
	unsigned counter = 0;
	char code = '1';
	for (unsigned i = 0; i < haplotype.size(); i++) {
		counter += 1;
		if (haplotype[i] == '2') code = '2';
		else code = (code == '?' || code == '2') ? code : haplotype[i];
		if (counter == wsize) {
			collapsed_haplotype.push_back(code);
			counter = 0;
			code = '1';
		}
	}
	// make it missing data if "?" is in the haplotype
	if (collapsed_haplotype.find("?") != std::string::npos) return "?";
	else return collapsed_haplotype;
}


