// $File: chp.cpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#include "chp.hpp"
#include <iostream>
#include <cstdlib>
void SEQLinco::PedigreeData::LoadVariants(std::vector<std::string> & names,
                                          std::vector<int> & positions, std::string chrom)
{
	for (unsigned i = 0; i < names.size(); ++i) {
		data.pd.columnHash.Push(data.GetMarkerID(names[i].c_str()));
		data.pd.columns.Push(1);
		data.pd.columnCount++;
		MarkerInfo * info = data.GetMarkerInfo(i);
		info->chromosome = (chrom == "X") ? 999 : atoi(chrom.c_str());
		info->position = (double)positions[i] * 0.01;
	}
}


void SEQLinco::PedigreeData::LoadSamples(std::vector< std::vector<std::string> > & samples)
{

	for (unsigned i = 0; i < samples.size(); ++i) {
		std::vector<std::string> fam_info(samples[i].begin(), samples[i].begin() + 5);
		std::vector<std::string> genotypes(samples[i].begin() + 5, samples[i].end());
		__AddPerson(fam_info, genotypes);
	}
	data.Sort();
	SortFamilies(data);
}


void SEQLinco::PedigreeData::__AddPerson(std::vector<std::string> & fam_info,
                                         std::vector<std::string> & genotypes)
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
	data.resize(0);
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
	errorCount = 0;
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


void SEQLinco::HaplotypeCoder::Apply(std::vector< std::vector< std::vector<std::string> > > & ghdata)
{
	data.resize(0);
	if (!ghdata.size()) return;

}


int SEQLinco::HaplotypeCoder::__AdjustSize(int n)
{
	if (__size < 0) return n;
	div_t divresult;
	divresult = div(n, __size);
	// reduce size by x such that rem + res * x = self.size - x
	return __size - (int)((__size - divresult.rem) / (double)(divresult.quot + 1));
}


