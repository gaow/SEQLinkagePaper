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

void SEQLinco::DataLoader::LoadVariants(Pedigree & ped,
                                        const VecString & names,
                                        const VecInt & positions,
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
		// adjust input position to make sure every position is unique
		int position = positions[i];
		while (std::find(vs.begin(), vs.end(), position) != vs.end()) position++;
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
		String c1 = genotypes[i].substr(0, 1).c_str();
		String c2 = genotypes[i].substr(1, 1).c_str();
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
		std::vector<unsigned> recombPositions(0);
		for (unsigned p = 0; p < haploVecs[f].size(); p++) {
			for (unsigned i = 2; i < haploVecs[f][p].size(); i++) {
				// recombination event detected
				if (!SEQLinco::hasEnding(haploVecs[f][p][i], ":") && !SEQLinco::hasEnding(haploVecs[f][p][i], "|")) {
					recombCount += 1;
					if (std::find(recombPositions.begin(), recombPositions.end(), i) == recombPositions.end()) {
						recombPositions.push_back(i);
					}
				}
				// use one of the likely haplotype configuration
				haploVecs[f][p][i] = (haploVecs[f][p][i].substr(0, 1) == "A") ? haploVecs[f][p][i].substr(1, 1) : haploVecs[f][p][i].substr(0, 1);
			}
		}
		std::sort(recombPositions.begin(), recombPositions.end());
		//
		// collapse haplotype vectors to haplotype Strings
		//
		VecVecString haploStrs(haploVecs[f].size());
		VecVecString patterns(recombPositions.size() + 1);
		for (unsigned p = 0; p < haploVecs[f].size(); p++) {
			for (unsigned r = 0; r <= recombPositions.size(); r++) {
				std::string haplotype = SEQLinco::collapse(haploVecs[f][p],
					(r > 0) ? recombPositions[r - 1] : 2,
					(r == recombPositions.size()) ? haploVecs[f][p].size() : recombPositions[r],
					__size);
				if (haplotype != "?" && std::find(patterns[r].begin(), patterns[r].end(), haplotype) == patterns[r].end()) {
					patterns[r].push_back(haplotype);
				}
				haploStrs[p].push_back(haplotype);
			}
		}
		//
		// convert haplotype super strings to haplotype patterns
		//
		for (unsigned p = 0; p < haploVecs[f].size(); p++) {
			// convert one haplotype's all recomb segments to a single super marker
			for (unsigned r = 0; r <= recombPositions.size(); r++) {
				std::sort(patterns[r].begin(), patterns[r].end());
				if (haploStrs[p][r] == "?") haploStrs[p][r] = "0";
				else haploStrs[p][r] = std::to_string(std::find(patterns[r].begin(), patterns[r].end(), haploStrs[p][r]) - patterns[r].begin() + 1);
			}
			// push combined haplotype to data
			if (!data.size() || (data.back()[0] != haploVecs[f][p][0] || data.back()[1] != haploVecs[f][p][1])) {
				VecString newperson(haploVecs[f][p].begin(), haploVecs[f][p].begin() + 2);
				data.push_back(newperson);
			}
			data[data.size() - 1].push_back(std::accumulate(haploStrs[p].begin(), haploStrs[p].end(), std::string("")));
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


