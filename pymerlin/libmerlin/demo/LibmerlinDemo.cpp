#include "Pedigree.h"
#include "MerlinFamily.h"
#include "MerlinHaplotype.h"
#include "MerlinSort.h"
#include <algorithm>
#include <vector>
#include <string>
#include <iterator>
#include <iostream>

void showPed(Pedigree & ped)
{
	printf("Loaded %d individuals\n", ped.count);
	for (int i = 0; i < std::min(ped.count, 10); i++) {
		printf("[%s]: %s, %s, %s, %d\t|\t",
			(const char *)ped[i].pid, (const char *)ped[i].famid,
			(const char *)ped[i].fatid, (const char *)ped[i].motid,
			ped[i].sex);
		for (int j = 0; j < ped.markerCount; ++j) {
			printf("%d%d\t", ped[i].markers[j].one, ped[i].markers[j].two);
		}
		printf("\n");
	}
	//
	printf("Loaded %d markers\n", ped.markerCount);
	// Estimate allele frequencies for all markers, verbose mode
	ped.EstimateFrequencies(1, false);
	// Get genotype statistics for markers
	for (int i = 0; i < ped.markerNames.Length(); i++) {
		printf("Statistics for marker [%s]\n", (const char *)ped.markerNames[i]);
		// Allele index starts with 1 not 0
		for (int j = 1; j <= ped.GetMarkerInfo(i)->CountAlleles(); j++) {
			printf("\tFrequency for allele %d: %f\n", j, ped.GetMarkerInfo(i)->freq[j]);
			printf("\tName for allele %d: %s\n", j, (const char *)ped.GetMarkerInfo(i)->GetAlleleLabel(j));
		}
	}
	return;
}


void readData(Pedigree & ped,
              const char * datfile, const char * pedfile, const char * mapfile)
{
	// The data file contains a description of the contents of the
	// pedigree file, including for example, a list of marker and
	// trait names
	ped.Prepare(datfile);
	// The pedigree file contains a list of individuals, stored one
	// per row, with specific information about each individual as
	// detailed in the data file.
	ped.Load(pedfile);
	SortFamilies(ped);
	ped.LoadMarkerMap(mapfile);
	return;
}


void loadVariants(Pedigree & ped, std::vector<std::string> & marker_ids,
                  std::vector<int> & marker_positions,
                  int chrom = 1)
{
	for (unsigned i = 0; i < marker_ids.size(); ++i) {
		ped.pd.columnHash.Push(ped.GetMarkerID(marker_ids[i].c_str()));
		ped.pd.columns.Push(1);
		ped.pd.columnCount++;
		MarkerInfo * info = ped.GetMarkerInfo(i);
		info->chromosome = chrom;
		info->position = (double)marker_positions[i] * 0.01;
	}
	return;
}


void addPerson(Pedigree & ped, std::vector<std::string> & fam_info,
               std::vector<std::string> & genotypes)
{
	// add person info
	bool failure = false;

	ped.AddPerson(fam_info[0].c_str(), fam_info[1].c_str(),
		fam_info[2].c_str(), fam_info[3].c_str(),
		ped.TranslateSexCode(fam_info[4].c_str(), failure));
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


void loadData(Pedigree & ped, int which = 1)
{
	if (which == 1) {
		//
		// haplo.dat
		//
		std::vector<std::string> marker_ids { "V1", "V2", "V3" };
		std::vector<int> marker_positions { 1, 2, 3 };
		loadVariants(ped, marker_ids, marker_positions);
		//
		// haplo.ped
		//
		std::vector<std::string> fam_info { "1", "1", "0", "0", "1" };
		std::vector<std::string> genotypes { "21", "21", "21" };
		addPerson(ped, fam_info, genotypes);
		fam_info = { "1", "2", "0", "0", "2" };
		genotypes = { "11", "11", "11" };
		addPerson(ped, fam_info, genotypes);
		fam_info = { "1", "3", "1", "2", "1" };
		genotypes = { "21", "21", "21" };
		addPerson(ped, fam_info, genotypes);
		fam_info = { "2", "1", "0", "0", "1" };
		genotypes = { "22", "21", "00" };
		addPerson(ped, fam_info, genotypes);
		fam_info = { "2", "2", "0", "0", "2" };
		genotypes = { "11", "11", "11" };
		addPerson(ped, fam_info, genotypes);
		fam_info = { "2", "3", "1", "2", "1" };
		genotypes = { "21", "21", "21" };
		addPerson(ped, fam_info, genotypes);
		fam_info = { "3", "1", "0", "0", "1" };
		genotypes = { "22", "21", "21" };
		addPerson(ped, fam_info, genotypes);
		fam_info = { "3", "2", "0", "0", "2" };
		genotypes = { "11", "11", "21" };
		addPerson(ped, fam_info, genotypes);
		fam_info = { "3", "3", "1", "2", "1" };
		genotypes = { "21", "21", "21" };
		addPerson(ped, fam_info, genotypes);
	}
	//
	// sort
	//
	ped.Sort();
	SortFamilies(ped);
}


void haplotyping(Pedigree & ped, String chrom)
{
	// activate these analysis options
	FamilyAnalysis::bestHaplotype = true;
	FamilyAnalysis::zeroRecombination = false;
	MerlinHaplotype::outputHorizontal = true;
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
	for (int i = 0; i < ped.familyCount; i++)
		if (engine.SelectFamily(ped.families[i])) {
			engine.Analyse();
			for (unsigned i = 0; i < engine.hapOutput.size(); ++i) {
				for (unsigned j = 0; j < engine.hapOutput[i].size(); ++j)
					std::cout << engine.hapOutput[i][j] << "\t";
				std::cout << std::endl;
			}

		}
	engine.CleanupGlobals();
}


int main(int argc, char ** argv)
{
	if (argc != 3) {
		printf("usage: %s <data source code: 1, 2, 3> <task: 1 or 2>\n", argv[0]);
		return 0;
	}

	Pedigree ped;
	if (atoi(argv[1]) == 1) readData(ped, "haplo.dat", "haplo.ped", "haplo.map");
	else if (atoi(argv[1]) == 2) readData(ped, "gene.dat", "gene.ped", "gene.map");
	else if (atoi(argv[1]) == 3) loadData(ped, 1);
	else ;
	if (atoi(argv[2]) == 1) showPed(ped);
	else if (atoi(argv[2]) == 2) haplotyping(ped, "1");
	else ;
}


