// $File: chp.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#ifndef _CHP_HPP_
#define _CHP_HPP_

#include <string>
#include "Core.hpp"
#include "Exception.hpp"
void reset_ped(Pedigree & ped)
{

	// FIXME: It does not make sense to me but I have to reset ped data object manually,
	// otherwise in a Python program that calls CHP::Apply multiple times in a loop,
	// ped.GetMarkerInfo(i)->freq.dim will not equal 0 after a few runs
	// complaining the same marker name has been previously used.
	// cannot yet figure out why as this is suppose to be a brand new ped object here!
	// UPDATE:
	// It might due to using it from Python. The swig generated wrapper did not delete it after use
	// anyways let me just manually clean up everything instead of wrestling with swig
	// UPDATE 2:
	// seems cleaning markerInfo would be sufficient

	// from PedigreeGlobal.h
	for (int i = 0; i < ped.markerInfoCount; i++) {
		ped.markerInfo[i]->freq.Clear();
		ped.markerInfo[i]->name.Clear();
		ped.markerInfo[i]->alleleLabels.Clear();
		ped.markerInfo[i]->alleleNumbers.Clear();
	}
	ped.markerNames.Clear();
	ped.markerLookup.Clear();
	ped.markerInfoByName.Clear();
	// from Pedigree.h
	for (int i = 0; i < ped.count; i++) {
		ped.persons[i]->famid.Clear();
		ped.persons[i]->pid.Clear();
		ped.persons[i]->motid.Clear();
		ped.persons[i]->fatid.Clear();
		ped.persons[i]->sex = ped.persons[i]->sibCount = ped.persons[i]->ngeno = ped.persons[i]->zygosity = 0;
	}
	for (int i = 0; i < ped.familyCount; i++) {
		ped.families[i]->famid.Clear();
		ped.families[i]->serial = ped.families[i]->first = ped.families[i]->last = ped.families[i]->count = ped.families[i]->founders = ped.families[i]->nonFounders = ped.families[i]->generations = 0;
	}
	ped.pd.columns.Clear();
	ped.pd.columnHash.Clear();
	ped.pd.columnCount = 0;
	ped.markerCount = ped.markerInfoCount = ped.markerInfoSize = 0;
	ped.count = ped.familyCount = ped.haveTwins = 0;
	ped.size = 10000;
}


namespace SEQLinco {

class CHP
{
public:
	CHP(const int wsize, const double padj = 0.01, const int verbose = 0) :
		__windowSize(wsize), __positionAdjustment(padj), __verbose(verbose) {};
	~CHP() {};
	CHP * clone() const { return new CHP(*this); }
	// input chrom must be 1 .. 22 and X
	// input samples follows PED convention, e.g., { "1", "1", "0", "0", "1", "21", "21", "21" }
	// first 5 cols are fid, sid, pid, mid and sex (0/1/2 coding), followed by genotypes
	// (1/2 coding, 0 for missing)
	VecVecString Apply(const std::string & chrom, const VecString & marker_names,
	                   const VecString & marker_positions, const VecVecString & samples)
	{
		Pedigree ped;

		reset_ped(ped);

		try {
			DataLoader dl;
			dl.LoadVariants(ped, marker_names, marker_positions, chrom, __positionAdjustment);
			dl.LoadSamples(ped, samples);
			MendelianErrorChecker mc;
			mc.Apply(ped);
			__mendelianErrorCount = mc.errorCount;
			GeneticHaplotyper gh(chrom);
			gh.Apply(ped);
			if (__verbose) gh.Print();
			HaplotypeCoder hc(__windowSize);
			hc.Apply(gh.data);
			__recombCount = hc.recombCount;
			if (__verbose) hc.Print();
			return hc.data;
		} catch (...) {
			const VecVecString nulldata(0);
			return nulldata;
		}
	}


	int countMendelianErrors() { return __mendelianErrorCount; }
	int countRecombs() { return __recombCount; }

private:
	int __mendelianErrorCount;
	int __recombCount;
	int __windowSize;
	double __positionAdjustment;
	int __verbose;
};
}
#endif
