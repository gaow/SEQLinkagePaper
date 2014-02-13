// $File: chp.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)
#ifndef _CHP_HPP_
#define _CHP_HPP_

#include <string>
#include <stdexcept>
#include "Exception.hpp"
#include "Core.hpp"
void resetPed(Pedigree & ped)
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
	// seems it's just the markerInfo part gives problems; ped.count and ped.familyCount are always 0

	// delete old pointer
	for (int i = 0; i < ped.markerInfoCount; i++)
		delete ped.markerInfo[i];
	delete [] ped.markerInfo;
	delete [] ped.markerInfoByInteger;
	// FIXME: Only Clear() is not going to delete the char * buffer pointer which unfortunately is protected ...
	// so there will still be memory leak here
	ped.markerNames.Clear();
	ped.markerLookup.Clear();
	ped.markerInfoByName.Clear();
	ped.markerCount = ped.markerInfoCount = ped.markerInfoSize = 0;
	// reset pointer
	ped.GrowMarkerInfo();
	// delete old pointer
	for (int i = 0; i < ped.count; i++)
		delete ped.persons[i];

	for (int i = 0; i < ped.familyCount; i++)
		delete ped.families[i];

	delete [] ped.families;
	delete [] ped.persons;
	ped.size = 10000;
	ped.count = ped.familyCount = ped.haveTwins = 0;
	// reset pointer
	ped.persons = new Person *[ped.size];
	ped.families = new Family * [1];
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

		resetPed(ped);

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
			// } catch (...) {
		} catch (std::exception e) {
			// std::clog << e.what() << std::endl;
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
