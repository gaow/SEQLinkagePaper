// $File: VCFStream.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

#ifndef _VCFSTM_HPP_
#define _VCFSTM_HPP_
#include <string>
#include <vector>
#include <stdexcept>
#include "Exception.hpp"
#include "VcfFileReader.h"

namespace SEQLinco {

typedef std::vector<int> VecInt;
typedef std::vector<std::string> VecString;

class VCFstream
{
public:
	/// Initialize VCFstream with VCF file
	/// \param vcf file
	/// \param vcf index
	VCFstream(const char * vcf) :
	{
		__reader.open(vcf, __header);
		try {
			__reader.readVcfIndex();
		} catch (std::exception & e) {
			throw RuntimeError("BAD VCF index");
		}
		__tabixPtr = __reader.getVcfIndex();
		if (__tabixPtr == NULL || tabixPtr->getFormat() != Tabix::FORMAT_VCF)
			throw RuntimeError("Failed to load VCF index");
	}


	~VCFstream()
	{
		__reader.close();
		delete [] __tabixPtr;
	}


	/// Copy function
	/// \return a copy of VCFstream
	VCFstream * clone() const
	{
		return new VCFstream(*this);
	}


	/// Get list of sample names
	/// \return list of sample names
	VecString GetSampleNames()
	{
		VecString samples(__header.getNumSamples());

		for (unsigned i = 0; i < samples.size(); ++i) {
			samples[i].assign(__header.getSampleName(i));
		}
		return samples;
	}


	/// Extract VCF region
	/// \param chrom
	/// \param start pos
	/// \param end pos
	void Extract(const char * chrom, int start, int end)
	{
		try {
			__reader.set1BasedReadSection(chrom, start, end);
		} catch (std::exception & e) {
			throw RuntimeError("Failed to extract VCF region");
		}
	}


	/// Point to next record with @line
	/// \return true if line is valid otherwise false
	bool Next()
	{
		return __reader.readRecord(line);
	}


	std::string GetChrom()
	{
		std::string chrom = line.getChromStr();

		return chrom;
	}


	int GetPosition()
	{
		return line.get1BasedPosition();
	}


	bool IsBiAllelic()
	{
		return line.getNumAlts() == 1;
	}


	/// Get sample genotype
	/// \param variant ID
	/// \param sample IDs
	/// \return list of sample genotypes
	Vecstring GetGenotypes(const VecInt & sid)
	{
		VecString genotypes(sid.size());

		for (unsigned i = 0; i < sid.size(); ++i) {
			int allele1 = line.getGT(sid[i], 0);
			int allele2 = line.getGT(sid[i], 1);
			allele1 = (allele1 == 0 || allele1 == 1) ? allele1 + 1 : 0;
			allele2 = (allele2 == 0 || allele2 == 1) ? allele2 + 1 : 0;
			genotypes[i] = std::to_string(allele1) + std::to_string(allele2);
		}
		return genotypes;
	}


	VcfRecord line;

private:
	VcfFileReader __reader;
	VcfHeader __header;
	const Tabix * __tabixPtr;
};
}
#endif
