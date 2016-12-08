#ifndef INCOMPATIBILITYCHECKER_CXX
#define INCOMPATIBILITYCHECKER_CXX

#include "IncompatibilityChecker.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderException.h"

namespace flashana {

  static IncompatibilityCheckerFactory __global_IncompatibilityChecker__;

  IncompatibilityChecker::IncompatibilityChecker(const std::string name)
    : BaseAlgorithm(kCustomAlgo, name)
    , _gap         ( 0.5    )
  {}

  void IncompatibilityChecker::_Configure_(const Config_t &pset)
  {
    _gap              = pset.get< double > ( "SegmentSize"      );
    _sigmaThreshold   = pset.get< double > ( "SigmaThreshold"   );
    _nBinsRequirement = pset.get< int    > ( "NBinsRequirement" );
  }

  bool IncompatibilityChecker::CheckIncompatibility(const Flash_t &flash, const Flash_t &flash_hypo) {

    std::cout << "IncompatibilityChecker::CheckIncompatibility starts" << std::endl;
    // Now we have two spectra: the 'true' one (flash), and the hypo one (flash_hypo)
    // We want to understand is they are incompatible
    std::cout << "pe for pmt 0 from hypo " << flash_hypo.pe_v[0] << std::endl;

    if (flash.pe_v.size() != flash_hypo.pe_v.size()) {
     throw OpT0FinderException("Flash and hypo flash pe vector size mismatch."); 
    }

    double totalPE_flash = 0.;
    double totalPE_hypo  = 0.;

    // Check bin by bin incompatibility
    int nIncompBins = 0;
    for (unsigned int pmt = 0; pmt < flash.pe_v.size(); pmt++) {
      double error  = std::sqrt(flash_hypo.pe_v[pmt]);
      double nsigma = (flash_hypo.pe_v[pmt] - flash.pe_v[pmt]) / error;
      if (nsigma > _sigmaThreshold) {
        nIncompBins ++;
        if (nIncompBins >= _nBinsRequirement){
          return true;
        }
      }
      totalPE_flash += flash.pe_v[pmt];
      totalPE_hypo  += flash_hypo.pe_v[pmt];
    }

    // Check overall incompatibility
    double error  = std::sqrt(totalPE_hypo);
    double nsigma = (totalPE_hypo - totalPE_flash) / error;
    if (nsigma > _sigmaThreshold) {
      return true;
    }

    std::cout << "HERE WE ARE!" << std::endl;
    return false;
  
  }
}


#endif
