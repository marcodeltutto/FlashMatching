////////////////////////////////////////////////////////////////////////
// Class:       CosmicFlashTagger
// Plugin Type: producer (art v2_05_00)
// File:        CosmicFlashTagger_module.cc
//
// Generated at Wed Nov 30 09:45:20 2016 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/IncompatibilityChecker.h"

#include "TString.h"
#include "TTree.h"

#include <memory>

class CosmicFlashTagger;


class CosmicFlashTagger : public art::EDProducer {
public:
  explicit CosmicFlashTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicFlashTagger(CosmicFlashTagger const &) = delete;
  CosmicFlashTagger(CosmicFlashTagger &&) = delete;
  CosmicFlashTagger & operator = (CosmicFlashTagger const &) = delete;
  CosmicFlashTagger & operator = (CosmicFlashTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  ::geoalgo::Trajectory GetTrajectory(art::Ptr<recob::Track> track, double xoffset);

  ::flashana::FlashMatchManager       _mgr;
  std::vector<::flashana::Flash_t>    beam_flashes;
  std::vector<art::Ptr<recob::Track>> track_v;
    
  std::string _track_producer;
  std::string _opflash_producer_beam;
  std::string _opflash_producer_cosmic;
  std::string _trigger_producer;
  double _flash_trange_start;
  double _flash_trange_end;
  double _min_track_length;
  bool _debug;

  TTree* _tree1;
  int _run, _subrun, _event, _matchid;
  int _n_beam_flashes, _n_tracks;
  std::vector<std::vector<double>> _beam_flash_spec, _track_hypo_spec;
};


CosmicFlashTagger::CosmicFlashTagger(fhicl::ParameterSet const & p)
{
  _track_producer          = p.get<std::string>("TrackProducer");
  _opflash_producer_beam   = p.get<std::string>("BeamOpFlashProducer");
  _opflash_producer_cosmic = p.get<std::string>("CosmicOpFlashProducer");
  _trigger_producer        = p.get<std::string>("TriggerProducer");
  _flash_trange_start      = p.get<double>     ("FlashVetoTimeStart");
  _flash_trange_end        = p.get<double>     ("FlashVetoTimeEnd");
  _min_track_length        = p.get<double>     ("MinimumTrackLength");
  _debug                   = p.get<bool>       ("DebugMode");

  _mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",&_run,"run/I");
  _tree1->Branch("subrun",&_subrun,"subrun/I");
  _tree1->Branch("event",&_event,"event/I");
  _tree1->Branch("n_beam_flashes",&_n_beam_flashes,"n_beam_flashes/I");
  _tree1->Branch("beam_flash_spec","std::vector<std::vector<double>>",&_beam_flash_spec);
  _tree1->Branch("ntracks",&_n_tracks,"n_tracks/I");
  _tree1->Branch("track_hypo_spec","std::vector<std::vector<double>>",&_track_hypo_spec);

  produces< std::vector<anab::CosmicTag>>();  
  //produces< art::Assns<anab::CosmicTag,   recob::Track>>();  
  //produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();
}

void CosmicFlashTagger::produce(art::Event & e)
{
  if(_debug) {
    std::cout << "CosmicFlashTagger::produce starts." << std::endl;
    std::cout << "This is run | subrun | event: " << e.id().run() << " | " << e.id().subRun() << " | " << e.id().event() << std::endl;
  }

  // Instantiate the output
  std::unique_ptr< std::vector< anab::CosmicTag > >                  cosmicTagTrackVector(       new std::vector<anab::CosmicTag>                  );
  //std::unique_ptr< art::Assns<anab::CosmicTag,   recob::Track > >    assnOutCosmicTagTrack(      new art::Assns<anab::CosmicTag,   recob::Track   >);
  //std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag > > assnOutCosmicTagPFParticle( new art::Assns<recob::PFParticle, anab::CosmicTag>);

  _mgr.Reset();

  ::art::ServiceHandle<geo::Geometry> geo;

  double Xoffset = -9999.;

  // Get Beam Flashes from the ART event
  ::art::Handle<std::vector<recob::OpFlash> > beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);
  if( !beamflash_h.isValid() || beamflash_h->empty() ) {
    std::cerr << "Don't have good flashes." << std::endl;
    return;
  }

  // Get Tracks from the ART event
  ::art::Handle<std::vector<recob::Track> > track_h;
  e.getByLabel(_track_producer,track_h);
  if( !track_h.isValid() || track_h->empty() )  {
    std::cerr << "Don't have tracks, or they are not valid." << std::endl;
    throw std::exception();
  }

  // Save beam flashes to file
  _n_beam_flashes = 0;
  for (size_t n = 0; n < beamflash_h->size(); n++) {

    auto const& flash = (*beamflash_h)[n];
    if (flash.Time() < 0. && flash.Time() > 50.) {
      continue;
    }
    if(flash.Time() < _flash_trange_start || _flash_trange_end < flash.Time()) {
      std::cout << "Flash is in veto region (flash time is " << flash.Time() << "). Continue." << std::endl;
      continue;
    } 
    
    _n_beam_flashes++;
    // resize vec
    _beam_flash_spec.resize(_n_beam_flashes);
    _beam_flash_spec[_n_beam_flashes-1].resize(geo->NOpDets());
    for (unsigned int i = 0; i < geo->NOpDets(); i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      _beam_flash_spec[_n_beam_flashes-1][opdet] = flash.PE(i);
    }

    // Construct a Flash_t
    ::flashana::Flash_t f;
    f.x = f.x_err = 0;
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    for (unsigned int i = 0; i < f.pe_v.size(); i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      f.pe_v[opdet] = flash.PE(i);
      f.pe_err_v[opdet] = sqrt(flash.PE(i));
    }
    f.time = flash.Time();
    //f.idx = flash_id;
    beam_flashes.resize(_n_beam_flashes);
    beam_flashes[_n_beam_flashes-1] = f;

  } // end of flash loop

  if(_debug) std::cerr << _n_beam_flashes << " flashes have been saved to file" << std::endl;


  _n_tracks = 0;
  for (size_t trk_idx=0; trk_idx<track_h->size(); trk_idx++) {

    if (_debug) std::cerr << "This is track " << trk_idx << std::endl;

    const art::Ptr<recob::Track> track_ptr(track_h,trk_idx);

    // Minimum required track length
    if (track_ptr->Length() < _min_track_length) {
      std::cerr << "Skipping this track, length is less than " << _min_track_length << " cm." << std::endl;
      continue;
    }
    _n_tracks++;
    _track_hypo_spec.resize(_n_tracks);

    track_v.resize(_n_tracks);
    track_v[_n_tracks-1] = track_ptr;

  } // end of track loop 



  


  std::cout << "beam_flashes.size() " << beam_flashes.size() << std::endl;
  std::cout << "track_v.size()      " << track_v.size()  << std::endl;

  // --- Loop over tracks ---
  for (unsigned int tt = 0; tt < track_v.size(); tt++) {
    std::cout << "    tt is " << tt << std::endl;

    // --- Loop over beam flashes ---
    for (unsigned int bf = 0; bf < beam_flashes.size(); bf++) {
      std::cout << "  bf is " << bf << std::endl;

      // Get the beam flash
      ::flashana::Flash_t flashBeam = beam_flashes[bf];

      // Calculate x offset, assuming this track caused this beam flash
      double Xoffset = flashBeam.time * 0.1114359;
      if(_debug) std::cerr << "Xoffset is " << Xoffset << std::endl;

      // Get track trajectory
      ::geoalgo::Trajectory trkTrj = this->GetTrajectory(track_v[tt], Xoffset);

      // From the trajectory construct a QCluster
      auto qcluster = ((flashana::LightPath*)(_mgr.GetCustomAlgo("LightPath")))->FlashHypothesis(trkTrj);

      // From the QCluster get the flash hypothesis using registered FlashHypothesis algorithm
      flashana::Flash_t flashHypo;
      flashHypo.pe_v.resize(geo->NOpDets());
      ((flashana::PhotonLibHypothesis*)(_mgr.GetAlgo(flashana::kFlashHypothesis)))->FillEstimate(qcluster,flashHypo);

      // CORE FUNCTION: Check if this beam flash and this flash hypothesis are incompatible
      bool areIncompatible = ((flashana::IncompatibilityChecker*)(_mgr.GetCustomAlgo("IncompatibilityChecker")))->CheckIncompatibility(flashBeam,flashHypo);
      std::cout << "for this track: " << areIncompatible << std::endl;

      if (areIncompatible == false) break;
      else if (areIncompatible && bf == beam_flashes.size() - 1) {
        // This track is not compatible with any of the beam flashes
        std::vector<float> endPt1;
        std::vector<float> endPt2;
        endPt1.resize(3);
        endPt1.push_back( -9999. );
        endPt2.resize(3);
        endPt2.push_back( -9999. );
        float cosmicScore = 1;
        anab::CosmicTagID_t tag_id = anab::CosmicTagID_t::kFlash_BeamIncompatible;
        cosmicTagTrackVector->emplace_back( endPt1, endPt2, cosmicScore, tag_id);
      }


    } // end of track trajectory loop
  } // end of beam flash loop





  _tree1->Fill();

  e.put( std::move(cosmicTagTrackVector)      );

  if(_debug) std::cout << "CosmicFlashTagger::produce ends." << std::endl;
}

::geoalgo::Trajectory CosmicFlashTagger::GetTrajectory(art::Ptr<recob::Track> track_ptr, double Xoffset) {

  ::geoalgo::Trajectory track_geotrj;
  track_geotrj.resize(track_ptr->NumberTrajectoryPoints(),::geoalgo::Vector(0.,0.,0.));
  for (size_t pt_idx=0; pt_idx < track_ptr->NumberTrajectoryPoints(); ++pt_idx) {
    auto const& pt = track_ptr->LocationAtPoint(pt_idx);
    track_geotrj[pt_idx][0] = pt[0] - Xoffset;
    track_geotrj[pt_idx][1] = pt[1];
    track_geotrj[pt_idx][2] = pt[2];
  }
  
  return track_geotrj;
}




DEFINE_ART_MODULE(CosmicFlashTagger)
