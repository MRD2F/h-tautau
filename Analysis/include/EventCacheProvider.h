/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/CacheTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Analysis/include/EventCandidate.h"
#include "SVfitAnaInterface.h"
#include "KinFitInterface.h"

namespace analysis {

class EventCacheProvider {
public:
    using LegPair = ntuple::LegPair;

    EventCacheProvider(){} //default constructor

    template<typename Event>
    EventCacheProvider(const Event& event){
        AddEvent(event);
    }

    template<typename Event>
    void AddEvent(const Event& event)
    {
        for(unsigned n = 0; n < event.kinFit_Higgs_index.size(); ++n){
            KinFitKey kinFitKey(ntuple::LegIndexToPair(event.kinFit_Higgs_index.at(n)),
                                ntuple::LegIndexToPair(event.kinFit_jetPairId.at(n)),
                                static_cast<UncertaintySource>(event.kinFit_unc_source.at(n)),
                                static_cast<UncertaintyScale>(event.kinFit_unc_scale.at(n)));
            kin_fit::FitResults kinFit_results(event.kinFit_m.at(n),event.kinFit_chi2.at(n),
                                                  0,event.kinFit_convergence.at(n));
            kinFit_map[kinFitKey] = kinFit_results;
        }

        for(unsigned n = 0; n < event.SVfit_Higgs_index.size(); ++n){
            SVFitKey SVFitKey(ntuple::LegIndexToPair(event.SVfit_Higgs_index.at(n)),
                              static_cast<UncertaintySource>(event.SVfit_unc_source.at(n)),
                              static_cast<UncertaintyScale>(event.SVfit_unc_scale.at(n)));
            sv_fit_ana::FitResults SVFit_results(event.SVfit_is_valid.at(n),
                                                    analysis::LorentzVectorM(event.SVfit_p4.at(n)),
                                                    analysis::LorentzVectorM(event.SVfit_p4_error.at(n)),
                                                    event.SVfit_mt.at(n),
                                                    event.SVfit_mt_error.at(n));
            SVFit_map[SVFitKey] = SVFit_results;
        }
    }

    bool TryGetKinFit(kin_fit::FitResults& kinfit_result, const LegPair& htt_pair, const LegPair& hbb_pair, UncertaintySource unc_source, UncertaintyScale unc_scale);
    bool TryGetSVFit(sv_fit_ana::FitResults& svfit_result, const LegPair& htt_pair, UncertaintySource unc_source, UncertaintyScale unc_scale);

    struct SVFitKey{
        LegPair htt_pair;
        UncertaintySource unc_source;
        UncertaintyScale unc_scale;

        SVFitKey();
        SVFitKey(LegPair _htt_pair,UncertaintySource _unc_source,UncertaintyScale _unc_scale);

        bool operator<(const SVFitKey& other) const;

        virtual ~SVFitKey(){} // virtual destructor
        SVFitKey( const SVFitKey & ) = default;
    };

    struct KinFitKey : SVFitKey {
        LegPair hbb_pair;

        KinFitKey();
        KinFitKey(LegPair _htt_pair, LegPair _hbb_pair,UncertaintySource _unc_source,UncertaintyScale _unc_scale);
        bool operator<(const KinFitKey& other) const;
    };



private:
    std::map<KinFitKey,kin_fit::FitResults> kinFit_map;
    std::map<SVFitKey,sv_fit_ana::FitResults> SVFit_map;

};

} // namespace analysis