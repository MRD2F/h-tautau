/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "AnalysisTools/Core/include/RootExt.h"

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TriggerResults.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Analysis/include/EventCandidate.h"
#include "h-tautau/Analysis/include/EventCacheProvider.h"
#include "HHTools/HH-btag/include/HH_BTag.h"

#include "SVfitAnaInterface.h"
#include "KinFitInterface.h"
#include "MT2.h"

#include "SignalObjectSelector.h"

#include <numeric>


namespace analysis {

class SummaryInfo {
public:
    using ProdSummary = ntuple::ProdSummary;

    explicit SummaryInfo(const ProdSummary& _summary, const Channel& _channel,
                         const std::string& _trigger_cfg = "");
    std::shared_ptr<const TriggerDescriptorCollection> GetTriggerDescriptors() const;
    const ProdSummary& operator*() const;
    const ProdSummary* operator->() const;
    const jec::JECUncertaintiesWrapper& GetJecUncertainties() const;

private:
    ProdSummary summary;
    std::shared_ptr<const TriggerDescriptorCollection> triggerDescriptors;
    std::shared_ptr<jec::JECUncertaintiesWrapper> jecUncertainties;

};

class EventInfo {
public:
    using Event = ntuple::Event;
    using LegPair = ntuple::LegPair;
    using HiggsBBCandidate = CompositeCandidate<JetCandidate, JetCandidate>;
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;
    using HiggsTTCandidate = CompositeCandidate<LepCandidate, LepCandidate>;


    std::array<size_t,2> GetSelectedBjetIndices() const;
    std::set<size_t> GetSelectedBjetIndicesSet() const;

    Channel GetChannel() const { return static_cast<Channel>(event_candidate.GetEvent().channelId); }

    EventInfo(EventCandidate&& _event_candidate, const SummaryInfo* _summaryInfo,
                  size_t _selected_htt_index, const SignalObjectSelector::SelectedSignalJets& _selected_signal_jets,
                  Period _period, JetOrdering _jet_ordering);


    EventInfo(const EventInfo& ) = default; //copy constructor
    virtual ~EventInfo(){} //destructor

    EventInfo& operator= ( const EventInfo& ) = default; //assignment


    const Event& operator*() const;
    const Event* operator->() const;

    const EventIdentifier& GetEventId() const;
    const TriggerResults& GetTriggerResults() const;
    const SummaryInfo& GetSummaryInfo() const;

    size_t GetNJets() const;
    size_t GetNFatJets() const;
    size_t GetHttIndex() const;
    const SignalObjectSelector::SelectedSignalJets& GetSelectedSignalJets() const;
    Period GetPeriod() const;
    JetOrdering GetJetOrdering() const;

    JetCollection SelectJets(double pt_cut = std::numeric_limits<double>::lowest(),
                             double eta_cut = std::numeric_limits<double>::max(),
                             bool applyPu = false, bool passBtag = false,
                             JetOrdering jet_ordering = JetOrdering::DeepCSV,
                             const std::set<size_t>& jet_to_exclude_indexes = {},
                             double low_eta_cut = 0,
                             analysis::UncertaintySource unc_source = analysis::UncertaintySource::None,
                             analysis::UncertaintyScale unc_scale = analysis::UncertaintyScale::Central);

    double GetHT(bool includeHbbJets, bool apply_pt_eta_cut);
    const FatJetCollection& GetFatJets();
    bool HasBjetPair() const;
    bool HasVBFjetPair() const;
    const JetCandidate& GetVBFJet(const size_t index);
    const JetCandidate& GetBJet(const size_t index);
    const HiggsBBCandidate& GetHiggsBB();
    size_t GetLegIndex(const size_t leg_id);
    static bool PassDefaultLegSelection(const ntuple::TupleLepton& lepton, Channel channel);

    const kin_fit::FitResults& GetKinFitResults(bool allow_calc = false, int verbosity = 0);
    const sv_fit_ana::FitResults& GetSVFitResults(bool allow_calc = false, int verbosity = 0);

    const std::vector<float>& GetJetScore(double pt_cut, double eta_cut, JetOrdering jet_ordering,
                                          bool apply_pu, bool pass_btag, double low_eta_cut = 0,
                                          analysis::UncertaintySource unc_source = analysis::UncertaintySource::None,
                                          analysis::UncertaintyScale unc_scale = analysis::UncertaintyScale::Central);

    LorentzVector GetResonanceMomentum(bool useSVfit, bool addMET, bool allow_calc = false);
    double GetMT2();
    const FatJetCandidate* SelectFatJet(double mass_cut, double deltaR_subjet_cut);
    void SetMvaScore(double _mva_score);
    double GetMvaScore() const;

    const LepCandidate& GetFirstLeg();
    const LepCandidate& GetSecondLeg();

    const LepCandidate& GetLeg(size_t leg_id)
    {
        if(leg_id == 1) return GetFirstLeg();
        if(leg_id == 2) return GetSecondLeg();
        throw exception("Invalid leg id = %1%.") % leg_id;
    }

    const HiggsTTCandidate& GetHiggsTT(bool useSVfit, bool allow_calc = false)
    {
        Lock lock(*mutex);
        if(useSVfit) {
            if(!higgs_tt_sv) {
                if(!GetSVFitResults(allow_calc).has_valid_momentum) throw exception("SVFit not converged");
                higgs_tt_sv = std::make_shared<HiggsTTCandidate>(GetFirstLeg(), GetSecondLeg(),
                                                                 GetSVFitResults(allow_calc).momentum);
            }
            return *higgs_tt_sv;
        }
        if(!higgs_tt)
            higgs_tt = std::make_shared<HiggsTTCandidate>(GetFirstLeg(), GetSecondLeg());
        return *higgs_tt;
    }

    LorentzVector GetHiggsTTMomentum(bool useSVfit, bool allow_calc = false)
    {
        return GetHiggsTT(useSVfit, allow_calc).GetMomentum();
    }

    EventCandidate& GetEventCandidate()
    {
        return event_candidate;
    }

    const JetCollection& GetJets();
    const MET& GetMET();

protected:
    EventCandidate event_candidate;
    EventCacheProvider eventCacheProvider;
    const SummaryInfo* summaryInfo;
    TriggerResults triggerResults;
    std::shared_ptr<Mutex> mutex;
    size_t selected_htt_index;

private:
    EventIdentifier eventIdentifier;
    SignalObjectSelector::SelectedSignalJets selected_signal_jets;
    Period period;
    JetOrdering jet_ordering;

    std::shared_ptr<HiggsBBCandidate> higgs_bb;
    std::shared_ptr<kin_fit::FitResults> kinfit_results;
    std::shared_ptr<sv_fit_ana::FitResults> svfit_results;
    std::vector<float> jet_score;
    boost::optional<double> mt2;
    double mva_score;
    std::shared_ptr<HiggsTTCandidate> higgs_tt, higgs_tt_sv;

};

//to be added isEmbedded flag
boost::optional<EventInfo> CreateEventInfo(const ntuple::Event& event,
                                               const SignalObjectSelector& signalObjectSelector,
                                               const SummaryInfo* summaryInfo = nullptr,
                                               Period period = analysis::Period::Run2017,
                                               JetOrdering jet_ordering = JetOrdering::DeepCSV,
                                               bool is_sync = false,
                                               UncertaintySource uncertainty_source = UncertaintySource::None,
                                               UncertaintyScale scale = UncertaintyScale::Central,
                                               bool debug = false);

} // namespace analysis
