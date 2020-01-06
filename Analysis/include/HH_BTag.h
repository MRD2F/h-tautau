#include <fstream>
#include <boost/algorithm/string.hpp>
#include <Math/VectorUtil.h>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

class HH_BTag {
public:
    HH_BTag(const std::string& model);
    // ~HH_BTag();

    float hello();
    // std::vector<float> GetScore();

    enum pdg {jet_valid = 0, jet_pt = 1, jet_eta = 2, rel_jet_M_pt = 3, rel_jet_E_pt = 4, jet_htt_deta = 5,
             jet_deepFlavour = 6, jet_htt_dphi = 7, sample_year= 8, channelId = 9, htt_pt = 10, htt_eta = 10,
             htt_met_dphi = 11, rel_met_pt_htt_pt = 12, htt_scalar_pt = 13};
private:
    int number_of_var;
    std::string model;
    //     std::shared_ptr<tensorflow::GraphDef> graph;
    //     std::shared_ptr<tensorflow::Session>* session;
    //     tensorflow::Tensor x;

};
