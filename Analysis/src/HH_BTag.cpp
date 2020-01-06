#include "h-tautau/Analysis/include/HH_BTag.h"

float HH_BTag::hello()
{
    std::cout << "nutella" << '\n';
    return 10;
}

    HH_BTag::HH_BTag(const std::string& _model): model(_model)
    {
        number_of_var = 12;
    }
    //     graph_par0(tensorflow::loadGraphDef(model), session_par0(tensorflow::createSession(graph.get()))/*,
    //     x(tensorflow::DT_FLOAT, {1, dnn_input_2017v1::NumberOfInputs})*/

    // float HH_BTag::GetScore()
    // {
    //
    // }
    // ~HH_BTag::HH_BTag()
    // {
    //     tensorflow::closeSession(session_par0);
    // }

//}
