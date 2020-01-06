#include <fstream>
#include <boost/algorithm/string.hpp>
#include "AnalysisTools/Run/include/program_main.h"

#include "h-tautau/Analysis/include/HH_BTag.h"

int main()
{
    HH_BTag test("model");
    auto x = test.hello();

}
