/*! b-jet selection recommended by b-tag POG.
If not specified otherwise, all definitions are taken from the b-tag POG 102X recommendation TWiki:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
namespace btag_2018 { 
    constexpr double DeepCSVL = 0.1241; // >
    constexpr double DeepCSVM = 0.4184; // >
    constexpr double DeepCSVT = 0.7527; // >

    constexpr double deepFlavourL = 0.0494; // >
    constexpr double deepFlavourM = 0.2770; // >
    constexpr double deepFlavourT = 0.7264; // >

    constexpr double pt = 20; // >
    constexpr double eta = 2.4; // <
}
} // namespace cuts