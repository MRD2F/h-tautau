/*! b-jet selection recommended by b-tag POG.
If not specified otherwise, all definitions are taken from the b-tag POG 80X recommendation TWiki:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
namespace btag_2016 {
    constexpr double CSVv2L = 0.460; // >
    constexpr double CSVv2M = 0.800; // >
    constexpr double CSVv2T = 0.935; // >

    constexpr double pt = 20; // >
    constexpr double eta = 2.4; // <
}
} // namespace cuts
