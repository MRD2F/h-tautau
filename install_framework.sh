#!/usr/bin/env bash
# Install hh-italian-group framework.
# This file is part of https://github.com/hh-italian-group/h-tautau.

declare -A INSTALL_MODES
INSTALL_MODES=( ["prod16"]="CMSSW_10_2_14 _amd64_gcc700" \
                ["prod17"]="CMSSW_10_2_14 _amd64_gcc700" \
                ["ana"]="CMSSW_10_2_14 _amd64_gcc700" \
                ["ana_osx"]="bbtautau None")
DEFAULT_N_JOBS=4

function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

if [ $# -lt 1 -o $# -gt 3 ] ; then
    echo "Usage: mode [n_jobs] [cmssw_release]"
    printf "\n\t%-20sinstallation mode. Supported modes: " "mode"
    join_by ", " "${!INSTALL_MODES[@]}"
    printf ".\n\t%-20sthe number of jobs to run simultaneously during the compilation." "n_jobs"
    printf " Default: $DEFAULT_N_JOBS.\n"
    printf "\t%-20sCMSSW release." "cmssw_release"
    printf " Defaults:\n"
    for mode in "${!INSTALL_MODES[@]}" ; do
        MODE_DESC=( ${INSTALL_MODES[$mode]} )
        printf "\t%-20s\t%-10s\t%s\n" "" "$mode" "${MODE_DESC[0]}"
    done
    exit 1
fi

MODE=$1
if [[ ! ${!INSTALL_MODES[@]} =~ $MODE ]] ; then
    echo "ERROR: unsupported installation mode '$MODE'."
    printf "Supported installation modes: "
    join_by ", " "${!INSTALL_MODES[@]}"
    printf ".\n"
    exit 1
fi

N_JOBS=$2
if [ "x$N_JOBS" = "x" ] ; then N_JOBS=$DEFAULT_N_JOBS ; fi
if ! [ $N_JOBS -eq $N_JOBS ] 2>/dev/null ; then
    echo "ERROR: invalid number of jobs '$N_JOBS'."
    exit 1
fi

RELEASE=$3
if [ "x$RELEASE" = "x" ] ; then
    MODE_DESC=( ${INSTALL_MODES[$MODE]} )
    RELEASE=${MODE_DESC[0]}
    if [ "$MODE" != "ana_osx" ] ; then
        export SCRAM_ARCH="slc$(cat /etc/redhat-release | sed -E 's/[^67]*([67])\..*/\1/')${MODE_DESC[1]}"
    fi
fi
if [ -e $RELEASE ] ; then
    echo "ERROR: Working area for $RELEASE already exists."
    exit 1
fi

function run_cmd {
    "$@"
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Error while rinning '$@'"
        exit 1
    fi
}

if [ "$MODE" != "ana_osx" ] ; then
    run_cmd scramv1 project CMSSW $RELEASE
    cd $RELEASE/src
    run_cmd eval `scramv1 runtime -sh`
else
    mkdir -p "$RELEASE"
    cd "$RELEASE"
fi

if [ $MODE = "prod16" ] ; then
    run_cmd git cms-init

    # MET filters
    #git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate #outdated
    #git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter

    # MET corrections
    #git cms-merge-topic cms-met:METRecipe_8020
    #run_cmd git cms-merge-topic cms-met:METRecipe_8020_for80Xintegration

    # Tau ID
    #run_cmd git cms-merge-topic -u cms-tau-pog:CMSSW_8_0_X_tau-pog_tauIDOnMiniAOD-legacy-backport-81Xv2
fi

if [ $MODE = "prod17" ] ; then
    run_cmd git cms-init

    # Electron MVA identification
#git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
#scram b -j 8

    # Add the area containing the MVA weights (from cms-data, to appear in “external”).
    # Note: the “external” area appears after “scram build” is run at least once, as above
#cd $CMSSW_BASE/external/$SCRAM_ARCH
#git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
#cd data/RecoEgamma/PhotonIdentification/data
#git checkout CMSSW_9_4_0_pre3_TnP
#cd $CMSSW_BASE/external/$SCRAM_ARCH
#git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
#    cd data/RecoEgamma/ElectronIdentification/data
#    git checkout CMSSW_9_4_0_pre3_TnP
    # Go back to the src/
    #run_cmd git cms-addpkg RecoMET/METFilters
    #run_cmd git cms-merge-topic cms-egamma:EgammaPostRecoTools
    #Add DeepTau code from Tau POG repository (note "-u" option preventing checkout of unnecessary stuff)
    #run_cmd git cms-merge-topic -u cms-tau-pog:CMSSW_10_2_X_tau-pog_DeepTau2017v2
    #Add 2017v2 training file
    #run_cmd wget https://github.com/cms-tau-pog/RecoTauTag-TrainingFiles/raw/DeepTau2017v2_alone/DeepTauId/deepTau_2017v2p6_e6.pb -P RecoTauTag/TrainingFiles/data/DeepTauId
    #cd $CMSSW_BASE/src
fi

# old SVfit packages
#git clone git@github.com:hh-italian-group/SVfit_standalone.git TauAnalysis/SVfitStandalone
#cd TauAnalysis/SVfitStandalone
#git checkout hh_italian
#cd ../..

# new SVfit packages
run_cmd git clone git@github.com:hh-italian-group/ClassicSVfit.git TauAnalysis/ClassicSVfit
cd TauAnalysis/ClassicSVfit
run_cmd git checkout hh-italian
cd ../..
run_cmd git clone git@github.com:hh-italian-group/SVfitTF.git TauAnalysis/SVfitTF


# HHKinFit2 packages
run_cmd git clone git@github.com:hh-italian-group/HHKinFit2.git HHKinFit2/HHKinFit2

# LeptonEfficiencies packages
run_cmd git clone git@github.com:hh-italian-group/LeptonEff-interface.git HTT-utilities
run_cmd git clone git@github.com:hh-italian-group/LeptonEfficiencies.git HTT-utilities/LepEffInterface/data

# Recoil Corrections
if [ $MODE = "prod16" -o $MODE = "prod17" ] ; then
    run_cmd git clone https://github.com/CMS-HTT/RecoilCorrections.git HTT-utilities/RecoilCorrections
fi

# Install analysis packages
declare -A ANA_PACKAGES
ANA_PACKAGES=( ["AnalysisTools"]="prod16:master prod17:master ana:master ana_osx:master" \
               ["h-tautau"]="prod16:prod_v4 prod17:prod_v5 ana:ana_v4 ana_osx:ana_v4" \
               ["hh-bbtautau"]="prod16:ana_v4 prod17:ana_v5 ana:ana_v4 ana_osx:ana_v4" )
GITHUB_USER=$(git config user.github)

for pkg in "${!ANA_PACKAGES[@]}" ; do
    pkg_descs="${ANA_PACKAGES[$pkg]}"
    branch="master"
    for desc in $pkg_descs ; do
        if [ "${desc%%:*}" = "$MODE" ] ; then
            branch=${desc##*:}
            break
        fi
    done

    run_cmd git clone git@github.com:hh-italian-group/${pkg}.git
    cd "$pkg"
    if [ "$branch" != "master" ] ; then
        run_cmd git checkout -b $branch origin/$branch
    fi
    git ls-remote git@github.com:$GITHUB_USER/${pkg}.git &> /dev/null
    RESULT=$?
    if [ $RESULT -eq 0 ] ; then
        git remote add $GITHUB_USER git@github.com:$GITHUB_USER/${pkg}.git
        git fetch $GITHUB_USER
    fi
    cd ..
done

if [ "$MODE" = "ana_osx" ] ; then
    BUILD_PATH=build
else
    BUILD_PATH=../build
fi

# Prepare analysis working area
run_cmd ./AnalysisTools/Run/install.sh "$BUILD_PATH" "${!ANA_PACKAGES[@]}"

if ! [[ $MODE =~ ana.* ]] ; then
    run_cmd scram b -j$N_JOBS
fi

echo "Framework has been successfully installed."
