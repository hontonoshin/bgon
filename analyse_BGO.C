#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"
#include "TAxis.h"
#include "TLine.h"
#include "TF1.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TString.h"
#include "TPad.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Plotting and histogram helpers
// -----------------------------------------------------------------------------
void SetBGOStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadRightMargin(0.06);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(kGray);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.033);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.046, "XYZ");
    gStyle->SetLabelSize(0.039, "XYZ");
    gStyle->SetEndErrorSize(3.0);
}

void FixLogAxis(TAxis* axis)
{
    if (!axis) return;
    axis->SetMoreLogLabels(false);
    axis->SetNoExponent(false);
    axis->SetNdivisions(510);
    axis->SetLabelOffset(0.008);
    axis->SetTitleOffset(1.35);
}

std::vector<double> LogEdges(int nBins, double xMin, double xMax)
{
    std::vector<double> edges(static_cast<std::size_t>(nBins) + 1U);
    const double logMin = std::log10(xMin);
    const double logMax = std::log10(xMax);
    for (int i = 0; i <= nBins; ++i) {
        edges[static_cast<std::size_t>(i)] =
            std::pow(10.0, logMin + (logMax - logMin) * i / nBins);
    }
    return edges;
}

TH1D* MakeLogH1(const char* name, int nBins, double xMin, double xMax)
{
    const auto edges = LogEdges(nBins, xMin, xMax);
    auto* h = new TH1D(name, "", nBins, edges.data());
    h->SetDirectory(nullptr);
    h->Sumw2();
    return h;
}

TH2D* MakeLogLogH2(const char* name,
                   int nXBins, double xMin, double xMax,
                   int nYBins, double yMin, double yMax)
{
    const auto xEdges = LogEdges(nXBins, xMin, xMax);
    const auto yEdges = LogEdges(nYBins, yMin, yMax);
    auto* h = new TH2D(name, "", nXBins, xEdges.data(), nYBins, yEdges.data());
    h->SetDirectory(nullptr);
    h->Sumw2();
    return h;
}

TH2D* MakeLogLinearH2(const char* name,
                      int nXBins, double xMin, double xMax,
                      int nYBins, double yMin, double yMax)
{
    const auto xEdges = LogEdges(nXBins, xMin, xMax);
    auto* h = new TH2D(name, "", nXBins, xEdges.data(), nYBins, yMin, yMax);
    h->SetDirectory(nullptr);
    h->Sumw2();
    return h;
}

void StyleLine(TH1* h, Color_t color, Style_t lineStyle = 1,
               int lineWidth = 2, Style_t markerStyle = 0)
{
    if (!h) return;
    h->SetLineColor(color);
    h->SetLineStyle(lineStyle);
    h->SetLineWidth(lineWidth);
    h->SetMarkerColor(color);
    if (markerStyle > 0) {
        h->SetMarkerStyle(markerStyle);
        h->SetMarkerSize(0.75);
    }
}

void SaveCanvas(TCanvas* canvas, const char* outDir, const char* stem)
{
    if (!canvas) return;
    const TString base = TString::Format("%s/%s", outDir, stem);
    canvas->SaveAs(base + ".png");
    canvas->SaveAs(base + ".pdf");
    std::cout << "  -> " << base << ".png/.pdf\n";
}

double SafeDivide(double numerator, double denominator)
{
    return denominator != 0.0 ? numerator / denominator : 0.0;
}

// -----------------------------------------------------------------------------
// Branch helpers. Aliases let the macro read both the preserved GAGG-style
// branch names and future BGO-specific branch names.
// -----------------------------------------------------------------------------
template <typename T>
bool SetFirstExistingBranch(TTree* tree,
                            std::initializer_list<const char*> names,
                            T* address,
                            bool required = false)
{
    if (!tree || !address) return false;
    for (const char* name : names) {
        if (tree->GetBranch(name)) {
            tree->SetBranchStatus(name, 1);
            tree->SetBranchAddress(name, address);
            return true;
        }
    }

    if (required) {
        std::cerr << "Required branch is missing. Tried:";
        for (const char* name : names) std::cerr << " " << name;
        std::cerr << "\n";
    }
    return false;
}

TH2* GetDetachedH2(TFile* file, const char* name)
{
    if (!file) return nullptr;
    auto* source = dynamic_cast<TH2*>(file->Get(name));
    if (!source) return nullptr;
    auto* clone = dynamic_cast<TH2*>(source->Clone(TString::Format("%s_analysis", name)));
    if (clone) clone->SetDirectory(nullptr);
    return clone;
}

// -----------------------------------------------------------------------------
// Weighted energy-binned statistics
// -----------------------------------------------------------------------------
struct BinnedStats
{
    explicit BinnedStats(const std::vector<double>& binEdges)
        : edges(binEdges),
          sumW(binEdges.size() - 1U, 0.0),
          sumW2(binEdges.size() - 1U, 0.0),
          sumWX(binEdges.size() - 1U, 0.0),
          sumWX2(binEdges.size() - 1U, 0.0),
          entries(binEdges.size() - 1U, 0)
    {}

    int FindBin(double x) const
    {
        if (x < edges.front() || x > edges.back()) return -1;
        if (x == edges.back()) return static_cast<int>(edges.size()) - 2;
        const auto it = std::upper_bound(edges.begin(), edges.end(), x);
        const int index = static_cast<int>(it - edges.begin()) - 1;
        return (index >= 0 && index < static_cast<int>(sumW.size())) ? index : -1;
    }

    void Fill(double energy, double value, double weight)
    {
        const int bin = FindBin(energy);
        if (bin < 0 || !std::isfinite(value) || !std::isfinite(weight)) return;
        const std::size_t i = static_cast<std::size_t>(bin);
        sumW[i] += weight;
        sumW2[i] += weight * weight;
        sumWX[i] += weight * value;
        sumWX2[i] += weight * value * value;
        ++entries[i];
    }

    double Mean(std::size_t i) const
    {
        return SafeDivide(sumWX[i], sumW[i]);
    }

    double Variance(std::size_t i) const
    {
        if (sumW[i] <= 0.0) return 0.0;
        const double mean = Mean(i);
        return std::max(0.0, sumWX2[i] / sumW[i] - mean * mean);
    }

    double RMS(std::size_t i) const
    {
        return std::sqrt(Variance(i));
    }

    double EffectiveEntries(std::size_t i) const
    {
        return sumW2[i] > 0.0 ? sumW[i] * sumW[i] / sumW2[i] : 0.0;
    }

    double MeanError(std::size_t i) const
    {
        const double nEff = EffectiveEntries(i);
        return nEff > 0.0 ? RMS(i) / std::sqrt(nEff) : 0.0;
    }

    std::vector<double> edges;
    std::vector<double> sumW;
    std::vector<double> sumW2;
    std::vector<double> sumWX;
    std::vector<double> sumWX2;
    std::vector<Long64_t> entries;
};

struct BinnedFraction
{
    explicit BinnedFraction(const std::vector<double>& binEdges)
        : edges(binEdges),
          denW(binEdges.size() - 1U, 0.0),
          denW2(binEdges.size() - 1U, 0.0),
          numW(binEdges.size() - 1U, 0.0)
    {}

    int FindBin(double x) const
    {
        if (x < edges.front() || x > edges.back()) return -1;
        if (x == edges.back()) return static_cast<int>(edges.size()) - 2;
        const auto it = std::upper_bound(edges.begin(), edges.end(), x);
        const int index = static_cast<int>(it - edges.begin()) - 1;
        return (index >= 0 && index < static_cast<int>(denW.size())) ? index : -1;
    }

    void Fill(double energy, bool condition, double weight)
    {
        const int bin = FindBin(energy);
        if (bin < 0 || !std::isfinite(weight)) return;
        const std::size_t i = static_cast<std::size_t>(bin);
        denW[i] += weight;
        denW2[i] += weight * weight;
        if (condition) numW[i] += weight;
    }

    double Value(std::size_t i) const
    {
        return SafeDivide(numW[i], denW[i]);
    }

    double Error(std::size_t i) const
    {
        if (denW[i] <= 0.0 || denW2[i] <= 0.0) return 0.0;
        const double p = Value(i);
        const double nEff = denW[i] * denW[i] / denW2[i];
        return nEff > 0.0 ? std::sqrt(std::max(0.0, p * (1.0 - p) / nEff)) : 0.0;
    }

    std::vector<double> edges;
    std::vector<double> denW;
    std::vector<double> denW2;
    std::vector<double> numW;
};

TH1D* MeanHistogram(const char* name, const BinnedStats& stats,
                    const char* xTitle, const char* yTitle)
{
    auto* h = new TH1D(name, "", static_cast<int>(stats.edges.size()) - 1,
                       stats.edges.data());
    h->SetDirectory(nullptr);
    for (std::size_t i = 0; i < stats.sumW.size(); ++i) {
        h->SetBinContent(static_cast<int>(i) + 1, stats.Mean(i));
        h->SetBinError(static_cast<int>(i) + 1, stats.MeanError(i));
    }
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle(yTitle);
    return h;
}

TH1D* RelativeRMSHistogram(const char* name, const BinnedStats& stats,
                           const char* xTitle, const char* yTitle)
{
    auto* h = new TH1D(name, "", static_cast<int>(stats.edges.size()) - 1,
                       stats.edges.data());
    h->SetDirectory(nullptr);
    for (std::size_t i = 0; i < stats.sumW.size(); ++i) {
        const double mean = stats.Mean(i);
        const double value = mean != 0.0 ? stats.RMS(i) / std::abs(mean) : 0.0;
        h->SetBinContent(static_cast<int>(i) + 1, value);
    }
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle(yTitle);
    return h;
}

TH1D* FractionHistogram(const char* name, const BinnedFraction& fraction,
                        const char* xTitle, const char* yTitle)
{
    auto* h = new TH1D(name, "", static_cast<int>(fraction.edges.size()) - 1,
                       fraction.edges.data());
    h->SetDirectory(nullptr);
    for (std::size_t i = 0; i < fraction.denW.size(); ++i) {
        h->SetBinContent(static_cast<int>(i) + 1, fraction.Value(i));
        h->SetBinError(static_cast<int>(i) + 1, fraction.Error(i));
    }
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle(yTitle);
    return h;
}

void PrepareLinePlot(TH1* h, const char* xTitle, const char* yTitle)
{
    if (!h) return;
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle(yTitle);
    h->GetYaxis()->SetTitleOffset(1.45);
    FixLogAxis(h->GetXaxis());
}

} // namespace

// =============================================================================
// Main analysis
// =============================================================================
void analyse_BGO(const char* fname = "bgo_run0.root",
                 const char* outDir = "bgo_plots")
{
    SetBGOStyle();
    TH1::AddDirectory(kFALSE);
    gSystem->mkdir(outDir, kTRUE);

    std::unique_ptr<TFile> file(TFile::Open(fname, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << fname << "\n";
        return;
    }

    auto* events = dynamic_cast<TTree*>(file->Get("NeutronStudy"));
    if (!events) {
        std::cerr << "NeutronStudy tree is missing in " << fname << "\n";
        return;
    }

    std::cout << "[analyse_BGO] input: " << fname << "\n";
    std::cout << "[analyse_BGO] output directory: " << outDir << "\n";
    std::cout << "[analyse_BGO] events: " << events->GetEntries() << "\n";

    // The detector has two 3-cm BGO layers along the beam direction.
    // The 0.1-mm inter-layer gap is air and is excluded from material depth.
    constexpr double kBGOActiveDepthCm = 6.0;
    constexpr double kBGOEnvelopeDepthMm = 60.1;
    // PDG atomic/nuclear properties table for BGO: nuclear interaction length.
    constexpr double kBGOInteractionLengthCm = 22.32;
    const double expectedInteractionProbability =
        1.0 - std::exp(-kBGOActiveDepthCm / kBGOInteractionLengthCm);

    const auto energyEdges = LogEdges(30, 0.1, 10.0);

    // -------------------------------------------------------------------------
    // Event tree variables
    // -------------------------------------------------------------------------
    Int_t eventID = 0;
    Double_t primaryE = 0.0;
    Double_t weight = 1.0;
    Double_t edepTotal = 0.0;
    Double_t edepNP = 0.0;
    Double_t edepHeavy = 0.0;
    Double_t edepLight = 0.0;
    Double_t edepEM = 0.0;
    Double_t edepPiK = 0.0;
    Double_t edepOther = 0.0;
    Double_t edepFraction = 0.0;
    Int_t anyHadronic = 0;

    Int_t nHits = 0;
    Int_t nUniqueTracks = 0;
    Int_t nSecondaries = 0;
    Int_t nSecN = 0;
    Int_t nSecP = 0;
    Int_t nSecGamma = 0;
    Int_t nSecElectron = 0;
    Int_t nSecPositron = 0;
    Int_t nSecPion = 0;
    Int_t nSecKaon = 0;
    Int_t nSecLightIon = 0;
    Int_t nSecHeavyIon = 0;
    Int_t nSecOther = 0;
    Double_t sumSecondaryKE = 0.0;

    Int_t primaryHadronic = 0;
    Int_t primaryInelastic = 0;
    Int_t nPrimaryInelastic = 0;
    Int_t firstProcessCode = 0;
    Double_t firstDepth = -1.0;
    Double_t firstRadius = -1.0;
    Double_t firstTime = -1.0;
    Double_t firstPreKE = -1.0;
    Int_t firstMultiplicity = 0;
    Double_t firstSecondaryKE = 0.0;

    Double_t totalTrackLength = 0.0;
    Double_t maxStepEdep = 0.0;
    Double_t depthCentroid = -1.0;
    Double_t depthRMS = -1.0;
    Double_t radialCentroid = -1.0;
    Double_t radialRMS = -1.0;
    Double_t edepFront = 0.0;
    Double_t edepMiddle = 0.0;
    Double_t edepRear = 0.0;

    Double_t firstHitTime = -1.0;
    Double_t lastHitTime = -1.0;
    Double_t timeSpan = -1.0;
    Double_t promptEdep = 0.0;
    Double_t delayedEdep = 0.0;

    Double_t leakTotal = 0.0;
    Double_t leakNP = 0.0;
    Double_t leakEM = 0.0;
    Double_t leakLight = 0.0;
    Double_t leakHeavy = 0.0;
    Double_t leakPiK = 0.0;
    Int_t primaryEscaped = 0;
    Double_t primaryExitKE = 0.0;

    Int_t maxCrystalID = -1;
    Double_t maxCrystalEdep = 0.0;
    Int_t nHitCrystals = 0;
    Double_t edepLayer0 = 0.0;
    Double_t edepLayer1 = 0.0;

    events->SetBranchStatus("*", 0);
    const bool havePrimary = SetFirstExistingBranch(
        events, {"primaryE_GeV"}, &primaryE, true);
    const bool haveTotal = SetFirstExistingBranch(
        events, {"edepBGO_MeV", "edepGAGG_MeV"}, &edepTotal, true);
    const bool haveWeight = SetFirstExistingBranch(
        events, {"weight"}, &weight, false);
    const bool haveAnyHadronic = SetFirstExistingBranch(
        events, {"hadronicInBGO", "hadronicInGAGG"}, &anyHadronic, true);

    SetFirstExistingBranch(events, {"eventID"}, &eventID);
    SetFirstExistingBranch(events, {"edepNP_BGO_MeV", "edepNP_GAGG_MeV"}, &edepNP);
    SetFirstExistingBranch(events, {"edepHeavyIon_BGO_MeV", "edepHeavyIon_GAGG_MeV"}, &edepHeavy);
    SetFirstExistingBranch(events, {"edepLightIon_BGO_MeV", "edepLightIon_GAGG_MeV"}, &edepLight);
    SetFirstExistingBranch(events, {"edepEM_BGO_MeV", "edepEM_GAGG_MeV"}, &edepEM);
    SetFirstExistingBranch(events, {"edepPiK_BGO_MeV", "edepPiK_GAGG_MeV"}, &edepPiK);
    SetFirstExistingBranch(events, {"edepOther_BGO_MeV"}, &edepOther);
    const bool haveStoredFraction = SetFirstExistingBranch(
        events, {"edepFraction"}, &edepFraction);

    SetFirstExistingBranch(events, {"nHits"}, &nHits);
    const bool haveExtended = SetFirstExistingBranch(
        events, {"nUniqueTracks"}, &nUniqueTracks);
    SetFirstExistingBranch(events, {"nSecondaries"}, &nSecondaries);
    SetFirstExistingBranch(events, {"nSecondaryNeutrons"}, &nSecN);
    SetFirstExistingBranch(events, {"nSecondaryProtons"}, &nSecP);
    SetFirstExistingBranch(events, {"nSecondaryGammas"}, &nSecGamma);
    SetFirstExistingBranch(events, {"nSecondaryElectrons"}, &nSecElectron);
    SetFirstExistingBranch(events, {"nSecondaryPositrons"}, &nSecPositron);
    SetFirstExistingBranch(events, {"nSecondaryPions"}, &nSecPion);
    SetFirstExistingBranch(events, {"nSecondaryKaons"}, &nSecKaon);
    SetFirstExistingBranch(events, {"nSecondaryLightIons"}, &nSecLightIon);
    SetFirstExistingBranch(events, {"nSecondaryHeavyIons"}, &nSecHeavyIon);
    SetFirstExistingBranch(events, {"nSecondaryOther"}, &nSecOther);
    SetFirstExistingBranch(events, {"sumSecondaryKE_MeV"}, &sumSecondaryKE);

    const bool havePrimaryHadronic = SetFirstExistingBranch(
        events, {"primaryNeutronHadronic"}, &primaryHadronic);
    const bool havePrimaryInelastic = SetFirstExistingBranch(
        events, {"primaryNeutronInelastic"}, &primaryInelastic);
    SetFirstExistingBranch(events, {"nPrimaryNeutronInelastic"}, &nPrimaryInelastic);
    SetFirstExistingBranch(events, {"firstInteractionProcessCode"}, &firstProcessCode);
    const bool haveFirstDepth = SetFirstExistingBranch(
        events, {"firstInteractionDepth_mm"}, &firstDepth);
    SetFirstExistingBranch(events, {"firstInteractionR_mm"}, &firstRadius);
    SetFirstExistingBranch(events, {"firstInteractionTime_ns"}, &firstTime);
    SetFirstExistingBranch(events, {"firstInteractionPreKE_MeV"}, &firstPreKE);
    SetFirstExistingBranch(events, {"firstInteractionMultiplicity"}, &firstMultiplicity);
    SetFirstExistingBranch(events, {"firstInteractionSecondaryKE_MeV"}, &firstSecondaryKE);

    SetFirstExistingBranch(events, {"totalTrackLength_mm"}, &totalTrackLength);
    SetFirstExistingBranch(events, {"maxStepEdep_MeV"}, &maxStepEdep);
    SetFirstExistingBranch(events, {"depthCentroid_mm"}, &depthCentroid);
    SetFirstExistingBranch(events, {"depthRMS_mm"}, &depthRMS);
    SetFirstExistingBranch(events, {"radialCentroid_mm"}, &radialCentroid);
    SetFirstExistingBranch(events, {"radialRMS_mm"}, &radialRMS);
    SetFirstExistingBranch(events, {"edepFrontThird_MeV"}, &edepFront);
    SetFirstExistingBranch(events, {"edepMiddleThird_MeV"}, &edepMiddle);
    SetFirstExistingBranch(events, {"edepRearThird_MeV"}, &edepRear);

    SetFirstExistingBranch(events, {"firstHitTime_ns"}, &firstHitTime);
    SetFirstExistingBranch(events, {"lastHitTime_ns"}, &lastHitTime);
    SetFirstExistingBranch(events, {"timeSpan_ns"}, &timeSpan);
    SetFirstExistingBranch(events, {"promptEdep10ns_MeV"}, &promptEdep);
    SetFirstExistingBranch(events, {"delayedEdepAfter10ns_MeV"}, &delayedEdep);

    SetFirstExistingBranch(events, {"leakKE_total_MeV"}, &leakTotal);
    SetFirstExistingBranch(events, {"leakKE_np_MeV"}, &leakNP);
    SetFirstExistingBranch(events, {"leakKE_em_MeV"}, &leakEM);
    SetFirstExistingBranch(events, {"leakKE_lightIon_MeV"}, &leakLight);
    SetFirstExistingBranch(events, {"leakKE_heavyIon_MeV"}, &leakHeavy);
    SetFirstExistingBranch(events, {"leakKE_piK_MeV"}, &leakPiK);
    const bool haveEscape = SetFirstExistingBranch(
        events, {"primaryNeutronEscaped"}, &primaryEscaped);
    SetFirstExistingBranch(events, {"primaryNeutronExitKE_MeV"}, &primaryExitKE);

    SetFirstExistingBranch(events, {"maxCrystalID"}, &maxCrystalID);
    SetFirstExistingBranch(events, {"maxCrystalEdep_MeV"}, &maxCrystalEdep);
    SetFirstExistingBranch(events, {"nHitCrystals"}, &nHitCrystals);
    SetFirstExistingBranch(events, {"edepLayer0_MeV"}, &edepLayer0);
    SetFirstExistingBranch(events, {"edepLayer1_MeV"}, &edepLayer1);

    if (!havePrimary || !haveTotal || !haveAnyHadronic) return;
    if (!haveWeight) weight = 1.0;

    events->SetCacheSize(64LL * 1024LL * 1024LL);
    events->AddBranchToCache("*", true);

    // -------------------------------------------------------------------------
    // Histograms and binned statistics filled from the event tree
    // -------------------------------------------------------------------------
    auto* hTotalEdep = MakeLogH1("hTotalEdep_analysis", 220, 1.0e-3, 1.0e4);
    auto* hEdepNP = MakeLogH1("hEdepNP_analysis", 220, 1.0e-3, 1.0e4);
    auto* hEdepHeavy = MakeLogH1("hEdepHeavy_analysis", 220, 1.0e-3, 1.0e4);
    auto* hEdepLight = MakeLogH1("hEdepLight_analysis", 220, 1.0e-3, 1.0e4);
    auto* hEdepEM = MakeLogH1("hEdepEM_analysis", 220, 1.0e-3, 1.0e4);
    auto* hEdepPiK = MakeLogH1("hEdepPiK_analysis", 220, 1.0e-3, 1.0e4);
    auto* hPrimarySpectrum = MakeLogH1("hPrimarySpectrum_analysis", 150, 0.1, 10.0);
    auto* hResponse = MakeLogH1("hResponse_analysis", 200, 1.0e-6, 2.0);
    auto* h2PrimaryEdep = MakeLogLogH2(
        "h2PrimaryEdep_analysis", 90, 0.1, 10.0, 120, 1.0e-3, 1.0e4);

    auto* hFirstDepth = new TH1D("hFirstDepth_analysis", "", 120, 0.0, kBGOEnvelopeDepthMm);
    hFirstDepth->SetDirectory(nullptr); hFirstDepth->Sumw2();
    auto* h2FirstDepth = MakeLogLinearH2(
        "h2FirstDepth_analysis", 80, 0.1, 10.0, 120, 0.0, kBGOEnvelopeDepthMm);

    auto* hNSecondaries = new TH1D("hNSecondaries_analysis", "", 400, -0.5, 399.5);
    auto* hNUniqueTracks = new TH1D("hNUniqueTracks_analysis", "", 500, -0.5, 499.5);
    auto* hNHitCrystals = new TH1D("hNHitCrystals_analysis", "", 33, -0.5, 32.5);
    hNSecondaries->SetDirectory(nullptr); hNSecondaries->Sumw2();
    hNUniqueTracks->SetDirectory(nullptr); hNUniqueTracks->Sumw2();
    hNHitCrystals->SetDirectory(nullptr); hNHitCrystals->Sumw2();

    auto* hTimeSpan = MakeLogH1("hTimeSpan_analysis", 180, 1.0e-4, 1.0e6);
    auto* hLeakage = MakeLogH1("hLeakage_analysis", 220, 1.0e-3, 1.0e4);
    auto* hPrimaryExitKE = MakeLogH1("hPrimaryExitKE_analysis", 220, 1.0e-3, 1.0e4);
    const auto edepCorrelationEdges = LogEdges(120, 1.0e-3, 1.0e4);
    auto* h2EdepVsMultiplicity = new TH2D(
        "h2EdepVsMultiplicity_analysis", "",
        120, edepCorrelationEdges.data(),
        300, -0.5, 599.5);
    h2EdepVsMultiplicity->SetDirectory(nullptr);
    h2EdepVsMultiplicity->Sumw2();

    BinnedStats statEdepAll(energyEdges);
    BinnedStats statEdepHadronic(energyEdges);
    BinnedStats statEdepNonHadronic(energyEdges);
    BinnedStats statResponse(energyEdges);

    BinnedStats statFracNP(energyEdges);
    BinnedStats statFracEM(energyEdges);
    BinnedStats statFracHeavy(energyEdges);
    BinnedStats statFracLight(energyEdges);
    BinnedStats statFracPiK(energyEdges);

    BinnedStats statMultN(energyEdges);
    BinnedStats statMultP(energyEdges);
    BinnedStats statMultEM(energyEdges);
    BinnedStats statMultPionKaon(energyEdges);
    BinnedStats statMultLight(energyEdges);
    BinnedStats statMultHeavy(energyEdges);
    BinnedStats statMultTotal(energyEdges);
    BinnedStats statUniqueTracks(energyEdges);
    BinnedStats statHitCrystals(energyEdges);

    BinnedStats statFirstDepth(energyEdges);
    BinnedStats statDepthCentroid(energyEdges);
    BinnedStats statDepthRMS(energyEdges);
    BinnedStats statRadialRMS(energyEdges);
    BinnedStats statFrontEdep(energyEdges);
    BinnedStats statMiddleEdep(energyEdges);
    BinnedStats statRearEdep(energyEdges);

    BinnedStats statPrompt(energyEdges);
    BinnedStats statDelayed(energyEdges);
    BinnedStats statPromptFraction(energyEdges);
    BinnedStats statTimeSpan(energyEdges);

    BinnedStats statLeakage(energyEdges);
    BinnedStats statLeakageFraction(energyEdges);
    BinnedStats statPrimaryExitFraction(energyEdges);
    BinnedStats statLayer0(energyEdges);
    BinnedStats statLayer1(energyEdges);

    BinnedFraction probAnyHadronic(energyEdges);
    BinnedFraction probPrimaryHadronic(energyEdges);
    BinnedFraction probPrimaryInelastic(energyEdges);
    BinnedFraction probPrimaryEscape(energyEdges);

    const std::array<double, 5> thresholdsMeV = {{1.0, 5.0, 10.0, 50.0, 100.0}};
    std::vector<BinnedFraction> detectionEfficiencies;
    detectionEfficiencies.reserve(thresholdsMeV.size());
    for (std::size_t i = 0; i < thresholdsMeV.size(); ++i)
        detectionEfficiencies.emplace_back(energyEdges);

    double sumW = 0.0;
    double sumW2 = 0.0;
    double sumEdepTotal = 0.0;
    double sumEdepNP = 0.0;
    double sumEdepEM = 0.0;
    double sumEdepHeavy = 0.0;
    double sumEdepLight = 0.0;
    double sumEdepPiK = 0.0;
    double sumEdepOther = 0.0;
    double sumLeakage = 0.0;
    double sumPrompt = 0.0;
    double sumDelayed = 0.0;
    double sumSecondaryCount = 0.0;
    double sumUniqueTrackCount = 0.0;

    double weightedAnyHadronic = 0.0;
    double weightedPrimaryHadronic = 0.0;
    double weightedPrimaryInelastic = 0.0;
    double weightedPrimaryEscape = 0.0;
    Long64_t rawAnyHadronic = 0;
    Long64_t rawPrimaryHadronic = 0;
    Long64_t rawPrimaryInelastic = 0;
    Long64_t rawPrimaryEscape = 0;

    const Long64_t nEntries = events->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        events->GetEntry(i);
        if (!haveWeight) weight = 1.0;
        if (!(weight > 0.0) || !std::isfinite(weight)) continue;
        if (!(primaryE > 0.0) || !std::isfinite(primaryE)) continue;

        const bool strictHadronic = havePrimaryHadronic ? (primaryHadronic != 0)
                                                        : (anyHadronic != 0);
        const bool strictInelastic = havePrimaryInelastic ? (primaryInelastic != 0)
                                                          : strictHadronic;
        const double response = haveStoredFraction
            ? edepFraction
            : edepTotal / (1000.0 * primaryE);

        sumW += weight;
        sumW2 += weight * weight;
        sumEdepTotal += edepTotal * weight;
        sumEdepNP += edepNP * weight;
        sumEdepEM += edepEM * weight;
        sumEdepHeavy += edepHeavy * weight;
        sumEdepLight += edepLight * weight;
        sumEdepPiK += edepPiK * weight;
        sumEdepOther += edepOther * weight;

        hPrimarySpectrum->Fill(primaryE, weight);
        if (edepTotal > 0.0) {
            hTotalEdep->Fill(edepTotal, weight);
            h2PrimaryEdep->Fill(primaryE, edepTotal, weight);
        }
        if (edepNP > 0.0) hEdepNP->Fill(edepNP, weight);
        if (edepEM > 0.0) hEdepEM->Fill(edepEM, weight);
        if (edepHeavy > 0.0) hEdepHeavy->Fill(edepHeavy, weight);
        if (edepLight > 0.0) hEdepLight->Fill(edepLight, weight);
        if (edepPiK > 0.0) hEdepPiK->Fill(edepPiK, weight);
        if (response > 0.0) hResponse->Fill(response, weight);

        statEdepAll.Fill(primaryE, edepTotal, weight);
        if (strictHadronic) statEdepHadronic.Fill(primaryE, edepTotal, weight);
        else statEdepNonHadronic.Fill(primaryE, edepTotal, weight);
        statResponse.Fill(primaryE, response, weight);

        if (edepTotal > 0.0) {
            statFracNP.Fill(primaryE, edepNP / edepTotal, weight);
            statFracEM.Fill(primaryE, edepEM / edepTotal, weight);
            statFracHeavy.Fill(primaryE, edepHeavy / edepTotal, weight);
            statFracLight.Fill(primaryE, edepLight / edepTotal, weight);
            statFracPiK.Fill(primaryE, edepPiK / edepTotal, weight);
        }

        probAnyHadronic.Fill(primaryE, anyHadronic != 0, weight);
        probPrimaryHadronic.Fill(primaryE, strictHadronic, weight);
        probPrimaryInelastic.Fill(primaryE, strictInelastic, weight);

        if (anyHadronic) {
            weightedAnyHadronic += weight;
            ++rawAnyHadronic;
        }
        if (strictHadronic) {
            weightedPrimaryHadronic += weight;
            ++rawPrimaryHadronic;
        }
        if (strictInelastic) {
            weightedPrimaryInelastic += weight;
            ++rawPrimaryInelastic;
        }

        for (std::size_t it = 0; it < thresholdsMeV.size(); ++it) {
            detectionEfficiencies[it].Fill(
                primaryE, edepTotal >= thresholdsMeV[it], weight);
        }

        if (haveExtended) {
            hNSecondaries->Fill(nSecondaries, weight);
            hNUniqueTracks->Fill(nUniqueTracks, weight);
            hNHitCrystals->Fill(nHitCrystals, weight);
            if (edepTotal > 0.0)
                h2EdepVsMultiplicity->Fill(edepTotal, nSecondaries, weight);

            statMultN.Fill(primaryE, nSecN, weight);
            statMultP.Fill(primaryE, nSecP, weight);
            statMultEM.Fill(primaryE, nSecGamma + nSecElectron + nSecPositron, weight);
            statMultPionKaon.Fill(primaryE, nSecPion + nSecKaon, weight);
            statMultLight.Fill(primaryE, nSecLightIon, weight);
            statMultHeavy.Fill(primaryE, nSecHeavyIon, weight);
            statMultTotal.Fill(primaryE, nSecondaries, weight);
            statUniqueTracks.Fill(primaryE, nUniqueTracks, weight);
            statHitCrystals.Fill(primaryE, nHitCrystals, weight);
            sumSecondaryCount += nSecondaries * weight;
            sumUniqueTrackCount += nUniqueTracks * weight;

            if (firstDepth >= 0.0) {
                hFirstDepth->Fill(firstDepth, weight);
                h2FirstDepth->Fill(primaryE, firstDepth, weight);
                statFirstDepth.Fill(primaryE, firstDepth, weight);
            }
            if (depthCentroid >= 0.0)
                statDepthCentroid.Fill(primaryE, depthCentroid, weight);
            if (depthRMS >= 0.0)
                statDepthRMS.Fill(primaryE, depthRMS, weight);
            if (radialRMS >= 0.0)
                statRadialRMS.Fill(primaryE, radialRMS, weight);

            statFrontEdep.Fill(primaryE, edepFront, weight);
            statMiddleEdep.Fill(primaryE, edepMiddle, weight);
            statRearEdep.Fill(primaryE, edepRear, weight);

            statPrompt.Fill(primaryE, promptEdep, weight);
            statDelayed.Fill(primaryE, delayedEdep, weight);
            statPromptFraction.Fill(
                primaryE, edepTotal > 0.0 ? promptEdep / edepTotal : 0.0, weight);
            if (timeSpan > 0.0) {
                hTimeSpan->Fill(timeSpan, weight);
                statTimeSpan.Fill(primaryE, timeSpan, weight);
            }
            sumPrompt += promptEdep * weight;
            sumDelayed += delayedEdep * weight;

            if (leakTotal > 0.0) hLeakage->Fill(leakTotal, weight);
            statLeakage.Fill(primaryE, leakTotal, weight);
            statLeakageFraction.Fill(
                primaryE, leakTotal / (1000.0 * primaryE), weight);
            sumLeakage += leakTotal * weight;

            probPrimaryEscape.Fill(primaryE, primaryEscaped != 0, weight);
            if (primaryEscaped) {
                weightedPrimaryEscape += weight;
                ++rawPrimaryEscape;
                if (primaryExitKE > 0.0) hPrimaryExitKE->Fill(primaryExitKE, weight);
            }
            statPrimaryExitFraction.Fill(
                primaryE, primaryExitKE / (1000.0 * primaryE), weight);

            statLayer0.Fill(primaryE, edepLayer0, weight);
            statLayer1.Fill(primaryE, edepLayer1, weight);
        }

        if ((i + 1) % 1000000 == 0 || i + 1 == nEntries) {
            std::cout << "  processed " << (i + 1) << " / " << nEntries << " events\r"
                      << std::flush;
        }
    }
    std::cout << "\n";
    events->ResetBranchAddresses();

    if (!(sumW > 0.0)) {
        std::cerr << "No valid weighted events were found.\n";
        return;
    }

    // -------------------------------------------------------------------------
    // First-interaction secondary tree
    // -------------------------------------------------------------------------
    auto* firstSecondaries = dynamic_cast<TTree*>(file->Get("FirstInteractionSecondaries"));
    std::array<TH1D*, 7> hSecondaryKE = {{nullptr, nullptr, nullptr, nullptr,
                                          nullptr, nullptr, nullptr}};
    for (int i = 0; i < 7; ++i) {
        hSecondaryKE[static_cast<std::size_t>(i)] = MakeLogH1(
            TString::Format("hSecondaryKE_class%d", i), 180, 1.0e-3, 1.0e4);
    }
    auto* hFragmentA = new TH1D("hFragmentA_analysis", "", 220, 0.5, 220.5);
    auto* hFragmentZ = new TH1D("hFragmentZ_analysis", "", 90, 0.5, 90.5);
    auto* hFragmentZA = new TH2D(
        "hFragmentZA_analysis", "", 90, 0.5, 90.5, 220, 0.5, 220.5);
    hFragmentA->SetDirectory(nullptr); hFragmentA->Sumw2();
    hFragmentZ->SetDirectory(nullptr); hFragmentZ->Sumw2();
    hFragmentZA->SetDirectory(nullptr); hFragmentZA->Sumw2();

    Long64_t nSecondaryRows = 0;
    if (firstSecondaries) {
        Double_t secWeight = 1.0;
        Double_t secKE = 0.0;
        Int_t secClass = 0;
        Int_t secZ = 0;
        Int_t secA = 0;
        Int_t isFragment = 0;

        firstSecondaries->SetBranchStatus("*", 0);
        SetFirstExistingBranch(firstSecondaries, {"weight"}, &secWeight);
        SetFirstExistingBranch(firstSecondaries, {"kineticEnergy_MeV"}, &secKE, true);
        SetFirstExistingBranch(firstSecondaries, {"particleClass"}, &secClass, true);
        SetFirstExistingBranch(firstSecondaries, {"Z"}, &secZ);
        SetFirstExistingBranch(firstSecondaries, {"A"}, &secA);
        SetFirstExistingBranch(firstSecondaries, {"isFragment"}, &isFragment);
        firstSecondaries->SetCacheSize(64LL * 1024LL * 1024LL);
        firstSecondaries->AddBranchToCache("*", true);

        nSecondaryRows = firstSecondaries->GetEntries();
        for (Long64_t i = 0; i < nSecondaryRows; ++i) {
            firstSecondaries->GetEntry(i);
            if (secClass >= 0 && secClass < 7 && secKE > 0.0)
                hSecondaryKE[static_cast<std::size_t>(secClass)]->Fill(secKE, secWeight);
            if (isFragment && secA > 0 && secZ > 0) {
                hFragmentA->Fill(secA, secWeight);
                hFragmentZ->Fill(secZ, secWeight);
                hFragmentZA->Fill(secZ, secA, secWeight);
            }
        }
        firstSecondaries->ResetBranchAddresses();
    }

    // Existing step-level energy-weighted histograms from the extended scorer.
    std::unique_ptr<TH2> h2DepthProfile(GetDetachedH2(file.get(), "PrimaryE_vs_DepthProfile"));
    std::unique_ptr<TH2> h2RadialProfile(GetDetachedH2(file.get(), "PrimaryE_vs_RadialProfile"));
    std::unique_ptr<TH2> h2CrystalEdep(GetDetachedH2(file.get(), "PrimaryE_vs_CrystalID"));

    // -------------------------------------------------------------------------
    // Derived histograms
    // -------------------------------------------------------------------------
    std::unique_ptr<TH1D> hProbAny(FractionHistogram(
        "hProbAny", probAnyHadronic,
        "Primary neutron energy (GeV)", "Interaction probability per incident neutron"));
    std::unique_ptr<TH1D> hProbPrimary(FractionHistogram(
        "hProbPrimary", probPrimaryHadronic,
        "Primary neutron energy (GeV)", "Interaction probability per incident neutron"));
    std::unique_ptr<TH1D> hProbInelastic(FractionHistogram(
        "hProbInelastic", probPrimaryInelastic,
        "Primary neutron energy (GeV)", "Interaction probability per incident neutron"));

    std::unique_ptr<TH1D> hMeanEdepAll(MeanHistogram(
        "hMeanEdepAll", statEdepAll,
        "Primary neutron energy (GeV)", "Mean BGO E_{dep} per incident neutron (MeV)"));
    std::unique_ptr<TH1D> hMeanEdepHad(MeanHistogram(
        "hMeanEdepHad", statEdepHadronic,
        "Primary neutron energy (GeV)", "Mean BGO E_{dep} per incident neutron (MeV)"));
    std::unique_ptr<TH1D> hMeanEdepNonHad(MeanHistogram(
        "hMeanEdepNonHad", statEdepNonHadronic,
        "Primary neutron energy (GeV)", "Mean BGO E_{dep} per incident neutron (MeV)"));
    std::unique_ptr<TH1D> hMeanResponse(MeanHistogram(
        "hMeanResponse", statResponse,
        "Primary neutron energy (GeV)", "Mean response fraction E_{dep}/E_{n}"));
    std::unique_ptr<TH1D> hResponseResolution(RelativeRMSHistogram(
        "hResponseResolution", statResponse,
        "Primary neutron energy (GeV)", "Relative response width RMS/#LTresponse#GT"));

    std::unique_ptr<TH1D> hFracNP(MeanHistogram(
        "hFracNP", statFracNP, "Primary neutron energy (GeV)",
        "Mean fraction of total BGO E_{dep}"));
    std::unique_ptr<TH1D> hFracEM(MeanHistogram(
        "hFracEM", statFracEM, "Primary neutron energy (GeV)",
        "Mean fraction of total BGO E_{dep}"));
    std::unique_ptr<TH1D> hFracHeavy(MeanHistogram(
        "hFracHeavy", statFracHeavy, "Primary neutron energy (GeV)",
        "Mean fraction of total BGO E_{dep}"));
    std::unique_ptr<TH1D> hFracLight(MeanHistogram(
        "hFracLight", statFracLight, "Primary neutron energy (GeV)",
        "Mean fraction of total BGO E_{dep}"));
    std::unique_ptr<TH1D> hFracPiK(MeanHistogram(
        "hFracPiK", statFracPiK, "Primary neutron energy (GeV)",
        "Mean fraction of total BGO E_{dep}"));

    std::array<std::unique_ptr<TH1D>, 5> hDetectionEfficiency;
    for (std::size_t i = 0; i < thresholdsMeV.size(); ++i) {
        hDetectionEfficiency[i].reset(FractionHistogram(
            TString::Format("hDetectionEff_%zu", i), detectionEfficiencies[i],
            "Primary neutron energy (GeV)", "Detection efficiency"));
    }

    // -------------------------------------------------------------------------
    // PLOT 01: Total BGO Edep and species components
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c01", "", 880, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();

        StyleLine(hTotalEdep, kBlack, 1, 3);
        PrepareLinePlot(hTotalEdep,
            "Total BGO energy deposition (MeV)",
            "Weighted events / logarithmic bin");
        hTotalEdep->SetMinimum(0.5);
        hTotalEdep->Draw("hist");

        struct Curve { TH1* h; Color_t color; const char* label; };
        const Curve curves[] = {
            {hEdepNP,    kBlue + 1,    "n + p"},
            {hEdepHeavy, kRed + 1,     "heavy ions (A > 4)"},
            {hEdepEM,    kGreen + 2,   "e^{#pm} + #gamma"},
            {hEdepPiK,   kOrange + 1,  "#pi / K"},
            {hEdepLight, kMagenta + 1, "light ions (A #leq 4)"}
        };
        auto* legend = new TLegend(0.52, 0.15, 0.93, 0.47);
        legend->AddEntry(hTotalEdep, "Total BGO", "l");
        for (const auto& curve : curves) {
            StyleLine(curve.h, curve.color, 2, 2);
            curve.h->Draw("hist same");
            legend->AddEntry(curve.h, curve.label, "l");
        }
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_01_total_edep");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 02: Primary energy vs BGO Edep
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c02", "", 860, 680);
        c->SetRightMargin(0.16);
        c->SetLogx(); c->SetLogy(); c->SetLogz();
        c->SetGrid();
        h2PrimaryEdep->GetXaxis()->SetTitle("Primary neutron energy (GeV)");
        h2PrimaryEdep->GetYaxis()->SetTitle("Total BGO E_{dep} (MeV)");
        h2PrimaryEdep->GetZaxis()->SetTitle("Weighted events");
        h2PrimaryEdep->GetYaxis()->SetTitleOffset(1.45);
        h2PrimaryEdep->GetZaxis()->SetTitleOffset(1.25);
        FixLogAxis(h2PrimaryEdep->GetXaxis());
        FixLogAxis(h2PrimaryEdep->GetYaxis());
        h2PrimaryEdep->Draw("colz");
        SaveCanvas(c, outDir, "bgo_02_primaryE_vs_edep");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 03: Hadronic probabilities. The analytical line now uses the actual
    // 6-cm BGO material depth and the BGO nuclear interaction length 22.32 cm.
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c03", "", 860, 680);
        c->SetLogx(); c->SetGrid();

        StyleLine(hProbAny.get(), kGray + 2, 2, 2, 20);
        StyleLine(hProbPrimary.get(), kBlue + 1, 1, 2, 21);
        StyleLine(hProbInelastic.get(), kRed + 1, 1, 2, 22);
        hProbAny->GetYaxis()->SetRangeUser(0.0, 0.50);
        hProbAny->GetYaxis()->SetTitleOffset(1.45);
        hProbAny->Draw("E1");
        hProbPrimary->Draw("E1 same");
        hProbInelastic->Draw("E1 same");

        auto* expected = new TLine(0.1, expectedInteractionProbability,
                                   10.0, expectedInteractionProbability);
        expected->SetLineColor(kBlue - 7);
        expected->SetLineStyle(3);
        expected->SetLineWidth(3);
        expected->Draw();

        auto* legend = new TLegend(0.17, 0.68, 0.92, 0.92);
        legend->AddEntry(hProbAny.get(), "Any selected hadronic process", "lp");
        legend->AddEntry(hProbPrimary.get(), "Primary-neutron hadronic process", "lp");
        legend->AddEntry(hProbInelastic.get(), "Primary-neutron inelastic process", "lp");
        legend->AddEntry(expected,
            TString::Format("1-exp(-%.1f/%.2f) = %.3f",
                            kBGOActiveDepthCm, kBGOInteractionLengthCm,
                            expectedInteractionProbability), "l");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_03_hadronic_probability");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 04: Primary spectrum
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c04", "", 860, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();
        StyleLine(hPrimarySpectrum, kBlue + 1, 1, 2);
        PrepareLinePlot(hPrimarySpectrum,
            "Primary neutron energy (GeV)",
            "Incident neutrons / logarithmic bin");
        hPrimarySpectrum->SetMinimum(0.5);
        hPrimarySpectrum->Draw("hist");
        SaveCanvas(c, outDir, "bgo_04_primary_spectrum");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 05: Species Edep spectra
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c05", "", 880, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();
        struct Curve { TH1* h; Color_t color; const char* label; };
        const Curve curves[] = {
            {hEdepNP,    kBlue + 1,    "n + p"},
            {hEdepEM,    kGreen + 2,   "e^{#pm} + #gamma"},
            {hEdepHeavy, kRed + 1,     "heavy ions (A > 4)"},
            {hEdepLight, kMagenta + 1, "light ions (A #leq 4)"},
            {hEdepPiK,   kOrange + 1,  "#pi / K"}
        };
        bool first = true;
        double maximum = 1.0;
        for (const auto& curve : curves) maximum = std::max(maximum, curve.h->GetMaximum());
        auto* legend = new TLegend(0.17, 0.17, 0.58, 0.47);
        for (const auto& curve : curves) {
            StyleLine(curve.h, curve.color, 1, 2);
            PrepareLinePlot(curve.h,
                "BGO E_{dep} by species (MeV)",
                "Weighted events / logarithmic bin");
            curve.h->GetYaxis()->SetRangeUser(0.5, maximum * 3.0);
            curve.h->Draw(first ? "hist" : "hist same");
            legend->AddEntry(curve.h, curve.label, "l");
            first = false;
        }
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_05_edep_by_species");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 06: Mean species fraction vs primary energy
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c06", "", 880, 680);
        c->SetLogx(); c->SetGrid();
        struct Curve { TH1* h; Color_t color; const char* label; };
        const Curve curves[] = {
            {hFracNP.get(),    kBlue + 1,    "n + p"},
            {hFracEM.get(),    kGreen + 2,   "e^{#pm} + #gamma"},
            {hFracHeavy.get(), kRed + 1,     "heavy ions (A > 4)"},
            {hFracLight.get(), kMagenta + 1, "light ions (A #leq 4)"},
            {hFracPiK.get(),   kOrange + 1,  "#pi / K"}
        };
        bool first = true;
        auto* legend = new TLegend(0.52, 0.55, 0.93, 0.86);
        for (const auto& curve : curves) {
            StyleLine(curve.h, curve.color, 1, 2);
            curve.h->GetYaxis()->SetRangeUser(0.0, 1.02);
            curve.h->GetYaxis()->SetTitleOffset(1.45);
            FixLogAxis(curve.h->GetXaxis());
            curve.h->Draw(first ? "hist" : "hist same");
            legend->AddEntry(curve.h, curve.label, "l");
            first = false;
        }
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_06_edep_fraction_vs_energy");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 07: Mean Edep per incident neutron by species
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c07", "", 940, 680);
        c->SetLogy(); c->SetGrid();
        c->SetBottomMargin(0.17);
        constexpr int nSpecies = 6;
        const char* labels[nSpecies] = {
            "Total", "n+p", "e^{#pm}/#gamma", "Heavy ions", "Light ions", "#pi/K"
        };
        const double means[nSpecies] = {
            sumEdepTotal / sumW,
            sumEdepNP / sumW,
            sumEdepEM / sumW,
            sumEdepHeavy / sumW,
            sumEdepLight / sumW,
            sumEdepPiK / sumW
        };
        auto* hBar = new TH1D("hMeanSpeciesBar", "", nSpecies, 0.0, nSpecies);
        hBar->SetDirectory(nullptr);
        for (int i = 0; i < nSpecies; ++i) {
            hBar->SetBinContent(i + 1, means[i]);
            hBar->GetXaxis()->SetBinLabel(i + 1, labels[i]);
        }
        hBar->GetYaxis()->SetTitle("Mean BGO E_{dep} per incident neutron (MeV)");
        hBar->GetYaxis()->SetTitleOffset(1.45);
        hBar->GetXaxis()->SetLabelSize(0.046);
        hBar->SetFillColor(kBlue + 1);
        hBar->SetLineColor(kBlue + 2);
        hBar->SetBarWidth(0.72);
        hBar->SetBarOffset(0.14);
        hBar->SetMinimum(1.0e-4);
        hBar->SetMaximum(std::max(1.0, means[0] * 8.0));
        hBar->Draw("bar2");

        TLatex latex;
        latex.SetTextAlign(22);
        latex.SetTextSize(0.028);
        for (int i = 0; i < nSpecies; ++i) {
            latex.DrawLatex(i + 0.5, means[i] * 1.35,
                            TString::Format("%.3g", means[i]));
        }
        SaveCanvas(c, outDir, "bgo_07_mean_edep_by_species");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 08: Conditional calorimeter response
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c08", "", 880, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();
        StyleLine(hMeanEdepAll.get(), kBlack, 1, 3, 20);
        StyleLine(hMeanEdepHad.get(), kRed + 1, 1, 2, 21);
        StyleLine(hMeanEdepNonHad.get(), kBlue + 1, 1, 2, 22);
        hMeanEdepAll->SetMinimum(1.0e-3);
        hMeanEdepAll->Draw("E1");
        hMeanEdepHad->Draw("E1 same");
        hMeanEdepNonHad->Draw("E1 same");
        auto* legend = new TLegend(0.17, 0.70, 0.64, 0.91);
        legend->AddEntry(hMeanEdepAll.get(), "All incident neutrons", "lp");
        legend->AddEntry(hMeanEdepHad.get(), "Primary-neutron hadronic events", "lp");
        legend->AddEntry(hMeanEdepNonHad.get(), "No primary-neutron hadronic event", "lp");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_08_conditional_mean_response");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 09: Mean response fraction and relative response width
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c09", "", 880, 680);
        c->SetLogx(); c->SetGrid();
        StyleLine(hMeanResponse.get(), kBlue + 1, 1, 2, 20);
        StyleLine(hResponseResolution.get(), kRed + 1, 1, 2, 21);
        hMeanResponse->GetYaxis()->SetTitle("Dimensionless response metric");
        hMeanResponse->GetYaxis()->SetTitleOffset(1.45);
        hMeanResponse->SetMinimum(0.0);
        const double ymax = std::max(hMeanResponse->GetMaximum(),
                                     hResponseResolution->GetMaximum());
        hMeanResponse->SetMaximum(std::max(0.1, 1.25 * ymax));
        hMeanResponse->Draw("E1");
        hResponseResolution->Draw("E1 same");
        auto* legend = new TLegend(0.50, 0.75, 0.93, 0.91);
        legend->AddEntry(hMeanResponse.get(), "#LT E_{dep}/E_{n} #GT", "lp");
        legend->AddEntry(hResponseResolution.get(), "RMS(response)/#LTresponse#GT", "lp");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_09_response_fraction_and_width");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 10: Detection efficiency for several deposited-energy thresholds
    // -------------------------------------------------------------------------
    {
        auto* c = new TCanvas("bgo_c10", "", 880, 680);
        c->SetLogx(); c->SetGrid();
        const std::array<Color_t, 5> colors = {{
            kBlue + 1, kGreen + 2, kOrange + 1, kMagenta + 1, kRed + 1
        }};
        bool first = true;
        auto* legend = new TLegend(0.59, 0.57, 0.93, 0.88);
        for (std::size_t i = 0; i < thresholdsMeV.size(); ++i) {
            auto* h = hDetectionEfficiency[i].get();
            StyleLine(h, colors[i], 1, 2, 20 + static_cast<int>(i));
            h->GetYaxis()->SetRangeUser(0.0, 1.02);
            h->GetYaxis()->SetTitleOffset(1.45);
            FixLogAxis(h->GetXaxis());
            h->Draw(first ? "E1" : "E1 same");
            legend->AddEntry(h,
                TString::Format("E_{dep} #geq %.0f MeV", thresholdsMeV[i]), "lp");
            first = false;
        }
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_10_detection_efficiency");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 11: First primary-neutron hadronic interaction depth
    // -------------------------------------------------------------------------
    double fittedLambdaMm = -1.0;
    double fittedLambdaErrorMm = -1.0;
    {
        auto* c = new TCanvas("bgo_c11", "", 880, 680);
        c->SetGrid();
        StyleLine(hFirstDepth, kBlue + 1, 1, 2);
        hFirstDepth->GetXaxis()->SetTitle("First interaction depth from front BGO face (mm)");
        hFirstDepth->GetYaxis()->SetTitle("Weighted primary-neutron interactions / bin");
        hFirstDepth->GetYaxis()->SetTitleOffset(1.45);
        hFirstDepth->SetMinimum(0.0);
        hFirstDepth->Draw("E1");

        TF1* fit = nullptr;
        if (haveFirstDepth && hFirstDepth->GetEntries() > 100.0) {
            fit = new TF1("bgo_depth_exponential", "[0]*exp(-x/[1])", 0.0, 60.0);
            fit->SetParameters(hFirstDepth->GetMaximum(),
                               10.0 * kBGOInteractionLengthCm);
            fit->SetParLimits(1, 10.0, 1000.0);
            fit->SetLineColor(kRed + 1);
            fit->SetLineWidth(3);
            hFirstDepth->Fit(fit, "QNR");
            fit->Draw("same");
            fittedLambdaMm = fit->GetParameter(1);
            fittedLambdaErrorMm = fit->GetParError(1);
        }

        auto* legend = new TLegend(0.48, 0.73, 0.93, 0.90);
        legend->AddEntry(hFirstDepth, "Simulation", "lep");
        if (fit) {
            legend->AddEntry(fit,
                TString::Format("Exponential fit: #lambda = %.1f #pm %.1f mm",
                                fittedLambdaMm, fittedLambdaErrorMm), "l");
        }
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_11_first_interaction_depth");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 12: First interaction depth vs primary energy
    // -------------------------------------------------------------------------
    if (haveFirstDepth) {
        auto* c = new TCanvas("bgo_c12", "", 860, 680);
        c->SetRightMargin(0.16);
        c->SetLogx(); c->SetLogz(); c->SetGrid();
        h2FirstDepth->GetXaxis()->SetTitle("Primary neutron energy (GeV)");
        h2FirstDepth->GetYaxis()->SetTitle("First interaction depth (mm)");
        h2FirstDepth->GetZaxis()->SetTitle("Weighted events");
        h2FirstDepth->GetYaxis()->SetTitleOffset(1.45);
        h2FirstDepth->GetZaxis()->SetTitleOffset(1.25);
        FixLogAxis(h2FirstDepth->GetXaxis());
        h2FirstDepth->Draw("colz");
        SaveCanvas(c, outDir, "bgo_12_first_depth_vs_primary_energy");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 13: Secondary and unique-track multiplicity distributions
    // -------------------------------------------------------------------------
    if (haveExtended) {
        auto* c = new TCanvas("bgo_c13", "", 880, 680);
        c->SetLogy(); c->SetGrid();
        StyleLine(hNSecondaries, kBlue + 1, 1, 2);
        StyleLine(hNUniqueTracks, kRed + 1, 1, 2);
        hNSecondaries->GetXaxis()->SetTitle("Multiplicity per incident neutron");
        hNSecondaries->GetYaxis()->SetTitle("Weighted events / bin");
        hNSecondaries->GetYaxis()->SetTitleOffset(1.45);
        hNSecondaries->SetMinimum(0.5);
        hNSecondaries->Draw("hist");
        hNUniqueTracks->Draw("hist same");
        auto* legend = new TLegend(0.56, 0.76, 0.93, 0.90);
        legend->AddEntry(hNSecondaries, "Produced secondaries", "l");
        legend->AddEntry(hNUniqueTracks, "Unique tracks with BGO hits", "l");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_13_multiplicity_distributions");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 14: Mean secondary multiplicity by particle group
    // -------------------------------------------------------------------------
    if (haveExtended) {
        std::array<std::unique_ptr<TH1D>, 6> multiplicityMeans = {{
            std::unique_ptr<TH1D>(MeanHistogram("hMultN", statMultN,
                "Primary neutron energy (GeV)", "Mean produced multiplicity per incident neutron")),
            std::unique_ptr<TH1D>(MeanHistogram("hMultP", statMultP,
                "Primary neutron energy (GeV)", "Mean produced multiplicity per incident neutron")),
            std::unique_ptr<TH1D>(MeanHistogram("hMultEM", statMultEM,
                "Primary neutron energy (GeV)", "Mean produced multiplicity per incident neutron")),
            std::unique_ptr<TH1D>(MeanHistogram("hMultPiK", statMultPionKaon,
                "Primary neutron energy (GeV)", "Mean produced multiplicity per incident neutron")),
            std::unique_ptr<TH1D>(MeanHistogram("hMultLI", statMultLight,
                "Primary neutron energy (GeV)", "Mean produced multiplicity per incident neutron")),
            std::unique_ptr<TH1D>(MeanHistogram("hMultHI", statMultHeavy,
                "Primary neutron energy (GeV)", "Mean produced multiplicity per incident neutron"))
        }};
        const std::array<Color_t, 6> colors = {{
            kBlue + 1, kCyan + 2, kGreen + 2, kOrange + 1, kMagenta + 1, kRed + 1
        }};
        const char* labels[6] = {
            "neutrons", "protons", "#gamma + e^{#pm}", "#pi / K", "light ions", "heavy ions"
        };

        auto* c = new TCanvas("bgo_c14", "", 880, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();
        bool first = true;
        auto* legend = new TLegend(0.57, 0.57, 0.93, 0.89);
        for (std::size_t i = 0; i < multiplicityMeans.size(); ++i) {
            StyleLine(multiplicityMeans[i].get(), colors[i], 1, 2);
            multiplicityMeans[i]->SetMinimum(1.0e-4);
            multiplicityMeans[i]->Draw(first ? "hist" : "hist same");
            legend->AddEntry(multiplicityMeans[i].get(), labels[i], "l");
            first = false;
        }
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_14_secondary_multiplicity_vs_energy");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 15: Kinetic-energy spectra of first-interaction secondaries
    // -------------------------------------------------------------------------
    if (firstSecondaries && nSecondaryRows > 0) {
        auto* c = new TCanvas("bgo_c15", "", 880, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();
        struct Curve { int cls; Color_t color; const char* label; };
        const Curve curves[] = {
            {1, kBlue + 1, "neutrons"},
            {2, kCyan + 2, "protons"},
            {5, kGreen + 2, "e^{#pm} + #gamma"},
            {6, kOrange + 1, "#pi / K"},
            {3, kMagenta + 1, "light ions"},
            {4, kRed + 1, "heavy ions"}
        };
        double maximum = 1.0;
        for (const auto& curve : curves)
            maximum = std::max(maximum, hSecondaryKE[curve.cls]->GetMaximum());
        bool first = true;
        auto* legend = new TLegend(0.58, 0.56, 0.93, 0.89);
        for (const auto& curve : curves) {
            auto* h = hSecondaryKE[curve.cls];
            StyleLine(h, curve.color, 1, 2);
            PrepareLinePlot(h,
                "Kinetic energy at first primary-neutron interaction (MeV)",
                "Weighted secondaries / logarithmic bin");
            h->GetYaxis()->SetRangeUser(0.5, maximum * 3.0);
            h->Draw(first ? "hist" : "hist same");
            legend->AddEntry(h, curve.label, "l");
            first = false;
        }
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_15_first_interaction_secondary_spectra");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 16: Fragment mass distribution
    // -------------------------------------------------------------------------
    if (hFragmentA->GetEntries() > 0.0) {
        auto* c = new TCanvas("bgo_c16", "", 880, 680);
        c->SetLogy(); c->SetGrid();
        StyleLine(hFragmentA, kBlue + 1, 1, 2);
        hFragmentA->GetXaxis()->SetTitle("Fragment mass number A");
        hFragmentA->GetYaxis()->SetTitle("Weighted fragments / mass bin");
        hFragmentA->GetYaxis()->SetTitleOffset(1.45);
        hFragmentA->SetMinimum(0.5);
        hFragmentA->Draw("hist");
        SaveCanvas(c, outDir, "bgo_16_fragment_mass_distribution");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 17: Fragment charge distribution
    // -------------------------------------------------------------------------
    if (hFragmentZ->GetEntries() > 0.0) {
        auto* c = new TCanvas("bgo_c17", "", 880, 680);
        c->SetLogy(); c->SetGrid();
        StyleLine(hFragmentZ, kRed + 1, 1, 2);
        hFragmentZ->GetXaxis()->SetTitle("Fragment atomic number Z");
        hFragmentZ->GetYaxis()->SetTitle("Weighted fragments / charge bin");
        hFragmentZ->GetYaxis()->SetTitleOffset(1.45);
        hFragmentZ->SetMinimum(0.5);
        hFragmentZ->Draw("hist");
        SaveCanvas(c, outDir, "bgo_17_fragment_charge_distribution");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 18: Fragment Z-A map
    // -------------------------------------------------------------------------
    if (hFragmentZA->GetEntries() > 0.0) {
        auto* c = new TCanvas("bgo_c18", "", 900, 700);
        c->SetRightMargin(0.16);
        c->SetLogz(); c->SetGrid();
        hFragmentZA->GetXaxis()->SetTitle("Fragment atomic number Z");
        hFragmentZA->GetYaxis()->SetTitle("Fragment mass number A");
        hFragmentZA->GetZaxis()->SetTitle("Weighted fragments");
        hFragmentZA->GetYaxis()->SetTitleOffset(1.45);
        hFragmentZA->GetZaxis()->SetTitleOffset(1.25);
        hFragmentZA->Draw("colz");
        SaveCanvas(c, outDir, "bgo_18_fragment_ZA_map");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 19: Mean longitudinal shower profile per incident neutron
    // -------------------------------------------------------------------------
    if (h2DepthProfile) {
        std::unique_ptr<TH1D> profile(
            h2DepthProfile->ProjectionY("hDepthProfileProjection", 1,
                                        h2DepthProfile->GetNbinsX(), "e"));
        profile->SetDirectory(nullptr);
        profile->Scale(1.0 / sumW, "width");
        auto* c = new TCanvas("bgo_c19", "", 880, 680);
        c->SetGrid();
        StyleLine(profile.get(), kBlue + 1, 1, 2);
        profile->GetXaxis()->SetTitle("Depth from front BGO face (mm)");
        profile->GetYaxis()->SetTitle("Mean dE_{dep}/dz per incident neutron (MeV/mm)");
        profile->GetYaxis()->SetTitleOffset(1.55);
        profile->SetMinimum(0.0);
        profile->Draw("hist");
        SaveCanvas(c, outDir, "bgo_19_longitudinal_shower_profile");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 20: Mean radial shower profile per incident neutron
    // -------------------------------------------------------------------------
    if (h2RadialProfile) {
        std::unique_ptr<TH1D> profile(
            h2RadialProfile->ProjectionY("hRadialProfileProjection", 1,
                                         h2RadialProfile->GetNbinsX(), "e"));
        profile->SetDirectory(nullptr);
        profile->Scale(1.0 / sumW, "width");
        auto* c = new TCanvas("bgo_c20", "", 880, 680);
        c->SetGrid();
        StyleLine(profile.get(), kRed + 1, 1, 2);
        profile->GetXaxis()->SetTitle("Radius from incident trajectory (mm)");
        profile->GetYaxis()->SetTitle("Mean dE_{dep}/dr per incident neutron (MeV/mm)");
        profile->GetYaxis()->SetTitleOffset(1.55);
        profile->SetMinimum(0.0);
        profile->Draw("hist");
        SaveCanvas(c, outDir, "bgo_20_radial_shower_profile");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 21: Shower depth centroid, depth RMS, and radial RMS vs energy
    // -------------------------------------------------------------------------
    if (haveExtended) {
        std::unique_ptr<TH1D> hDepthCent(MeanHistogram(
            "hDepthCent", statDepthCentroid,
            "Primary neutron energy (GeV)", "Mean shower spatial scale (mm)"));
        std::unique_ptr<TH1D> hDepthSpread(MeanHistogram(
            "hDepthSpread", statDepthRMS,
            "Primary neutron energy (GeV)", "Mean shower spatial scale (mm)"));
        std::unique_ptr<TH1D> hRadialSpread(MeanHistogram(
            "hRadialSpread", statRadialRMS,
            "Primary neutron energy (GeV)", "Mean shower spatial scale (mm)"));
        auto* c = new TCanvas("bgo_c21", "", 880, 680);
        c->SetLogx(); c->SetGrid();
        StyleLine(hDepthCent.get(), kBlue + 1, 1, 2, 20);
        StyleLine(hDepthSpread.get(), kRed + 1, 1, 2, 21);
        StyleLine(hRadialSpread.get(), kGreen + 2, 1, 2, 22);
        hDepthCent->SetMinimum(0.0);
        hDepthCent->Draw("E1");
        hDepthSpread->Draw("E1 same");
        hRadialSpread->Draw("E1 same");
        auto* legend = new TLegend(0.56, 0.73, 0.93, 0.90);
        legend->AddEntry(hDepthCent.get(), "Depth centroid", "lp");
        legend->AddEntry(hDepthSpread.get(), "Longitudinal RMS", "lp");
        legend->AddEntry(hRadialSpread.get(), "Radial RMS", "lp");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_21_shower_moments_vs_energy");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 22: Front/middle/rear thirds of the calorimeter
    // -------------------------------------------------------------------------
    if (haveExtended) {
        std::unique_ptr<TH1D> hFront(MeanHistogram(
            "hFront", statFrontEdep,
            "Primary neutron energy (GeV)", "Mean BGO E_{dep} per incident neutron (MeV)"));
        std::unique_ptr<TH1D> hMiddle(MeanHistogram(
            "hMiddle", statMiddleEdep,
            "Primary neutron energy (GeV)", "Mean BGO E_{dep} per incident neutron (MeV)"));
        std::unique_ptr<TH1D> hRear(MeanHistogram(
            "hRear", statRearEdep,
            "Primary neutron energy (GeV)", "Mean BGO E_{dep} per incident neutron (MeV)"));
        auto* c = new TCanvas("bgo_c22", "", 880, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();
        StyleLine(hFront.get(), kBlue + 1, 1, 2);
        StyleLine(hMiddle.get(), kGreen + 2, 1, 2);
        StyleLine(hRear.get(), kRed + 1, 1, 2);
        hFront->SetMinimum(1.0e-4);
        hFront->Draw("hist");
        hMiddle->Draw("hist same");
        hRear->Draw("hist same");
        auto* legend = new TLegend(0.60, 0.74, 0.93, 0.90);
        legend->AddEntry(hFront.get(), "Front third", "l");
        legend->AddEntry(hMiddle.get(), "Middle third", "l");
        legend->AddEntry(hRear.get(), "Rear third", "l");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_22_longitudinal_partition_vs_energy");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 23: Prompt and delayed energy deposition
    // -------------------------------------------------------------------------
    if (haveExtended) {
        std::unique_ptr<TH1D> hPromptMean(MeanHistogram(
            "hPromptMean", statPrompt,
            "Primary neutron energy (GeV)", "Mean BGO E_{dep} per incident neutron (MeV)"));
        std::unique_ptr<TH1D> hDelayedMean(MeanHistogram(
            "hDelayedMean", statDelayed,
            "Primary neutron energy (GeV)", "Mean BGO E_{dep} per incident neutron (MeV)"));
        auto* c = new TCanvas("bgo_c23", "", 880, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();
        StyleLine(hPromptMean.get(), kBlue + 1, 1, 2);
        StyleLine(hDelayedMean.get(), kRed + 1, 1, 2);
        hPromptMean->SetMinimum(1.0e-6);
        hPromptMean->Draw("hist");
        hDelayedMean->Draw("hist same");
        auto* legend = new TLegend(0.57, 0.76, 0.93, 0.90);
        legend->AddEntry(hPromptMean.get(), "Prompt: #Deltat #leq 10 ns", "l");
        legend->AddEntry(hDelayedMean.get(), "Delayed: #Deltat > 10 ns", "l");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_23_prompt_delayed_edep");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 24: Approximate kinetic-energy leakage and deposited fraction
    // -------------------------------------------------------------------------
    if (haveExtended) {
        std::unique_ptr<TH1D> hLeakFrac(MeanHistogram(
            "hLeakFrac", statLeakageFraction,
            "Primary neutron energy (GeV)", "Mean energy fraction per incident neutron"));
        auto* c = new TCanvas("bgo_c24", "", 880, 680);
        c->SetLogx(); c->SetGrid();
        StyleLine(hMeanResponse.get(), kBlue + 1, 1, 2, 20);
        StyleLine(hLeakFrac.get(), kRed + 1, 1, 2, 21);
        hMeanResponse->GetYaxis()->SetTitle("Mean energy fraction per incident neutron");
        hMeanResponse->SetMinimum(0.0);
        const double ymax = std::max(hMeanResponse->GetMaximum(), hLeakFrac->GetMaximum());
        hMeanResponse->SetMaximum(std::max(0.1, 1.25 * ymax));
        hMeanResponse->Draw("E1");
        hLeakFrac->Draw("E1 same");
        auto* legend = new TLegend(0.56, 0.76, 0.93, 0.90);
        legend->AddEntry(hMeanResponse.get(), "Deposited: E_{dep}/E_{n}", "lp");
        legend->AddEntry(hLeakFrac.get(), "Approximate escaping kinetic energy/E_{n}", "lp");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_24_deposition_and_leakage_fraction");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 25: Primary-neutron escape probability
    // -------------------------------------------------------------------------
    if (haveEscape) {
        std::unique_ptr<TH1D> hEscapeProbability(FractionHistogram(
            "hEscapeProbability", probPrimaryEscape,
            "Primary neutron energy (GeV)", "Primary-neutron escape probability"));
        auto* c = new TCanvas("bgo_c25", "", 880, 680);
        c->SetLogx(); c->SetGrid();
        StyleLine(hEscapeProbability.get(), kBlue + 1, 1, 2, 20);
        hEscapeProbability->GetYaxis()->SetRangeUser(0.0, 1.02);
        hEscapeProbability->GetYaxis()->SetTitleOffset(1.45);
        hEscapeProbability->Draw("E1");
        SaveCanvas(c, outDir, "bgo_25_primary_neutron_escape_probability");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOTS 26-27: Mean crystal energy maps for each layer
    // -------------------------------------------------------------------------
    if (h2CrystalEdep) {
        for (int layer = 0; layer < 2; ++layer) {
            auto* map = new TH2D(
                TString::Format("hCrystalMap_layer%d", layer), "",
                4, -0.5, 3.5, 4, -0.5, 3.5);
            map->SetDirectory(nullptr);
            for (int iy = 0; iy < 4; ++iy) {
                for (int ix = 0; ix < 4; ++ix) {
                    const int crystalID = layer * 16 + iy * 4 + ix;
                    const int yBin = h2CrystalEdep->GetYaxis()->FindBin(crystalID);
                    const double energySum = h2CrystalEdep->Integral(
                        1, h2CrystalEdep->GetNbinsX(), yBin, yBin);
                    map->SetBinContent(ix + 1, iy + 1, energySum / sumW);
                }
            }
            auto* c = new TCanvas(TString::Format("bgo_c%d", 26 + layer), "", 820, 700);
            c->SetRightMargin(0.17);
            c->SetGrid();
            map->GetXaxis()->SetTitle("Crystal x index");
            map->GetYaxis()->SetTitle("Crystal y index");
            map->GetZaxis()->SetTitle("Mean E_{dep} per incident neutron (MeV)");
            map->GetYaxis()->SetTitleOffset(1.35);
            map->GetZaxis()->SetTitleOffset(1.30);
            map->Draw("colz text");
            SaveCanvas(c, outDir,
                TString::Format("bgo_%02d_crystal_energy_map_layer%d", 26 + layer, layer));
            delete c;
            delete map;
        }
    }

    // -------------------------------------------------------------------------
    // PLOT 28: Mean number of hit crystals vs primary energy
    // -------------------------------------------------------------------------
    if (haveExtended) {
        std::unique_ptr<TH1D> hHitCrystalMean(MeanHistogram(
            "hHitCrystalMean", statHitCrystals,
            "Primary neutron energy (GeV)",
            "Mean number of crystals with E_{dep} > 0 per incident neutron"));
        auto* c = new TCanvas("bgo_c28", "", 880, 680);
        c->SetLogx(); c->SetGrid();
        StyleLine(hHitCrystalMean.get(), kBlue + 1, 1, 2, 20);
        hHitCrystalMean->SetMinimum(0.0);
        hHitCrystalMean->Draw("E1");
        SaveCanvas(c, outDir, "bgo_28_hit_crystal_multiplicity_vs_energy");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 29: Mean Edep in the two BGO layers
    // -------------------------------------------------------------------------
    if (haveExtended) {
        std::unique_ptr<TH1D> hLayer0Mean(MeanHistogram(
            "hLayer0Mean", statLayer0,
            "Primary neutron energy (GeV)",
            "Mean BGO E_{dep} per incident neutron (MeV)"));
        std::unique_ptr<TH1D> hLayer1Mean(MeanHistogram(
            "hLayer1Mean", statLayer1,
            "Primary neutron energy (GeV)",
            "Mean BGO E_{dep} per incident neutron (MeV)"));
        auto* c = new TCanvas("bgo_c29", "", 880, 680);
        c->SetLogx(); c->SetLogy(); c->SetGrid();
        StyleLine(hLayer0Mean.get(), kBlue + 1, 1, 2);
        StyleLine(hLayer1Mean.get(), kRed + 1, 1, 2);
        hLayer0Mean->SetMinimum(1.0e-4);
        hLayer0Mean->Draw("hist");
        hLayer1Mean->Draw("hist same");
        auto* legend = new TLegend(0.61, 0.77, 0.93, 0.90);
        legend->AddEntry(hLayer0Mean.get(), "Front BGO layer", "l");
        legend->AddEntry(hLayer1Mean.get(), "Rear BGO layer", "l");
        legend->Draw();
        SaveCanvas(c, outDir, "bgo_29_layer_edep_vs_energy");
        delete c;
    }

    // -------------------------------------------------------------------------
    // PLOT 30: Edep vs produced-secondary multiplicity
    // -------------------------------------------------------------------------
    if (haveExtended) {
        auto* c = new TCanvas("bgo_c30", "", 860, 680);
        c->SetRightMargin(0.16);
        c->SetLogx(); c->SetLogz(); c->SetGrid();
        h2EdepVsMultiplicity->GetXaxis()->SetTitle("Total BGO E_{dep} (MeV)");
        h2EdepVsMultiplicity->GetYaxis()->SetTitle("Produced-secondary multiplicity");
        h2EdepVsMultiplicity->GetZaxis()->SetTitle("Weighted events");
        h2EdepVsMultiplicity->GetYaxis()->SetTitleOffset(1.45);
        h2EdepVsMultiplicity->GetZaxis()->SetTitleOffset(1.25);
        FixLogAxis(h2EdepVsMultiplicity->GetXaxis());
        h2EdepVsMultiplicity->Draw("colz");
        SaveCanvas(c, outDir, "bgo_30_edep_vs_secondary_multiplicity");
        delete c;
    }

    // -------------------------------------------------------------------------
    // Energy-binned CSV summary
    // -------------------------------------------------------------------------
    {
        std::ofstream csv(TString::Format("%s/bgo_energy_binned_summary.csv", outDir).Data());
        csv << std::setprecision(10);
        csv << "E_low_GeV,E_high_GeV,E_center_GeV,sum_weight,mean_edep_MeV,"
               "rms_edep_MeV,mean_response,response_rel_rms,"
               "prob_any_hadronic,prob_primary_hadronic,prob_primary_inelastic,"
               "eff_edep_ge_1MeV,eff_edep_ge_5MeV,eff_edep_ge_10MeV,"
               "eff_edep_ge_50MeV,eff_edep_ge_100MeV,"
               "mean_secondary_multiplicity,mean_unique_tracks,"
               "mean_first_depth_mm,mean_depth_centroid_mm,mean_radial_rms_mm,"
               "mean_leakage_MeV,prob_primary_escape\n";

        for (std::size_t i = 0; i < statEdepAll.sumW.size(); ++i) {
            const double eLow = energyEdges[i];
            const double eHigh = energyEdges[i + 1U];
            const double eCenter = std::sqrt(eLow * eHigh);
            const double meanResponse = statResponse.Mean(i);
            const double relWidth = meanResponse != 0.0
                ? statResponse.RMS(i) / std::abs(meanResponse) : 0.0;
            csv << eLow << ',' << eHigh << ',' << eCenter << ','
                << statEdepAll.sumW[i] << ','
                << statEdepAll.Mean(i) << ',' << statEdepAll.RMS(i) << ','
                << meanResponse << ',' << relWidth << ','
                << probAnyHadronic.Value(i) << ','
                << probPrimaryHadronic.Value(i) << ','
                << probPrimaryInelastic.Value(i) << ',';
            for (std::size_t it = 0; it < detectionEfficiencies.size(); ++it)
                csv << detectionEfficiencies[it].Value(i) << ',';
            csv << statMultTotal.Mean(i) << ','
                << statUniqueTracks.Mean(i) << ','
                << statFirstDepth.Mean(i) << ','
                << statDepthCentroid.Mean(i) << ','
                << statRadialRMS.Mean(i) << ','
                << statLeakage.Mean(i) << ','
                << probPrimaryEscape.Value(i) << '\n';
        }
    }

    // -------------------------------------------------------------------------
    // Human-readable summary
    // -------------------------------------------------------------------------
    const double rawN = static_cast<double>(nEntries);
    const double meanTotal = sumEdepTotal / sumW;
    const double probabilityAny = weightedAnyHadronic / sumW;
    const double probabilityPrimary = weightedPrimaryHadronic / sumW;
    const double probabilityInelastic = weightedPrimaryInelastic / sumW;
    const double probabilityEscape = weightedPrimaryEscape / sumW;

    std::ofstream summary(TString::Format("%s/bgo_summary.txt", outDir).Data());
    summary << std::fixed << std::setprecision(6);
    summary << "========== EXTENDED BGO SUMMARY ==========\n";
    summary << "Input file: " << fname << "\n";
    summary << "Raw event entries: " << nEntries << "\n";
    summary << "Sum of event weights: " << sumW << "\n";
    summary << "Effective weighted entries: "
            << (sumW2 > 0.0 ? sumW * sumW / sumW2 : 0.0) << "\n\n";

    summary << "Geometry and analytical interaction estimate:\n";
    summary << "  Active BGO material depth: " << kBGOActiveDepthCm << " cm\n";
    summary << "  BGO nuclear interaction length: " << kBGOInteractionLengthCm << " cm\n";
    summary << "  Expected 1-exp(-L/lambda): "
            << 100.0 * expectedInteractionProbability << " %\n";
    if (fittedLambdaMm > 0.0) {
        summary << "  Fitted interaction length from depth distribution: "
                << fittedLambdaMm / 10.0 << " +/- "
                << fittedLambdaErrorMm / 10.0 << " cm\n";
    }

    summary << "\nInteraction probabilities:\n";
    summary << "  Any selected hadronic process: " << 100.0 * probabilityAny
            << " % (raw " << rawAnyHadronic << "/" << nEntries << ")\n";
    summary << "  Primary-neutron hadronic: " << 100.0 * probabilityPrimary
            << " % (raw " << rawPrimaryHadronic << "/" << nEntries << ")\n";
    summary << "  Primary-neutron inelastic: " << 100.0 * probabilityInelastic
            << " % (raw " << rawPrimaryInelastic << "/" << nEntries << ")\n";

    summary << "\nMean BGO Edep per incident neutron:\n";
    summary << "  Total: " << meanTotal << " MeV\n";
    summary << "  n+p: " << sumEdepNP / sumW << " MeV\n";
    summary << "  e/gamma: " << sumEdepEM / sumW << " MeV\n";
    summary << "  heavy ions: " << sumEdepHeavy / sumW << " MeV\n";
    summary << "  light ions: " << sumEdepLight / sumW << " MeV\n";
    summary << "  pi/K: " << sumEdepPiK / sumW << " MeV\n";
    summary << "  other: " << sumEdepOther / sumW << " MeV\n";

    summary << "\nIntegrated detection efficiencies:\n";
    for (std::size_t it = 0; it < thresholdsMeV.size(); ++it) {
        double numerator = 0.0;
        double denominator = 0.0;
        for (std::size_t ib = 0; ib < detectionEfficiencies[it].denW.size(); ++ib) {
            numerator += detectionEfficiencies[it].numW[ib];
            denominator += detectionEfficiencies[it].denW[ib];
        }
        summary << "  Edep >= " << thresholdsMeV[it] << " MeV: "
                << 100.0 * SafeDivide(numerator, denominator) << " %\n";
    }

    if (haveExtended) {
        summary << "\nExtended event observables:\n";
        summary << "  Mean produced-secondary multiplicity: "
                << sumSecondaryCount / sumW << "\n";
        summary << "  Mean unique tracks with BGO hits: "
                << sumUniqueTrackCount / sumW << "\n";
        summary << "  Mean prompt Edep (<=10 ns): " << sumPrompt / sumW << " MeV\n";
        summary << "  Mean delayed Edep (>10 ns): " << sumDelayed / sumW << " MeV\n";
        summary << "  Mean approximate escaping kinetic energy: "
                << sumLeakage / sumW << " MeV\n";
        summary << "  Primary-neutron escape probability: "
                << 100.0 * probabilityEscape << " % (raw "
                << rawPrimaryEscape << "/" << nEntries << ")\n";
        summary << "  First-interaction secondary rows above 0.1 MeV: "
                << nSecondaryRows << "\n";
    }
    summary << "==========================================\n";
    summary.close();

    std::cout << "\n========== EXTENDED BGO SUMMARY ==========\n";
    std::cout << "  Total events:                     " << nEntries << "\n";
    std::cout << "  Sum of weights:                   " << sumW << "\n";
    std::cout << "  Any hadronic process:             "
              << TString::Format("%.3f%%", 100.0 * probabilityAny) << "\n";
    std::cout << "  Primary-neutron hadronic:         "
              << TString::Format("%.3f%%", 100.0 * probabilityPrimary) << "\n";
    std::cout << "  Primary-neutron inelastic:        "
              << TString::Format("%.3f%%", 100.0 * probabilityInelastic) << "\n";
    std::cout << "  Analytical expectation (6/22.32): "
              << TString::Format("%.3f%%", 100.0 * expectedInteractionProbability) << "\n";
    if (fittedLambdaMm > 0.0) {
        std::cout << "  Fitted interaction length:        "
                  << TString::Format("%.3f +/- %.3f cm",
                                     fittedLambdaMm / 10.0,
                                     fittedLambdaErrorMm / 10.0) << "\n";
    }
    std::cout << "  Mean total BGO Edep:              "
              << TString::Format("%.6f MeV/incident neutron", meanTotal) << "\n";
    std::cout << "  Mean approximate leakage KE:      "
              << TString::Format("%.6f MeV/incident neutron", sumLeakage / sumW) << "\n";
    std::cout << "  Primary-neutron escape:           "
              << TString::Format("%.3f%%", 100.0 * probabilityEscape) << "\n";
    std::cout << "  Summary text:                     " << outDir << "/bgo_summary.txt\n";
    std::cout << "  Energy-binned CSV:                "
              << outDir << "/bgo_energy_binned_summary.csv\n";
    std::cout << "==========================================\n\n";
}
