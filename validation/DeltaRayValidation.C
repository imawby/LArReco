/**
 *  @file   LArReco/validation/Validation.C
 *
 *  @brief  Implementation of validation functionality
 *
 *  $Log: $
 */
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream> 

#include "DeltaRayValidation.h"

void DeltaRayValidation(const std::string &inputFileName, const Parameters &parameters)
{
    CosmicRayVector cosmicRayVector; DeltaRayVector deltaRayVector;

    ReadTree(inputFileName, cosmicRayVector, deltaRayVector);

    if (parameters.m_printOverallRecoMetrics)
        DisplayOverallRecoMetrics(deltaRayVector);

    if (parameters.m_histogramOutput)
    {
        CosmicRayMCHistogramCollection cosmicRayMCHistogramCollection;
        FillCosmicRayMCHistogramCollection(cosmicRayVector, cosmicRayMCHistogramCollection);

        DeltaRayMCHistogramCollection deltaRayMCHistogramCollection;
        FillDeltaRayMCHistogramCollection(deltaRayVector, deltaRayMCHistogramCollection);

        DeltaRayRecoHistogramCollection deltaRayRecoHistogramCollection;
        FillDeltaRayRecoHistogramCollection(deltaRayVector, deltaRayRecoHistogramCollection);

        ProcessHistograms(deltaRayMCHistogramCollection, deltaRayRecoHistogramCollection);

        WriteHistograms(cosmicRayMCHistogramCollection, deltaRayMCHistogramCollection, deltaRayRecoHistogramCollection);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ReadTree(const std::string &inputFileName, CosmicRayVector &cosmicRayVector, DeltaRayVector &deltaRayVector)
{
    TFile * validationFile = new TFile(inputFileName.c_str(), "READ"); 
    TTree * validationTree = (TTree*)validationFile->Get("Validation");
    
    CosmicRay cosmicRay; 

    validationTree->SetBranchAddress("mcE_CR", &cosmicRay.m_energy);
    validationTree->SetBranchAddress("mcPX_CR", &cosmicRay.m_momentum.m_x);
    validationTree->SetBranchAddress("mcPY_CR", &cosmicRay.m_momentum.m_y);
    validationTree->SetBranchAddress("mcPZ_CR", &cosmicRay.m_momentum.m_z);
    validationTree->SetBranchAddress("nMCHitsTotal_CR", &cosmicRay.m_nMCHitsTotal);
    validationTree->SetBranchAddress("nMCHitsU_CR", &cosmicRay.m_nMCHitsU);  
    validationTree->SetBranchAddress("nMCHitsV_CR", &cosmicRay.m_nMCHitsV);
    validationTree->SetBranchAddress("nMCHitsW_CR", &cosmicRay.m_nMCHitsW);
    validationTree->SetBranchAddress("nReconstructableChildDRs", &cosmicRay.m_nReconstructableChildDRs);
    validationTree->SetBranchAddress("nCorrectChildDRs", &cosmicRay.m_nCorrectChildDRs);

    FloatVector *mcE_DR(nullptr), *mcPX_DR(nullptr), *mcPY_DR(nullptr), *mcPZ_DR(nullptr);
    IntVector *nMCHitsTotal_DR(nullptr), *nMCHitsU_DR(nullptr), *nMCHitsV_DR(nullptr), *nMCHitsW_DR(nullptr);
    FloatVector *mcVertexX_DR(nullptr), *mcVertexY_DR(nullptr), *mcVertexZ_DR(nullptr), *mcEndX_DR(nullptr), *mcEndY_DR(nullptr), *mcEndZ_DR(nullptr);
    IntVector *nAboveThresholdMatches_DR(nullptr), *isCorrect_DR(nullptr), *isCorrectParentLink_DR(nullptr);
    IntVector *bestMatchNHitsTotal_DR(nullptr), *bestMatchNHitsU_DR(nullptr), *bestMatchNHitsV_DR(nullptr), *bestMatchNHitsW_DR(nullptr);
    IntVector *bestMatchNSharedHitsTotal_DR(nullptr), *bestMatchNSharedHitsU_DR(nullptr), *bestMatchNSharedHitsV_DR(nullptr), *bestMatchNSharedHitsW_DR(nullptr);
    IntVector *bestMatchNParentTrackHitsTotal_DR(nullptr), *bestMatchNParentTrackHitsU_DR(nullptr), *bestMatchNParentTrackHitsV_DR(nullptr), *bestMatchNParentTrackHitsW_DR(nullptr);
    IntVector *bestMatchNOtherTrackHitsTotal_DR(nullptr), *bestMatchNOtherTrackHitsU_DR(nullptr), *bestMatchNOtherTrackHitsV_DR(nullptr), *bestMatchNOtherTrackHitsW_DR(nullptr);
    IntVector *bestMatchNOtherShowerHitsTotal_DR(nullptr), *bestMatchNOtherShowerHitsU_DR(nullptr), *bestMatchNOtherShowerHitsV_DR(nullptr), *bestMatchNOtherShowerHitsW_DR(nullptr);
    IntVector *totalDRHitsInBestMatchParentCR_DR(nullptr), *uDRHitsInBestMatchParentCR_DR(nullptr), *vDRHitsInBestMatchParentCR_DR(nullptr), *wDRHitsInBestMatchParentCR_DR(nullptr);

    validationTree->SetBranchAddress("mcE_DR", &mcE_DR);
    validationTree->SetBranchAddress("mcPX_DR", &mcPX_DR);
    validationTree->SetBranchAddress("mcPY_DR", &mcPY_DR);
    validationTree->SetBranchAddress("mcPZ_DR", &mcPZ_DR);
    validationTree->SetBranchAddress("nMCHitsTotal_DR", &nMCHitsTotal_DR);
    validationTree->SetBranchAddress("nMCHitsU_DR", &nMCHitsU_DR);  
    validationTree->SetBranchAddress("nMCHitsV_DR", &nMCHitsV_DR);
    validationTree->SetBranchAddress("nMCHitsW_DR", &nMCHitsW_DR);
    validationTree->SetBranchAddress("nAboveThresholdMatches_DR", &nAboveThresholdMatches_DR);
    validationTree->SetBranchAddress("isCorrect_DR", &isCorrect_DR);
    validationTree->SetBranchAddress("isCorrectParentLink_DR", &isCorrectParentLink_DR);
    validationTree->SetBranchAddress("bestMatchNHitsTotal_DR", &bestMatchNHitsTotal_DR);    
    validationTree->SetBranchAddress("bestMatchNHitsU_DR", &bestMatchNHitsU_DR);  
    validationTree->SetBranchAddress("bestMatchNHitsV_DR", &bestMatchNHitsV_DR);
    validationTree->SetBranchAddress("bestMatchNHitsW_DR", &bestMatchNHitsW_DR);
    validationTree->SetBranchAddress("bestMatchNSharedHitsTotal_DR", &bestMatchNSharedHitsTotal_DR);    
    validationTree->SetBranchAddress("bestMatchNSharedHitsU_DR", &bestMatchNSharedHitsU_DR);        
    validationTree->SetBranchAddress("bestMatchNSharedHitsV_DR", &bestMatchNSharedHitsV_DR);
    validationTree->SetBranchAddress("bestMatchNSharedHitsW_DR", &bestMatchNSharedHitsW_DR);
    validationTree->SetBranchAddress("bestMatchNParentTrackHitsTotal_DR", &bestMatchNParentTrackHitsTotal_DR);    
    validationTree->SetBranchAddress("bestMatchNParentTrackHitsU_DR", &bestMatchNParentTrackHitsU_DR);
    validationTree->SetBranchAddress("bestMatchNParentTrackHitsV_DR", &bestMatchNParentTrackHitsV_DR);
    validationTree->SetBranchAddress("bestMatchNParentTrackHitsW_DR", &bestMatchNParentTrackHitsW_DR);
    validationTree->SetBranchAddress("bestMatchNOtherTrackHitsTotal_DR", &bestMatchNOtherTrackHitsTotal_DR);
    validationTree->SetBranchAddress("bestMatchNOtherTrackHitsU_DR", &bestMatchNOtherTrackHitsU_DR);
    validationTree->SetBranchAddress("bestMatchNOtherTrackHitsV_DR", &bestMatchNOtherTrackHitsV_DR);
    validationTree->SetBranchAddress("bestMatchNOtherTrackHitsW_DR", &bestMatchNOtherTrackHitsW_DR);
    validationTree->SetBranchAddress("bestMatchNOtherShowerHitsTotal_DR", &bestMatchNOtherShowerHitsTotal_DR);    
    validationTree->SetBranchAddress("bestMatchNOtherShowerHitsU_DR", &bestMatchNOtherShowerHitsU_DR);    
    validationTree->SetBranchAddress("bestMatchNOtherShowerHitsV_DR", &bestMatchNOtherShowerHitsV_DR);
    validationTree->SetBranchAddress("bestMatchNOtherShowerHitsW_DR", &bestMatchNOtherShowerHitsW_DR);
    validationTree->SetBranchAddress("totalDRHitsInBestMatchParentCR_DR", &totalDRHitsInBestMatchParentCR_DR);    
    validationTree->SetBranchAddress("uDRHitsInBestMatchParentCR_DR", &uDRHitsInBestMatchParentCR_DR);
    validationTree->SetBranchAddress("vDRHitsInBestMatchParentCR_DR", &vDRHitsInBestMatchParentCR_DR);
    validationTree->SetBranchAddress("wDRHitsInBestMatchParentCR_DR", &wDRHitsInBestMatchParentCR_DR);

    DeltaRay deltaRay;
    for (Int_t i(0); i < (Int_t)validationTree->GetEntries(); ++i)
    {
        validationTree->GetEntry(i);

        //std::cout << "COSMIC RAY" << std::endl;
        //cosmicRay.Print();

        for (int i = 0; i < cosmicRay.m_nReconstructableChildDRs; ++i)
        {
            deltaRay.m_energy = mcE_DR->at(i);
            deltaRay.m_momentum = SimpleThreeVector(mcPX_DR->at(i), mcPY_DR->at(i), mcPZ_DR->at(i));
            deltaRay.m_nMCHitsTotal = nMCHitsTotal_DR->at(i);
            deltaRay.m_nMCHitsU = nMCHitsU_DR->at(i);
            deltaRay.m_nMCHitsV = nMCHitsV_DR->at(i);
            deltaRay.m_nMCHitsW = nMCHitsW_DR->at(i);

            if ((cosmicRay.m_momentum.GetMagnitude() < std::numeric_limits<float>::epsilon()) || (deltaRay.m_momentum.GetMagnitude() < std::numeric_limits<float>::epsilon()))
            {
                deltaRay.m_openingAngleFromMuon = -1.f;
            }
            else
            {
                deltaRay.m_openingAngleFromMuon = cosmicRay.m_momentum.GetOpeningAngle(deltaRay.m_momentum);
            }
    
            deltaRay.m_nAboveThresholdMatches = nAboveThresholdMatches_DR->at(i);
            deltaRay.m_isCorrect = isCorrect_DR->at(i);
            deltaRay.m_isCorrectParentLink = isCorrectParentLink_DR->at(i);

            deltaRay.m_bestMatchNHitsTotal = bestMatchNHitsTotal_DR->at(i);
            deltaRay.m_bestMatchNHitsU = bestMatchNHitsU_DR->at(i);
            deltaRay.m_bestMatchNHitsV = bestMatchNHitsV_DR->at(i);
            deltaRay.m_bestMatchNHitsW = bestMatchNHitsW_DR->at(i);

            deltaRay.m_bestMatchNSharedHitsTotal = bestMatchNSharedHitsTotal_DR->at(i);
            deltaRay.m_bestMatchNSharedHitsU = bestMatchNSharedHitsU_DR->at(i);
            deltaRay.m_bestMatchNSharedHitsV = bestMatchNSharedHitsV_DR->at(i);
            deltaRay.m_bestMatchNSharedHitsW = bestMatchNSharedHitsW_DR->at(i);
    
            deltaRay.m_bestMatchNParentTrackHitsTotal = bestMatchNParentTrackHitsTotal_DR->at(i);
            deltaRay.m_bestMatchNParentTrackHitsU = bestMatchNParentTrackHitsU_DR->at(i);
            deltaRay.m_bestMatchNParentTrackHitsV = bestMatchNParentTrackHitsV_DR->at(i);
            deltaRay.m_bestMatchNParentTrackHitsW = bestMatchNParentTrackHitsW_DR->at(i);

            deltaRay.m_bestMatchNOtherTrackHitsTotal = bestMatchNOtherTrackHitsTotal_DR->at(i);
            deltaRay.m_bestMatchNOtherTrackHitsU = bestMatchNOtherTrackHitsU_DR->at(i);
            deltaRay.m_bestMatchNOtherTrackHitsV = bestMatchNOtherTrackHitsV_DR->at(i);
            deltaRay.m_bestMatchNOtherTrackHitsW = bestMatchNOtherTrackHitsW_DR->at(i);

            deltaRay.m_bestMatchNOtherShowerHitsTotal = bestMatchNOtherShowerHitsTotal_DR->at(i);
            deltaRay.m_bestMatchNOtherShowerHitsU = bestMatchNOtherShowerHitsU_DR->at(i);
            deltaRay.m_bestMatchNOtherShowerHitsV = bestMatchNOtherShowerHitsV_DR->at(i);
            deltaRay.m_bestMatchNOtherShowerHitsW = bestMatchNOtherShowerHitsW_DR->at(i);

            deltaRay.m_totalDRHitsInBestMatchParentCR = totalDRHitsInBestMatchParentCR_DR->at(i);
            deltaRay.m_uDRHitsInBestMatchParentCR = uDRHitsInBestMatchParentCR_DR->at(i);
            deltaRay.m_vDRHitsInBestMatchParentCR = vDRHitsInBestMatchParentCR_DR->at(i);
            deltaRay.m_wDRHitsInBestMatchParentCR = wDRHitsInBestMatchParentCR_DR->at(i);

            //std::cout << "DELTA RAY" << std::endl;
            //deltaRay.Print();
            deltaRayVector.push_back(deltaRay);
        }

        cosmicRayVector.push_back(cosmicRay);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DisplayOverallRecoMetrics(const DeltaRayVector &deltaRayVector)
{
    unsigned int nCorrectParentLinks(0), nCorrectlyReconstructedDRs(0);
    unsigned int nZeroMatches(0), nOneMatches(0), nTwoMatches(0), nThreePlusMatches(0);
    
    for (const DeltaRay &deltaRay : deltaRayVector)
    {
        if (deltaRay.m_isCorrectParentLink)
            ++nCorrectParentLinks;

        if (deltaRay.m_isCorrect)
            ++nCorrectlyReconstructedDRs;
        
        switch (deltaRay.m_nAboveThresholdMatches)
        {
        case 0:
            ++nZeroMatches;
            break;
        case 1:
            ++nOneMatches;
            break;
        case 2:
            ++nTwoMatches;
            break;
        default:
            ++nThreePlusMatches;
        }
    }

    const unsigned int nReconstructableDRs(deltaRayVector.size());
    const float fCorrectParentLinks(nReconstructableDRs == 0 ? 0.f : static_cast<float>(nCorrectParentLinks) / static_cast<float>(nReconstructableDRs));
    const float fCorrectlyReconstructedDRs(nReconstructableDRs == 0 ? 0.f : static_cast<float>(nCorrectlyReconstructedDRs) / static_cast<float>(nReconstructableDRs));
    const float fZeroMatches(nReconstructableDRs == 0 ? 0.f : static_cast<float>(nZeroMatches) / static_cast<float>(nReconstructableDRs));
    const float fOneMatches(nReconstructableDRs == 0 ? 0.f : static_cast<float>(nOneMatches) / static_cast<float>(nReconstructableDRs));
    const float fTwoMatches(nReconstructableDRs == 0 ? 0.f : static_cast<float>(nTwoMatches) / static_cast<float>(nReconstructableDRs));
    const float fThreePlusMatches(nReconstructableDRs == 0 ? 0.f : static_cast<float>(nThreePlusMatches) / static_cast<float>(nReconstructableDRs));

    std::cout << "nReconstructableDeltaRays: " << nReconstructableDRs << std::endl
              << "nCorrectlyReconstructedDeltaRays: " << nCorrectlyReconstructedDRs << " (" << fCorrectlyReconstructedDRs * 100.f << "%)" << std::endl
              << "nCorrectParentLinks: " << nCorrectParentLinks << " (" << fCorrectParentLinks * 100.f << "%)" << std::endl
              << "Above threshold pfo match distribution: "
              << "|0: " << fZeroMatches * 100.f << "%|, |1: " << fOneMatches * 100.f << "%|, |2: " << fTwoMatches * 100.f << "%|, |3+: " << fThreePlusMatches * 100.f <<  "%|" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FillCosmicRayMCHistogramCollection(const CosmicRayVector &cosmicRayVector, CosmicRayMCHistogramCollection &cosmicRayMCHistogramCollection)
{
   if (!cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays)
   {
        cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays = new TH1F("hReconstructableChildDeltaRays_CR", "hReconstructableChildDeltaRays_CR", 40000, -100., 1900.);
        cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays->SetTitle(";nReconstructableChildDeltaRays;Occurance");
        cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays->GetXaxis()->SetRangeUser(0.f, 50.f);
    }
    
    for (const CosmicRay &cosmicRay : cosmicRayVector)
    {
        cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays->Fill(cosmicRay.m_nReconstructableChildDRs);
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FillDeltaRayMCHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection)
{
   if (!deltaRayMCHistogramCollection.m_hEnergyDistribution)
   {
        deltaRayMCHistogramCollection.m_hEnergyDistribution = new TH1F("hEnergyDistribution_DR", "hEnergyDistribution_DR", 40000, -100., 1900.);
        deltaRayMCHistogramCollection.m_hEnergyDistribution->SetTitle(";TrueDREnergy [GeV];Occurance");
        //deltaRayMCHistogramCollection.m_hEnergyDistribution->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayMCHistogramCollection.m_hTotalHitDistribution)
   {
        deltaRayMCHistogramCollection.m_hTotalHitDistribution = new TH1F("hTotalHitDistribution_DR", "hTotalHitDistribution_DR", 40000, -100., 1900.);
        deltaRayMCHistogramCollection.m_hTotalHitDistribution->SetTitle(";nMCTotalHits;Occurance");
        //deltaRayMCHistogramCollection.m_hTotalHitDistribution->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayMCHistogramCollection.m_hOpeningAngleDistribution)
   {
        deltaRayMCHistogramCollection.m_hOpeningAngleDistribution = new TH1F("hOpeningAngleDistribution_DR", "hOpeningAngleDistribution_DR", 40000, -100., 1900.);
        deltaRayMCHistogramCollection.m_hOpeningAngleDistribution->SetTitle(";TrueOpeningAngle [degrees];Occurance");
        //deltaRayMCHistogramCollection.m_hOpeningAngleDistribution->GetXaxis()->SetRangeUser(0.f, 50.f);
    }      
    
    for (const DeltaRay &deltaRay : deltaRayVector)
    {
        deltaRayMCHistogramCollection.m_hEnergyDistribution->Fill(deltaRay.m_energy);
        deltaRayMCHistogramCollection.m_hTotalHitDistribution->Fill(deltaRay.m_nMCHitsTotal);
        deltaRayMCHistogramCollection.m_hOpeningAngleDistribution->Fill(deltaRay.m_openingAngleFromMuon);        
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FillDeltaRayRecoHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection)
{
   if (!deltaRayRecoHistogramCollection.m_hCompleteness)
   {
        deltaRayRecoHistogramCollection.m_hCompleteness = new TH1F("hCompleteness_DR", "hCompleteness_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hCompleteness->SetTitle(";Completeness;Occurance");
        //deltaRayRecoHistogramCollection.m_hCompleteness->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hPurity)
   {
        deltaRayRecoHistogramCollection.m_hPurity = new TH1F("hPurity_DR", "hPurity_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hPurity->SetTitle(";Purity;Occurance");
        //deltaRayRecoHistogramCollection.m_hPurity->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hAboveThresholdMatches)
   {
        deltaRayRecoHistogramCollection.m_hAboveThresholdMatches = new TH1F("hAboveThresholdMatches_DR", "hAboveThresholdMatches_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->SetTitle(";nAboveThresholdMatches;Occurance");
        //deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal)
   {
        deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal = new TH1F("hParentTrackHitsTotal_DR", "hParentTrackHitsTotal_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->SetTitle(";nParentTrackHitsTotal;Occurance");
        //deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal)
   {
        deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal = new TH1F("hOtherTrackHitsTotal_DR", "hOtherTrackHitsTotal_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->SetTitle(";nOtherTrackHitsTotal;Occurance");
        //deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal)
   {
        deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal = new TH1F("hOtherShowerHitsTotal_DR", "hOtherShowerHitsTotal_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->SetTitle(";nOtherShowerHitsTotal;Occurance");
        //deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay)
   {
        deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay = new TH1F("hTotalHitsTakenByCosmicRay_DR", "hTotalHitsTakenByCosmicRay_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->SetTitle(";hTotalHitsTakenByCosmicRay;Occurance");
        //deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hEfficiency_Energy)
   {
        deltaRayRecoHistogramCollection.m_hEfficiency_Energy = new TH1F("hEfficiency_Energy_DR", "hEfficiency_Energy_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hEfficiency_Energy->SetTitle(";TrueDREnergy [GeV];Efficiency");
        //deltaRayRecoHistogramCollection.m_hEfficiency_Energy->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits)
   {
        deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits = new TH1F("hEfficiency_TotalHits_DR", "hEfficiency_TotalHits_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->SetTitle(";nMCTotalHits;Efficiency");
        //deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle)
   {
        deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle = new TH1F("hEfficiency_OpeningAngle_DR", "_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->SetTitle(";TrueOpeningAngle [degrees];Efficiency");
        //deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy)
   {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy = new TH1F("hCorrectParentLink_Energy_DR", "hCorrectParentLink_Energy_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy->SetTitle(";TrueDREnergy [GeV];CorrectParentLinkFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits)
   {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits = new TH1F("hCorrectParentLink_TotalHits_DR", "hCorrectParentLink_TotalHits_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits->SetTitle(";nMCTotalHits;CorrectParentLinkFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle)
   {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle = new TH1F("hCorrectParentLink_OpeningAngle_DR", "hCorrectParentLink_OpeningAngle_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle->SetTitle(";TrueOpeningAngle [degrees];CorrectParentLinkFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle->GetXaxis()->SetRangeUser(0.f, 50.f);
    }      

   if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy)
   {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy = new TH1F("hCorrectEvent_Energy_DR", "hCorrectEvent_Energy_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy->SetTitle(";TrueDREnergy [GeV];CorrectlyReconstructedFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits)
   {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits = new TH1F("hCorrectEvent_TotalHits_DR", "hCorrectEvent_TotalHits_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits->SetTitle(";nMCTotalHits;CorrectlyReconstructedFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle)
   {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle = new TH1F("hCorrectEvent_OpeningAngle_DR", "hCorrectEvent_OpeningAngle_DR", 40000, -100., 1900.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle->SetTitle(";TrueOpeningAngle [degrees];CorrectlyReconstructedFraction");
        //deltaRayRecoHistogramCollection.m_-hCorrectEvent_OpeningAngle>GetXaxis()->SetRangeUser(0.f, 50.f);
    }   
   //ISOBEL - CHANGE CORRECT EVENT TO SOMETHING ELSE (NOT ACTUALLY AN EVENT)
   

    for (const DeltaRay &deltaRay : deltaRayVector)
    {
        deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->Fill(deltaRay.m_nAboveThresholdMatches);

        if (deltaRay.m_nAboveThresholdMatches > 0)
        {
            deltaRayRecoHistogramCollection.m_hCompleteness->Fill(static_cast<float>(deltaRay.m_bestMatchNSharedHitsTotal) / static_cast<float>(deltaRay.m_nMCHitsTotal));
            deltaRayRecoHistogramCollection.m_hPurity->Fill(static_cast<float>(deltaRay.m_bestMatchNSharedHitsTotal) / static_cast<float>(deltaRay.m_bestMatchNHitsTotal));
            deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->Fill(deltaRay.m_bestMatchNParentTrackHitsTotal);
            deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->Fill(deltaRay.m_bestMatchNOtherTrackHitsTotal);
            deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->Fill(deltaRay.m_bestMatchNOtherShowerHitsTotal);
            deltaRayRecoHistogramCollection.m_hEfficiency_Energy->Fill(deltaRay.m_energy);
            deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->Fill(deltaRay.m_nMCHitsTotal);            
            deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->Fill(deltaRay.m_openingAngleFromMuon);              
        }

        if (deltaRay.m_isCorrectParentLink)
        {
            deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->Fill(deltaRay.m_totalDRHitsInBestMatchParentCR);
            deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy->Fill(deltaRay.m_energy);
            deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits->Fill(deltaRay.m_nMCHitsTotal);
            deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle->Fill(deltaRay.m_openingAngleFromMuon);
        }

        if (deltaRay.m_isCorrect)
        {
            deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy->Fill(deltaRay.m_energy);
            deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits->Fill(deltaRay.m_nMCHitsTotal);
            deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle->Fill(deltaRay.m_openingAngleFromMuon);             
        }
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessHistograms(DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection, DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection)
{
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_Energy, deltaRayMCHistogramCollection.m_hEnergyDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits, deltaRayMCHistogramCollection.m_hTotalHitDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle, deltaRayMCHistogramCollection.m_hOpeningAngleDistribution);

    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy, deltaRayMCHistogramCollection.m_hEnergyDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits, deltaRayMCHistogramCollection.m_hTotalHitDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle, deltaRayMCHistogramCollection.m_hOpeningAngleDistribution);

    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy, deltaRayMCHistogramCollection.m_hEnergyDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits, deltaRayMCHistogramCollection.m_hTotalHitDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle, deltaRayMCHistogramCollection.m_hOpeningAngleDistribution);

    deltaRayRecoHistogramCollection.m_hCompleteness->Scale(1. / static_cast<double>(deltaRayRecoHistogramCollection.m_hCompleteness->GetEntries()));
    deltaRayRecoHistogramCollection.m_hPurity->Scale(1. / static_cast<double>(deltaRayRecoHistogramCollection.m_hPurity->GetEntries()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DivideHistogram(TH1F *&recoHistogram, TH1F *&mcDistribution)
{
    for (int n = -1; n <= recoHistogram->GetXaxis()->GetNbins(); ++n)
    {
        const float found = recoHistogram->GetBinContent(n + 1);
        const float all = mcDistribution->GetBinContent(n + 1);
        const float efficiency = (all > 0.f) ? found / all : 0.f;
        const float error = (all > found) ? std::sqrt(efficiency * (1. - efficiency) / all) : 0.f;
        recoHistogram->SetBinContent(n + 1, efficiency);
        recoHistogram->SetBinError(n + 1, error);
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void WriteHistograms(CosmicRayMCHistogramCollection &cosmicRayMCHistogramCollection, DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection,
    DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection)
{
    TFile *histogramFile = new TFile("ValidationHistograms.root", "CREATE");

    cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays->Write("hReconstructableChildDeltaRays");

    deltaRayMCHistogramCollection.m_hEnergyDistribution->Write("hEnergyDistribution");
    deltaRayMCHistogramCollection.m_hTotalHitDistribution->Write("hTotalHitDistribution");
    deltaRayMCHistogramCollection.m_hOpeningAngleDistribution->Write("hOpeningAngleDistribution");

    deltaRayRecoHistogramCollection.m_hCompleteness->Write("hCompleteness");
    deltaRayRecoHistogramCollection.m_hPurity->Write("hPurity");
    deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->Write("hAboveThresholdMatches");
    deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->Write("hParentTrackHitsTotal");
    deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->Write("hOtherTrackHitsTotal");
    deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->Write("hOtherShowerHitsTotal");
    deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->Write("hTotalHitsTakenByCosmicRay");
    deltaRayRecoHistogramCollection.m_hEfficiency_Energy->Write("hEfficiency_Energy");
    deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->Write("hEfficiency_TotalHits");
    deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->Write("hEfficiency_OpeningAngle");
    deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy->Write("hCorrectParentLink_Energy");
    deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits->Write("hCorrectParentLink_TotalHits");
    deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle->Write("hCorrectParentLink_OpeningAngle");
    deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy->Write("hCorrectEvent_Energy");
    deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits->Write("hCorrectEvent_TotalHits");
    deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle->Write("hCorrectEvent_OpeningAngle");
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

float SimpleThreeVector::GetMagnitude()
{
    return std::sqrt((this->m_x * this->m_x) + (this->m_y * this->m_y) + (this->m_z * this->m_z));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float SimpleThreeVector::GetOpeningAngle(const SimpleThreeVector &otherVector)
{
    const float dotProduct = (this->m_x * otherVector.m_x) + (this->m_y * otherVector.m_y) + (this->m_z * otherVector.m_z);
    const float thisMagnitude = this->GetMagnitude();
    const float otherMagnitude = otherVector->GetMagnitude();

    if ((thisMagnitude < std::numeric_limits<float>::epsilon()) || (otherMagnitude < std::numeric_limits<float>::epsilon()))
    {
        std::cout << "SimpleThreeVector::GetOpeningAngle - Vector has zero length" << std::endl;
        throw;
    }

    const float cosTheta(dotProduct / (thisMagnitude * otherMagnitude));

    if (std::fabs(cosTheta > 1.f))
    {
        std::cout << "SimpleThreeVector::GetOpeningAngle - CosTheta unphysical" << std::endl;
        throw;
    }

    const float openingAngle(std::acos(cosTheta) * 180.f / M_PI);

    return openingAngle;
}
