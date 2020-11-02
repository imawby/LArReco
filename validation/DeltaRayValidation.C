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

        DeltaRayContaminationHistogramCollection deltaRayContaminationHistogramCollection;
        FillDeltaRayContaminationHistogramCollection(deltaRayVector, deltaRayContaminationHistogramCollection);

        WriteHistograms(cosmicRayMCHistogramCollection, deltaRayMCHistogramCollection, deltaRayRecoHistogramCollection, deltaRayContaminationHistogramCollection);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ReadTree(const std::string &inputFileName, CosmicRayVector &cosmicRayVector, DeltaRayVector &deltaRayVector)
{
    TFile * validationFile = new TFile(inputFileName.c_str(), "READ");
    
    TTree * validationTree = (TTree*)validationFile->Get("Validation");
    TTree * otherShowerContaminationTree = (TTree*)validationFile->Get("OtherShowerContamination");
    TTree * otherTrackContaminationTree = (TTree*)validationFile->Get("OtherTrackContamination");
    TTree * parentTrackContaminationTree = (TTree*)validationFile->Get("ParentTrackContamination");
    TTree * CRLHitsInCRTree = (TTree*)validationFile->Get("CRLHitsInCR");

    int eventNumber;    
    CosmicRay cosmicRay;
    validationTree->SetBranchAddress("eventNumber", &eventNumber);
    validationTree->SetBranchAddress("mcE_CR", &cosmicRay.m_energy);
    validationTree->SetBranchAddress("mcPX_CR", &cosmicRay.m_momentum.m_x);
    validationTree->SetBranchAddress("mcPY_CR", &cosmicRay.m_momentum.m_y);
    validationTree->SetBranchAddress("mcPZ_CR", &cosmicRay.m_momentum.m_z);
    validationTree->SetBranchAddress("nMCHitsTotal_CR", &cosmicRay.m_nMCHitsTotal);
    validationTree->SetBranchAddress("nMCHitsU_CR", &cosmicRay.m_nMCHitsU);  
    validationTree->SetBranchAddress("nMCHitsV_CR", &cosmicRay.m_nMCHitsV);
    validationTree->SetBranchAddress("nMCHitsW_CR", &cosmicRay.m_nMCHitsW);
    validationTree->SetBranchAddress("nReconstructableChildCRLs", &cosmicRay.m_nReconstructableChildCRLs);
    validationTree->SetBranchAddress("nCorrectChildCRLs", &cosmicRay.m_nCorrectChildCRLs);

    FloatVector *mcE_CRL(nullptr), *mcPX_CRL(nullptr), *mcPY_CRL(nullptr), *mcPZ_CRL(nullptr);
    IntVector *nMCHitsTotal_CRL(nullptr), *nMCHitsU_CRL(nullptr), *nMCHitsV_CRL(nullptr), *nMCHitsW_CRL(nullptr);
    IntVector *nAboveThresholdMatches_CRL(nullptr), *isCorrect_CRL(nullptr), *isCorrectParentLink_CRL(nullptr), *isBestMatchedCorrectParentLink_CRL(nullptr);
    IntVector *bestMatchNHitsTotal_CRL(nullptr), *bestMatchNHitsU_CRL(nullptr), *bestMatchNHitsV_CRL(nullptr), *bestMatchNHitsW_CRL(nullptr);
    IntVector *bestMatchNSharedHitsTotal_CRL(nullptr), *bestMatchNSharedHitsU_CRL(nullptr), *bestMatchNSharedHitsV_CRL(nullptr), *bestMatchNSharedHitsW_CRL(nullptr);
    IntVector *bestMatchNParentTrackHitsTotal_CRL(nullptr), *bestMatchNParentTrackHitsU_CRL(nullptr), *bestMatchNParentTrackHitsV_CRL(nullptr), *bestMatchNParentTrackHitsW_CRL(nullptr);
    IntVector *bestMatchNOtherTrackHitsTotal_CRL(nullptr), *bestMatchNOtherTrackHitsU_CRL(nullptr), *bestMatchNOtherTrackHitsV_CRL(nullptr), *bestMatchNOtherTrackHitsW_CRL(nullptr);
    IntVector *bestMatchNOtherShowerHitsTotal_CRL(nullptr), *bestMatchNOtherShowerHitsU_CRL(nullptr), *bestMatchNOtherShowerHitsV_CRL(nullptr), *bestMatchNOtherShowerHitsW_CRL(nullptr);
    IntVector *totalCRLHitsInBestMatchParentCR_CRL(nullptr), *uCRLHitsInBestMatchParentCR_CRL(nullptr), *vCRLHitsInBestMatchParentCR_CRL(nullptr), *wCRLHitsInBestMatchParentCR_CRL(nullptr);
    
    validationTree->SetBranchAddress("mcE_CRL", &mcE_CRL);
    validationTree->SetBranchAddress("mcPX_CRL", &mcPX_CRL);
    validationTree->SetBranchAddress("mcPY_CRL", &mcPY_CRL);
    validationTree->SetBranchAddress("mcPZ_CRL", &mcPZ_CRL);
    validationTree->SetBranchAddress("nMCHitsTotal_CRL", &nMCHitsTotal_CRL);
    validationTree->SetBranchAddress("nMCHitsU_CRL", &nMCHitsU_CRL);  
    validationTree->SetBranchAddress("nMCHitsV_CRL", &nMCHitsV_CRL);
    validationTree->SetBranchAddress("nMCHitsW_CRL", &nMCHitsW_CRL);
    validationTree->SetBranchAddress("nAboveThresholdMatches_CRL", &nAboveThresholdMatches_CRL);
    validationTree->SetBranchAddress("isCorrect_CRL", &isCorrect_CRL);
    validationTree->SetBranchAddress("isCorrectParentLink_CRL", &isCorrectParentLink_CRL);
    validationTree->SetBranchAddress("isBestMatchedCorrectParentLink_CRL", &isBestMatchedCorrectParentLink_CRL);
    validationTree->SetBranchAddress("bestMatchNHitsTotal_CRL", &bestMatchNHitsTotal_CRL);    
    validationTree->SetBranchAddress("bestMatchNHitsU_CRL", &bestMatchNHitsU_CRL);  
    validationTree->SetBranchAddress("bestMatchNHitsV_CRL", &bestMatchNHitsV_CRL);
    validationTree->SetBranchAddress("bestMatchNHitsW_CRL", &bestMatchNHitsW_CRL);
    validationTree->SetBranchAddress("bestMatchNSharedHitsTotal_CRL", &bestMatchNSharedHitsTotal_CRL);    
    validationTree->SetBranchAddress("bestMatchNSharedHitsU_CRL", &bestMatchNSharedHitsU_CRL);        
    validationTree->SetBranchAddress("bestMatchNSharedHitsV_CRL", &bestMatchNSharedHitsV_CRL);
    validationTree->SetBranchAddress("bestMatchNSharedHitsW_CRL", &bestMatchNSharedHitsW_CRL);
    validationTree->SetBranchAddress("bestMatchNParentTrackHitsTotal_CRL", &bestMatchNParentTrackHitsTotal_CRL);    
    validationTree->SetBranchAddress("bestMatchNParentTrackHitsU_CRL", &bestMatchNParentTrackHitsU_CRL);
    validationTree->SetBranchAddress("bestMatchNParentTrackHitsV_CRL", &bestMatchNParentTrackHitsV_CRL);
    validationTree->SetBranchAddress("bestMatchNParentTrackHitsW_CRL", &bestMatchNParentTrackHitsW_CRL);
    validationTree->SetBranchAddress("bestMatchNOtherTrackHitsTotal_CRL", &bestMatchNOtherTrackHitsTotal_CRL);
    validationTree->SetBranchAddress("bestMatchNOtherTrackHitsU_CRL", &bestMatchNOtherTrackHitsU_CRL);
    validationTree->SetBranchAddress("bestMatchNOtherTrackHitsV_CRL", &bestMatchNOtherTrackHitsV_CRL);
    validationTree->SetBranchAddress("bestMatchNOtherTrackHitsW_CRL", &bestMatchNOtherTrackHitsW_CRL);
    validationTree->SetBranchAddress("bestMatchNOtherShowerHitsTotal_CRL", &bestMatchNOtherShowerHitsTotal_CRL);    
    validationTree->SetBranchAddress("bestMatchNOtherShowerHitsU_CRL", &bestMatchNOtherShowerHitsU_CRL);    
    validationTree->SetBranchAddress("bestMatchNOtherShowerHitsV_CRL", &bestMatchNOtherShowerHitsV_CRL);
    validationTree->SetBranchAddress("bestMatchNOtherShowerHitsW_CRL", &bestMatchNOtherShowerHitsW_CRL);
    validationTree->SetBranchAddress("totalCRLHitsInBestMatchParentCR_CRL", &totalCRLHitsInBestMatchParentCR_CRL);    
    validationTree->SetBranchAddress("uCRLHitsInBestMatchParentCR_CRL", &uCRLHitsInBestMatchParentCR_CRL);
    validationTree->SetBranchAddress("vCRLHitsInBestMatchParentCR_CRL", &vCRLHitsInBestMatchParentCR_CRL);
    validationTree->SetBranchAddress("wCRLHitsInBestMatchParentCR_CRL", &wCRLHitsInBestMatchParentCR_CRL);


    CosmicRayOwnershipMap cosmicRayOwnershipMap;
    //int count(0);
    for (Int_t i(0); i < (Int_t)validationTree->GetEntries(); ++i)
    {
        validationTree->GetEntry(i);

        //++count;

        for (Int_t j = 0; j < cosmicRay.m_nReconstructableChildCRLs; ++j)
        {
            DeltaRay deltaRay;
            deltaRay.m_energy = mcE_CRL->at(j);
            deltaRay.m_momentum = SimpleThreeVector(mcPX_CRL->at(j), mcPY_CRL->at(j), mcPZ_CRL->at(j));
            deltaRay.m_nMCHitsTotal = nMCHitsTotal_CRL->at(j);
            deltaRay.m_nMCHitsU = nMCHitsU_CRL->at(j);
            deltaRay.m_nMCHitsV = nMCHitsV_CRL->at(j);
            deltaRay.m_nMCHitsW = nMCHitsW_CRL->at(j);
            
            if ((cosmicRay.m_momentum.GetMagnitude() < std::numeric_limits<float>::epsilon()) || (deltaRay.m_momentum.GetMagnitude() < std::numeric_limits<float>::epsilon()))
            {
                deltaRay.m_openingAngleFromMuon = -1.f;
            }
            else
            {
                deltaRay.m_openingAngleFromMuon = cosmicRay.m_momentum.GetOpeningAngle(deltaRay.m_momentum);
            }
            
            deltaRay.m_nAboveThresholdMatches = nAboveThresholdMatches_CRL->at(j);
            deltaRay.m_isCorrect = isCorrect_CRL->at(j);
            deltaRay.m_isCorrectParentLink = isCorrectParentLink_CRL->at(j);
            deltaRay.m_isBestMatchedCorrectParentLink = isBestMatchedCorrectParentLink_CRL->at(j);
            
            deltaRay.m_bestMatchNHitsTotal = bestMatchNHitsTotal_CRL->at(j);
            deltaRay.m_bestMatchNHitsU = bestMatchNHitsU_CRL->at(j);
            deltaRay.m_bestMatchNHitsV = bestMatchNHitsV_CRL->at(j);
            deltaRay.m_bestMatchNHitsW = bestMatchNHitsW_CRL->at(j);

            deltaRay.m_bestMatchNSharedHitsTotal = bestMatchNSharedHitsTotal_CRL->at(j);
            deltaRay.m_bestMatchNSharedHitsU = bestMatchNSharedHitsU_CRL->at(j);
            deltaRay.m_bestMatchNSharedHitsV = bestMatchNSharedHitsV_CRL->at(j);
            deltaRay.m_bestMatchNSharedHitsW = bestMatchNSharedHitsW_CRL->at(j);
    
            deltaRay.m_bestMatchNParentTrackHitsTotal = bestMatchNParentTrackHitsTotal_CRL->at(j);
            deltaRay.m_bestMatchNParentTrackHitsU = bestMatchNParentTrackHitsU_CRL->at(j);
            deltaRay.m_bestMatchNParentTrackHitsV = bestMatchNParentTrackHitsV_CRL->at(j);
            deltaRay.m_bestMatchNParentTrackHitsW = bestMatchNParentTrackHitsW_CRL->at(j);

            deltaRay.m_bestMatchNOtherTrackHitsTotal = bestMatchNOtherTrackHitsTotal_CRL->at(j);
            deltaRay.m_bestMatchNOtherTrackHitsU = bestMatchNOtherTrackHitsU_CRL->at(j);
            deltaRay.m_bestMatchNOtherTrackHitsV = bestMatchNOtherTrackHitsV_CRL->at(j);
            deltaRay.m_bestMatchNOtherTrackHitsW = bestMatchNOtherTrackHitsW_CRL->at(j);

            deltaRay.m_bestMatchNOtherShowerHitsTotal = bestMatchNOtherShowerHitsTotal_CRL->at(j);
            deltaRay.m_bestMatchNOtherShowerHitsU = bestMatchNOtherShowerHitsU_CRL->at(j);
            deltaRay.m_bestMatchNOtherShowerHitsV = bestMatchNOtherShowerHitsV_CRL->at(j);
            deltaRay.m_bestMatchNOtherShowerHitsW = bestMatchNOtherShowerHitsW_CRL->at(j);

            deltaRay.m_totalCRLHitsInBestMatchParentCR = totalCRLHitsInBestMatchParentCR_CRL->at(j);
            deltaRay.m_uCRLHitsInBestMatchParentCR = uCRLHitsInBestMatchParentCR_CRL->at(j);
            deltaRay.m_vCRLHitsInBestMatchParentCR = vCRLHitsInBestMatchParentCR_CRL->at(j);
            deltaRay.m_wCRLHitsInBestMatchParentCR = wCRLHitsInBestMatchParentCR_CRL->at(j);
            
            deltaRayVector.push_back(deltaRay);
            cosmicRayOwnershipMap[eventNumber][i + 1].push_back(deltaRay);
        }

        cosmicRayVector.push_back(cosmicRay);
        /*
        if (count == 2)
            break;
        */
    }


    /*
    int bestMatchOtherShowerHitsEventNumber, bestMatchOtherShowerHitsID_CR;
    IntVector *bestMatchOtherShowerHitsID_CRL(nullptr);
    FloatVector *bestMatchOtherShowerHitsDistance_CRL(nullptr);

    otherShowerContaminationTree->SetBranchAddress("bestMatchOtherShowerHitsEventNumber", &bestMatchOtherShowerHitsEventNumber);
    otherShowerContaminationTree->SetBranchAddress("bestMatchOtherShowerHitsID_CR", &bestMatchOtherShowerHitsID_CR);
    otherShowerContaminationTree->SetBranchAddress("bestMatchOtherShowerHitsID_CRL", &bestMatchOtherShowerHitsID_CRL);
    otherShowerContaminationTree->SetBranchAddress("bestMatchOtherShowerHitsDistance_CRL", &bestMatchOtherShowerHitsDistance_CRL);
    
    for (Int_t i(0); i < (Int_t)otherShowerContaminationTree->GetEntries(); ++i)
    {
        otherShowerContaminationTree->GetEntry(i);

        const DeltaRayVector childDeltaRays(cosmicRayOwnershipMap.at(bestMatchOtherShowerHitsEventNumber).at(bestMatchOtherShowerHitsID_CR));
        
        for (Int_t j(0); j < (Int_t)bestMatchOtherShowerHitsID_CRL->size(); ++j)
        {
            const DeltaRay &deltaRayCopy(childDeltaRays.at(bestMatchOtherShowerHitsID_CRL->at(j) - 1));
            auto iter(std::find(deltaRayVector.begin(), deltaRayVector.end(), deltaRayCopy));
            (*iter).m_bestMatchOtherShowerHitsDistance.push_back(bestMatchOtherShowerHitsDistance_CRL->at(j));
        }
        
    }
    
    int bestMatchOtherTrackHitsEventNumber, bestMatchOtherTrackHitsID_CR;
    IntVector *bestMatchOtherTrackHitsID_CRL(nullptr);
    FloatVector *bestMatchOtherTrackHitsDistance_CRL(nullptr);

    otherTrackContaminationTree->SetBranchAddress("bestMatchOtherTrackHitsEventNumber", &bestMatchOtherTrackHitsEventNumber);
    otherTrackContaminationTree->SetBranchAddress("bestMatchOtherTrackHitsID_CR", &bestMatchOtherTrackHitsID_CR);
    otherTrackContaminationTree->SetBranchAddress("bestMatchOtherTrackHitsID_CRL", &bestMatchOtherTrackHitsID_CRL);
    otherTrackContaminationTree->SetBranchAddress("bestMatchOtherTrackHitsDistance_CRL", &bestMatchOtherTrackHitsDistance_CRL);

    for (Int_t i(0); i < (Int_t)otherTrackContaminationTree->GetEntries(); ++i)
    {
        otherTrackContaminationTree->GetEntry(i);

        const DeltaRayVector childDeltaRays(cosmicRayOwnershipMap.at(bestMatchOtherTrackHitsEventNumber).at(bestMatchOtherTrackHitsID_CR));

        for (Int_t j(0); j < (Int_t)bestMatchOtherTrackHitsID_CRL->size(); ++j)
        {
            const DeltaRay &deltaRayCopy(childDeltaRays.at(bestMatchOtherTrackHitsID_CRL->at(j) - 1));
            auto iter(std::find(deltaRayVector.begin(), deltaRayVector.end(), deltaRayCopy));
            (*iter).m_bestMatchOtherTrackHitsDistance.push_back(bestMatchOtherTrackHitsDistance_CRL->at(j));
        }
    }    

    int bestMatchParentTrackHitsEventNumber, bestMatchParentTrackHitsID_CR;
    IntVector *bestMatchParentTrackHitsID_CRL(nullptr);
    FloatVector *bestMatchParentTrackHitsDistance_CRL(nullptr);

    parentTrackContaminationTree->SetBranchAddress("bestMatchParentTrackHitsEventNumber", &bestMatchParentTrackHitsEventNumber);    
    parentTrackContaminationTree->SetBranchAddress("bestMatchParentTrackHitsID_CR", &bestMatchParentTrackHitsID_CR);
    parentTrackContaminationTree->SetBranchAddress("bestMatchParentTrackHitsID_CRL", &bestMatchParentTrackHitsID_CRL);
    parentTrackContaminationTree->SetBranchAddress("bestMatchParentTrackHitsDistance_CRL", &bestMatchParentTrackHitsDistance_CRL);

    for (Int_t i(0); i < (Int_t)parentTrackContaminationTree->GetEntries(); ++i)
    {
        parentTrackContaminationTree->GetEntry(i);

        if (bestMatchParentTrackHitsID_CR > 2)
            break;
        
        const DeltaRayVector childDeltaRays(cosmicRayOwnershipMap.at(bestMatchParentTrackHitsEventNumber).at(bestMatchParentTrackHitsID_CR));

        for (Int_t j(0); j < (Int_t)bestMatchParentTrackHitsID_CRL->size(); ++j)
        {
            const DeltaRay &deltaRayCopy(childDeltaRays.at(bestMatchParentTrackHitsID_CRL->at(j) - 1));
            auto iter(std::find(deltaRayVector.begin(), deltaRayVector.end(), deltaRayCopy));
            (*iter).m_bestMatchParentTrackHitsDistance.push_back(bestMatchParentTrackHitsDistance_CRL->at(j));
        }
    }       
    
    int bestMatchCRLHitsInCREventNumber, bestMatchCRLHitsInCRID_CR;
    IntVector *bestMatchCRLHitsInCRID_CRL(nullptr);
    FloatVector *bestMatchCRLHitsInCRDistance_CRL(nullptr);

    CRLHitsInCRTree->SetBranchAddress("bestMatchCRLHitsInCREventNumber", &bestMatchCRLHitsInCREventNumber);
    CRLHitsInCRTree->SetBranchAddress("bestMatchCRLHitsInCRID_CR", &bestMatchCRLHitsInCRID_CR);
    CRLHitsInCRTree->SetBranchAddress("bestMatchCRLHitsInCRID_CRL", &bestMatchCRLHitsInCRID_CRL);
    CRLHitsInCRTree->SetBranchAddress("bestMatchCRLHitsInCRDistance_CRL", &bestMatchCRLHitsInCRDistance_CRL);

    for (Int_t i(0); i < (Int_t)CRLHitsInCRTree->GetEntries(); ++i)
    {
        CRLHitsInCRTree->GetEntry(i);
        
        if (bestMatchCRLHitsInCRID_CR > 2)
            break;
        
        const DeltaRayVector childDeltaRays(cosmicRayOwnershipMap.at(bestMatchCRLHitsInCREventNumber).at(bestMatchCRLHitsInCRID_CR));

        for (Int_t j(0); j < (Int_t)bestMatchCRLHitsInCRID_CRL->size(); ++j)
        {
            const DeltaRay &deltaRayCopy(childDeltaRays.at(bestMatchCRLHitsInCRID_CRL->at(j) - 1));
            auto iter(std::find(deltaRayVector.begin(), deltaRayVector.end(), deltaRayCopy));            
            (*iter).m_bestMatchCRLHitsInCRDistance.push_back(bestMatchCRLHitsInCRDistance_CRL->at(j));
        }
    }
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DisplayOverallRecoMetrics(const DeltaRayVector &deltaRayVector)
{
    unsigned int nCorrectParentLinks(0), nCorrectlyReconstructedCRLs(0);
    unsigned int nZeroMatches(0), nOneMatches(0), nTwoMatches(0), nThreePlusMatches(0);
    
    for (const DeltaRay &deltaRay : deltaRayVector)
    {
        if (deltaRay.m_isCorrectParentLink)
            ++nCorrectParentLinks;

        if (deltaRay.m_isCorrect)
            ++nCorrectlyReconstructedCRLs;
        
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

    const unsigned int nReconstructableCRLs(deltaRayVector.size());
    const float fCorrectParentLinks(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nCorrectParentLinks) / static_cast<float>(nReconstructableCRLs));
    const float fCorrectlyReconstructedCRLs(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nCorrectlyReconstructedCRLs) / static_cast<float>(nReconstructableCRLs));
    const float fZeroMatches(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nZeroMatches) / static_cast<float>(nReconstructableCRLs));
    const float fOneMatches(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nOneMatches) / static_cast<float>(nReconstructableCRLs));
    const float fTwoMatches(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nTwoMatches) / static_cast<float>(nReconstructableCRLs));
    const float fThreePlusMatches(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nThreePlusMatches) / static_cast<float>(nReconstructableCRLs));

    std::cout << "nReconstructableDeltaRays: " << nReconstructableCRLs << std::endl
              << "nCorrectlyReconstructedDeltaRays: " << nCorrectlyReconstructedCRLs << " (" << fCorrectlyReconstructedCRLs * 100.f << "%)" << std::endl
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
        cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays->Fill(cosmicRay.m_nReconstructableChildCRLs);
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FillDeltaRayMCHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection)
{
   if (!deltaRayMCHistogramCollection.m_hEnergyDistribution)
   {
        deltaRayMCHistogramCollection.m_hEnergyDistribution = new TH1F("hEnergyDistribution_CRL", "hEnergyDistribution_CRL", 100, 0., 1.);
        deltaRayMCHistogramCollection.m_hEnergyDistribution->SetTitle(";TrueCRLEnergy [GeV];Occurance");
        //deltaRayMCHistogramCollection.m_hEnergyDistribution->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayMCHistogramCollection.m_hTotalHitDistribution)
   {
        deltaRayMCHistogramCollection.m_hTotalHitDistribution = new TH1F("hTotalHitDistribution_CRL", "hTotalHitDistribution_CRL", 40000, 0., 100.);
        deltaRayMCHistogramCollection.m_hTotalHitDistribution->SetTitle(";nMCTotalHits;Occurance");
        //deltaRayMCHistogramCollection.m_hTotalHitDistribution->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayMCHistogramCollection.m_hOpeningAngleDistribution)
   {
        deltaRayMCHistogramCollection.m_hOpeningAngleDistribution = new TH1F("hOpeningAngleDistribution_CRL", "hOpeningAngleDistribution_CRL", 40000, 0., 180.);
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
        deltaRayRecoHistogramCollection.m_hCompleteness = new TH1F("hCompleteness_CRL", "hCompleteness_CRL", 100, 0., 1.);
        deltaRayRecoHistogramCollection.m_hCompleteness->SetTitle(";Completeness;Occurance");
        //deltaRayRecoHistogramCollection.m_hCompleteness->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hPurity)
   {
        deltaRayRecoHistogramCollection.m_hPurity = new TH1F("hPurity_CRL", "hPurity_CRL", 100, 0., 1.);
        deltaRayRecoHistogramCollection.m_hPurity->SetTitle(";Purity;Occurance");
        //deltaRayRecoHistogramCollection.m_hPurity->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hAboveThresholdMatches)
   {
        deltaRayRecoHistogramCollection.m_hAboveThresholdMatches = new TH1F("hAboveThresholdMatches_CRL", "hAboveThresholdMatches_CRL", 40000, 0., 4.);
        deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->SetTitle(";nAboveThresholdMatches;Occurance");
        //deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal)
   {
        deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal = new TH1F("hParentTrackHitsTotal_CRL", "hParentTrackHitsTotal_CRL", 40000, 0., 50.);
        deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->SetTitle(";nParentTrackHitsTotal;Occurance");
        //deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal)
   {
        deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal = new TH1F("hOtherTrackHitsTotal_CRL", "hOtherTrackHitsTotal_CRL", 40000, 0., 20.);
        deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->SetTitle(";nOtherTrackHitsTotal;Occurance");
        //deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal)
   {
        deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal = new TH1F("hOtherShowerHitsTotal_CRL", "hOtherShowerHitsTotal_CRL", 40000, 0., 20.);
        deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->SetTitle(";nOtherShowerHitsTotal;Occurance");
        //deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay)
   {
        deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay = new TH1F("hTotalHitsTakenByCosmicRay_CRL", "hTotalHitsTakenByCosmicRay_CRL", 40000, 0., 50.);
        deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->SetTitle(";hTotalHitsTakenByCosmicRay;Occurance");
        //deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hEfficiency_Energy)
   {
        deltaRayRecoHistogramCollection.m_hEfficiency_Energy = new TH1F("hEfficiency_Energy_CRL", "hEfficiency_Energy_CRL", 100, 0., 1.);
        deltaRayRecoHistogramCollection.m_hEfficiency_Energy->SetTitle(";TrueCRLEnergy [GeV];Efficiency");
        //deltaRayRecoHistogramCollection.m_hEfficiency_Energy->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits)
   {
        deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits = new TH1F("hEfficiency_TotalHits_CRL", "hEfficiency_TotalHits_CRL", 40000, 0., 100.);
        deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->SetTitle(";nMCTotalHits;Efficiency");
        //deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle)
   {
        deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle = new TH1F("hEfficiency_OpeningAngle_CRL", "_CRL", 40000, 0., 180.);
        deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->SetTitle(";TrueOpeningAngle [degrees];Efficiency");
        //deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy)
   {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy = new TH1F("hCorrectParentLink_Energy_CRL", "hCorrectParentLink_Energy_CRL", 100, 0., 1.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy->SetTitle(";TrueCRLEnergy [GeV];CorrectParentLinkFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits)
   {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits = new TH1F("hCorrectParentLink_TotalHits_CRL", "hCorrectParentLink_TotalHits_CRL", 40000, 0., 100.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits->SetTitle(";nMCTotalHits;CorrectParentLinkFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle)
   {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle = new TH1F("hCorrectParentLink_OpeningAngle_CRL", "hCorrectParentLink_OpeningAngle_CRL", 40000, 0., 180.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle->SetTitle(";TrueOpeningAngle [degrees];CorrectParentLinkFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle->GetXaxis()->SetRangeUser(0.f, 50.f);
    }      

   if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy)
   {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy = new TH1F("hCorrectEvent_Energy_CRL", "hCorrectEvent_Energy_CRL", 100, 0., 1.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy->SetTitle(";TrueCRLEnergy [GeV];CorrectlyReconstructedFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits)
   {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits = new TH1F("hCorrectEvent_TotalHits_CRL", "hCorrectEvent_TotalHits_CRL", 40000, 0., 100.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits->SetTitle(";nMCTotalHits;CorrectlyReconstructedFraction");
        //deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits->GetXaxis()->SetRangeUser(0.f, 50.f);
    }

   if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle)
   {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle = new TH1F("hCorrectEvent_OpeningAngle_CRL", "hCorrectEvent_OpeningAngle_CRL", 40000, 0., 180.);
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
            deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->Fill(deltaRay.m_totalCRLHitsInBestMatchParentCR);
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

void FillDeltaRayContaminationHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayContaminationHistogramCollection &deltaRayContaminationHistogramCollection)
{
   if (!deltaRayContaminationHistogramCollection.m_hOtherShowerHitsDistance)
   {
        deltaRayContaminationHistogramCollection.m_hOtherShowerHitsDistance = new TH1F("hOtherShowerHitsDistance_CRL", "hOtherShowerHitsDistance_CRL", 40000, 0., 500.);
        deltaRayContaminationHistogramCollection.m_hOtherShowerHitsDistance->SetTitle(";Distance From MC CRL [cm];Occurance");
   }

   if (!deltaRayContaminationHistogramCollection.m_hOtherTrackHitsDistance)
   {
        deltaRayContaminationHistogramCollection.m_hOtherTrackHitsDistance = new TH1F("hOtherTrackHitsDistance_CRL", "hOtherTrackHitsDistance_CRL", 40000, 0., 500.);
        deltaRayContaminationHistogramCollection.m_hOtherTrackHitsDistance->SetTitle(";Distance From MC CRL [cm];Occurance");
   }

   if (!deltaRayContaminationHistogramCollection.m_hParentTrackHitsDistance)
   {
        deltaRayContaminationHistogramCollection.m_hParentTrackHitsDistance = new TH1F("hParentTrackHitsDistance_CRL", "hParentTrackHitsDistance_CRL", 40000, 0., 500.);
        deltaRayContaminationHistogramCollection.m_hParentTrackHitsDistance->SetTitle(";Distance From MC CRL [cm];Occurance");
   }  

   if (!deltaRayContaminationHistogramCollection.m_hCRLHitsInCRDistance)
   {
        deltaRayContaminationHistogramCollection.m_hCRLHitsInCRDistance = new TH1F("hCRLHitsInCRDistance_CRL", "hCRLHitsInCRDistance_CRL", 40000, 0., 500.);
        deltaRayContaminationHistogramCollection.m_hCRLHitsInCRDistance->SetTitle(";Distance From MC CR [cm];Occurance");
   }     

   for (const DeltaRay &deltaRay : deltaRayVector)
   {
       for (const float distance : deltaRay.m_bestMatchOtherShowerHitsDistance)
           deltaRayContaminationHistogramCollection.m_hOtherShowerHitsDistance->Fill(distance);

       for (const float distance : deltaRay.m_bestMatchOtherTrackHitsDistance)
           deltaRayContaminationHistogramCollection.m_hOtherTrackHitsDistance->Fill(distance);

       for (const float distance : deltaRay.m_bestMatchParentTrackHitsDistance)
           deltaRayContaminationHistogramCollection.m_hParentTrackHitsDistance->Fill(distance);       
       
       for (const float distance : deltaRay.m_bestMatchCRLHitsInCRDistance)
           deltaRayContaminationHistogramCollection.m_hCRLHitsInCRDistance->Fill(distance);
   }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteHistograms(CosmicRayMCHistogramCollection &cosmicRayMCHistogramCollection, DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection,
    DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection, DeltaRayContaminationHistogramCollection &deltaRayContaminationHistogramCollection)
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

    deltaRayContaminationHistogramCollection.m_hOtherShowerHitsDistance->Write("hOtherShowerHitsDistance");
    deltaRayContaminationHistogramCollection.m_hOtherTrackHitsDistance->Write("hOtherTrackHitsDistance");
    deltaRayContaminationHistogramCollection.m_hParentTrackHitsDistance->Write("hParentTrackHitsDistance");
    deltaRayContaminationHistogramCollection.m_hCRLHitsInCRDistance->Write("hCRLHitsInCRDistance");
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

float SimpleThreeVector::GetMagnitude() const
{
    return std::sqrt((this->m_x * this->m_x) + (this->m_y * this->m_y) + (this->m_z * this->m_z));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float SimpleThreeVector::GetOpeningAngle(const SimpleThreeVector &otherVector)
{
    const float dotProduct = (this->m_x * otherVector.m_x) + (this->m_y * otherVector.m_y) + (this->m_z * otherVector.m_z);
    const float thisMagnitude = this->GetMagnitude();
    const float otherMagnitude = otherVector.GetMagnitude();

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
