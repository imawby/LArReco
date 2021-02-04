/**
 *  @file   LArReco/validation/Validation.C
 *
 *  @brief  Implementation of validation functionality
 *
 *  $Log: $
 */
#include "TH1F.h"
#include "TH2F.h"
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
        CreateCosmicRayMCHistogramCollection(cosmicRayMCHistogramCollection);
        FillCosmicRayMCHistogramCollection(cosmicRayVector, cosmicRayMCHistogramCollection);

        DeltaRayMCHistogramCollection deltaRayMCHistogramCollection;
        CreateDeltaRayMCHistogramCollection(deltaRayMCHistogramCollection);
        FillDeltaRayMCHistogramCollection(deltaRayVector, deltaRayMCHistogramCollection);

        DeltaRayRecoHistogramCollection deltaRayRecoHistogramCollection;
        CreateDeltaRayRecoHistogramCollection(deltaRayRecoHistogramCollection);
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

    CosmicRay cosmicRay;
    validationTree->SetBranchAddress("eventNumber", &cosmicRay.m_eventNumber);
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
    IntVector *nAboveThresholdMatches_CRL(nullptr), *isCorrect_CRL(nullptr), *isCorrectParentLink_CRL(nullptr);
    IntVector *bestMatchNHitsTotal_CRL(nullptr), *bestMatchNHitsU_CRL(nullptr), *bestMatchNHitsV_CRL(nullptr), *bestMatchNHitsW_CRL(nullptr);
    IntVector *bestMatchNSharedHitsTotal_CRL(nullptr), *bestMatchNSharedHitsU_CRL(nullptr), *bestMatchNSharedHitsV_CRL(nullptr), *bestMatchNSharedHitsW_CRL(nullptr);
    IntVector *bestMatchNParentTrackHitsTotal_CRL(nullptr), *bestMatchNParentTrackHitsU_CRL(nullptr), *bestMatchNParentTrackHitsV_CRL(nullptr), *bestMatchNParentTrackHitsW_CRL(nullptr);
    IntVector *bestMatchNOtherTrackHitsTotal_CRL(nullptr), *bestMatchNOtherTrackHitsU_CRL(nullptr), *bestMatchNOtherTrackHitsV_CRL(nullptr), *bestMatchNOtherTrackHitsW_CRL(nullptr);
    IntVector *bestMatchNOtherShowerHitsTotal_CRL(nullptr), *bestMatchNOtherShowerHitsU_CRL(nullptr), *bestMatchNOtherShowerHitsV_CRL(nullptr), *bestMatchNOtherShowerHitsW_CRL(nullptr);
    IntVector *totalCRLHitsInBestMatchParentCR_CRL(nullptr), *uCRLHitsInBestMatchParentCR_CRL(nullptr), *vCRLHitsInBestMatchParentCR_CRL(nullptr), *wCRLHitsInBestMatchParentCR_CRL(nullptr);

    IntVector *bestMatchOtherShowerHitsID_CRL(nullptr), *bestMatchOtherTrackHitsID_CRL(nullptr), *bestMatchParentTrackHitsID_CRL(nullptr), *bestMatchCRLHitsInCRID_CRL(nullptr);
    FloatVector *bestMatchOtherShowerHitsDistance_CRL(nullptr), *bestMatchOtherTrackHitsDistance_CRL(nullptr), *bestMatchParentTrackHitsDistance_CRL(nullptr);
    FloatVector *bestMatchCRLHitsInCRDistance_CRL(nullptr);
    
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

    validationTree->SetBranchAddress("bestMatchOtherShowerHitsID_CRL", &bestMatchOtherShowerHitsID_CRL);
    validationTree->SetBranchAddress("bestMatchOtherShowerHitsDistance_CRL", &bestMatchOtherShowerHitsDistance_CRL);
    validationTree->SetBranchAddress("bestMatchOtherTrackHitsID_CRL", &bestMatchOtherTrackHitsID_CRL);
    validationTree->SetBranchAddress("bestMatchOtherTrackHitsDistance_CRL", &bestMatchOtherTrackHitsDistance_CRL);
    validationTree->SetBranchAddress("bestMatchParentTrackHitsID_CRL", &bestMatchParentTrackHitsID_CRL);
    validationTree->SetBranchAddress("bestMatchParentTrackHitsDistance_CRL", &bestMatchParentTrackHitsDistance_CRL);
    validationTree->SetBranchAddress("bestMatchCRLHitsInCRID_CRL", &bestMatchCRLHitsInCRID_CRL);    
    validationTree->SetBranchAddress("bestMatchCRLHitsInCRDistance_CRL", &bestMatchCRLHitsInCRDistance_CRL);
    
    for (Int_t i(0); i < (Int_t)validationTree->GetEntries(); ++i)
    {
        validationTree->GetEntry(i);

        const float pTot(std::sqrt(cosmicRay.m_momentum.m_x * cosmicRay.m_momentum.m_x + cosmicRay.m_momentum.m_y * cosmicRay.m_momentum.m_y +
            cosmicRay.m_momentum.m_z * cosmicRay.m_momentum.m_z));

        float theta0XZ = std::atan2(cosmicRay.m_momentum.m_x, cosmicRay.m_momentum.m_z);
        theta0XZ *= (180.f / M_PI);
        float theta0YZ = std::asin(cosmicRay.m_momentum.m_y / pTot);
        theta0YZ *= (180.f / M_PI);
        
        for (Int_t j = 0; j < cosmicRay.m_nReconstructableChildCRLs; ++j)
        {
            DeltaRay deltaRay;
            deltaRay.m_energy = mcE_CRL->at(j);
            deltaRay.m_momentum = SimpleThreeVector(mcPX_CRL->at(j), mcPY_CRL->at(j), mcPZ_CRL->at(j));
            deltaRay.m_nMCHitsTotal = nMCHitsTotal_CRL->at(j);
            deltaRay.m_nMCHitsU = nMCHitsU_CRL->at(j);
            deltaRay.m_nMCHitsV = nMCHitsV_CRL->at(j);
            deltaRay.m_nMCHitsW = nMCHitsW_CRL->at(j);

            deltaRay.m_parentMuonTheta0XZ = theta0XZ;
            deltaRay.m_parentMuonTheta0YZ = theta0YZ;
            
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

            if (bestMatchOtherShowerHitsID_CRL->size() != bestMatchOtherShowerHitsDistance_CRL->size())
            {
                std::cout << "OTHER SHOWER HITS VECTORS ARE UNEQUAL LENGHTS" << std::endl;
                throw;
            }

            for (unsigned int k = 0; k < bestMatchOtherShowerHitsID_CRL->size(); ++k)
            {
                if (bestMatchOtherShowerHitsID_CRL->at(k) == (j+1))
                    deltaRay.m_bestMatchOtherShowerHitsDistance.push_back(bestMatchOtherShowerHitsDistance_CRL->at(k));
            }

            if (bestMatchOtherTrackHitsID_CRL->size() != bestMatchOtherTrackHitsDistance_CRL->size())
            {
                std::cout << "OTHER TRACK HITS VECTORS ARE UNEQUAL LENGHTS" << std::endl;
                throw;
            }
            
            for (unsigned int k = 0; k < bestMatchOtherTrackHitsID_CRL->size(); ++k)
            {
                if (bestMatchOtherTrackHitsID_CRL->at(k) == (j+1))
                    deltaRay.m_bestMatchOtherTrackHitsDistance.push_back(bestMatchOtherTrackHitsDistance_CRL->at(k));
            }            

            if (bestMatchParentTrackHitsID_CRL->size() != bestMatchParentTrackHitsDistance_CRL->size())
            {
                std::cout << "PARENT TRACK HITS VECTORS ARE UNEQUAL LENGHTS" << std::endl;
                throw;
            }
            
            for (unsigned int k = 0; k < bestMatchParentTrackHitsID_CRL->size(); ++k)
            {
                if (bestMatchParentTrackHitsID_CRL->at(k) == (j+1))
                    deltaRay.m_bestMatchParentTrackHitsDistance.push_back(bestMatchParentTrackHitsDistance_CRL->at(k));
            }
            
            if (bestMatchCRLHitsInCRID_CRL->size() != bestMatchCRLHitsInCRDistance_CRL->size())
            {
                std::cout << "HITS IN PARENT CR VECTORS ARE UNEQUAL LENGHTS" << std::endl;
                throw;
            }
            
            for (unsigned int k = 0; k < bestMatchCRLHitsInCRID_CRL->size(); ++k)
            {
                if (bestMatchCRLHitsInCRID_CRL->at(k) == (j+1))
                    deltaRay.m_bestMatchCRLHitsInCRDistance.push_back(bestMatchCRLHitsInCRDistance_CRL->at(k));
            }            
            
            deltaRayVector.push_back(deltaRay);
        }

        cosmicRayVector.push_back(cosmicRay);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DisplayOverallRecoMetrics(const DeltaRayVector &deltaRayVector)
{
    unsigned int nCorrectParentLinks(0), nCorrectlyReconstructedCRLs(0);
    unsigned int nZeroMatches(0), nOneMatches(0), nTwoMatches(0), nThreePlusMatches(0), nAboveThresholdMatches(0);

    unsigned int CRHitsInDRTotal(0), DRHitsInCRTotal(0);

    float completenessSum(0), puritySum(0);
    
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

        if (deltaRay.m_nAboveThresholdMatches > 0)
        {
            ++nAboveThresholdMatches;
            completenessSum += static_cast<float>(deltaRay.m_bestMatchNSharedHitsTotal) / static_cast<float>(deltaRay.m_nMCHitsTotal);
            puritySum += static_cast<float>(deltaRay.m_bestMatchNSharedHitsTotal) / static_cast<float>(deltaRay.m_bestMatchNHitsTotal);
            CRHitsInDRTotal += deltaRay.m_bestMatchNParentTrackHitsTotal;
        }

        if (deltaRay.m_isCorrectParentLink)
            DRHitsInCRTotal += deltaRay.m_totalCRLHitsInBestMatchParentCR;
    }

    const float averageCompleteness(completenessSum / static_cast<float>(nAboveThresholdMatches));
    const float averagePurity(puritySum / static_cast<float>(nAboveThresholdMatches));
    
    const unsigned int nReconstructableCRLs(deltaRayVector.size());
    const float fCorrectParentLinks(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nCorrectParentLinks) / static_cast<float>(nReconstructableCRLs));
    const float fCorrectlyReconstructedCRLs(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nCorrectlyReconstructedCRLs) / static_cast<float>(nReconstructableCRLs));
    const float fZeroMatches(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nZeroMatches) / static_cast<float>(nReconstructableCRLs));
    const float fOneMatches(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nOneMatches) / static_cast<float>(nReconstructableCRLs));
    const float fTwoMatches(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nTwoMatches) / static_cast<float>(nReconstructableCRLs));
    const float fThreePlusMatches(nReconstructableCRLs == 0 ? 0.f : static_cast<float>(nThreePlusMatches) / static_cast<float>(nReconstructableCRLs));

    std::cout << "averageCompleteness: " << averageCompleteness << std::endl;
    std::cout << "averagePurity: " << averagePurity << std::endl;
    
    std::cout << "CRHitsInDRTotal: " << CRHitsInDRTotal << std::endl;
    std::cout << "DRHitsInCRTotal: " << DRHitsInCRTotal << std::endl;
    
    std::cout << "nReconstructableDeltaRays: " << nReconstructableCRLs << std::endl
              << "nCorrectlyReconstructedDeltaRays: " << nCorrectlyReconstructedCRLs << " (" << fCorrectlyReconstructedCRLs * 100.f << "%)" << std::endl
              << "nCorrectParentLinks: " << nCorrectParentLinks << " (" << fCorrectParentLinks * 100.f << "%)" << std::endl
              << "Above threshold pfo match distribution: "
              << "|0: " << fZeroMatches * 100.f << "%|, |1: " << fOneMatches * 100.f << "%|, |2: " << fTwoMatches * 100.f << "%|, |3+: " << fThreePlusMatches * 100.f <<  "%|" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FillCosmicRayMCHistogramCollection(const CosmicRayVector &cosmicRayVector, CosmicRayMCHistogramCollection &cosmicRayMCHistogramCollection)
{
   int eventNumberCounter(0), nMuonsInReadoutWindow(0), nCRLsInReadoutWindow(0);
    
   for (const CosmicRay &cosmicRay : cosmicRayVector)
   {
       if (cosmicRay.m_eventNumber != eventNumberCounter)
       {
           cosmicRayMCHistogramCollection.m_hTotalCRsWithReconstructableCRLs->Fill(nMuonsInReadoutWindow);
           cosmicRayMCHistogramCollection.m_hTotalReconstructableCRLs->Fill(nCRLsInReadoutWindow);
           
           nMuonsInReadoutWindow = 0;
           nCRLsInReadoutWindow = 0;
           
           if (eventNumberCounter == 9)
           {
               eventNumberCounter = 0;
           }
           else
           {
               ++eventNumberCounter;
           }
       }

       ++nMuonsInReadoutWindow;
       nCRLsInReadoutWindow += cosmicRay.m_nReconstructableChildCRLs;
       
       cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays->Fill(cosmicRay.m_nReconstructableChildCRLs);
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FillDeltaRayMCHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection)
{
    for (const DeltaRay &deltaRay : deltaRayVector)
    {
        int lowestViewNHits(0);
        if (deltaRay.m_nMCHitsU < deltaRay.m_nMCHitsV)
        {
            lowestViewNHits = (deltaRay.m_nMCHitsU < deltaRay.m_nMCHitsW ? deltaRay.m_nMCHitsU : deltaRay.m_nMCHitsW);
        }
        else 
        {
            lowestViewNHits = (deltaRay.m_nMCHitsV < deltaRay.m_nMCHitsW ? deltaRay.m_nMCHitsV : deltaRay.m_nMCHitsW);
        }
        
        deltaRayMCHistogramCollection.m_hEnergyDistribution->Fill(deltaRay.m_energy);
        deltaRayMCHistogramCollection.m_hTotalHitDistribution->Fill(deltaRay.m_nMCHitsTotal);
        deltaRayMCHistogramCollection.m_hOpeningAngleDistribution->Fill(deltaRay.m_openingAngleFromMuon);
        deltaRayMCHistogramCollection.m_hParentMuonTheta0XZDistribution->Fill(deltaRay.m_parentMuonTheta0XZ);  
        deltaRayMCHistogramCollection.m_hParentMuonTheta0YZDistribution->Fill(deltaRay.m_parentMuonTheta0YZ);
        deltaRayMCHistogramCollection.m_hLowestViewNHitsDistribution->Fill(lowestViewNHits);
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FillDeltaRayRecoHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection)
{
    for (const DeltaRay &deltaRay : deltaRayVector)
    {
        const float completeness(static_cast<float>(deltaRay.m_bestMatchNSharedHitsTotal) / static_cast<float>(deltaRay.m_nMCHitsTotal));
        const float purity(static_cast<float>(deltaRay.m_bestMatchNSharedHitsTotal) / static_cast<float>(deltaRay.m_bestMatchNHitsTotal));

        const float completenessU(deltaRay.m_nMCHitsU == 0 ? 1.f : static_cast<float>(deltaRay.m_bestMatchNSharedHitsU) / static_cast<float>(deltaRay.m_nMCHitsU));
        const float completenessV(deltaRay.m_nMCHitsV == 0 ? 1.f : static_cast<float>(deltaRay.m_bestMatchNSharedHitsV) / static_cast<float>(deltaRay.m_nMCHitsV));
        const float completenessW(deltaRay.m_nMCHitsW == 0 ? 1.f : static_cast<float>(deltaRay.m_bestMatchNSharedHitsW) / static_cast<float>(deltaRay.m_nMCHitsW));

        float lowestCompletenessView(0.f);
        if (completenessU < completenessV)
        {
            lowestCompletenessView = (completenessU < completenessW ? completenessU : completenessW);
        }
        else 
        {
            lowestCompletenessView = (completenessV < completenessW ? completenessV : completenessW);
        }

        int lowestViewNHits(0);
        if (deltaRay.m_nMCHitsU < deltaRay.m_nMCHitsV)
        {
            lowestViewNHits = (deltaRay.m_nMCHitsU < deltaRay.m_nMCHitsW ? deltaRay.m_nMCHitsU : deltaRay.m_nMCHitsW);
        }
        else 
        {
            lowestViewNHits = (deltaRay.m_nMCHitsV < deltaRay.m_nMCHitsW ? deltaRay.m_nMCHitsV : deltaRay.m_nMCHitsW);
        }
        
        deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->Fill(deltaRay.m_nAboveThresholdMatches);

        if (deltaRay.m_nAboveThresholdMatches == 0)
        {
            deltaRayRecoHistogramCollection.m_hMatches0_TotalHits->Fill(deltaRay.m_nMCHitsTotal);
            deltaRayRecoHistogramCollection.m_hMatches0_OpeningAngle->Fill(deltaRay.m_openingAngleFromMuon);
            deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0XZ->Fill(deltaRay.m_parentMuonTheta0XZ);
            deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0YZ->Fill(deltaRay.m_parentMuonTheta0YZ);
            deltaRayRecoHistogramCollection.m_hMatches0_LowestViewNHits->Fill(lowestViewNHits);
        }
        else if (deltaRay.m_nAboveThresholdMatches == 1)
        {
            deltaRayRecoHistogramCollection.m_hMatches1_TotalHits->Fill(deltaRay.m_nMCHitsTotal);
            deltaRayRecoHistogramCollection.m_hMatches1_OpeningAngle->Fill(deltaRay.m_openingAngleFromMuon);
            deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0XZ->Fill(deltaRay.m_parentMuonTheta0XZ);
            deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0YZ->Fill(deltaRay.m_parentMuonTheta0YZ);
            deltaRayRecoHistogramCollection.m_hMatches1_LowestViewNHits->Fill(lowestViewNHits);
        }
        else
        {
            deltaRayRecoHistogramCollection.m_hMatchesMultiple_TotalHits->Fill(deltaRay.m_nMCHitsTotal);
            deltaRayRecoHistogramCollection.m_hMatchesMultiple_OpeningAngle->Fill(deltaRay.m_openingAngleFromMuon);
            deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0XZ->Fill(deltaRay.m_parentMuonTheta0XZ);
            deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0YZ->Fill(deltaRay.m_parentMuonTheta0YZ);
            deltaRayRecoHistogramCollection.m_hMatchesMultiple_LowestViewNHits->Fill(lowestViewNHits);
        }

        if (deltaRay.m_nAboveThresholdMatches > 0)
        {
            deltaRayRecoHistogramCollection.m_hCompleteness->Fill(completeness);
            deltaRayRecoHistogramCollection.m_hLowestCompletenessView->Fill(lowestCompletenessView);
            deltaRayRecoHistogramCollection.m_hPurity->Fill(purity);
            
            deltaRayRecoHistogramCollection.m_hCompletenessVsHits->Fill(completeness, deltaRay.m_nMCHitsTotal, 1.f);
            
            deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->Fill(deltaRay.m_bestMatchNParentTrackHitsTotal);
            deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->Fill(deltaRay.m_bestMatchNOtherTrackHitsTotal);
            deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->Fill(deltaRay.m_bestMatchNOtherShowerHitsTotal);
            
            deltaRayRecoHistogramCollection.m_hEfficiency_Energy->Fill(deltaRay.m_energy);
            deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->Fill(deltaRay.m_nMCHitsTotal);            
            deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->Fill(deltaRay.m_openingAngleFromMuon);
            deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0XZ->Fill(deltaRay.m_parentMuonTheta0XZ);
            deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0YZ->Fill(deltaRay.m_parentMuonTheta0YZ);
            deltaRayRecoHistogramCollection.m_hEfficiency_LowestViewNHits->Fill(lowestViewNHits);   
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
            deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0XZ->Fill(deltaRay.m_parentMuonTheta0XZ);
            deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0YZ->Fill(deltaRay.m_parentMuonTheta0YZ);
            deltaRayRecoHistogramCollection.m_hCorrectEvent_LowestViewNHits->Fill(lowestViewNHits);   
        }
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessHistograms(DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection, DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection)
{
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_Energy, deltaRayMCHistogramCollection.m_hEnergyDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits, deltaRayMCHistogramCollection.m_hTotalHitDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle, deltaRayMCHistogramCollection.m_hOpeningAngleDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0XZ, deltaRayMCHistogramCollection.m_hParentMuonTheta0XZDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0YZ, deltaRayMCHistogramCollection.m_hParentMuonTheta0YZDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hEfficiency_LowestViewNHits, deltaRayMCHistogramCollection.m_hLowestViewNHitsDistribution);
    
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy, deltaRayMCHistogramCollection.m_hEnergyDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits, deltaRayMCHistogramCollection.m_hTotalHitDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle, deltaRayMCHistogramCollection.m_hOpeningAngleDistribution);

    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy, deltaRayMCHistogramCollection.m_hEnergyDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits, deltaRayMCHistogramCollection.m_hTotalHitDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle, deltaRayMCHistogramCollection.m_hOpeningAngleDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0XZ, deltaRayMCHistogramCollection.m_hParentMuonTheta0XZDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0YZ, deltaRayMCHistogramCollection.m_hParentMuonTheta0YZDistribution);
    DivideHistogram(deltaRayRecoHistogramCollection.m_hCorrectEvent_LowestViewNHits, deltaRayMCHistogramCollection.m_hLowestViewNHitsDistribution);    

    deltaRayRecoHistogramCollection.m_hCompleteness->Scale(1. / static_cast<double>(deltaRayRecoHistogramCollection.m_hCompleteness->GetEntries()));
    deltaRayRecoHistogramCollection.m_hLowestCompletenessView->Scale(1. / static_cast<double>(deltaRayRecoHistogramCollection.m_hLowestCompletenessView->GetEntries()));
    deltaRayRecoHistogramCollection.m_hPurity->Scale(1. / static_cast<double>(deltaRayRecoHistogramCollection.m_hPurity->GetEntries()));
    deltaRayRecoHistogramCollection.m_hCompletenessVsHits->Scale(1. / static_cast<double>(deltaRayRecoHistogramCollection.m_hCompletenessVsHits->GetEntries()));

    deltaRayRecoHistogramCollection.m_hMatches0_TotalHits->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches0_TotalHits->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatches0_OpeningAngle->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches0_OpeningAngle->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0XZ->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0XZ->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0YZ->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0YZ->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatches0_LowestViewNHits->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches0_LowestViewNHits->GetEntries()));

    deltaRayRecoHistogramCollection.m_hMatches1_TotalHits->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches1_TotalHits->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatches1_OpeningAngle->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches1_OpeningAngle->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0XZ->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0XZ->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0YZ->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0YZ->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatches1_LowestViewNHits->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatches1_LowestViewNHits->GetEntries()));

    deltaRayRecoHistogramCollection.m_hMatchesMultiple_TotalHits->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatchesMultiple_TotalHits->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatchesMultiple_OpeningAngle->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatchesMultiple_OpeningAngle->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0XZ->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0XZ->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0YZ->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0YZ->GetEntries()));
    deltaRayRecoHistogramCollection.m_hMatchesMultiple_LowestViewNHits->Scale(1.f / static_cast<float>(deltaRayRecoHistogramCollection.m_hMatchesMultiple_LowestViewNHits->GetEntries()));    
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

    cosmicRayMCHistogramCollection.m_hTotalCRsWithReconstructableCRLs->Write("hTotalCRsWithReconstructableCRLs");     
    cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays->Write("hReconstructableChildDeltaRays");
    cosmicRayMCHistogramCollection.m_hTotalReconstructableCRLs->Write("hTotalReconstructableCRLs");         

    deltaRayMCHistogramCollection.m_hEnergyDistribution->Write("hEnergyDistribution");
    deltaRayMCHistogramCollection.m_hTotalHitDistribution->Write("hTotalHitsDistribution");
    deltaRayMCHistogramCollection.m_hOpeningAngleDistribution->Write("hOpeningAngleDistribution");
    deltaRayMCHistogramCollection.m_hParentMuonTheta0XZDistribution->Write("hParentMuonTheta0XZDistribution");
    deltaRayMCHistogramCollection.m_hParentMuonTheta0YZDistribution->Write("hParentMuonTheta0YZDistribution");
    deltaRayMCHistogramCollection.m_hLowestViewNHitsDistribution->Write("hLowestViewNHitsDistribution");
    
    deltaRayRecoHistogramCollection.m_hCompleteness->Write("hCompleteness");
    deltaRayRecoHistogramCollection.m_hLowestCompletenessView->Write("hLowestCompletenessView");
    deltaRayRecoHistogramCollection.m_hPurity->Write("hPurity");

    deltaRayRecoHistogramCollection.m_hCompletenessVsHits->Write("hCompletenessVsHits");
    deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->Write("hAboveThresholdMatches");
    
    deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->Write("hParentTrackHitsTotal");
    deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->Write("hOtherTrackHitsTotal");
    deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->Write("hOtherShowerHitsTotal");
    deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->Write("hTotalHitsTakenByCosmicRay");
    
    deltaRayRecoHistogramCollection.m_hEfficiency_Energy->Write("hEfficiency_Energy");
    deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->Write("hEfficiency_TotalHits");
    deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->Write("hEfficiency_OpeningAngle");
    deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0XZ->Write("hEfficiency_ParentMuonTheta0XZ");
    deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0YZ->Write("hEfficiency_ParentMuonTheta0YZ");
    deltaRayRecoHistogramCollection.m_hEfficiency_LowestViewNHits->Write("hEfficiency_LowestViewNHits");
    
    deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy->Write("hCorrectParentLink_Energy");
    deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits->Write("hCorrectParentLink_TotalHits");
    deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle->Write("hCorrectParentLink_OpeningAngle");
    
    deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy->Write("hCorrectEvent_Energy");
    deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits->Write("hCorrectEvent_TotalHits");
    deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle->Write("hCorrectEvent_OpeningAngle");
    deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0XZ->Write("hCorrectEvent_ParentMuonTheta0XZ");
    deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0YZ->Write("hCorrectEvent_ParentMuonTheta0YZ");
    deltaRayRecoHistogramCollection.m_hCorrectEvent_LowestViewNHits->Write("hCorrectEvent_LowestViewNHits");

    deltaRayContaminationHistogramCollection.m_hOtherShowerHitsDistance->Write("hOtherShowerHitsDistance");
    deltaRayContaminationHistogramCollection.m_hOtherTrackHitsDistance->Write("hOtherTrackHitsDistance");
    deltaRayContaminationHistogramCollection.m_hParentTrackHitsDistance->Write("hParentTrackHitsDistance");
    deltaRayContaminationHistogramCollection.m_hCRLHitsInCRDistance->Write("hCRLHitsInCRDistance");

    deltaRayRecoHistogramCollection.m_hMatches0_TotalHits->Write("hMatches0_TotalHits");
    deltaRayRecoHistogramCollection.m_hMatches0_OpeningAngle->Write("hMatches0_OpeningAngle");
    deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0XZ->Write("hMatches0_ParentMuonTheta0XZ");
    deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0YZ->Write("hMatches0_ParentMuonTheta0YZ");
    deltaRayRecoHistogramCollection.m_hMatches0_LowestViewNHits->Write("hMatches0_LowestViewNHits");

    deltaRayRecoHistogramCollection.m_hMatches1_TotalHits->Write("hMatches1_TotalHits");
    deltaRayRecoHistogramCollection.m_hMatches1_OpeningAngle->Write("hMatches1_OpeningAngle");
    deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0XZ->Write("hMatches1_ParentMuonTheta0XZ");
    deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0YZ->Write("hMatches1_ParentMuonTheta0YZ");
    deltaRayRecoHistogramCollection.m_hMatches1_LowestViewNHits->Write("hMatches1_LowestViewNHits");

    deltaRayRecoHistogramCollection.m_hMatchesMultiple_TotalHits->Write("hMatchesMultiple_TotalHits");
    deltaRayRecoHistogramCollection.m_hMatchesMultiple_OpeningAngle->Write("hMatchesMultiple_OpeningAngle");
    deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0XZ->Write("hMatchesMultiple_ParentMuonTheta0XZ");
    deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0YZ->Write("hMatchesMultiple_ParentMuonTheta0YZ");
    deltaRayRecoHistogramCollection.m_hMatchesMultiple_LowestViewNHits->Write("hMatchesMultiple_LowestViewNHits");    
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void CreateCosmicRayMCHistogramCollection(CosmicRayMCHistogramCollection &cosmicRayMCHistogramCollection)
{
   if (!cosmicRayMCHistogramCollection.m_hTotalCRsWithReconstructableCRLs)
   {
        cosmicRayMCHistogramCollection.m_hTotalCRsWithReconstructableCRLs = new TH1F("hTotalCRsWithReconstructableCRLs_CR", "hTotalCRsWithReconstructableCRLs_CR", 50, 0., 100.);
        cosmicRayMCHistogramCollection.m_hTotalCRsWithReconstructableCRLs->SetTitle(";nCRsWithReconstructableCRLsInReadoutWindow;Occurance");
   }
    
   if (!cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays)
   {
        cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays = new TH1F("hReconstructableChildDeltaRays_CR", "hReconstructableChildDeltaRays_CR", 40000, 0., 10.);
        cosmicRayMCHistogramCollection.m_hReconstructableChildDeltaRays->SetTitle(";nReconstructableChildDeltaRays;Occurance");
   }

   if (!cosmicRayMCHistogramCollection.m_hTotalReconstructableCRLs)
   {
        cosmicRayMCHistogramCollection.m_hTotalReconstructableCRLs = new TH1F("hTotalReconstructableCRLs_CR", "hTotalReconstructableCRLs_CR", 50, 0., 100.);
        cosmicRayMCHistogramCollection.m_hTotalReconstructableCRLs->SetTitle(";nReconstructableCRLsInReadoutWindow;Occurance");
   }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateDeltaRayMCHistogramCollection(DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection)
{
    if (!deltaRayMCHistogramCollection.m_hEnergyDistribution)
    {
        deltaRayMCHistogramCollection.m_hEnergyDistribution = new TH1F("hEnergyDistribution_CRL", "hEnergyDistribution_CRL", 100, 0., 1.);
        deltaRayMCHistogramCollection.m_hEnergyDistribution->SetTitle(";TrueCRLEnergy [GeV];Occurance");
    }

    if (!deltaRayMCHistogramCollection.m_hTotalHitDistribution)
    {
        deltaRayMCHistogramCollection.m_hTotalHitDistribution = new TH1F("hTotalHitDistribution_CRL", "hTotalHitDistribution_CRL", 100, 0., 100.);
        deltaRayMCHistogramCollection.m_hTotalHitDistribution->SetTitle(";nMCTotalHits;Occurance");
    }

    if (!deltaRayMCHistogramCollection.m_hOpeningAngleDistribution)
    {
        deltaRayMCHistogramCollection.m_hOpeningAngleDistribution = new TH1F("hOpeningAngleDistribution_CRL", "hOpeningAngleDistribution_CRL", 100, 0., 180.);
        deltaRayMCHistogramCollection.m_hOpeningAngleDistribution->SetTitle(";TrueOpeningAngle [degrees];Occurance");
    }      

    if (!deltaRayMCHistogramCollection.m_hParentMuonTheta0XZDistribution)
    {
        deltaRayMCHistogramCollection.m_hParentMuonTheta0XZDistribution = new TH1F("hParentMuonTheta0XZDistribution_CRL", "hParentMuonTheta0XZDistribution_CRL", 100, -180., 180.);
        deltaRayMCHistogramCollection.m_hParentMuonTheta0XZDistribution->SetTitle(";ParentMuonTheta0XZ [degrees];Occurance");
    }

    if (!deltaRayMCHistogramCollection.m_hParentMuonTheta0YZDistribution)
    {
        deltaRayMCHistogramCollection.m_hParentMuonTheta0YZDistribution = new TH1F("hParentMuonTheta0YZDistribution_CRL", "hParentMuonTheta0YZDistribution_CRL", 200, -180., 180.);
        deltaRayMCHistogramCollection.m_hParentMuonTheta0YZDistribution->SetTitle(";ParentMuonTheta0YZ [degrees];Occurance");
    }

    if (!deltaRayMCHistogramCollection.m_hLowestViewNHitsDistribution)
    {
        deltaRayMCHistogramCollection.m_hLowestViewNHitsDistribution = new TH1F("hLowestViewNHitsDistribution_CRL", "hLowestViewNHitsDistribution_CRL", 100, 0., 100.);
        deltaRayMCHistogramCollection.m_hLowestViewNHitsDistribution->SetTitle("; nMCHits;Occurance");
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateDeltaRayRecoHistogramCollection(DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection)
{
    if (!deltaRayRecoHistogramCollection.m_hCompleteness)
    {
        deltaRayRecoHistogramCollection.m_hCompleteness = new TH1F("hCompleteness_CRL", "hCompleteness_CRL", 100, 0., 1.1);
        deltaRayRecoHistogramCollection.m_hCompleteness->SetTitle(";Completeness;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hLowestCompletenessView)
    {
        deltaRayRecoHistogramCollection.m_hLowestCompletenessView = new TH1F("hLowestCompletenessView_CRL", "hLowestCompletenessView_CRL", 100, -0.1, 1.1);
        deltaRayRecoHistogramCollection.m_hLowestCompletenessView->SetTitle(";Completeness;Occurance");
    }
    
    if (!deltaRayRecoHistogramCollection.m_hPurity)
    {
        deltaRayRecoHistogramCollection.m_hPurity = new TH1F("hPurity_CRL", "hPurity_CRL", 100, 0., 1.1);
        deltaRayRecoHistogramCollection.m_hPurity->SetTitle(";Purity;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hCompletenessVsHits)
    {
         deltaRayRecoHistogramCollection.m_hCompletenessVsHits = new TH2F("hCompletenessVsHits_CRL", "hCompletenessVsHits_CRL", 100, 0., 1.1, 50, 0., 200.);
         deltaRayRecoHistogramCollection.m_hCompletenessVsHits->SetTitle(";Completeness;Hits");
    }

    if (!deltaRayRecoHistogramCollection.m_hAboveThresholdMatches)
    {
        deltaRayRecoHistogramCollection.m_hAboveThresholdMatches = new TH1F("hAboveThresholdMatches_CRL", "hAboveThresholdMatches_CRL", 40000, 0., 4.);
        deltaRayRecoHistogramCollection.m_hAboveThresholdMatches->SetTitle(";nAboveThresholdMatches;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal)
    {
        deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal = new TH1F("hParentTrackHitsTotal_CRL", "hParentTrackHitsTotal_CRL", 40000, 0., 50.);
        deltaRayRecoHistogramCollection.m_hParentTrackHitsTotal->SetTitle(";nParentTrackHitsTotal;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal)
    {
        deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal = new TH1F("hOtherTrackHitsTotal_CRL", "hOtherTrackHitsTotal_CRL", 40000, 0., 20.);
        deltaRayRecoHistogramCollection.m_hOtherTrackHitsTotal->SetTitle(";nOtherTrackHitsTotal;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal)
    {
        deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal = new TH1F("hOtherShowerHitsTotal_CRL", "hOtherShowerHitsTotal_CRL", 40000, 0., 20.);
        deltaRayRecoHistogramCollection.m_hOtherShowerHitsTotal->SetTitle(";nOtherShowerHitsTotal;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay)
    {
        deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay = new TH1F("hTotalHitsTakenByCosmicRay_CRL", "hTotalHitsTakenByCosmicRay_CRL", 40000, 0., 50.);
        deltaRayRecoHistogramCollection.m_hTotalHitsTakenByCosmicRay->SetTitle(";hTotalHitsTakenByCosmicRay;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hEfficiency_Energy)
    {
        deltaRayRecoHistogramCollection.m_hEfficiency_Energy = new TH1F("hEfficiency_Energy_CRL", "hEfficiency_Energy_CRL", 100, 0., 1.);
        deltaRayRecoHistogramCollection.m_hEfficiency_Energy->SetTitle(";TrueCRLEnergy [GeV];Efficiency");
    }

    if (!deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits)
    {
        deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits = new TH1F("hEfficiency_TotalHits_CRL", "hEfficiency_TotalHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hEfficiency_TotalHits->SetTitle(";nMCTotalHits;Efficiency");
    }

    if (!deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle)
    {
        deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle = new TH1F("hEfficiency_OpeningAngle_CRL", "_CRL", 100, 0., 180.);
        deltaRayRecoHistogramCollection.m_hEfficiency_OpeningAngle->SetTitle(";TrueOpeningAngle [degrees];Efficiency");
    }

    if (!deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0XZ)
    {
        deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0XZ = new TH1F("hEfficiency_ParentMuonTheta0XZ_CRL", "_CRL", 100, -180., 180.);
        deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0XZ->SetTitle(";ParentMuonTheta0XZ [degrees];Efficiency");
    }

    if (!deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0YZ)
    {
        deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0YZ = new TH1F("hEfficiency_ParentMuonTheta0YZ_CRL", "_CRL", 200, -180., 180.);
        deltaRayRecoHistogramCollection.m_hEfficiency_ParentMuonTheta0YZ->SetTitle(";ParentMuonTheta0YZ [degrees];Efficiency");
    }

    if (!deltaRayRecoHistogramCollection.m_hEfficiency_LowestViewNHits)
    {
        deltaRayRecoHistogramCollection.m_hEfficiency_LowestViewNHits = new TH1F("hEfficiency_LowestViewNHits_CRL", "_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hEfficiency_LowestViewNHits->SetTitle("; nMCHits;Efficiency");
    }     

    if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy)
    {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy = new TH1F("hCorrectParentLink_Energy_CRL", "hCorrectParentLink_Energy_CRL", 100, 0., 1.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_Energy->SetTitle(";TrueCRLEnergy [GeV];CorrectParentLinkFraction");
    }

    if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits)
    {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits = new TH1F("hCorrectParentLink_TotalHits_CRL", "hCorrectParentLink_TotalHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_TotalHits->SetTitle(";nMCTotalHits;CorrectParentLinkFraction");
    }

    if (!deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle)
    {
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle = new TH1F("hCorrectParentLink_OpeningAngle_CRL", "hCorrectParentLink_OpeningAngle_CRL", 100, 0., 180.);
        deltaRayRecoHistogramCollection.m_hCorrectParentLink_OpeningAngle->SetTitle(";TrueOpeningAngle [degrees];CorrectParentLinkFraction");
    }      

    if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy)
    {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy = new TH1F("hCorrectEvent_Energy_CRL", "hCorrectEvent_Energy_CRL", 100, 0., 1.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_Energy->SetTitle(";True Delta Ray Energy [GeV];Correct Reconstruction Fraction");
    }

    if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits)
    {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits = new TH1F("hCorrectEvent_TotalHits_CRL", "hCorrectEvent_TotalHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_TotalHits->SetTitle(";nMCTotalHits;Correct Reconstruction Fraction");
    }

    if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle)
    {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle = new TH1F("hCorrectEvent_OpeningAngle_CRL", "hCorrectEvent_OpeningAngle_CRL", 100, 0., 180.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_OpeningAngle->SetTitle(";True Opening Angle [degrees];Correct Reconstruction Fraction");
    }

    if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0XZ)
    {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0XZ = new TH1F("hCorrectEvent_ParentMuonTheta0XZ_CRL", "hCorrectEvent_ParentMuonTheta0XZ_CRL", 100, -180., 180.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0XZ->SetTitle(";ParentMuonTheta0XZ [degrees];Correct Reconstruction Fraction");
    }

    if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0YZ)
    {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0YZ = new TH1F("hCorrectEvent_ParentMuonTheta0YZ_CRL", "hCorrectEvent_ParentMuonTheta0YZ_CRL", 200, -180., 180.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_ParentMuonTheta0YZ->SetTitle(";ParentMuonTheta0YZ [degrees];Correct Reconstruction Fraction");
    }

    if (!deltaRayRecoHistogramCollection.m_hCorrectEvent_LowestViewNHits)
    {
        deltaRayRecoHistogramCollection.m_hCorrectEvent_LowestViewNHits = new TH1F("hCorrectEvent_LowestViewNHits_CRL", "hCorrectEvent_LowestViewNHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hCorrectEvent_LowestViewNHits->SetTitle("; nMCHits;Correct Reconstruction Fraction");
    }         

    if (!deltaRayRecoHistogramCollection.m_hMatches0_TotalHits)
    {
        deltaRayRecoHistogramCollection.m_hMatches0_TotalHits = new TH1F("hMatches0_TotalHits_CRL", "hMatches0_TotalHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hMatches0_TotalHits->SetTitle(";TotalHits;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches0_OpeningAngle)
    {
        deltaRayRecoHistogramCollection.m_hMatches0_OpeningAngle = new TH1F("hMatches0_OpeningAngle_CRL", "hMatches0_OpeningAngle_CRL", 100, 0., 180.);
        deltaRayRecoHistogramCollection.m_hMatches0_OpeningAngle->SetTitle(";OpeningAngle [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0XZ)
    {
        deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0XZ = new TH1F("hMatches0_ParentMuonTheta0XZ_CRL", "hMatches0_ParentMuonTheta0XZ_CRL", 100, -180., 180.);
        deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0XZ->SetTitle(";ParentMuonTheta0XZ [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0YZ)
    {
        deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0YZ = new TH1F("hMatches0_ParentMuonTheta0YZ_CRL", "hMatches0_ParentMuonTheta0YZ_CRL", 200, -180., 180.);
        deltaRayRecoHistogramCollection.m_hMatches0_ParentMuonTheta0YZ->SetTitle(";ParentMuonTheta0YZ [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches0_LowestViewNHits)
    {
        deltaRayRecoHistogramCollection.m_hMatches0_LowestViewNHits = new TH1F("hMatches0_LowestViewNHits_CRL", "hMatches0_LowestViewNHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hMatches0_LowestViewNHits->SetTitle("; LowestViewNHits;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches1_TotalHits)
    {
        deltaRayRecoHistogramCollection.m_hMatches1_TotalHits = new TH1F("hMatches1_TotalHits_CRL", "hMatches1_TotalHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hMatches1_TotalHits->SetTitle(";TotalHits;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches1_OpeningAngle)
    {
        deltaRayRecoHistogramCollection.m_hMatches1_OpeningAngle = new TH1F("hMatches1_OpeningAngle_CRL", "hMatches1_OpeningAngle_CRL", 100, 0., 180.);
        deltaRayRecoHistogramCollection.m_hMatches1_OpeningAngle->SetTitle(";OpeningAngle [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0XZ)
    {
        deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0XZ = new TH1F("hMatches1_ParentMuonTheta0XZ_CRL", "hMatches1_ParentMuonTheta0XZ_CRL", 100, -180., 180.);
        deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0XZ->SetTitle(";ParentMuonTheta0XZ [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0YZ)
    {
        deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0YZ = new TH1F("hMatches1_ParentMuonTheta0YZ_CRL", "hMatches1_ParentMuonTheta0YZ_CRL", 200, -180., 180.);
        deltaRayRecoHistogramCollection.m_hMatches1_ParentMuonTheta0YZ->SetTitle(";ParentMuonTheta0YZ [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatches1_LowestViewNHits)
    {
        deltaRayRecoHistogramCollection.m_hMatches1_LowestViewNHits = new TH1F("hMatches1_LowestViewNHits_CRL", "hMatches1_LowestViewNHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hMatches1_LowestViewNHits->SetTitle("; LowestViewNHits;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatchesMultiple_TotalHits)
    {
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_TotalHits = new TH1F("hMatchesMultiple_TotalHits_CRL", "hMatchesMultiple_TotalHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_TotalHits->SetTitle(";TotalHits;Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatchesMultiple_OpeningAngle)
    {
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_OpeningAngle = new TH1F("hMatchesMultiple_OpeningAngle_CRL", "hMatchesMultiple_OpeningAngle_CRL", 100, 0., 180.);
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_OpeningAngle->SetTitle(";OpeningAngle [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0XZ)
    {
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0XZ = new TH1F("hMatchesMultiple_ParentMuonTheta0XZ_CRL", "hMatchesMultiple_ParentMuonTheta0XZ_CRL", 100, -180., 180.);
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0XZ->SetTitle(";ParentMuonTheta0XZ [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0YZ)
    {
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0YZ = new TH1F("hMatchesMultiple_ParentMuonTheta0YZ_CRL", "hMatchesMultiple_ParentMuonTheta0YZ_CRL", 200, -180., 180.);
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_ParentMuonTheta0YZ->SetTitle(";ParentMuonTheta0YZ [degrees];Occurance");
    }

    if (!deltaRayRecoHistogramCollection.m_hMatchesMultiple_LowestViewNHits)
    {
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_LowestViewNHits = new TH1F("hMatchesMultiple_LowestViewNHits_CRL", "hMatchesMultiple_LowestViewNHits_CRL", 100, 0., 100.);
        deltaRayRecoHistogramCollection.m_hMatchesMultiple_LowestViewNHits->SetTitle("; LowestViewNHits;Occurance");
    }      

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
