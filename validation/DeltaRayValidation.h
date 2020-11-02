/**
 *  @file   LArReco/validation/DeltaRayValidation.h
 *
 *  @brief  Header file for delta ray validation functionality
 *
 *  $Log: $
 */
#ifndef DELTA_RAY_VALIDATION_H
#define DELTA_RAY_VALIDATION_H 1

typedef std::vector<int> IntVector;
typedef std::vector<float> FloatVector;

/**
 * @brief   Parameters class
 */
class Parameters
{
public:
    /**
     *  @brief  Default constructor
     */
    Parameters();

    bool                    m_printOverallRecoMetrics;
    bool                    m_histogramOutput;
    std::string             m_histogramFileName;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SimpleThreeVector class
 */
class SimpleThreeVector
{
public:
    /**
     *  @brief  Default constructor
     */
    SimpleThreeVector();

    /**
     *  @brief  Constructor
     *
     *  @param  x the x value
     *  @param  y the y value
     *  @param  z the z value
     */
    SimpleThreeVector(const float x, const float y, const float z); 

    float GetMagnitude() const; 
    
    /**
     *  @brief  Calculate the opening angle of the vector with another vector
     *
     *  @param  otherVector the other vector
     *
     *  @return  the opening angle in degrees
     */
    float GetOpeningAngle(const SimpleThreeVector &otherVector); 

    
    float               m_x;                            ///< The x value
    float               m_y;                            ///< The y value
    float               m_z;                            ///< The z value
};

typedef std::vector<SimpleThreeVector> SimpleThreeVectorList;

/**
 *  @brief  Simple three vector subtraction operator
 *
 *  @param  lhs first vector, from which the second is subtracted
 *  @param  rhs second vector, which is subtracted from the first
 */
SimpleThreeVector operator-(const SimpleThreeVector &lhs, const SimpleThreeVector &rhs);

/**
 *  @brief  Simple three vector addition operator
 *
 *  @param  lhs first vector, from which the second is added
 *  @param  rhs second vector, which is added to the first
 */
SimpleThreeVector operator+(const SimpleThreeVector &lhs, const SimpleThreeVector &rhs);

//------------------------------------------------------------------------------------------------------------------------------------------

/*
 *  @brief CosmicRay class
 */
class CosmicRay
{
public:
    /**
     *  @brief  Constructor
     */
    CosmicRay();

    void Print();
    
    float               m_energy;                       ///< The energy
    SimpleThreeVector   m_momentum;                     ///< The momentum
    int                 m_nMCHitsTotal;                 ///< The total number of mc hits
    int                 m_nMCHitsU;                     ///< The number of u mc hits
    int                 m_nMCHitsV;                     ///< The number of v mc hits
    int                 m_nMCHitsW;                     ///< The number of w mc hits
    int                 m_nReconstructableChildCRLs;     ///< The number of reconstructable child delta rays 
    
    int                 m_nCorrectChildCRLs;             ///< The number of correctly reconstructed child delta rays
};

typedef std::vector<CosmicRay> CosmicRayVector;

   
//------------------------------------------------------------------------------------------------------------------------------------------

/*
 *  @brief DeltaRay class
 */
class DeltaRay
{
public:
    /**
     *  @brief  Constructor
     */
    DeltaRay();

    bool operator== (const DeltaRay &deltaRay)
    {
        return (std::fabs(this->m_energy - deltaRay.m_energy) < std::numeric_limits<float>::epsilon() &&
                std::fabs(this->m_momentum.m_x - deltaRay.m_momentum.m_x) < std::numeric_limits<float>::epsilon());
    }
 

    void Print();

    float               m_energy;                       ///< The energy
    SimpleThreeVector   m_momentum;                     ///< The momentum
    int                 m_nMCHitsTotal;                 ///< The total number of mc hits
    int                 m_nMCHitsU;                     ///< The number of u mc hits
    int                 m_nMCHitsV;                     ///< The number of v mc hits
    int                 m_nMCHitsW;                     ///< The number of w mc hits
    float               m_openingAngleFromMuon;         ///< The opening angle wrt MC muon momentum vector
    
    int m_nAboveThresholdMatches;
    int m_isCorrect;
    int m_isCorrectParentLink;
    int m_isBestMatchedCorrectParentLink;

    int m_bestMatchNHitsTotal;
    int m_bestMatchNHitsU;
    int m_bestMatchNHitsV;
    int m_bestMatchNHitsW;

    int m_bestMatchNSharedHitsTotal;
    int m_bestMatchNSharedHitsU;
    int m_bestMatchNSharedHitsV;
    int m_bestMatchNSharedHitsW;
    
    int m_bestMatchNParentTrackHitsTotal;
    int m_bestMatchNParentTrackHitsU;
    int m_bestMatchNParentTrackHitsV;
    int m_bestMatchNParentTrackHitsW;

    int m_bestMatchNOtherTrackHitsTotal;
    int m_bestMatchNOtherTrackHitsU;
    int m_bestMatchNOtherTrackHitsV;
    int m_bestMatchNOtherTrackHitsW;

    int m_bestMatchNOtherShowerHitsTotal;
    int m_bestMatchNOtherShowerHitsU;
    int m_bestMatchNOtherShowerHitsV;
    int m_bestMatchNOtherShowerHitsW;

    int m_totalCRLHitsInBestMatchParentCR;
    int m_uCRLHitsInBestMatchParentCR;
    int m_vCRLHitsInBestMatchParentCR;
    int m_wCRLHitsInBestMatchParentCR;

    FloatVector m_bestMatchOtherShowerHitsDistance;
    FloatVector m_bestMatchOtherTrackHitsDistance;
    FloatVector m_bestMatchParentTrackHitsDistance;
    FloatVector m_bestMatchCRLHitsInCRDistance;  
};

typedef std::vector<DeltaRay> DeltaRayVector;
typedef std::vector<DeltaRay*> DeltaRayPointerVector;

typedef std::map<int, std::map<int, DeltaRayVector>> CosmicRayOwnershipMap;

//------------------------------------------------------------------------------------------------------------------------------------------

class TH1F;

/**
 *  @brief  CosmicRayMCHistogramCollection class
 */
class CosmicRayMCHistogramCollection
{
public:
    /**
     *  @brief  Default constructor
     */
    CosmicRayMCHistogramCollection();

    TH1F                   *m_hReconstructableChildDeltaRays;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DeltaRayMCHistogramCollection class
 */
class DeltaRayMCHistogramCollection
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayMCHistogramCollection();

    TH1F                   *m_hEnergyDistribution;
    TH1F                   *m_hTotalHitDistribution;       
    TH1F                   *m_hOpeningAngleDistribution;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DeltaRayRecoHistogramCollection class
 */
class DeltaRayRecoHistogramCollection
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayRecoHistogramCollection();

    TH1F                   *m_hCompleteness;
    TH1F                   *m_hPurity;

    TH1F                   *m_hAboveThresholdMatches;

    TH1F                   *m_hParentTrackHitsTotal;
    TH1F                   *m_hOtherTrackHitsTotal;
    TH1F                   *m_hOtherShowerHitsTotal;
    TH1F                   *m_hTotalHitsTakenByCosmicRay;

    TH1F                   *m_hEfficiency_Energy;
    TH1F                   *m_hEfficiency_TotalHits;
    TH1F                   *m_hEfficiency_OpeningAngle;  
    
    TH1F                   *m_hCorrectParentLink_Energy;
    TH1F                   *m_hCorrectParentLink_TotalHits;    
    TH1F                   *m_hCorrectParentLink_OpeningAngle;

    TH1F                   *m_hCorrectEvent_Energy;
    TH1F                   *m_hCorrectEvent_TotalHits;
    TH1F                   *m_hCorrectEvent_OpeningAngle;  
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DeltaRayRecoHistogramCollection class
 */
class DeltaRayContaminationHistogramCollection
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayContaminationHistogramCollection();

    TH1F                   *m_hOtherShowerHitsDistance;
    TH1F                   *m_hOtherTrackHitsDistance;
    TH1F                   *m_hParentTrackHitsDistance;
    TH1F                   *m_hCRLHitsInCRDistance;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Delta Ray Validation - Main entry point for analysis
 *
 *  @param  inputFiles the regex identifying the input root files
 *  @param  parameters the parameters
 */
void DeltaRayValidation(const std::string &inputFileName, const Parameters &parameters);

void ReadTree(const std::string &inputFileName, CosmicRayVector &cosmicRayVector, DeltaRayVector &deltaRayVector);

void DisplayOverallRecoMetrics(const DeltaRayVector &deltaRayVector);

void FillCosmicRayMCHistogramCollection(const CosmicRayVector &cosmicRayVector, CosmicRayMCHistogramCollection &cosmicRayMCHistogramCollection);

void FillDeltaRayMCHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection);

void FillDeltaRayRecoHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection);

void ProcessHistograms(DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection, DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection);

void FillDeltaRayContaminationHistogramCollection(const DeltaRayVector &deltaRayVector, DeltaRayContaminationHistogramCollection &deltaRayContaminationHistogramCollection);

void DivideHistogram(TH1F *&recoHistogram, TH1F *&mcDistribution);

void WriteHistograms(CosmicRayMCHistogramCollection &cosmicRayMCHistogramCollection, DeltaRayMCHistogramCollection &deltaRayMCHistogramCollection,
    DeltaRayRecoHistogramCollection &deltaRayRecoHistogramCollection, DeltaRayContaminationHistogramCollection &deltaRayContaminationHistogramCollection);

//------------------------------------------------------------------------------------------------------------------------------------------

Parameters::Parameters() :
    m_printOverallRecoMetrics(true),
    m_histogramOutput(false),
    m_histogramFileName("")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

SimpleThreeVector::SimpleThreeVector() :
    m_x(-std::numeric_limits<float>::max()),
    m_y(-std::numeric_limits<float>::max()),
    m_z(-std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SimpleThreeVector::SimpleThreeVector(const float x, const float y, const float z) :
    m_x(x),
    m_y(y),
    m_z(z)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SimpleThreeVector operator-(const SimpleThreeVector &lhs, const SimpleThreeVector &rhs)
{
    return SimpleThreeVector(lhs.m_x - rhs.m_x, lhs.m_y - rhs.m_y, lhs.m_z - rhs.m_z);
}

//------------------------------------------------------------------------------------------------------------------------------------------

SimpleThreeVector operator+(const SimpleThreeVector &lhs, const SimpleThreeVector &rhs)
{
    return SimpleThreeVector(lhs.m_x + rhs.m_x, lhs.m_y + rhs.m_y, lhs.m_z + rhs.m_z);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRay::CosmicRay() :
    m_energy(0.f),
    m_momentum(0.f, 0.f, 0.f),
    m_nMCHitsTotal(0),
    m_nMCHitsU(0),
    m_nMCHitsV(0),
    m_nMCHitsW(0),
    m_nReconstructableChildCRLs(0),
    m_nCorrectChildCRLs(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DeltaRay::DeltaRay() :
    m_energy(0.f),
    m_momentum(0.f, 0.f, 0.f),
    m_nMCHitsTotal(0),
    m_nMCHitsU(0),
    m_nMCHitsV(0),
    m_nMCHitsW(0),
    m_openingAngleFromMuon(0.f),
    m_nAboveThresholdMatches(0),
    m_isCorrect(false),
    m_isCorrectParentLink(false),
    m_bestMatchNHitsTotal(0),
    m_bestMatchNHitsU(0),
    m_bestMatchNHitsV(0),
    m_bestMatchNHitsW(0),
    m_bestMatchNSharedHitsTotal(0),
    m_bestMatchNSharedHitsU(0),
    m_bestMatchNSharedHitsV(0),
    m_bestMatchNSharedHitsW(0),
    m_bestMatchNParentTrackHitsTotal(0),
    m_bestMatchNParentTrackHitsU(0),
    m_bestMatchNParentTrackHitsV(0),
    m_bestMatchNParentTrackHitsW(0),
    m_bestMatchNOtherTrackHitsTotal(0),
    m_bestMatchNOtherTrackHitsU(0),
    m_bestMatchNOtherTrackHitsV(0),
    m_bestMatchNOtherTrackHitsW(0),
    m_bestMatchNOtherShowerHitsTotal(0),
    m_bestMatchNOtherShowerHitsU(0),
    m_bestMatchNOtherShowerHitsV(0),
    m_bestMatchNOtherShowerHitsW(0),
    m_totalCRLHitsInBestMatchParentCR(0),
    m_uCRLHitsInBestMatchParentCR(0),
    m_vCRLHitsInBestMatchParentCR(0),
    m_wCRLHitsInBestMatchParentCR(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayMCHistogramCollection::CosmicRayMCHistogramCollection() :
    m_hReconstructableChildDeltaRays(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DeltaRayMCHistogramCollection::DeltaRayMCHistogramCollection() :
    m_hEnergyDistribution(nullptr),
    m_hTotalHitDistribution(nullptr),
    m_hOpeningAngleDistribution(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DeltaRayRecoHistogramCollection::DeltaRayRecoHistogramCollection() : 
    m_hCompleteness(nullptr),
    m_hPurity(nullptr),
    m_hAboveThresholdMatches(nullptr),
    m_hParentTrackHitsTotal(nullptr),
    m_hOtherTrackHitsTotal(nullptr),
    m_hOtherShowerHitsTotal(nullptr),
    m_hTotalHitsTakenByCosmicRay(nullptr),
    m_hEfficiency_Energy(nullptr),
    m_hEfficiency_TotalHits(nullptr),
    m_hEfficiency_OpeningAngle(nullptr),
    m_hCorrectParentLink_Energy(nullptr),
    m_hCorrectParentLink_TotalHits(nullptr),
    m_hCorrectParentLink_OpeningAngle(nullptr),
    m_hCorrectEvent_Energy(nullptr),
    m_hCorrectEvent_TotalHits(nullptr),
    m_hCorrectEvent_OpeningAngle(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DeltaRayContaminationHistogramCollection::DeltaRayContaminationHistogramCollection() : 
    m_hOtherShowerHitsDistance(nullptr),
    m_hOtherTrackHitsDistance(nullptr),
    m_hParentTrackHitsDistance(nullptr),
    m_hCRLHitsInCRDistance(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline void CosmicRay::Print()
{
    std::cout << "energy: " << this->m_energy << std::endl;
    std::cout << "momentum: " << "(" << this->m_momentum.m_x << ", " << this->m_momentum.m_y << ", " << this->m_momentum.m_z << ")" << std::endl;
    std::cout << "nMCHits: " << this->m_nMCHitsTotal << " (" << this->m_nMCHitsU << ", " << this->m_nMCHitsV << ", " << this->m_nMCHitsW << ")" << std::endl;
    std::cout << "nReconstructableChildCRLs: " << this->m_nReconstructableChildCRLs << std::endl;
    std::cout << "nCorrectChildCRLs: " << this->m_nCorrectChildCRLs << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DeltaRay::Print()
{
    std::cout << "energy" << this->m_energy << std::endl;
    std::cout << "momentum: " << "(" << this->m_momentum.m_x << ", " << this->m_momentum.m_y << ", " << this->m_momentum.m_z << ")" << std::endl;
    std::cout << "nMCHits: " << this->m_nMCHitsTotal << " (" << this->m_nMCHitsU << ", " << this->m_nMCHitsV << ", " << this->m_nMCHitsW << ")" << std::endl;
    std::cout << "nAboveThresholdMatches: " << this->m_nAboveThresholdMatches << std::endl;
    std::cout << "isCorrect: " << this->m_isCorrect << std::endl;
    std::cout << "isCorrectParentLink: " << this->m_isCorrectParentLink << std::endl;

    std::cout << std::endl;

    std::cout << "bestMatchNHitsTotal: " << this->m_bestMatchNHitsTotal << " (" << this->m_bestMatchNHitsU << ", " << this->m_bestMatchNHitsV << ", " << this->m_bestMatchNHitsW << ") " << std::endl;
    std::cout << "bestMatchNSharedHitsTotal: " << this->m_bestMatchNSharedHitsTotal << " (" << this->m_bestMatchNSharedHitsU << ", " << this->m_bestMatchNSharedHitsV << ", "
              << this->m_bestMatchNSharedHitsW << std::endl;
    std::cout << "bestMatchNParentTrackHitsTotal: " << this->m_bestMatchNParentTrackHitsTotal << " (" << this->m_bestMatchNParentTrackHitsU << ", " << this->m_bestMatchNParentTrackHitsV
              << ", " << m_bestMatchNParentTrackHitsW << ")" << std::endl;
    for (unsigned long i = 0; i < m_bestMatchParentTrackHitsDistance.size(); ++i)
    {
        std::cout << m_bestMatchParentTrackHitsDistance.at(i) << (i == m_bestMatchParentTrackHitsDistance.size() - 1 ? "" : ", ");
    }
    std::cout << std::endl;
    std::cout << "bestMatchNOtherTrackHitsTotal: " << this->m_bestMatchNOtherTrackHitsTotal << " (" << this->m_bestMatchNOtherTrackHitsU << ", " << this->m_bestMatchNOtherTrackHitsV
              << ", " << this->m_bestMatchNOtherTrackHitsW << ")" << std::endl;
    for (unsigned long i = 0; i < m_bestMatchOtherTrackHitsDistance.size(); ++i)
    {
        std::cout << m_bestMatchOtherTrackHitsDistance.at(i) << (i == m_bestMatchOtherTrackHitsDistance.size() - 1 ? "" : ", ");
    }
    std::cout << std::endl;    
    std::cout << "bestMatchNOtherShowerHitsTotal: " << this->m_bestMatchNOtherShowerHitsTotal << " (" << this->m_bestMatchNOtherShowerHitsU << ", " << this->m_bestMatchNOtherShowerHitsV
              << ", " << this->m_bestMatchNOtherShowerHitsW << ")" << std::endl;
    for (unsigned long i = 0; i < m_bestMatchOtherShowerHitsDistance.size(); ++i)
    {
        std::cout << m_bestMatchOtherShowerHitsDistance.at(i) << (i == m_bestMatchOtherShowerHitsDistance.size() - 1 ? "" : ", ");
    }
    std::cout << std::endl;     
    std::cout << "totalCRLHitsInBestMatchParentCR: " << this->m_totalCRLHitsInBestMatchParentCR << " (" << this->m_uCRLHitsInBestMatchParentCR << ", " << this->m_vCRLHitsInBestMatchParentCR
              << ", " << this->m_wCRLHitsInBestMatchParentCR << std::endl;
    for (unsigned long i = 0; i < m_bestMatchCRLHitsInCRDistance.size(); ++i)
    {
        std::cout << m_bestMatchCRLHitsInCRDistance.at(i) << (i == m_bestMatchCRLHitsInCRDistance.size() - 1 ? "" : ", ");
    }
    std::cout << std::endl; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

#endif // #ifndef DELTA_RAY_VALIDATION_H
