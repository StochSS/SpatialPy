#ifndef NRM_CONSTANT_V5_HPP
#define NRM_CONSTANT_V5_HPP

// needs proper error handling/logging
// does NOT handle change in number of reaction channels
// would be nice to update to remove propensity 0 channels from unordered_maps
// (then activeChannelCounter would be unordered_map.size but more logic needed in update (and rebuild?)

// class to select next reaction for Next Reaction Method using hash table

// assumes (best performance when) many (order of thousands) reaction channels,

// a lot could be modified! talk to Kevin Sanft for improvements.
// needs a test suite

//change from v3: use active (nonzero) channels instead of total channels
//change from v4: adapted from Sanft's test code to SpatialPy
// removes vectors and replace with hashmap (unordered_map) to allow adding or removing "channels"
// removed PROFILE code (look at old version for old profiling data)

#include <iostream>
#include <vector>
#include <random>
//#include <queue>
#include <vector>
#include <unordered_map>
#include <utility>
#include <ctime>
#include <chrono>
#include <limits>

using namespace Spatialpy ;
class NRMConstant_v5 {
    typedef double firingTimeT;
    typedef std::size_t reactionIndexT;//assumes channels have unique nonnegative ID (within range of size_t)
    typedef std::pair<firingTimeT,reactionIndexT> timeIndexPair;

    double previousFiringTime;//keep the time of the last reaction that fired (from last call of selectReaction)
    
    // these two could be combined
//    std::vector<firingTimeT> nextFiringTime;//for each reaction channel
//    std::vector<std::pair<int,int> > binIndexAndPositionInBin;//for each reaction channel
    std::unordered_map<reactionIndexT,firingTimeT> nextFiringTime;
    std::unordered_map<reactionIndexT,std::pair<int,int> > binIndexAndPositionInBin;//for each reaction channel
    
    std::exponential_distribution<> exponential;

    std::size_t minBin;//corresponds to bin of latest time step
    std::size_t numberOfBins;

    //range of current hash table
    double currentLowerBound;//corresponds to left endpoint of bin 0 of current hash table
    double currentUpperBound;//corresponds to right endpoint of last bin in current hash table
    
    std::size_t rxnCountThisBuildOrRebuild;
    std::vector<std::size_t> rxnsBetweenBuildOrRebuild;//keep entry for each build/rebuild
    std::vector<double> lowerBoundsBuildOrRebuild;//keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
    std::vector<double> upperBoundsBuildOrRebuild;
    
    std::vector<std::vector<std::pair<double, std::size_t> > > theHashTable;//pair = firing time, reaction index
    
    double endTime;//simulation end time
    
    std::size_t strategyBins;
    double strategyWidth;
    
 public:
	NRMConstant_v5();
	
    // activeChannels is count of nonzero propensities
    std::size_t activeChannelCounter;
    template <typename generatorType>
    void build(std::vector<double>& propensities, generatorType& generator,
               double propensitySum, std::size_t activeChannels,
               double timeOffset=0.0, 
               double simulationEndTime=std::numeric_limits<double>::max()) ;
    
    // activeChannels is count of nonzero propensities
    template <typename generatorType>
    void build(std::vector<Spatialpy::Particle> particles, generatorType& generator, 
               double propensitySum, std::size_t activeChannels, 
               double timeOffset=0.0, 
               double simulationEndTime=std::numeric_limits<double>::max()) ;
    
    void printHashTable();//for debug/dev

    timeIndexPair selectReaction(); // returns pair of firing time, rxn index

    // reactionIndex is particle id
    template <typename generatorType>
    void update(std::size_t reactionIndex, double newPropensity, double currentTime, generatorType& generator) ;
    
	virtual ~NRMConstant_v5();
    
private:
    int computeBinIndex(double firingTime);
    //std::size_t insertIntoBin(std::size_t binIndex, double firingTime, std::size_t reactionIndex);
    bool rebuild();//returns false if all propensities are 0
    void setBinNumberAndBounds(double newLowerBound, double propensitySum, int activeChannels);
    
};

#endif
