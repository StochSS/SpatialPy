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

#include "particle.hpp"

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
    template <typename generatorType>
    void build(std::vector<double>& propensities, generatorType& generator, double propensitySum, std::size_t activeChannels, double timeOffset=0.0, double simulationEndTime=std::numeric_limits<double>::max()) {
//        std::cout << "in build..." << std::endl;

        activeChannelCounter=activeChannels;
        endTime=simulationEndTime;
        previousFiringTime=timeOffset;//initialize (could set to nan or something instead)
        
//        nextFiringTime.resize(propensities.size());
//        binIndexAndPositionInBin.resize(propensities.size());

        if (propensitySum==0.0) {
            std::cerr << "ERROR: propensity sum is zero on initial build. Terminating." << std::endl;
            exit(1);
        }
        
        setBinNumberAndBounds(timeOffset,propensitySum,activeChannels);//choose lowerBound, upperBound, and numberOfBins
        minBin=0;

        std::vector<std::pair<double, std::size_t> > emptyVector;
        
        theHashTable.resize(numberOfBins,emptyVector);
        
        std::size_t active_channel_count=0; // for verifying build input
        
        for (std::size_t i=0; i!=propensities.size(); ++i) {
            double firingTime;
            if (propensities[i]==0.0) {
                firingTime=std::numeric_limits<double>::infinity();
            }
            else {
                firingTime=exponential(generator)/propensities[i]+timeOffset;
                ++active_channel_count;
            }

            nextFiringTime[i]=firingTime;//.insert(std::make_pair<reactionIndexT,firingTimeT>(i,firingTime));
//            std::cout << "nextFiringTime["<<i<<"]="<<nextFiringTime[i]<<"\n";
            //insert into hashTable
            int bin=computeBinIndex(firingTime);
            if (bin>=0) {
                theHashTable[bin].push_back(std::make_pair(firingTime,i));//place this rxn at back of bin
                binIndexAndPositionInBin[i]=std::make_pair(bin,theHashTable[bin].size()-1);
            } else {
                binIndexAndPositionInBin[i]=std::make_pair<int,int>(-1,-1);// bin (and index within bin) is -1 if not in the hash table
            }
        }

        //set rxn counter to 0
        rxnCountThisBuildOrRebuild=0;

        //record info from build
        rxnsBetweenBuildOrRebuild.push_back(0);
        lowerBoundsBuildOrRebuild.push_back(currentLowerBound);//keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
        upperBoundsBuildOrRebuild.push_back(currentUpperBound);
        
//        printHashTable();
        
        if (active_channel_count!=activeChannelCounter) {
            std::cout << "ERROR: active channel count is inconsistent.\n";
            exit(1);
        }
        
//        std::cout << "exiting build...\n";
    }
    
    // activeChannels is count of nonzero propensities
    template <typename generatorType>
    void build(std::vector<Spatialpy::Particle> particles, generatorType& generator, double propensitySum, std::size_t activeChannels, double timeOffset=0.0, double simulationEndTime=std::numeric_limits<double>::max()) {
        //        std::cout << "in build..." << std::endl;
        
        activeChannelCounter=activeChannels;
        endTime=simulationEndTime;
        previousFiringTime=timeOffset;//initialize (could set to nan or something instead)
        
        //        nextFiringTime.resize(propensities.size());
        //        binIndexAndPositionInBin.resize(propensities.size());
        
        if (propensitySum==0.0) {
            std::cerr << "ERROR: propensity sum is zero on initial build. Terminating." << std::endl;
            exit(1);
        }
        
        setBinNumberAndBounds(timeOffset,propensitySum,activeChannels);//choose lowerBound, upperBound, and numberOfBins
        minBin=0;
        
        std::vector<std::pair<double, std::size_t> > emptyVector;
        
        theHashTable.resize(numberOfBins,emptyVector);
        
        std::size_t active_channel_count=0; // for verifying build input
        
        std::cout << "?? is index into particles vector a valid unique id?\n";
        for (std::size_t i=0; i!=particles.size(); ++i) {
            double firingTime;
            double this_propensity = particles[i].srrate+particles[i].sdrate;
            if (this_propensity==0.0) {
                firingTime=std::numeric_limits<double>::infinity();
            }
            else {
                firingTime=exponential(generator)/this_propensity+timeOffset;
                ++active_channel_count;
            }
            
            nextFiringTime[i]=firingTime;//.insert(std::make_pair<reactionIndexT,firingTimeT>(i,firingTime));
            //            std::cout << "nextFiringTime["<<i<<"]="<<nextFiringTime[i]<<"\n";
            //insert into hashTable
            int bin=computeBinIndex(firingTime);
            if (bin>=0) {
                theHashTable[bin].push_back(std::make_pair(firingTime,i));//place this rxn at back of bin
                binIndexAndPositionInBin[i]=std::make_pair(bin,theHashTable[bin].size()-1);
            } else {
                binIndexAndPositionInBin[i]=std::make_pair<int,int>(-1,-1);// bin (and index within bin) is -1 if not in the hash table
            }
        }
        
        //set rxn counter to 0
        rxnCountThisBuildOrRebuild=0;
        
        //record info from build
        rxnsBetweenBuildOrRebuild.push_back(0);
        lowerBoundsBuildOrRebuild.push_back(currentLowerBound);//keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
        upperBoundsBuildOrRebuild.push_back(currentUpperBound);
        
        //        printHashTable();
        
        if (active_channel_count!=activeChannelCounter) {
            std::cout << "ERROR: active channel count is inconsistent.\n";
            exit(1);
        }
        
        //        std::cout << "exiting build...\n";

    }
    
    void printHashTable();//for debug/dev

    timeIndexPair selectReaction(); // returns pair of firing time, rxn index

    template <typename generatorType>
    void update(std::size_t reactionIndex, double newPropensity, double currentTime, generatorType& generator) {
//        std::cout << "updating reactionIndex=" << reactionIndex << ", newPropensity=" << newPropensity << ", currentTime=" << currentTime << std::endl;

        double firingTime;
        if (newPropensity==0) {
            firingTime=std::numeric_limits<double>::infinity();
            if (nextFiringTime.at(reactionIndex)!=std::numeric_limits<double>::infinity()) {
                activeChannelCounter--;
            }
        }
        else {
            firingTime=exponential(generator)/newPropensity+currentTime;
            if (nextFiringTime.at(reactionIndex)==std::numeric_limits<double>::infinity()) {
                activeChannelCounter++;
            }
        }
//        std::cout << "new firingTime=" << firingTime << std::endl;
        nextFiringTime[reactionIndex]=firingTime;
        int newBin=computeBinIndex(firingTime);
        
//        std::cout << "newBin=" << newBin << std::endl;
        int oldBin=binIndexAndPositionInBin[reactionIndex].first;
//        std::cout << "oldBin=" << oldBin << std::endl;
        int oldPositionInBin=binIndexAndPositionInBin[reactionIndex].second;
        if (newBin!=oldBin) {
            if (oldBin>=0) {
                //remove from old bin
                if (theHashTable[oldBin].size()>1) {
                    
                    //take last element in old bin and place in this element's spot
                    theHashTable[oldBin][oldPositionInBin]=theHashTable[oldBin].back();
                    std::size_t movedElementIndex=theHashTable[oldBin][oldPositionInBin].second;
                    //update old last element's ...within bin index
                    binIndexAndPositionInBin[movedElementIndex].second=oldPositionInBin;
                }
                theHashTable[oldBin].pop_back();
            }
            binIndexAndPositionInBin[reactionIndex].first=newBin;
            if (newBin>=0) {
                theHashTable[newBin].push_back(std::pair<double,std::size_t>(firingTime,reactionIndex));
                binIndexAndPositionInBin[reactionIndex].second=theHashTable[newBin].size()-1;
            }
        }
        else {
            //just update firing time
            if (newBin>=0) {
                theHashTable[newBin][binIndexAndPositionInBin[reactionIndex].second].first=firingTime;
            }
        }
        
    }
    
	virtual ~NRMConstant_v5();
    
private:
    int computeBinIndex(double firingTime);
    //std::size_t insertIntoBin(std::size_t binIndex, double firingTime, std::size_t reactionIndex);
    bool rebuild();//returns false if all propensities are 0
    void setBinNumberAndBounds(double newLowerBound, double propensitySum, int activeChannels);
    std::size_t activeChannelCounter;
    
};

#endif
