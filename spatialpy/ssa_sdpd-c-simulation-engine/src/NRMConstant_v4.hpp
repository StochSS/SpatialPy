#ifndef NRM_CONSTANT_V4_HPP
#define NRM_CONSTANT_V4_HPP

// class to select next reaction for Next Reaction Method using hash table

// assumes (best performance when) many (order of thousands) reaction channels,
// and

//change from v3: use active (nonzero) channels instead of total channels

#include <iostream>
#include <vector>
#include <random>
//#include <queue>
#include <vector>
#include <utility>
#include <ctime>
#include <chrono>
#include <limits>

class NRMConstant_v4 {
    typedef double firingTimeT;
    typedef int reactionIndexT;//negative are error or similar
    typedef std::pair<firingTimeT,reactionIndexT> timeIndexPair;

    double previousFiringTime;//keep the time of the last reaction that fired (from last call of selectReaction)
    
    std::vector<firingTimeT> nextFiringTime;//for each reaction channel
    std::vector<std::pair<int,int> > binIndexAndPositionInBin;//for each reaction channel

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
    
//#if NRM_STRAT==1
//    std::size_t strategyBins;
//    double strategyUpperBound;
//#endif
//#if NRM_STRAT==2
    std::size_t strategyBins;
    double strategyWidth;
//#endif
    
 public:
	NRMConstant_v4();
	
    template <typename generatorType>
    void build(std::vector<double>& propensities, generatorType& generator, double propensitySum, std::size_t activeChannels, double timeOffset=0.0, double simulationEndTime=std::numeric_limits<double>::max()) {
//        std::cout << "in build..." << std::endl;

        activeChannelCounter=activeChannels;
        endTime=simulationEndTime;
#ifdef PROFILE
        std::chrono::time_point<std::chrono::system_clock> buildTimeStart, buildTimeEnd;
        std::chrono::duration<double> elapsed;
        buildTimeStart=std::chrono::system_clock::now();
#endif
        previousFiringTime=timeOffset;//initialize (could set to nan or something instead)
        
        nextFiringTime.resize(propensities.size());
        binIndexAndPositionInBin.resize(propensities.size());
//        double propensitySum=std::accumulate(propensities.cbegin(),propensities.cend(),0.0);

        if (propensitySum==0.0) {
            std::cerr << "ERROR: propensity sum is zero on initial build. Terminating." << std::endl;
            exit(1);
        }
        
        setBinNumberAndBounds(timeOffset,propensitySum);//choose lowerBound, upperBound, and numberOfBins
        minBin=0;

        std::vector<std::pair<double, std::size_t> > emptyVector;
        
        theHashTable.resize(numberOfBins,emptyVector);
        
        for (std::size_t i=0; i!=propensities.size(); ++i) {
            double firingTime;
            if (propensities[i]==0.0) {
                firingTime=std::numeric_limits<double>::infinity();
            }
            else {
                firingTime=exponential(generator)/propensities[i]+timeOffset;
            }
            nextFiringTime[i]=firingTime;
//            std::cout << "nextFiringTime["<<i<<"]="<<nextFiringTime[i]<<"\n";
            //insert into hashTable
            int bin=computeBinIndex(firingTime);
            if (bin>=0) {
                theHashTable[bin].push_back(std::make_pair(firingTime,i));//place this rxn at back of bin
                binIndexAndPositionInBin[i]=std::make_pair(bin,theHashTable[bin].size()-1);
            } else {
                binIndexAndPositionInBin[i]=std::make_pair<int,int>(-1,-1);
            }
        }

        //set rxn counter to 0
        rxnCountThisBuildOrRebuild=0;

        //record info from build
        rxnsBetweenBuildOrRebuild.push_back(0);
        lowerBoundsBuildOrRebuild.push_back(currentLowerBound);//keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
        upperBoundsBuildOrRebuild.push_back(currentUpperBound);
        
#ifdef PROFILE
        buildTimeEnd=std::chrono::system_clock::now();
        elapsed=buildTimeEnd-buildTimeStart;
        buildCost=elapsed.count();
#endif

//        std::cout << "exiting build...\n";
    }
    
    void printHashTable();//for debug/dev

    timeIndexPair selectReaction(); // returns pair of firing time, rxn index

    template <typename generatorType>
    void update(std::size_t reactionIndex, double newPropensity, double currentTime, generatorType& generator) {
//        std::cout << "updating reactionIndex=" << reactionIndex << ", newPropensity=" << newPropensity << ", currentTime=" << currentTime << std::endl;
#ifdef PROFILE
        std::chrono::time_point<std::chrono::system_clock> updateTimeStart, updateTimeEnd;
        std::chrono::duration<double> elapsed;
        updateTimeStart=std::chrono::system_clock::now();
#endif

        double firingTime;
        if (newPropensity==0) {
            firingTime=std::numeric_limits<double>::infinity();
            if (nextFiringTime[reactionIndex]!=std::numeric_limits<double>::infinity()) {
                --activeChannelCounter;
            }
        }
        else {
            firingTime=exponential(generator)/newPropensity+currentTime;
            if (nextFiringTime[reactionIndex]==std::numeric_limits<double>::infinity()) {
                ++activeChannelCounter;
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
        
#ifdef PROFILE
        updateTimeEnd=std::chrono::system_clock::now();
        elapsed=updateTimeEnd-updateTimeStart;
        updateCost+=elapsed.count();
#endif
   
    }
    
	virtual ~NRMConstant_v4();
//#if NRM_STRAT==1
//    void setStrategyVariables(std::size_t strategyBins, double strategyUpperBound);
//#endif
//#if NRM_STRAT==2
//    void setStrategyVariables(std::size_t strategyBins, double strategyWidth);
//#endif
    
private:
    int computeBinIndex(double firingTime);
    //std::size_t insertIntoBin(std::size_t binIndex, double firingTime, std::size_t reactionIndex);
    bool rebuild();//returns false if all propensities are 0
    void setBinNumberAndBounds(double newLowerBound, double propensitySum);
    std::size_t activeChannelCounter;
    
#ifdef PROFILE
    double buildCost;
    double rebuildCost;
    double updateCost;
    double selectCost;
    std::size_t callsToRebuild;
    std::vector<std::size_t> searchBins;//how many bins searched before nonempty for each call to select reaction
    std::vector<std::size_t> searchDepth;//how many elements were in the first nonempty bin (ie the bin where the reaction occurred)
public:
    double getBuildCost();
    double getRebuildCost();
    std::size_t getCallsToRebuild();
    double getUpdateCost();
    double getSelectCost();
    std::vector<std::size_t>& getSearchBinsRef();
    std::vector<std::size_t>& getSearchDepthRef();
    
#endif
};

#endif
