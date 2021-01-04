#include "NRMConstant_v4.hpp"
#include <cmath>

NRMConstant_v4::NRMConstant_v4() : previousFiringTime(0.0), exponential(1) {
#ifdef PROFILE
    callsToRebuild=0;
    rebuildCost=0.0;
#endif

}


void NRMConstant_v4::printHashTable() {
    std::cout << "minBin=" << minBin << std::endl;
    std::cout << "lower bound="<<currentLowerBound << std::endl;
    std::cout << "upper bound="<<currentUpperBound << std::endl;
    std::cout << "approx bin width=" << (currentUpperBound-currentLowerBound)/numberOfBins << std::endl;
    for (std::size_t i=0; i!=theHashTable.size(); ++i) {
        std::cout << "[i]=" << i << ": (size=" << theHashTable[i].size() << ")\n";
        for (std::size_t j=0; j!=theHashTable[i].size(); ++j) {
            std::cout << theHashTable[i][j].second << "(" << theHashTable[i][j].first << "), ";
        }
        std::cout << "\n";
    }
}

NRMConstant_v4::timeIndexPair NRMConstant_v4::selectReaction() {
//std::cout << "in selectReaction...\n";
    //while front of queue is empty, pop it
#ifdef PROFILE
    std::chrono::time_point<std::chrono::system_clock> selectTimeStart, selectTimeEnd;
    std::chrono::duration<double> elapsed;
    selectTimeStart=std::chrono::system_clock::now();
    std::size_t emptySearchBinCount=1;//start counting at 1
#endif
    while (theHashTable[minBin].size()==0) {
        ++minBin;
#ifdef PROFILE
        ++emptySearchBinCount;
#endif
        if (minBin==theHashTable.size()) {
//            std::cout << "rebuilding...\n";
            if (!rebuild()) {
                return std::make_pair(std::numeric_limits<double>::max(),-1);
            }
        }
    }
    //now that we have a bin with elements, find minimum
    int minTimeRxnIndex=0;
    for (std::size_t i=1;i<theHashTable[minBin].size(); ++i) {
        if (theHashTable[minBin][i].first<theHashTable[minBin][minTimeRxnIndex].first) {
            minTimeRxnIndex=i;
        }
    }

#ifdef PROFILE
//    std::cout << "pushing emptySearchBinCount=" << emptySearchBinCount << std::endl;
//    std::cout << "pushing searchDepth=" << theHashTable[minBin].size() << std::endl;
    searchBins.push_back(emptySearchBinCount);
    searchDepth.push_back(theHashTable[minBin].size());
#endif

    
    previousFiringTime=theHashTable[minBin][minTimeRxnIndex].first;
    ++rxnCountThisBuildOrRebuild;
//    std::cout << "returning time=" <<theHashTable[minBin][minTimeRxnIndex].first<<", reaction index=" << theHashTable[minBin][minTimeRxnIndex].second << " in selectReaction\n";
#ifdef PROFILE
    selectTimeEnd=std::chrono::system_clock::now();
    elapsed=selectTimeEnd-selectTimeStart;
    selectCost+=elapsed.count();
#endif

    return theHashTable[minBin][minTimeRxnIndex];
}

NRMConstant_v4::~NRMConstant_v4() {
}

int NRMConstant_v4::computeBinIndex(double firingTime) {
    if (firingTime>currentUpperBound) {
//        std::cout << "computeBinIndex("<<firingTime<<") returning binIndex=" << -1 << std::endl;
        return -1;
    }
    int binIndex=(int)((firingTime-currentLowerBound)/(currentUpperBound-currentLowerBound)*numberOfBins);
    int returnVal=std::max(binIndex,(int)minBin);
//    std::cout << "computeBinIndex("<<firingTime<<") returning binIndex=" << binIndex << std::endl;
    return returnVal;
}

//returns false if all propensities are 0
bool NRMConstant_v4::rebuild() {
//    std::cout << "in rebuild...\n";

#ifdef PROFILE
    std::chrono::time_point<std::chrono::system_clock> rebuildTimeStart, rebuildTimeEnd;
    std::chrono::duration<double> elapsed;
    rebuildTimeStart=std::chrono::system_clock::now();
#endif

    //estimate propensitySum based on number of steps since last build or rebuild
    double propensitySumEstimate;
    if (rxnCountThisBuildOrRebuild>0) {
//        std::cout << "fired " << rxnCountThisBuildOrRebuild << " since last build or rebuild\n";
        propensitySumEstimate=(double)rxnCountThisBuildOrRebuild/(currentUpperBound-currentLowerBound);
//        std::cout << "propensity sum estimate = " << propensitySumEstimate << "\n";
    }
    else {
        //this shouldn't happen...if we got here, 0 reactions fired since last rebuild!
        //will only arise in toy problems
        //but need a strategy...
        double previousBinWidth=(currentUpperBound-currentLowerBound)/(double)numberOfBins;
        propensitySumEstimate=1.0/previousBinWidth;
        std::cout << "WARNING: 0 reactions fired before rebuild.\n";
//        std::cout << "previousBinWidth was " << previousBinWidth <<"\n";
    }
    
    setBinNumberAndBounds(currentUpperBound,propensitySumEstimate);//
    minBin=0;
    
    std::vector<std::pair<double, std::size_t> > emptyVector;

    //theHashTable.clear();
    theHashTable.resize(numberOfBins,emptyVector);
    
    //std::cout << "inserting into hash table...\n";
    bool allPropensitiesZero=true;
    for (std::size_t i=0; i!=nextFiringTime.size(); ++i) {
        //insert into hashTable
        if (!std::isinf(nextFiringTime[i])) {
            allPropensitiesZero=false;
        }
        int bin=computeBinIndex(nextFiringTime[i]);
        if (bin>=0) {
//            std::cout << "inserting into bin " << bin << std::endl;
            theHashTable[bin].push_back(std::make_pair(nextFiringTime[i],i));//place this rxn at back of bin
            binIndexAndPositionInBin[i]=std::make_pair(bin,theHashTable[bin].size()-1);
        }
        else {
            binIndexAndPositionInBin[i]=std::make_pair<int,int>(-1,-1);
        }
    }
    
    //record info from build
    rxnsBetweenBuildOrRebuild.push_back(rxnCountThisBuildOrRebuild);
    lowerBoundsBuildOrRebuild.push_back(currentLowerBound);//keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
    upperBoundsBuildOrRebuild.push_back(currentUpperBound);

    //set rxn counter to 0
    rxnCountThisBuildOrRebuild=0;

#ifdef PROFILE
    ++callsToRebuild;
    rebuildTimeEnd=std::chrono::system_clock::now();
    elapsed=rebuildTimeEnd-rebuildTimeStart;
    rebuildCost+=elapsed.count();
#endif

//    std::cout << "after rebuild, hash table:\n";
//    printHashTable();
    
    return !allPropensitiesZero;
}

//fixed number of bins
//currentUpperBound=
void NRMConstant_v4::setBinNumberAndBounds(double newLowerBound, double propensitySum) {
//    std::cout << "...in setBinNumberAndBounds...newLowerBound=" << newLowerBound << std::endl;
    currentLowerBound=newLowerBound;
    if (currentLowerBound>endTime) {
        std::cerr << "ERROR: calling rebuild when simulation end time exceeded. Terminating.\n";
        exit(1);
    }
    std::size_t NumberOfChannels=nextFiringTime.size();//
    
//    std::cout << "propensitySum is " << propensitySum << ", so step size is " << 1.0/propensitySum << "\n";
    double binWidth=16.0/propensitySum;
//    std::cout << "binWidth=" << binWidth << "\n";
  
    //sanft: temporary change for testing, edit: permanent change for v4
//    std::size_t maxNumberOfBins=(std::size_t)20.0*sqrt(NumberOfChannels);
//    std::size_t maxNumberOfBins=(std::size_t)3.0*sqrt(NumberOfChannels);
    //use active channels instead of total channels
    std::size_t maxNumberOfBins=(std::size_t)20.0*sqrt(activeChannelCounter);
    
//    std::cout << "activeChannels=" << activeChannelCounter << "\n";
//    std::cout << "maxNumberOfBins=" << maxNumberOfBins << "\n";
    
    //how many bins to get beyond end time?
    std::size_t endTimeBins=std::numeric_limits<std::size_t>::max();
    //if we know the end time:
    if (endTime<std::numeric_limits<double>::max()) {
//        std::cout << "endTime is set, how many bins needed to exceed...\n";
        endTimeBins=(std::size_t)ceil( (endTime-currentLowerBound)/binWidth);
//        std::cout << "endTimeBins="<<endTimeBins << "\n";
        //use number of bins that takes us past endTime, or maxNumberOfBins, whichever is smaller
        numberOfBins=std::min(endTimeBins,maxNumberOfBins);
    }
    else {
        numberOfBins=maxNumberOfBins;
    }
    currentUpperBound=currentLowerBound+binWidth*(double)numberOfBins;

    //if max bins gets us very close to the simulation end time, bump it up beyond max to avoid
    //an additional rebuild
    if (endTime<std::numeric_limits<double>::max()) {
        std::cout << "checking to see if we can use fewer than max bins...\n";
        if (currentUpperBound<endTime && ((double)(endTimeBins-maxNumberOfBins)/((double)numberOfBins))<.2) {
            numberOfBins=endTimeBins;
            currentUpperBound=currentLowerBound+binWidth*(double)numberOfBins;
        }
    }
    
//    std::cout << "currentLowerBound=" << currentLowerBound << ", currentUpperBound=" << currentUpperBound << std::endl;
//    std::cout << "binWidth=" << binWidth << std::endl;
//    std::cout << "numberOfBins=" << numberOfBins << std::endl;
}


#ifdef PROFILE
double NRMConstant_v4::getBuildCost() {
    return buildCost;
}

double NRMConstant_v4::getRebuildCost() {
    return rebuildCost;
}

std::size_t NRMConstant_v4::getCallsToRebuild() {
    return callsToRebuild;
}

double NRMConstant_v4::getUpdateCost() {
    return updateCost;
}

double NRMConstant_v4::getSelectCost() {
    return selectCost;
}

std::vector<std::size_t>& NRMConstant_v4::getSearchBinsRef() {
    return searchBins;
}
std::vector<std::size_t>& NRMConstant_v4::getSearchDepthRef() {
    return searchDepth;
}


}
#endif

