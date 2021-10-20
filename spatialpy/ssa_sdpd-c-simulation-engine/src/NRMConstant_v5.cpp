/**
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <cmath>

#include "NRMConstant_v5.hpp"

namespace Spatialpy{

    NRMConstant_v5::NRMConstant_v5() : previousFiringTime(0.0), exponential(1) {}

    bool NRMConstant_v5::build(std::vector<double>& propensities,
               double propensitySum, std::size_t activeChannels,
               std::mt19937_64& rng, double timeOffset,
               double simulationEndTime) {
        printf("NRMConstant_v5::build(");
        printf("propensitySum=%e,",propensitySum);
        printf("activeChannels=%lu,",activeChannels);
        printf("timeOffset=%e,",timeOffset);
        printf("simulationEndTime=%e,",simulationEndTime);
        printf(")\n");

        activeChannelCounter=activeChannels;
        endTime=simulationEndTime;
        previousFiringTime=timeOffset;
        // clear out maps - reset all data
        nextFiringTime.clear();
        binIndexAndPositionInBin.clear();
        if (propensitySum==0.0) {
            return false;
        }
        if(!setBinNumberAndBounds(timeOffset,propensitySum,activeChannels)){
            //choose lowerBound, upperBound, and numberOfBins
            return false;
        }
        minBin=0;
        std::vector<std::pair<double, std::size_t> > emptyVector;
        theHashTable.clear();
        theHashTable.resize(numberOfBins,emptyVector);
        std::size_t active_channel_count=0; //for verifying build input
        int numNonZero=0;

        for (std::size_t i=0; i!=propensities.size(); ++i) {
            double firingTime;
            if (propensities[i]==0.0) {
                firingTime=std::numeric_limits<double>::infinity();
            } else {
                firingTime=exponential(rng)/propensities[i]+timeOffset;
                ++active_channel_count;
            }

            nextFiringTime[i]=firingTime; 
            //insert into hashTable
            int bin=computeBinIndex(firingTime);
            if (bin>=0) {
                theHashTable[bin].push_back(std::make_pair(firingTime,i)); //place this rxn at back of bin
                binIndexAndPositionInBin[i]=std::make_pair(bin,theHashTable[bin].size()-1);
                numNonZero++;
            } else {
                binIndexAndPositionInBin[i]=std::make_pair<int,int>(-1,-1); //bin (and index within bin) is -1 if not in the hash table
            }
        }
        std::cout << " numNonZero="<<numNonZero<<"\n";
        //set rxn counter to 0
        rxnCountThisBuildOrRebuild=0;
        //record info from build
        rxnsBetweenBuildOrRebuild.push_back(0);
        lowerBoundsBuildOrRebuild.push_back(currentLowerBound); //keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
        upperBoundsBuildOrRebuild.push_back(currentUpperBound);

        if (active_channel_count!=activeChannelCounter) {
            std::cout << "ERROR: active channel count is inconsistent.\n";
            return false;
        }

        printHashTable();
        return true;
    }


    // reactionIndex is particle id
    void NRMConstant_v5::update(std::size_t reactionIndex, double newPropensity, double currentTime, std::mt19937_64& rng) {

        printf("NRMConstant_v5::update(");
        printf("reactionIndex=%lu,",reactionIndex);
        printf("newPropensity=%e,",newPropensity);
        printf("currentTime=%e,",currentTime);
        printf(")\n");
        double firingTime;
        if (newPropensity <= 0.0) {
            firingTime=std::numeric_limits<double>::infinity();
            if (nextFiringTime.at(reactionIndex)!=std::numeric_limits<double>::infinity()) {
                activeChannelCounter--;
            }
        } else {
            firingTime=exponential(rng)/newPropensity+currentTime;
            if (nextFiringTime.at(reactionIndex)==std::numeric_limits<double>::infinity()) {
                activeChannelCounter++;
            }
        }
        nextFiringTime[reactionIndex]=firingTime;
        int newBin=computeBinIndex(firingTime);

        int oldBin=binIndexAndPositionInBin[reactionIndex].first;
        int oldPositionInBin=binIndexAndPositionInBin[reactionIndex].second;
        std::cout << "\tnewBin="<<newBin<<" oldBin="<<oldBin<<"  oldPositionInBin="<<oldPositionInBin<<"\n";
        if (newBin!=oldBin) {
            std::cout << "\tnewBin!=oldBin\n";
            if (oldBin>=0) {
                //remove from old bin
                if (theHashTable[oldBin].size()>1) {
                    std::cout << "\tdoing swap\n";
                    //take last element in old bin and place in this element's spot
                    theHashTable[oldBin][oldPositionInBin]=theHashTable[oldBin].back();
                    std::size_t movedElementIndex=theHashTable[oldBin][oldPositionInBin].second;
                    //update old last element's ...within bin index
                    binIndexAndPositionInBin[movedElementIndex].second=oldPositionInBin;
                }
                std::cout << "\ttheHashTable[oldBin].size()="<<theHashTable[oldBin].size()<<" (before)\n";
                theHashTable[oldBin].pop_back();
                std::cout << "\ttheHashTable[oldBin].size()="<<theHashTable[oldBin].size()<<" (after)\n";
            }
            binIndexAndPositionInBin[reactionIndex].first=newBin;
            if (newBin>=0) {
                theHashTable[newBin].push_back(std::pair<double,std::size_t>(firingTime,reactionIndex));
                binIndexAndPositionInBin[reactionIndex].second=theHashTable[newBin].size()-1;
            }
        } else {
            std::cout << "\tnewBin==oldBin\n";
            //just update firing time
            if (newBin>=0) {
                theHashTable[newBin][binIndexAndPositionInBin[reactionIndex].second].first=firingTime;
            }
        }

    }

    void NRMConstant_v5::printHashTable() {
        std::cout << "NRMConstant_v5::printHashTable() ";
        std::cout << "###########################################\n";
        std::cout << " minBin=" << minBin; 
        std::cout << " lower="<<currentLowerBound;
        std::cout << " upper="<<currentUpperBound;
        std::cout << " width=" << (currentUpperBound-currentLowerBound)/numberOfBins << std::endl;
        for (std::size_t i=0; i!=theHashTable.size(); ++i) {
            std::cout << "[" << i << "] (sz=" << theHashTable[i].size() << "): ";
            double mint = std::numeric_limits<double>::infinity();
            double maxt = -1 * std::numeric_limits<double>::infinity();
            for (std::size_t j=0; j!=theHashTable[i].size(); ++j) {
                if(theHashTable[i][j].first < mint){ mint = theHashTable[i][j].first; }
                if(theHashTable[i][j].first > maxt){ maxt = theHashTable[i][j].first; }
                //std::cout << theHashTable[i][j].second << "(" << theHashTable[i][j].first << "), ";
                printf("%lu(%.3e) ",theHashTable[i][j].second, theHashTable[i][j].first);
            }
            std::cout << " ["<<mint<< " "<<maxt<<"] ["<<(maxt-mint)<<"]\n";
        }
        std::cout << "###########################################\n";
    }

    NRMConstant_v5::timeIndexPair NRMConstant_v5::selectReaction() {
        std::cout << "NRMConstant_v5::selectReaction()\n";
        printHashTable();
        //while front of queue is empty, pop it
        std::cout << " minBin="<<minBin;
        std::cout << " theHashTable[minBin].size()="<<theHashTable[minBin].size()<<"\n";
        while (theHashTable[minBin].size()==0) {
            ++minBin;
            std::cout << " ++minBin="<<minBin;
            std::cout << " theHashTable[minBin].size()="<<theHashTable[minBin].size()<<"\n";
            if (minBin==theHashTable.size()) {
                std::cout << " rebuild\n";
                if (!rebuild()) {
                    std::cout << " rebuild returned False\n";
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
        std::cout << " minTimeRxnIndex="<<minTimeRxnIndex<<"\n";
        previousFiringTime=theHashTable[minBin][minTimeRxnIndex].first;
        std::cout << " previousFiringTime="<<previousFiringTime<<"\n";
        ++rxnCountThisBuildOrRebuild;
        std::cout << " rxnCountThisBuildOrRebuild=" << rxnCountThisBuildOrRebuild << "\n";
        std::cout << "=theHashTable[minBin][minTimeRxnIndex]=(";
        std::cout <<theHashTable[minBin][minTimeRxnIndex].first<<" , ";
        std::cout <<theHashTable[minBin][minTimeRxnIndex].second<<")\n";

        return theHashTable[minBin][minTimeRxnIndex];
    }

    NRMConstant_v5::~NRMConstant_v5() {
    }

    int NRMConstant_v5::computeBinIndex(double firingTime) {
        if (firingTime>currentUpperBound) {
            //std::cout << "computeBinIndex("<<firingTime<<")==-1\n";
            return -1;
        }
        int binIndex=(int)((firingTime-currentLowerBound)/(currentUpperBound-currentLowerBound)*numberOfBins);
        int returnVal=std::max(binIndex,(int)minBin);
        int ret = std::min<int>(returnVal, numberOfBins-1);
        std::cout << "computeBinIndex("<<firingTime<<")=="<<ret<<"\n";
        return ret;
    }

    //returns false if all propensities are 0
    bool NRMConstant_v5::rebuild() {
        std::cout << "NRMConstant_v5::rebuild()\n";

        //estimate propensitySum based on number of steps since last build or rebuild
        double propensitySumEstimate;
        if (rxnCountThisBuildOrRebuild>0) {
            propensitySumEstimate=(double)rxnCountThisBuildOrRebuild/(currentUpperBound-currentLowerBound);
        } else {
            //this shouldn't happen...if we got here, 0 reactions fired since last
            //rebuild!  This should only arise in toy problems, but need a strategy.
            double previousBinWidth=(currentUpperBound-currentLowerBound)/(double)numberOfBins;
            propensitySumEstimate=1.0/previousBinWidth;
            std::cout << "WARNING: 0 reactions fired before rebuild.\n";
        }

        if(!setBinNumberAndBounds(currentUpperBound,propensitySumEstimate,activeChannelCounter)){//
            minBin=0;
            return false;
        }
        minBin=0;

        std::vector<std::pair<double, std::size_t> > emptyVector;

        theHashTable.clear();
        theHashTable.resize(numberOfBins,emptyVector);

        int numNonZero=0;
        for (auto it = nextFiringTime.cbegin(); it != nextFiringTime.end(); it++) {
            //insert into hashTable
            int bin=computeBinIndex(it->second); // it->second is firing time
            if (bin>=0) {
                theHashTable[bin].push_back(std::make_pair(it->second,it->first)); //place this rxn at back of bin
                binIndexAndPositionInBin[it->first]=std::make_pair(bin,theHashTable[bin].size()-1);
                numNonZero++;
            } else {
                binIndexAndPositionInBin[it->first]=std::make_pair<int,int>(-1,-1);
            }
        }
        std::cout << " numNonZero="<<numNonZero<<"\n";


        //record info from build
        rxnsBetweenBuildOrRebuild.push_back(rxnCountThisBuildOrRebuild);
        lowerBoundsBuildOrRebuild.push_back(currentLowerBound); //keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
        upperBoundsBuildOrRebuild.push_back(currentUpperBound);

        //set rxn counter to 0
        rxnCountThisBuildOrRebuild=0;

        printHashTable();
        return numNonZero>0;
    }

    //fixed number of bins
    bool NRMConstant_v5::setBinNumberAndBounds(double newLowerBound, double propensitySum, int activeChannels) {
        if (activeChannels==0) {
            return false;
        }

        currentLowerBound=newLowerBound;
        if (currentLowerBound>endTime) {
            return false;
        }
        double binWidth=16.0/propensitySum;
        // heuristic
        std::size_t maxNumberOfBins=(std::size_t)20.0*sqrt(activeChannels);

        //how many bins to get beyond end time?
        std::size_t endTimeBins=std::numeric_limits<std::size_t>::max();
        //if we know the end time:
        if (endTime<std::numeric_limits<double>::max()) {
            endTimeBins=(std::size_t)ceil( (endTime-currentLowerBound)/binWidth);
            //use number of bins that takes us past endTime, or maxNumberOfBins,
            //whichever is smaller
            numberOfBins=std::min(endTimeBins,maxNumberOfBins);
        } else {
            numberOfBins=maxNumberOfBins;
        }
        currentUpperBound=currentLowerBound+binWidth*(double)numberOfBins;

        //if max bins gets us very close to the simulation end time, bump it up
        //beyond max to avoid an additional rebuild
        if (endTime<std::numeric_limits<double>::max()) {
            if(currentUpperBound<endTime &&
               ((double)(endTimeBins-maxNumberOfBins)/((double)numberOfBins))<.2){
                numberOfBins=endTimeBins;
                currentUpperBound=currentLowerBound+binWidth*(double)numberOfBins;
            }
        }
        std::cout << "NRMConstant_v5::setBinNumberAndBounds()";
        std::cout << " currentLowerBound="<<currentLowerBound;
        std::cout << " currentUpperBound="<<currentUpperBound;
        std::cout << " numberOfBins="<<numberOfBins<<"\n";
        return true;
    }
}
