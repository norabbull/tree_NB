/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mathstuff;

import java.util.Iterator;
import java.util.Map;
import org.apache.commons.math3.stat.Frequency;

/**
 * Estimate the equivalent of the FWHM for the specified distribution. While this
 * may not be entirely appropriate here, all we are trying to do is get a measure
 * of the dispersion of reads. 
 * In most cases, we would expect (hope?) to find all the reads within a single
 * region on the pre-miRNA. However, there are occasions where the reads are
 * 'scattered' along the sequence. In this case it isn't possible to calculate
 * a FWHM. To check for the second scenario, we begin by checking the distance
 * between the first (closest to 5' end) and last (closest to 3´end) reads. If 
 * it is greater than half the length of the sequence, then we aren't going to 
 * be able to calculate a 'standard' FWHM. In this case we simply report this 
 * 5' - 3' distance
 * 
 * @author simonray
 */
public class FWHM {
    
    /**
     * add baseline cutoff?
     * @param distribution - distribution to analyze
     * @param distLength - true width of the distribution (i.e., distribution 
     *                     may be a fixed size so that different FWHM estimates
     *                     can be compared, but the true width of meaningful 
     *                     values may be less. In this case ,this parameter can 
     *                     be used
     * @return double   the estimated FWHM 
     */
    public static double getFWHMfromFrequency(Frequency distribution, int distLength){
        
        long max = 0;
        int maxIndex = 0;
        long min = 0;
        int minIndex = 0;
        int firstRead = -1;
        int lastRead = -1;
        
        
        
        Iterator itES = distribution.entrySetIterator();
        int v=0;
        while (itES.hasNext()){
            Map.Entry pair = (Map.Entry) itES.next();
            if((long) pair.getValue() > 0 && firstRead < 0)
                firstRead = v;
            if((long)pair.getValue() > 0)
                lastRead = v;
            if((long) pair.getValue() > max) {
                max = (long) pair.getValue();
                maxIndex = v;
            }
            
            if((long) pair.getValue() < min) {
                min = (long) pair.getValue();
                minIndex = v;
            }
            
            v++;
        }
        if(lastRead-firstRead > distLength/2)
            return (double)lastRead-firstRead;
        
        double midValue = ((double) max - (double)min)/2.0;
        
        if(midValue == 0) 
            return 0.0;
        
        
        Map.Entry pair = null;
        Map.Entry lastPair = null;
        v = 0;
        itES = distribution.entrySetIterator();        
        while (itES.hasNext()){
            pair = (Map.Entry) itES.next();
            if((long) pair.getValue() > midValue) 
                break;
            lastPair = pair;
            v++;
        }
        

        double yLow = distribution.getCount(v);
        double yHigh = distribution.getCount(v+1);

        double m = yHigh - yLow;
        double xLow;
        if(yHigh==yLow)
            xLow = (long)lastPair.getKey() + ((long)pair.getKey() - (long)lastPair.getKey())/2.0;
        else
            xLow = (long)lastPair.getKey() + (midValue - yLow)/(yHigh-yLow);


        int vResume = v++;
        int vv = 0;
        itES = distribution.entrySetIterator();        
        while (itES.hasNext()){
            pair = (Map.Entry) itES.next();
            if(vv<vResume){
                lastPair = pair;
                vv++;
                continue;
            }
            if((long) pair.getValue() < midValue) 
                break;
            lastPair = pair;
            vv++;
        }

        
        yHigh = distribution.getCount(vv-1);
        yLow = distribution.getCount(vv);
        double xHigh;
        if(yHigh==yLow)
            xHigh = (long)lastPair.getKey() + ((long)pair.getKey() - (long)lastPair.getKey())/2.0;
        else
            xHigh= (long)lastPair.getKey() + (midValue - yLow)/(yHigh-yLow);
        
        return xHigh-xLow;
        
    }
    
    
    
    /**
     * a simpler measure of the spread of reads
     * 
     * @param distribution the distribution to be characterized
     * @param distLength - true width of the distribution (i.e., distribution 
     *                     may be a fixed size so that different FWHM estimates
     *                     can be compared, but the true width of meaningful 
     *                     values may be less. In this case ,this parameter can 
     *                     be used
     * @return double the distance between the first and last position of the reads
     */
    public static double firstMinusLastReadDistance(Frequency distribution, int distLength){
        int firstRead = -1;
        int lastRead = -1;

        Iterator itES = distribution.entrySetIterator();
        int v=0;
        while (itES.hasNext()){
            Map.Entry pair = (Map.Entry) itES.next();
            if((long) pair.getValue() > 0 && firstRead < 0)
                firstRead = v;
            if((long)pair.getValue() > 0)
                lastRead = v;
            
            v++;
        }
        
        
        
        return (double)lastRead-firstRead;

        

        
    }
}
