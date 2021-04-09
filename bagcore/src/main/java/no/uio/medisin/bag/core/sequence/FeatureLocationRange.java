/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.sequence;

/**
 * specifies the range within which a feature is found
 * by a feature, this can be a physical entity such as a miRNA relative to
 * a reference sequence (such as a pre-miRNA or genome) or the range of a 
 * value (such as the location of mapped reads on a pre-miRNA)
 * @author simonray
 */
public class FeatureLocationRange {
    
    private long plusMinStart           = -1;
    private long plusMaxStart           = -1;
    private long plusMinEnd             = -1;
    private long plusMaxEnd             = -1;
    private long minusMinStart          = -1;
    private long minusMaxStart          = -1;
    private long minusMinEnd            = -1;
    private long minusMaxEnd            = -1;
    private double  minusAverageStart   = 0.0;
    private double  plusAverageStart    = 0.0;

    
    
    
    
    
    @Override
    public String toString(){
        return 
                minusMinStart + "\t"
                + minusMaxStart + "\t"
                + minusMinEnd + "\t"
                + minusMaxEnd + "\t"    
                + minusAverageStart + "\t"
                + plusMinStart + "\t" 
                + plusMaxStart + "\t"
                + plusMaxStart + "\t"
                + plusMinEnd + "\t"
                + minusAverageStart + "\n";
    }
    
    
    
    
    /**
     * @return the plusMinStart
     */
    public long getPlusMinStart() {
        return plusMinStart;
    }

    /**
     * @param plusMinStart the plusMinStart to set
     */
    public void setPlusMinStart(long plusMinStart) {
        this.plusMinStart = plusMinStart;
    }

    /**
     * @return the plusMaxStart
     */
    public long getPlusMaxStart() {
        return plusMaxStart;
    }

    /**
     * @param plusMaxStart the plusMaxStart to set
     */
    public void setPlusMaxStart(long plusMaxStart) {
        this.plusMaxStart = plusMaxStart;
    }

    /**
     * @return the plusMinEnd
     */
    public long getPlusMinEnd() {
        return plusMinEnd;
    }

    /**
     * @param plusMinEnd the plusMinEnd to set
     */
    public void setPlusMinEnd(long plusMinEnd) {
        this.plusMinEnd = plusMinEnd;
    }

    /**
     * @return the plusMaxEnd
     */
    public long getPlusMaxEnd() {
        return plusMaxEnd;
    }

    /**
     * @param plusMaxEnd the plusMaxEnd to set
     */
    public void setPlusMaxEnd(long plusMaxEnd) {
        this.plusMaxEnd = plusMaxEnd;
    }

    /**
     * @return the minusMinStart
     */
    public long getMinusMinStart() {
        return minusMinStart;
    }

    /**
     * @param minusMinStart the minusMinStart to set
     */
    public void setMinusMinStart(long minusMinStart) {
        this.minusMinStart = minusMinStart;
    }

    /**
     * @return the minusMaxStart
     */
    public long getMinusMaxStart() {
        return minusMaxStart;
    }

    /**
     * @param minusMaxStart the minusMaxStart to set
     */
    public void setMinusMaxStart(long minusMaxStart) {
        this.minusMaxStart = minusMaxStart;
    }

    /**
     * @return the minusMinEnd
     */
    public long getMinusMinEnd() {
        return minusMinEnd;
    }

    /**
     * @param minusMinEnd the minusMinEnd to set
     */
    public void setMinusMinEnd(long minusMinEnd) {
        this.minusMinEnd = minusMinEnd;
    }

    /**
     * @return the minusMaxEnd
     */
    public long getMinusMaxEnd() {
        return minusMaxEnd;
    }

    /**
     * @param minusMaxEnd the minusMaxEnd to set
     */
    public void setMinusMaxEnd(long minusMaxEnd) {
        this.minusMaxEnd = minusMaxEnd;
    }

    /**
     * @return the minusAverageStart
     */
    public double getMinusAverageStart() {
        return minusAverageStart;
    }

    /**
     * @param minusAverageStart the minusAverageStart to set
     */
    public void setMinusAverageStart(double minusAverageStart) {
        this.minusAverageStart = minusAverageStart;
    }

    /**
     * @return the plusAverageStart
     */
    public double getPlusAverageStart() {
        return plusAverageStart;
    }

    /**
     * @param plusAverageStart the plusAverageStart to set
     */
    public void setPlusAverageStart(double plusAverageStart) {
        this.plusAverageStart = plusAverageStart;
    }
}
