/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.mirna;

import java.util.HashMap;
import no.uio.medisin.bag.core.sequence.SimpleRNASequence;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * This class is simply for storing the features of the miRNA.
 * calculations and predictions are performed in the class {@link no.uio.medisin.bag.core.CharacterizedPriMiRNA}
 * because predictions are based on a pri-miRNA/pre-miRNA/miRNA triplet.
 * i.e., one pri-miRNA will generally have many pre-miRNA/miRNA candidate pairs
 * and the prediction is based on the characteristics of all three structures.
 * Because there are so many features, and because this feature list may change,
 * they are stored in a HashMap

 * @author weibo / Simon Rayner
 */
public class deprecatedMiRNA extends SimpleRNASequence{
    
    static Logger logger = LogManager.getRootLogger();    

    public deprecatedMiRNA(){
        super();
    }
    
    public deprecatedMiRNA(String id,String seq){
        super(id,seq);
    }

    private String              host;
    private int                 largestInternalLoop;
    private int                 numOfInternalLoops;
    private int                 numOfUnpairedBases;
    private double              fractOfUnpairedBases;

    private char                firstBase;
    private double              stability;
    private char                dangleBaseOne;
    private char                dangleBaseTwo;

    private int                 miStart;
    private int                 miEnd;

    private String              note="";

    private HashMap             featureSet;

    
    
    
    
    /**
     * because the list of features we might calculate can change, we also store
     * in a hash map to make it simpler to collectively pass information
     */    
    public void buildFeatureSet(){
        featureSet=new HashMap();
        
        featureSet.put("miRNA_id",                  this.getId());
        featureSet.put("miRNA_sequence",            this.getSeq());
        featureSet.put("miRNA_structure",           this.getStructureStr());
        featureSet.put("miRNA_energy",              this.getEnergy());
        featureSet.put("miRNA_size",                this.getLength());
        featureSet.put("miRNA_GC_content",          this.GCfraction() );
        featureSet.put("miRNA_A_content",           this.Afraction());
        featureSet.put("miRNA_U_content",           this.Ufraction());
        featureSet.put("miRNA_G_content",           this.Gfraction());
        featureSet.put("miRNA_C_content",           this.Cfraction());
        featureSet.put("strand",                    this.getStrand());

        featureSet.put("miStart",                   this.getStartPos());
        featureSet.put("miEnd",                     this.getEndPos());
    }

    
    
    
    public HashMap getFeatureSet(){
        return featureSet;
    }

    /**
     * @return the maxInternalLoopSize
     */
    public int getInternalLoopSize() {
        return largestInternalLoop;
    }

    /**
     * @param maxInternalLoopSize the maxInternalLoopSize to set
     */
    public void setLargestInternalLoop(int maxInternalLoopSize) {
        this.largestInternalLoop = maxInternalLoopSize;
    }

    /**
     * @return the internalLoop_num
     */
    public int getInternalLoop_num() {
        return numOfInternalLoops;
    }

    /**
     * @param internalLoop_num the internalLoop_num to set
     */
    public void setNumOfInternalLoops(int internalLoop_num) {
        this.numOfInternalLoops = internalLoop_num;
    }

    /**
     * @return the unpairedBase_num
     */
    public int getUnpairedBase_num() {
        return numOfUnpairedBases;
    }

    /**
     * @param unpairedBase_num the unpairedBase_num to set
     */
    public void setNumberOfUnpairedBases(int unpairedBase_num) {
        this.numOfUnpairedBases = unpairedBase_num;
    }

    /**
     * @return the unpairedBase_rate
     */
    public double getUnpairedBase_rate() {
        return fractOfUnpairedBases;
    }

    /**
     * @param unpairedBase_rate the unpairedBase_rate to set
     */
    public void setFractOfUnpairedBases(double unpairedBase_rate) {
        this.fractOfUnpairedBases = unpairedBase_rate;
    }

    /**
     * @return the firstBase
     */
    public char getFirstBase() {
        return firstBase;
    }

    /**
     * @param firstBase the firstBase to set
     */
    public void setFirstBase(char firstBase) {
        this.firstBase = firstBase;
    }

    /**
     * @return the stability
     */
    public double getStability() {
        return stability;
    }

    /**
     * @param stability the stability to set
     */
    public void setStability(double stability) {
        this.stability = stability;
    }

    /**
     * @return the dangleBaseOne
     */
    public char getDangleBaseOne() {
        return dangleBaseOne;
    }

    /**
     * @param dangleBaseOne the dangleBaseOne to set
     */
    public void setDangleBaseOne(char dangleBaseOne) {
        this.dangleBaseOne = dangleBaseOne;
    }

    /**
     * @return the dangleBaseTwo
     */
    public char getDangleBaseTwo() {
        return dangleBaseTwo;
    }

    /**
     * @param dangleBaseTwo the dangleBaseTwo to set
     */
    public void setDangleBaseTwo(char dangleBaseTwo) {
        this.dangleBaseTwo = dangleBaseTwo;
    }

    /**
     * @return the miEnd
     */
    public int getMiEnd() {
        return miEnd;
    }

    /**
     * @param miEnd the miEnd to set
     */
    public void setMiRNAEndPos(int miEnd) {
        this.miEnd = miEnd;
    }

    /**
     * @return the miStart
     */
    public int getMiStart() {
        return miStart;
    }

    /**
     * @param miStart the miStart to set
     */
    public void setMiRNAStartPos(int miStart) {
        this.miStart = miStart;
    }

    /**
     * @return the note
     */
    public String getNote() {
        return note;
    }

    /**
     * @param note the note to set
     */
    public void setNote(String note) {
        this.note = note;
    }

    /**
     * @return the host
     */
    public String getHost() {
        return host;
    }

    /**
     * @param host the host to set
     */
    public void setHost(String host) {
        this.host = host;
    }

}
