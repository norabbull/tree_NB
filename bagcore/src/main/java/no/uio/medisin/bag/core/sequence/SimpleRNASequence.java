
package no.uio.medisin.bag.core.sequence;

import java.util.HashMap;


/**
 * simple RNA class
 * stores basic information about an RNASequence
 * 
 * extends SimSeq 
 * extended by PriMiRNA, PreMiRNA and MiRNA classes
 * @author weibo
 */
public class SimpleRNASequence extends SimpleSeq {
    
    private String              structureString="";         // in Vienna Bracket Notation
    private float               energy=0;
    private float               GC_content=0;
    private int                 GU_num=0;
    private int                 pair_num=0;

    protected String            upperStructLine1;
    protected String            upperStructLine2;
    protected String            lowerStructLine4;
    protected String            lowerStructLine5;
    protected String            middleStructLine; // The middle line showing base pairing

    protected HashMap           featureSet = new HashMap();     // stores additional characteristics of this miRNA
    
    
    /**
     * Basic Class Constructor
     * 
     */
    public SimpleRNASequence(){
        super();
        featureSet = new HashMap();
    }
    
    
    
    
    /**
     * Class Constructor from Sequence specified as Strings
     * 
     * @param id String
     * @param seq String
     */
    public SimpleRNASequence(String id, String seq){
        super(id,seq);
    }
    
    
    
    /**
     * Class Constructor from SimpleSeq
     * 
     * @param seq SimpleSeq
     */
    public SimpleRNASequence(SimpleSeq seq){
        super.setID(seq.getId());
        super.setName(seq.getName());
        super.setSeq(seq.getSeq());
        super.setAbsStartInQuerySeq(seq.getStartPos());
        super.setAbsEndInQuerySeq(seq.getEndPos());
        super.setChromosome(seq.getChromosome());
        super.setStartPos(seq.getStartPos());
        super.setEndPos(seq.getEndPos());
        super.setStrand(seq.getStrand());
    }
    
    
    
    public SimpleRNASequence(SimpleRNASequence srnas){

        super(srnas);
        structureString     = srnas.structureString;
        energy              = srnas.energy;
        GC_content          = srnas.GC_content;
        GU_num              = srnas.GU_num;
        pair_num            = srnas.pair_num;
        
    }
    
    
    
    
    
    /**
     * get the FeatureSet
     * 
     * @return HashMap FeatureSet
     */
    public HashMap getFeatureSet(){
        return featureSet;
    }

    /**
     * get the value for the query Key from the Feature table
     * 
     * @param qKey String
     * @return String value for this qKey
     */
    public String getFeature(String qKey){
        return (String) this.getFeatureSet().get(qKey);
    }
    
    /**
     * add key, value pair to feature set
     * 
     * @param pKey
     * @param pValue 
     */
    public void putFeature(String pKey, String pValue){
        this.getFeatureSet().put(pKey, pValue);
    }
    
    
    /**
     * @return str
     */
    public String getStructureStr() {
        return structureString;
    }

    /**
     * @param str Structure string
     */
    public void setStructureString(String str) {
        this.structureString = str;
    }

    /**
     * @return energy
     */
    public float getEnergy() {
        return energy;
    }

    /**
     * @param  energy float
     */
    public void setEnergy(float energy) {
        this.energy = energy;
    }

    /**
     * @return String miRStructLine1
     */
    public String getStructLine1() {
        return upperStructLine1;
    }

    /**
     * @param miRStructLine1 String  
     */
    public void setStructLine1(String miRStructLine1) {
        this.upperStructLine1 = miRStructLine1;
    }

    /**
     * @return String miRStructLine2
     */
    public String getStructLine2() {
        return upperStructLine2;
    }

    /**
     * @param miRStructLine2 String
     */
    public void setStructLine2(String miRStructLine2) {
        this.upperStructLine2 = miRStructLine2;
    }

    /**
     * @return String miRStructLine3
     */
    public String getStructLine3() {
        return lowerStructLine4;
    }

    /**
     * @param miRStructLine3 String
     */
    public void setStructLine3(String miRStructLine3) {
        this.lowerStructLine4 = miRStructLine3;
    }

    /**
     * @return String miRStructLine4
     */
    public String getStructLine4() {
        return lowerStructLine5;
    }

    /**
     * @param miRStructLine4 String
     */
    public void setStructLine4(String miRStructLine4) {
        this.lowerStructLine5 = miRStructLine4;
    }

}