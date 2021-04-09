
package no.uio.medisin.bag.core.mirna;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import no.uio.medisin.bag.core.sequence.SimpleRNASequence;

/**
 * Contains additional parameters which serve to characterize the pri-miRNA
 * and which can be later used in miRNA prediction
 * 
 * Since a pri-miRNA can contain multiple pre-miRNAs, it's not clear how a
 * pri-miRNA should be used. i.e., in the original implementation of the class
 * the pri-miRNA was considered to encapsulate a single pre-miRNA, but this isn't
 * an accurate representation of what a pri-miRNA really is.
 * 
 * @author weibo and Simon Rayner
 */
public class PriMiRNA extends SimpleRNASequence{

    public PriMiRNA() {
        super();
    }

    public PriMiRNA(String id, String seq) {
        super(id, seq);
    }

    // positions count from 1
    private int         basalSegEnd = 0;
    private int         basalSegSize = 0;
    private String      terminalLoopSeq;
    private int         terminalLoopStart = 0;
    private int         terminalLoopEnd = 0;
    private int         terminalLoopSize = 0;
    private int         internalLoopSize = 0;
    private int         internalLoop_num = 0;
    private int         numberOfUnpairedBasesInStem = 0;
    private double      unpairedBaseRate = 0;
    private String      seq5;                       // 5' seq of pri-miRNA
    private String      str5;                       // 5' structure of pri-miRNA
    private String      seq3;                       // 3' seq of pri-miRNA
    private String      str3;                       // 3' structure of pri-miRNA
    private String      middleBase;                 // mid-base of terminal loop (can be empty)
    private String      priLine1;                   // unpaired nt on 5' arm in the pretty plot
    private String      priLine2;                   // paired nt on 5' arm in the pretty plot
    private String      priLine3;                   // unpaired nt on 3' arm in the pretty plot
    private String      priLine4;                   // paired nt on 3' arm in the pretty plot
    private String      priPlot;                    // <-- is this used ?
    private int[]       strIndex;                   // this tells us which nt in the 5' arm maps to which nt in the 3' arm

    


    //private HashMap     featureSet;


    private ArrayList<PreMiRNA> pres = new ArrayList<>();
    private int index=0;

    /**
     * add a pre-miRNA that is positioned within this pri-miRNA
     * 
     * @param preRNA PreMiRNA. The pre-miRNA to be added to the set
     */
    public void addProduct(PreMiRNA preRNA){
        pres.add(preRNA);
        index+=1;
    }

    /**
     * 
     * @return number of entries in the product list
     */
    public int sizeOfProduct(){
        return pres.size();
    }

    
    /**
     * There are too many possible characteristics that could be used to define
     * a pri-miRNA so we store them in a HashMap. We start by adding the parameters
     * that are defined as part of the pri-miRNA class
     */
    public void buildFeatureSet(){
        
        //featureSet=new HashMap();
        featureSet.put("priRNA_id", this.getId());
        featureSet.put("priRNA_start", this.getStartPos());
        featureSet.put("priRNA_end", this.getEndPos());
        featureSet.put("priRNA_sequence", this.getSeq());
        featureSet.put("priRNA_structure", this.getStructureStr());
        featureSet.put("priRNA_energy", this.getEnergy());
        featureSet.put("priRNA_size", this.getLength());
        featureSet.put("priRNA_plot", CharPriMiRNA.prettyPlotToString(this));
        featureSet.put("priRNA_GC_content", this.GCfraction() );
        featureSet.put("priRNA_A_content", this.Afraction());
        featureSet.put("priRNA_U_content", this.Ufraction());
        featureSet.put("priRNA_G_content", this.Gfraction());
        featureSet.put("priRNA_C_content", this.Cfraction());
    }
    
    
    
    
    /**
     * get the key values as a tab delimited string
     * 
     * @return String FeatureSet
     */
    public String printFeatureSetKeysAsTSV(){
        String featureKeyStr = "NAME|ID:\t";
        Set set = featureSet.entrySet();
        Iterator itSet = set.iterator();
        
        while(itSet.hasNext()){
            Map.Entry me = (Map.Entry)itSet.next();
            featureKeyStr = featureKeyStr.concat((String)me.getKey() + "\t");
            
        }
        return featureKeyStr.trim();
        
    }
    
    
    
    
    /**
     * print this pri-miRNA's features as a tab delimited string
     * 
     * @return String tab-delimited features
     */
    public String printFeatureSetValuesAsTSV(){
        
        String featureValueStr = this.getId() + "|" + this.getName() + "\t";
        Set set = featureSet.entrySet();
        Iterator itSet = set.iterator();
        
        while(itSet.hasNext()){
            Map.Entry me = (Map.Entry)itSet.next();
            featureValueStr = featureValueStr.concat(me.getValue() + "\t");
            
        }
        return featureValueStr.trim();
    }
    
    
    
    
    
    /**
     * @return the priLine1
     */
    public String getPriLine1() {
        return priLine1;
    }

    /**
     * @param priLine1 String
     */
    public void setPriLine1(String priLine1) {
        this.priLine1 = priLine1;
    }

    /**
     * @return String priLine2
     */
    public String getPriLine2() {
        return priLine2;
    }

    /**
     * @param priLine2 String
     */
    public void setPriLine2(String priLine2) {
        this.priLine2 = priLine2;
    }

    /**
     * @return String priLine3
     */
    public String getPriLine3() {
        return priLine3;
    }

    /**
     * @param priLine3 String
     */
    public void setPriLine3(String priLine3) {
        this.priLine3 = priLine3;
    }

    /**
     * @return String priLine4
     */
    public String getPriLine4() {
        return priLine4;
    }

    /**
     * @param priLine4 String
     */
    public void setPriLine4(String priLine4) {
        this.priLine4 = priLine4;
    }

    

    

    /**
     * @return int[] strIndex
     */
    public int[] getStrIndex() {
        return strIndex;
    }

    /**
     * @param    strIndex inrt[]
     */
    public void setStrIndex(int[] strIndex) {
        this.strIndex = strIndex;
    }

    /**
     * @return the seq5
     */
    public String getSeq5() {
        return seq5;
    }

    /**
     * @param seq5 the seq5 to set
     */
    public void setSeq5(String seq5) {
        this.seq5 = seq5;
    }

    /**
     * @return the str5
     */
    public String getStr5() {
        return str5;
    }

    /**
     * @param str5 the str5 to set
     */
    public void setStr5(String str5) {
        this.str5 = str5;
    }

    /**
     * @return the seq3
     */
    public String getSeq3() {
        return seq3;
    }

    /**
     * @param seq3 the seq3 to set
     */
    public void setSeq3(String seq3) {
        this.seq3 = seq3;
    }

    /**
     * @return the str3
     */
    public String getStr3() {
        return str3;
    }

    /**
     * @param str3 the str3 to set
     */
    public void setStr3(String str3) {
        this.str3 = str3;
    }

    /**
     * @return the middleBase
     */
    public String getMiddleBase() {
        return middleBase;
    }

    /**
     * @param middleBase the middleBase to set
     */
    public void setMiddleBase(String middleBase) {
        this.middleBase = middleBase;
    }

    
}
