
package no.uio.medisin.bag.core.sequence;

import no.uio.medisin.bag.core.mirna.MfeFoldRNA;
import no.uio.medisin.bag.core.mirna.PriMiRNA;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * this class finds hairpin structures within an input string in Vienna Bracket Notation.
 * 
 * @see <a href=http://rna.tbi.univie.ac.at/help.html#A6>  http://rna.tbi.univie.ac.at/help.html</a>
 * 
 * 
 * the order of calling methods after instantiation:
 * slidingWindow();
 * scanStemLoop();
 * noRepeat();
 * @author weibo and simon rayner
 */
public class FoldableRNASequence extends SimpleRNASequence{
    
    static Logger logger = LogManager.getRootLogger();    
    
    private static final int        MIN_FRAGMENT_SIZE   = 19;

    //private int                     window=500;
    //private int                     step=250;
    //private int                     start=1;
    private int                     minHairpinLength = 20;

    //private SimpleSeq               simpleSeq;
    
    private ArrayList<SimpleSeq>    fragmentList;
    private ArrayList<PriMiRNA>     primiRNAList;

    /*
      private static String regSL="([\\(\\.]*\\()(\\.+)(\\)[\\)\\.]*)";
      
      This can be broken down into
     
        {'('} + {zero or more '.' or '('} + {stem loop} + {zero or more '.' or ')'} + {')'}
    */
    private static String regSL="(\\.*\\([\\.\\(]*\\()(\\.+)(\\)[\\.\\)]*\\)\\.*)";
    private static Pattern pattern = Pattern.compile(regSL);
    private Matcher m;

    
    /**
     * Empty Class Constructor
     * 
     */
    public FoldableRNASequence(){
        
    }
    
    
    
    /**
     * Class constructor from input sequence
     * 
     * @param seq  SimpleSeq
     */
    public FoldableRNASequence(SimpleSeq seq){
        super(seq);
        
    }
    
    
    
    

    /**
     * Fold the sequence and scan resulting structure for a hairpin
     * @return ArrayList of PriMiRNAs
     */
    public ArrayList<PriMiRNA> foldAndScanSequence(){
        
        //SimpleRNASequence rna = new SimpleRNASequence(this.getSimpleSequence());
        foldRNASequence();

        return extractHairpinsFromBracketNotation();
    }

    
    
    
    /**
     * 
     * fold an RNA sequence using RNAFold and return predicted structure
     * in Vienna Bracket Notation
     * 
     * @param rnaSequence : the sequence to be folded
     */
    public static void foldRNASequence(SimpleRNASequence rnaSequence){
                
        rnaSequence.setEnergy(MfeFoldRNA.foldSequence(rnaSequence.getSeq()));
        rnaSequence.setStructureString(MfeFoldRNA.getStructure());
        logger.info("Energy" + rnaSequence.getEnergy());

    }

    
    
    
    /**
     * 
     * fold an RNA sequence using RNAFold and return predicted structure
     * in Vienna Bracket Notation
     * 
     */
    public void foldRNASequence(){
                
        this.setEnergy(MfeFoldRNA.foldSequence(this.getSeq()));
        this.setStructureString(MfeFoldRNA.getStructure());
        logger.info("Energy" + this.getEnergy());

    }

    
    
    
    /**
     * 
     * extract hairpins from a rna secondary structure defined using the 
     * Vienna Bracket Notation
     * 
     * 
     * @return ArrayList of PrimiRNA sequences
     * 
     */
    public ArrayList<PriMiRNA> extractHairpinsFromBracketNotation(){
        
        ArrayList<PriMiRNA> thisPrimiRNAList=new ArrayList<>();
        
        String str = this.getStructureStr(); 
        m = pattern.matcher(str);

        while(m.find()){

            // strip off the hanging ends from the structure
            // ..(((..((((...)))).)))... 
            //          becomes
            //   (((..((((...)))).)))   
            int start1      = str.lastIndexOf(")", m.start(1)) + 1; //replace m.start(1)
            String left     = str.substring(start1, m.end(1));//the 5' side of the stemloop
            String right    = m.group(3);//the 3' side of the stemloop
            
            int n = b2Num(left)+b2Num(right);//compare the two sides of the stem
            int slStart=0, slEnd=0, l=0, r=0;
            
            //if str is like (..((...)).. return ..((...))..
            if(n>0){
                l=bIndex(left,"(",n)+1; //new start of left
                left=left.substring(l); //new left
                slStart=start1+l; //start of stemloop, count from 0
                slEnd=m.end(3);//count from 1
            }
            
            //if str is like ..((...))..) return ..((...))..
            else if(n<0){
                r=bIndex(right,")",n)-1;
                right=right.substring(0,r+1);
                slStart=start1;
                slEnd=m.start(3)+r+1;//count from 1
            }
            
            //if str is like ..((...)).. return ..((...))..
            else{
                slStart=start1;
                slEnd=m.end(3);//count from 1
            }
            
            String subId    = this.getName()+"_mir_"; //the id of the stemloop
            String subSeq   = this.getSeq().substring(slStart, slEnd); //seq of the stemloop
            String subStr   = left + m.group(2) + right; //structure of the stemloop


            if(subStr.length() < minHairpinLength)
                continue;

            //create a new pri-miRNA
            PriMiRNA pri = new PriMiRNA(subId, subSeq);

            //if the stemloop has two or more loops
            if(hasMultipleLoops(pri)) continue;
            if(pattern.matcher(pri.getStructureStr()).matches() == false) continue;
            
            pri.setName(this.getName());
            pri.setAbsStartInQuerySeq(slStart+this.getStartPos());
            pri.setAbsEndInQuerySeq(slEnd-1+this.getStartPos());
            pri.setID(pri.getId()+pri.getStartPos()+"-"+pri.getEndPos());
            
            thisPrimiRNAList.add(pri);

            
        }
        return thisPrimiRNAList;
    }

    
    
    
    /**
     * 
     * remove repeat hairpins caused by overlapping fragments
     * 
     */
    public void removeDuplicatePrimiRNAs(){
        HashMap map=new HashMap();
        for(PriMiRNA pri:primiRNAList)
            map.put(pri.getId(), pri);

        primiRNAList = new ArrayList<PriMiRNA>(map.values());
    }

    
    
    
    /**
     * transform bracket-dot string to number and return the sum
     * @param String str: structure with bracket-dot notation string
     * @return int: the sum of the str( each '(' is 1, ')' is -1, '.' is 0
     */
    private int b2Num(String str){
        int num=0;
        for(int i=0;i<str.length();i++){
            if(str.charAt(i)=='(')
                num+=1;
            else if(str.charAt(i)==')')
                num-=1;

        }
        return num;
    }
    
    
    
    
    /**
     * find the index of the nth '(' or ')'
     * @param String p: the structure string
     * @param String s: "(" or ")"
     * @param int n: the '(' or ')' number
     * @return
     */
    private int bIndex(String p,String s,int n){
        int m=Math.abs(n);
        int c=0;
        if(s.equals("(")){
            for(int i=0;i<m;i++){
                c=p.indexOf(s,c)+1;
            }
            c=c-1;
        }
        if(s.equals(")")){
            c=p.length()-1;
            for(int i=0;i<m;i++){
                c=p.lastIndexOf(s, c)-1;
            }
            c=c+1;
        }

        return c;
    }

    
    
    
    /**
     * check if the structure of a sequence has two or more loops
     * 
     * @param  rna SimpleRNASequence seq: sequence to be tested
     * @return boolean; if have two or more return true, or false
     * 
     */
    public static boolean hasMultipleLoops(SimpleRNASequence rna){
        
        foldRNASequence(rna);
        int end5=rna.getStructureStr().lastIndexOf("(");
        int start3=rna.getStructureStr().indexOf(")");
        if(end5>=start3){
            return true;
        }
        return false;
    }

    /**
     * @return the cutoff
     */
    public int getDistance() {
        return minHairpinLength;
    }

    /**
     * @param minHairpinLen: shortest hairpin to pass
     */
    public void setMinHairpinLength(int minHairpinLen) {
        this.minHairpinLength = minHairpinLen;
    }

    /**
     * @return the priList
     */
    public ArrayList<PriMiRNA> getPriList() {
        return primiRNAList;
    }

    /**
     * @param priList the priList to set
     */
    public void setPriList(ArrayList<PriMiRNA> priList) {
        this.primiRNAList = priList;
    }

    /**
     * @return the segList
     */
    public ArrayList<SimpleSeq> getSegList() {
        return fragmentList;
    }

    /**
     * @param segList the segList to set
     */
    public void setSegList(ArrayList<SimpleSeq> segList) {
        this.fragmentList = segList;
    }


}
