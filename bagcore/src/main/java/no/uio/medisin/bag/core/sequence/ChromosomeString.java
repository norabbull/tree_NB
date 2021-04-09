/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.sequence;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * a simple utility class to extract chromosome information from a string
 * This class is primarily intended for providing consistent representation of
 * chromosome names specified in a genome FASTA file. This may vary between entries
 * with nomenclature such as 'chr1', 'CHR1', 'Chr1' and '1' all possible. Additionally
 * we want a constant numerical representation for the X, Y and M chromosomes to allow
 * us to sort chromosomes, or data defined in terms of chromosomes (such as a mapped read
 * which can have a start position, a strand and a chromosome)
 * Thus, we assign a value of -1, -2 and -3 to the M, X & Y chromosomes respectively
 * 
 * @author simonray
 */
public class ChromosomeString implements Comparable<ChromosomeString>{
    
    static Logger logger = LogManager.getLogger();
    
    static final int CHRMINT    = -1;
    static final int CHRXINT    = -2;
    static final int CHRYINT    = -3;
    static final int CHRUNKINT  = -99;
    
    private static String   chromosome;
    private static int      chrNo;
    
    /*
    add special nomenclature for FruitFly (Drosophila melanogaster) chromosomes
    */
    static final int CHR_DME2RHET    = -1121;
    static final int CHR_DME2R       = -112;
    static final int CHR_DME2L       = -113;
    static final int CHR_DME3R       = -114;
    static final int CHR_DME3L       = -115;
    
    public ChromosomeString(String chr){
        chromosome = chr;
        chrNo = ChromosomeString.GetChromNumber(chr);
    }
    
    
    
    /**
     * return an integer representing the chromosome
     * @param chrString - the string representing the chromosome
     * @return a 
     */
    public static int GetChromNumber(String chrString){
        chromosome = chrString;
        chrString = chrString.toUpperCase();
        if(chrString.contains("CHR")){
            chrString = chrString.replace("CHR", "");
        }
        
        if(chrString.contains("M") || chrString.contains("MT"))
            return CHRMINT;
        
        if(chrString.contains("X"))
            return CHRXINT;
        
        if(chrString.contains("Y"))
            return CHRYINT;
        /*
        since Heterochromatin is important in expression regulation
        so still addd a if loop to capture that eventhough only three in genome
        check here https://en.wikipedia.org/wiki/Heterochromatin
        by Joey
        */
        if (chrString.contains("2R")){
            if (chrString.endsWith("Het")){
                return CHR_DME2RHET;
            }
            else{
                return CHR_DME2R;
            }
        }
        if (chrString.contains("2L"))
            return CHR_DME2L;
        if (chrString.contains("3L"))
            return CHR_DME3L;
        if (chrString.contains("3R"))
            return CHR_DME3R;
        /*
        at this point we should be left with a number
        */
        try{
            return Integer.parseInt(chrString);           
        }
        catch(NumberFormatException exNM){
            logger.debug("couldn't get chromosome for string <"+ chromosome + ">");
            return CHRUNKINT;
        }
        
    }

    
    
    /**
     * 
     * @param chrString String
     */
    public void setChromNumber(String chrString){
        
        chromosome = chrString.trim();
        
        chrString = chrString.toUpperCase();
        if(chrString.contains("CHR")){
            chrString = chrString.replace("CHR", "");
        }
        
        if(chrString.contains("M") || chrString.contains("MT"))
            chrNo = CHRMINT;
        
        if(chrString.contains("X"))
            chrNo = CHRXINT;
        
        if(chrString.contains("Y"))
            chrNo = CHRYINT;
        /*
        since Heterochromatin is important in expression regulation
        so still addd a if loop to capture that eventhough only three in genome
        check here https://en.wikipedia.org/wiki/Heterochromatin
        by Joey
        */
        if (chrString.contains("2R")){
            if (chrString.endsWith("Het")){
                chrNo = CHR_DME2RHET;
            }
            else{
                chrNo = CHR_DME2R;
            }
        }
        if (chrString.contains("2L"))
            chrNo = CHR_DME2L;
        if (chrString.contains("3L"))
            chrNo = CHR_DME3L;
        if (chrString.contains("3R"))
            chrNo = CHR_DME3R;
        /*
        at this point we should be left with a number
        */
        try{
            chrNo = Integer.parseInt(chrString);           
        }
        catch(NumberFormatException exNM){
            logger.info("couldn't get chromsome type for string <"+ chromosome + ">");
            logger.error("couldn't get chromsome type for string <"+ chromosome + ">");
            chrNo = CHRUNKINT;
        }
        
    }
    /**
     * @return the chrNo
     */
    public int getChrNo() {
        return chrNo;
    }

    /**
     * @param aChrNo the chrNo to set
     */
    public static void setChrNo(int aChrNo) {
        chrNo = aChrNo;
    }
    
    
    
    
    
    
    /**
     * compare by numerical representation of chromosome
     * 
     * @param cStr ChromosomeString
     * @return int outcome of comparison
     */
    @Override 
    public int compareTo(ChromosomeString cStr){
        return this.getChrNo() - cStr.getChrNo();
    }
    
    
    
    
    
    
    /**
     * return the chromosome in user friendly format, i.e. convert M, X and Y from
     * integer back to their letter representation
     * 
     * @param chr int
     * @return String chromosome as String
     * 
     */
    public static String getChromosomeAsString(int chr){
        switch(chr){
            case -1:
                return "M";
                
            case -2:
                return "X";
                
            case -3:
                return "Y";
            case  -1121:
                return "2RHet";
            case -112:
                return "2R";
            case -113:
                return "2L";
            case -114:
                return "3R";
            case -115:
                return "3L";    
           
            default:
                try{
                    return Integer.toString(chr);
                }
                catch(NumberFormatException exNF){
                    logger.warn("found a non standard chromosome name <" + chr + ">");
                    logger.warn("this usually occurs when the reference is a scafford sequence");
                    logger.warn("setting the chromosome name to <?>");
                    return "?";
                }
        }
    }
}
