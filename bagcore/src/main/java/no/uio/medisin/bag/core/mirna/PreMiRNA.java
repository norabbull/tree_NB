/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.mirna;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.core.sequence.ChromosomeString;
import no.uio.medisin.bag.core.ShortPubMedEntry;
import no.uio.medisin.bag.core.sequence.FeatureLocationRange;
import no.uio.medisin.bag.core.sequence.SimpleRNASequence;
import no.uio.medisin.bag.core.sequence.Strand;
import no.uio.medisin.bag.core.sequence.StrandString;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.Frequency;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * stores sequence and structure information for a pre-miRNA and its child miRNAs.
 * this is a more general representation of the 'classical' representation of a
 * pre-miRNA which contains no more than two miRNAs, one on each strand. This
 * class can store all the reads that were mapped to the pre-miRNA and generate
 * summaries in terms of the miRNA distribution along the pre-miRNA
 * 
 * @author weibo / simon rayner
 */
public class PreMiRNA extends SimpleRNASequence{

    static Logger               logger                      = LogManager.getLogger();
    
    public final static     String        INDENT1 = "--";
    public final static     String        INDENT2 = "----";
    

    private String              host                            ="";
    private String              host3code                       ="";
    private String              miID;
    private String              parent;                         // this will be the enclosing pri-miRNA

    private String              note                            ="";
    private String              dbxrefs                         ="";
    
    private String              mirGeneDBGeneName;
    private String              mirGeneDBFamily;
    private String              mirGeneDBNodeOfOriginFamily;
    private String              mirGeneDBNodeOfOriginGene;
    
    private String      				preStrUpperLine1;               // unpaired nt on 5' arm in the pretty plot
    private String      				preStrUpperLine2;               // paired nt on 5' arm in the pretty plot
    private String      				preStrMiddleLine3;              // pairing in the pretty plot
    private String      				preStrLowerLine5;               // unpaired nt on 3' arm in the pretty plot
    private String      				preStrLowerLine4;               // paired nt on 3' arm in the pretty plot

    private HashMap             featureSet;
    /*
     * the following can probably be moved to the featureSet hashMap
     */
    private int                 maxInternalLoopSize             =0;
    private int                 numOfInternalLoops              =0;
    private int                 numOfUnpairedBases              =0;
    private double              fractOfUnpairedBases            =0;

    private int                 lowerStemSize                   =0;
    private int                 upperStemSize                   =0;
    private int                 topStemSize                     =0;
    private int                 upperStart                      =0;
    private int                 upperEnd                        =0;
    
    private int                 terminalLoopStart               =0;
    private int                 terminalLoopEnd                 =0;
    

    
    private long                readCount;

    private ArrayList<MiRNA>                miRNAList;
    private ArrayList<ShortPubMedEntry>     pubmedRefList;
    
    private int index=0;
    
    private static final String             miIDPattern             = "MI[\\d()\\.]{7,}";
    private static final String             miNamePattern           = "[a-z]{3,4}-[a-z]{3}-[\\w()-]{1,}"; //"[a-z]{3,4}-mir-\\S{1,}\\b";

    private static Pattern                  patternMI;
    private static Pattern                  patternName;

    public static final int                 COL_CHR                 = 0;
    public static final int                 COL_SOURCE              = 1;
    public static final int                 COL_FEATURETYPE         = 2;
    public static final int                 COL_STARTPOS            = 3;
    public static final int                 COL_ENDPOS              = 4;
    public static final int                 COL_SCORE               = 5;
    public static final int                 COL_STRAND              = 6;
    public static final int                 COL_FRAME               = 7;
    public static final int                 COL_ATTRIBUTES          = 8;
    
    private static final int                MIRGENEDB_GENE_NAME     = 0;
    private static final int                MIRGENEDB_MIRBASE_ID    = 1;
    private static final int                MIRGENEDB_FAMILY        = 2;
    private static final int                MIRGENEDB_ACC3P         = 3;
    private static final int                MIRGENEDB_ACC5P         = 4;
    private static final int                MIRGENEDB_CHROM         = 5;
    private static final int                MIRGENEDB_START         = 6;
    private static final int                MIRGENEDB_END           = 7;
    private static final int                MIRGENEDB_STRAND        = 8;
    private static final int                MIRGENEDB_GENE_ORIGIN   = 9;
    private static final int                MIRGENEDB_FAMILY_ORIGIN = 10;
    
    public PreMiRNA(){
        super();
        miRNAList           = new ArrayList<>();
        pubmedRefList       = new ArrayList<>();
        featureSet          = new HashMap();
        patternMI           = Pattern.compile(miIDPattern);
        patternName         = Pattern.compile(miNamePattern);
        
    }
    
    
    public PreMiRNA(String id,String seq){
        super(id,seq);
        miRNAList           = new ArrayList<>();
        pubmedRefList       = new ArrayList<>();
        featureSet          = new HashMap();

    }
    
    
    public PreMiRNA(PreMiRNA premirna){
        super(premirna);
        
        this.host = premirna.host;
        this.host3code = premirna.host3code;
        this.miID = premirna.miID;
        this.parent = premirna.parent;
        
        super.setStartPos(premirna.getStartPos());
        super.setEndPos(premirna.getEndPos());
        this.note = premirna.note;
        this.dbxrefs = premirna.dbxrefs;

        this.maxInternalLoopSize = premirna.maxInternalLoopSize;
        this.numOfInternalLoops = premirna.numOfInternalLoops;
        this.numOfUnpairedBases = premirna.numOfUnpairedBases;
        this.fractOfUnpairedBases = premirna.fractOfUnpairedBases;

        this.lowerStemSize = premirna.lowerStemSize;
        this.upperStemSize = premirna.upperStemSize;
        this.topStemSize = premirna.topStemSize;
        this.upperStart = premirna.upperStart;
        this.upperEnd = premirna.upperEnd;
        
        miRNAList           = new ArrayList<>();
        pubmedRefList       = new ArrayList<>();
        this.miRNAList.addAll(premirna.miRNAList);
        this.pubmedRefList.addAll(premirna.pubmedRefList);
        
        featureSet          = new HashMap();
        featureSet = deepCopyHash(premirna.featureSet);    

        this.readCount = premirna.readCount;
        this.index = premirna.index;
                        
    }

    
    
    
    /**
     * copy operation for generic HashMap. This is beyond my skill level to understand, i found
     * the code online. 
     * 
     * @param <K1> HashMap Key
     * @param <K2> HashMap Key
     * @param <V>  Value
     * @param original the HashMap to be copied (I think)
     * @return  the copy of the Hashmap
     */
    public static <K1, K2, V> HashMap<K1, HashMap<K2, V>> deepCopyHash(HashMap<K1, HashMap<K2, V>> original){

        HashMap<K1, HashMap<K2, V>> copy = new HashMap<>();
        for(Entry<K1, HashMap<K2, V>> entry : original.entrySet()){
            copy.put(entry.getKey(), new HashMap<>(entry.getValue()));
        }
        return copy;
    }    
    
    
    
    
    /**
     * parse out an miRBase entry for both the pre-miRNA and the child miRNAs
     * 
     * @param lines ArrayList of Strings containing the relevant GFF entries
     */
    public void parseMiRBaseGFFEntry(ArrayList<String> lines){
        String line = "";
        try{
            line = lines.get(0);
            /*
                chr1            chromosome
                source          n/a here
                miRNA           feature type (n/a)
                start pos
                end pos
                score           n/a here               
                strand          (+/-)
                frame           n/a here
                attributes      e.g. ID=MIMAT0027619;Alias=MIMAT0027619;Name=hsa-miR-6859-3p

            */
            /*
            the original line was encoding the chromosome as a number. I'm not sure why this was necessary
            and whether it breaks anything by changing to a simple string representation
            */
//            this.setChromosome(ChromosomeString.getChromosomeAsString(ChromosomeString.GetChromNumber(line.split("\t")[COL_CHR].trim())));
            this.setChromosome(line.split("\t")[COL_CHR].trim());
            setStartPos(Integer.parseInt(line.split("\t")[COL_STARTPOS].trim()));
            setEndPos(Integer.parseInt(line.split("\t")[COL_ENDPOS].trim()));
            this.setLength(getEndPos()-getStartPos() +1);
            this.setStrand(StrandString.guessStrand(line.split("\t")[COL_STRAND].trim()).toString());

            String attribs[] = line.split("\t")[COL_ATTRIBUTES].split(";");

            for (String attribStr: attribs){
                String attribType = attribStr.split("=")[0].trim();
                String attribValue = attribStr.split("=")[1].trim();
                switch (attribType){
                    case "ID":
                        miID = attribValue;
                        break;

                    case "Alias":                            
                        break;

                    case "Name":
                        this.setName(attribValue);
                        break;

                    default:
                        logger.warn("unknown attribute in parsing pre-miRNA entry \n" + line);
                        logger.info("unknown attribute in parsing pre-miRNA entry \n" + line);
                        break;
                }       
            }
            
            // .. and the child miRNAs
            for(int i=1;i<lines.size();i++){
                MiRNA miRNA = new MiRNA();
                miRNA.parseMiRBaseGFFEntry(lines.get(i));
                this.getMiRNAList().add(miRNA);
            }
        }
        catch(RuntimeException ex){
            logger.error("error parsing miRBase gff reference file pre-miRNA entry \n" + line);
            logger.info("error parsing miRBase gff reference file pre-miRNA entry \n" + line);
            logger.info(ex.toString());
            throw new RuntimeException("error parsing miRBase gff reference file pre-miRNA entry \n" + line);
        }
        
    }
    
     /**
     * parse out an miRBase entry for both the pre-miRNA and the child miRNAs
     * 
     * @param line GFF entry 
     */
    public void parseMiRBaseEntry(String line){
//        String line = "";
        try{
            
            /*
                chr1            chromosome
                source          n/a here
                miRNA           feature type (n/a)
                start pos
                end pos
                score           n/a here               
                strand          (+/-)
                frame           n/a here
                attributes      e.g. ID=MIMAT0027619;Alias=MIMAT0027619;Name=hsa-miR-6859-3p

            */

            this.setChromosome(ChromosomeString.getChromosomeAsString(ChromosomeString.GetChromNumber(line.split("\t")[COL_CHR].trim())));

            setStartPos(Integer.parseInt(line.split("\t")[COL_STARTPOS].trim()));
            setEndPos(Integer.parseInt(line.split("\t")[COL_ENDPOS].trim()));
            this.setLength(getEndPos()-getStartPos() +1);
            this.setStrand(StrandString.guessStrand(line.split("\t")[COL_STRAND].trim()).toString());

            String attribs[] = line.split("\t")[COL_ATTRIBUTES].split(";");

            for (String attribStr: attribs){
                String attribType = attribStr.split("=")[0].trim();
                String attribValue = attribStr.split("=")[1].trim();
                switch (attribType){
                    case "ID":
                        miID = attribValue;
                        break;

                    case "Alias":                            
                        break;

                    case "Name":
                        this.setName(attribValue);
                        break;

                    default:
                        logger.warn("unknown attribute in parsing pre-miRNA entry \n" + line);
                        logger.info("unknown attribute in parsing pre-miRNA entry \n" + line);
                        break;
                }       
            }
            
//            // .. and the child miRNAs
//            for(int i=1;i<lines.size();i++){
//                MiRNA miRNA = new MiRNA();
//                miRNA.parseMiRBaseEntry(lines.get(i));
//                this.getMiRNAList().add(miRNA);
//            }
        }
        catch(RuntimeException ex){
            logger.error("error parsing miRBase gff reference file pre-miRNA entry \n" + line);
            logger.info("error parsing miRBase gff reference file pre-miRNA entry \n" + line);
            logger.info(ex.toString());
            throw new RuntimeException("error parsing miRBase gff reference file pre-miRNA entry \n" + line);
        }
        
    }


    /**
     * Parse out a secondary structure entry from the miRBase STR file
     *
     * There appears to be 8 lines / entry. For example:
     * 
     *  1.  >cel-let-7 (-42.90)   [cel-let-7-5p:17-38] [cel-let-7-3p:60-81]
     *  2.  
     *  3.  ------uaca    gga             U              ---  aaua 
     *  4.            cugu   uccggUGAGGUAG AGGUUGUAUAGUUu   gg    u
     *  5.            ||||   ||||||||||||| ||||||||||||||   ||     
     *  6.            gaca   aggCCAUUCCAUC UUUAACGUAUCaag   cc    u
     *  7.  agcuucucaa    --g             U              ugg  acca 
     *  8.
     * @param lines
     */
    public void parseMiRBaseSStructureEntry(ArrayList<String> lines) {
    	String headerLine = lines.get(0); 

    	logger.debug("\n\n");
    	logger.debug("+" + StringUtils.repeat("*", 78) + "+");
    	logger.debug(headerLine);
    	logger.debug("+" + StringUtils.repeat("*", 78) + "+");

    	
    	String line = "";
      Boolean hasMiRs = this.getMiRNAList()!= null && this.getMiRNAList().size()>0;
      MiRNA miR1;
      MiRNA miR2;
      
      try{
      	this.setEnergy(Float.parseFloat(StringUtils.chop(headerLine.split(" ")[1]).substring(1)));
      	
      	this.preStrUpperLine1  = lines.get(2);
      	this.preStrUpperLine2  = lines.get(3);
      	this.preStrMiddleLine3 = lines.get(4);
      	this.preStrLowerLine4  = lines.get(5);
      	this.preStrLowerLine5  = lines.get(6);
      	
      	String seq = "";
      	
      	// Upper Arm
      	for(int nt=0; nt<preStrUpperLine1.length();nt++) {
      		if(StringUtils.substring(preStrUpperLine1, nt, nt+1).equals("-"))
      			continue;
      		if(StringUtils.substring(preStrUpperLine1, nt, nt+1).equals(" "))
      			seq = seq.concat(preStrUpperLine2.substring(nt, nt+1));
      		else
      			seq = seq.concat(preStrUpperLine1.substring(nt, nt+1));
      	}
      	
        // middle line
      	if(!preStrMiddleLine3.substring(preStrMiddleLine3.length() - 1).equals(" ")) {
      		seq = seq.concat(preStrMiddleLine3.substring(preStrMiddleLine3.length() - 1));
      	}
        // Lower Arm
      	for(int nt=preStrUpperLine1.length()-1;nt>=0;nt--) {
      		if(StringUtils.substring(preStrLowerLine5, nt, nt+1).equals("-"))
      			continue;
      		if(StringUtils.substring(preStrLowerLine5, nt, nt+1).equals(" "))
      			seq = seq.concat(preStrLowerLine4.substring(nt, nt+1));
      		else
      			seq = seq.concat(preStrLowerLine5.substring(nt, nt+1));
      	}
      	
      	this.setSeq(seq);
      	
      	/*
      	 *  parse miRNA information, if needed
      	 *  There can be one or two miRNAs, and if there is one, 
      	 *  it can be on the upper or lower arm
      	 */      	
        if(!hasMiRs) {
         	miR1 = new MiRNA();
        	String n1 = StringUtils.chop(headerLine.split("\\s+")[2]).substring(1).split(":")[0];
        	int start1 = Integer.parseInt(StringUtils.chop(headerLine.split("\\s+")[2]).substring(1).split(":")[1].split("-")[0]);
        	int end1 = Integer.parseInt(StringUtils.chop(headerLine.split("\\s+")[2]).substring(1).split(":")[1].split("-")[1]);
        	miR1.setName(n1);
        	miR1.setSeq(seq.substring(start1-1, end1));
        	this.addProduct(miR1);
        	
        	if(StringUtils.countMatches(headerLine, "\\[")==2) {
	        	miR2 = new MiRNA();
	        	String n2 = StringUtils.chop(headerLine.split("\\s+")[3]).substring(1).split(":")[0];
	        	int start2 = Integer.parseInt(StringUtils.chop(headerLine.split("\\s+")[3]).substring(1).split(":")[1].split("-")[0]);
	        	int end2 = Integer.parseInt(StringUtils.chop(headerLine.split("\\s+")[3]).substring(1).split(":")[1].split("-")[1]);
	        	miR2.setName(n2);
	        	miR2.setSeq(seq.substring(start2-1, end2));
	        	this.addProduct(miR2);
        	}        	
        	
        }
        
      	
      	
      }
      catch(RuntimeException ex){
          logger.error("error parsing miRBase SStructure reference file pre-miRNA entry \n" + line);
          logger.info("error parsing miRBase SStructure reference file pre-miRNA entry \n" + line);
          logger.info(ex.toString());
          throw new RuntimeException("error parsing miRBase SStructure reference file pre-miRNA entry \n" + line);
      }
    	
    }
    
    /**
     * in many cases the hairpin sequence associated with the pre-miRNA is longer than 
     * what is most likely the real pre-miRNA sequence. In some cases this may be a 
     * problem, in which case we want to be more conservative and use a shorter sequence
     * (for example, for searching for duplicate pre-miRNA like sequences within the genome)
     * 
     * In the following hairpin, it would be sufficient to select the sequence 
     * after the first bulge from the miRNA
     * 
     * 
     *                    \/  
     *  1.  ------uaca    gga             U              ---  aaua 
     *  2.            cugu   uccggUGAGGUAG AGGUUGUAUAGUUu   gg    u
     *  3.            ||||   ||||||||||||| ||||||||||||||   ||     
     *  4.            gaca   aggCCAUUCCAUC UUUAACGUAUCaag   cc    u
     *  5.  agcuucucaa    --g             U              ugg  acca 
     *  
     */
    public void trimBackHairpinStructure() {
    	/*
    	 * find miRNA on upper strand (by presence of upper case letters)
    	 * find first bulge upstream
    	 */
    	Boolean miROnUpper = true;
    	Boolean miROnLower = true;
    	if(this.getName().equals("xxxxx"))
    		logger.debug("take a break ");
    	/*
    	 * Upper Arm
    	 */
    	int ntUpper = -1;
    	while( ntUpper<preStrUpperLine1.length()) {
    		ntUpper++;
    		if(ntUpper==preStrUpperLine1.length()) {
    			miROnUpper = false;
    			break;
    		}
    		if(StringUtils.substring(preStrUpperLine1, ntUpper, ntUpper+1).equals("-"))
    			continue;
    		if(StringUtils.substring(preStrUpperLine1, ntUpper, ntUpper+1).equals(" ")) {
    			Character upper2NT = preStrUpperLine2.charAt(ntUpper);
    			if(Character.isUpperCase(upper2NT))
    				break;
    		}else {
    			Character upper1NT = preStrUpperLine1.charAt(ntUpper);
    			if(Character.isUpperCase(upper1NT))
    				break;
    		}
    	}
    	int ntUpperBreak = ntUpper-1;
    	
    	
    	logger.debug("Upper Arm - at NT <" + ntUpperBreak + ">");
    	
    	/* now trim back along preStrUpperLine1 or preStrUpperLine2
    	 * if we didn't find an miRNA on the 5' arm, then we don't do anything for now
    	 * if we are on preStrUpperLine1, we trim back until we hit a " " or "-"
    	 * if we are on preStrUpperLine2, we trim back until we hit a " ", then switch to preStrUpperLine1
    	 * so, we always start from preStrUpperLine2
    	 */
    	if(miROnUpper) {
	    	Character upper2NT = preStrUpperLine2.charAt(ntUpper);
	    	Boolean miRStartsOnUpper2 = Character.isUpperCase(upper2NT);
	
	    	// find first unpaired nt by searching for first " " in 5' direction 
				if(miRStartsOnUpper2) {
					while(ntUpperBreak>=0) {
		    		if(StringUtils.substring(preStrUpperLine2, ntUpperBreak, ntUpperBreak+1).equals(" "))
		    			break;
		    		ntUpperBreak--;
		    	}				
				}
    	
				// Switch to  Upper1 and look for first " " or "-" 
				while(ntUpperBreak>=0) {
	    		if(StringUtils.substring(preStrUpperLine1, ntUpperBreak, ntUpperBreak+1).equals(" "))
	    			break;;
	    		if(StringUtils.substring(preStrUpperLine1, ntUpperBreak, ntUpperBreak+1).equals("-")) {
	    			ntUpperBreak--;
	    			break;
	    		}
	  			ntUpperBreak--;
	    	}				
	    	// finally, need to count how many "-" there are to the end of the sequence
				if(ntUpperBreak>0)
					ntUpperBreak -= StringUtils.countMatches(preStrUpperLine1.substring(0, ntUpperBreak), "-");
    	}
			

    	/*
    	 * Lower Arm
    	 */
    	int ntLower = -1;
    	while( ntLower<preStrLowerLine5.length()) {
    		ntLower++;
    		if(ntLower==preStrLowerLine5.length()) {
    			miROnLower = false;
    			break;
    		}
    		if(StringUtils.substring(preStrLowerLine5, ntLower, ntLower+1).equals("-"))
    			continue;
    		if(StringUtils.substring(preStrLowerLine5, ntLower, ntLower+1).equals(" ")) {
    			Character lower4NT = preStrLowerLine4.charAt(ntLower);
    			if(Character.isUpperCase(lower4NT))
    				break;
    		}else {
    			Character lower5NT = preStrLowerLine5.charAt(ntLower);
    			if(Character.isUpperCase(lower5NT))
    				break;
    		}
    	}

    	/* now trim back along preStrUpperLine1 or preStrUpperLine2
    	 * if we are on preStrUpperLine1, we trim back until we hit a " " or "-"
    	 * if we are on preStrUpperLine2, we trim back until we hit a " ", then switch to preStrUpperLine1
    	 * so, we always start from preStrUpperLine2
    	 */
    	int ntLowerBreak = ntLower-1;
    	if(miROnLower) {
	    	Character lower4NT = preStrLowerLine4.charAt(ntLower);
				if(Character.isUpperCase(lower4NT)) {
					while(ntLowerBreak>=0) {
		    		if(StringUtils.substring(preStrLowerLine4, ntLowerBreak, ntLowerBreak+1).equals(" "))
		    			break;;
		    		ntLowerBreak--;
		    	}				
				}
    	
				while(ntLowerBreak>=0) {
	    		if(StringUtils.substring(preStrLowerLine5, ntLowerBreak, ntLowerBreak+1).equals(" "))
	    			break;;
	    		if(StringUtils.substring(preStrLowerLine5, ntLowerBreak, ntLowerBreak+1).equals("-")) {
	    			break;
	    		}
	    		ntLowerBreak--;
	    	}				
	  		ntLowerBreak++;
				ntLowerBreak -= StringUtils.countMatches(preStrLowerLine5.substring(0, ntLowerBreak), "-");
    	}
			
			if(!miROnUpper) {
				ntUpperBreak = ntLowerBreak-1;
			}
			if(!miROnLower) {
				ntLowerBreak = ntUpperBreak-1;
				if(ntLowerBreak<0) ntLowerBreak=0;
			}
			// take subsequence based on ntUpperBreak and ntLowerBreak
    	logger.debug(preStrUpperLine1);
    	logger.debug(preStrUpperLine2);
    	logger.debug(preStrMiddleLine3);
    	logger.debug(preStrLowerLine4);
    	logger.debug(preStrLowerLine5);
    	logger.debug(this.getSeq());
    	logger.debug(StringUtils.repeat(" ", ntUpperBreak+1) + this.getSeq().substring(ntUpperBreak+1, this.getSeq().length()-ntLowerBreak));
    	logger.debug("done"); 
    }
    

    
    
    /**
     * parses out an pre-miRNA entry line from MirGeneDB (http://mirgenedb.org/)
     * The entries don't appear to be exportable, so need to view a specific species
     * and then select and copy. Based on this, an entry would be tab delimited 
     * with the following fields
     * 
     *   Gene name                  +	
     *   MiRBase ID	
     *   Family                     +   miRNA family
     *   5p accession                   miRNA on the 5' strand
     *   3p accession                   miRNA on the 3' strand
     *   Chromosome	
     *   Start	
     *   End	
     *   Strand	
     *   Node of origin (gene)      +   
     *   Node of origin (family)    +
     * 
     * The entries marked with '+' are miRGeneDB specific and defined at 
     * the pre-miRNA level
     * 
     * e.g. 
     *  Dre-Let-7-P1a	dre-let-7a-2	LET-7	MIMAT0001759	
     *   // None	chr15	20399528	20399594	+	Teleostei	Bilateria
     * 
     * thus, for each line, we will add the child miRNA(s)
     * @param mirLine GFF entry to parse
     * 
     */
    public void parseMirGeneDBEntry(String mirLine){
        
        this.setMirGeneDBGeneName(mirLine.split("\t")[PreMiRNA.MIRGENEDB_GENE_NAME].trim());
        this.setName(mirLine.split("\t")[PreMiRNA.MIRGENEDB_MIRBASE_ID].trim());
        this.setMirGeneDBFamily(mirLine.split("\t")[PreMiRNA.MIRGENEDB_FAMILY].trim());
        String miRNA5pID = mirLine.split("\t")[PreMiRNA.MIRGENEDB_ACC5P].trim();
        String miRNA3pID = mirLine.split("\t")[PreMiRNA.MIRGENEDB_ACC3P].trim();
        this.setChromosome(mirLine.split("\t")[PreMiRNA.MIRGENEDB_CHROM].trim());
        try{
            if(mirLine.split("\t")[PreMiRNA.MIRGENEDB_START].trim().equals("None"))
                this.setStartPos(-1);
            else
                this.setStartPos(Integer.parseInt(mirLine.split("\t")[PreMiRNA.MIRGENEDB_START].trim()));
            
            if(mirLine.split("\t")[PreMiRNA.MIRGENEDB_START].trim().equals("None"))
                this.setEndPos(-1);
            else
                this.setEndPos(Integer.parseInt(mirLine.split("\t")[PreMiRNA.MIRGENEDB_END].trim()));
            
        }
        catch(NumberFormatException ioNB){
            logger.info("error parsing start/stop position in MirGeneDBEntry");
            logger.info(mirLine.split("\t")[PreMiRNA.MIRGENEDB_START].trim() + "/" 
                    + mirLine.split("\t")[PreMiRNA.MIRGENEDB_END].trim());
            
            throw new NumberFormatException("error parsing start/stop position in MirGeneDBEntry");
        }
        this.setStrand(mirLine.split("\t")[PreMiRNA.MIRGENEDB_STRAND].trim());
        this.setMirGeneDBNodeOfOriginFamily(mirLine.split("\t")[PreMiRNA.MIRGENEDB_FAMILY_ORIGIN].trim());
        this.setMirGeneDBNodeOfOriginGene(mirLine.split("\t")[PreMiRNA.MIRGENEDB_GENE_ORIGIN].trim());
        
        if(miRNA5pID.equals("None")==false){
            MiRNA miRNA5p = new MiRNA(miRNA5pID, this.getChromosome(), 0, 0, Strand.PLUS.toString());
            miRNA5p.setMimatID(miRNA5pID);
            this.addProduct(miRNA5p);            
        }
        
        if(miRNA3pID.equals("None")==false){
            MiRNA miRNA3p = new MiRNA(miRNA3pID, this.getChromosome(), 0, 0, Strand.MINUS.toString());
            miRNA3p.setMimatID(miRNA3pID);
            this.addProduct(miRNA3p);
        }
    }
    
    
    /**
     * check there is structure information for this pre-miRNA
     * the structure information should be specified for all 
     * five structure lines, so we only check the first.
     * 
     * @return
     */
    public Boolean hasStructureInfo() {
    	return preStrUpperLine1!=null;
    }
    /**
     * Format is
     * >cel-miR-86-3p MIMAT0020776 Caenorhabditis elegans miR-86-3p
     * 
     * i.e. long miR name|MIMAT ID|full organism name|short miR name
     * 
     * Full organism name is abbreviated to the first two words
     * e.g. 'Kaposi sarcoma-associated herpesvirus' --> 'Kaposi sarcoma-associated'
     * @param headerLine
     */
    public void parseMiRBaseFAHeaderLine(String headerLine) {
      patternName           = Pattern.compile(miNamePattern);
      patternMI             = Pattern.compile(miIDPattern);
      logger.debug(INDENT1 + headerLine);
      
      Matcher matcherName = patternName.matcher(headerLine);

      if(matcherName.find()==false){
        logger.info("oops");
      }
    	this.setName(matcherName.group(0));


      Matcher matcherMI = patternMI.matcher(headerLine);
      matcherMI.find();      
      if(matcherName.group(0)==null){
        logger.info("oops");
      }
    	this.setMimatID(matcherMI.group(0));
      logger.debug(INDENT2 + "-->" +  matcherMI.group(0));
      
   }
    
    /**
     * parse out miRBase name from header line of FastA entry and return
     * @param headerLine
     * @return
     *     private static final String             miIDPattern             = "MI\\d{7}";
    private static final String             miNamePattern           = "[a-z]{3,4}-mir-\\S{1,}\\b";
    //private static Pattern                  patternMI;
    //private static Matcher                  matcherMI; 
    private static Pattern  patternMI;

     */
    public static String parseNameFromMiRBaseFAHeaderLine(String headerLine) {
      patternName           = Pattern.compile(miNamePattern);
      Matcher matcherName = patternMI.matcher(headerLine);
      matcherName.find(); 
      return matcherName.group(0);
      
    	//return headerLine.split(" ")[0].substring(1).trim();
      // private static final String             miIDPattern             = "MI\\d{7}";
      // private static final String             miNamePattern           = "[a-z]{3,4}-mir-\\S{1,}\\b";

    }
    
    
    /**
     * parse out MIMATID from header line of FastA entry and return
     * @param headerLine
     * @return
     */
    public static String parseIDFromMiRBaseFAHeaderLine(String headerLine) {
      patternMI           = Pattern.compile(miIDPattern);
      Matcher matcherMI = patternMI.matcher(headerLine);
      matcherMI.find();      
      return matcherMI.group(0);
    }
    
    
    /**
     * parse out miRBase name from header line of Secondary Structure entry and return
     * @param headerLine
     * @return
     */
    public static String parseNameFromMiRBaseSSHeaderLine(String headerLine) {
    	return headerLine.split(" ")[0].substring(1).trim();
    }
    
    /**
     * count all reads mapped to this miRNA. These will always be mapped to an miRNA
     * 
     * @return total number of counts
     * 
     */
    public long countReads(){
        setReadCount(0);
        this.getMiRNAList().stream().forEach((miRNA) -> {
            setReadCount(getReadCount() + miRNA.getReadCount());
        });
        return getReadCount();
        
    }
    
    
    /*
    return sequence in FASTA format
    */
    @Override
    public String toFastA(){
        return ">" + this.getName() + "|" + this.getMiID() + System.lineSeparator()
                + this.getSeq() + System.lineSeparator();
    }
    
    
    /**
     * provides a summary of how the mapped reads are distributed across the
     * pre-miRNA
     * 
     * We want a list of start, stop, read count
     * this is something that we can plot in a pretty format
     * 
     * @return String pretty plot of the pre-miRNA
     */
    public String prettyPlotReadsNoSeq(){
        String plus = "+" + StringUtils.repeat("-", getEndPos()-this.getStartPos()+1) + "+\n";
        String minus = "";

        for(MiRNA miRNA: this.getMiRNAList()){
            if(miRNA.getStrand().equals(Strand.PLUS.toString())){
                String thisMiRLine = StringUtils.repeat(" ", miRNA.getStartPos()-this.getStartPos()) 
                        + StringUtils.repeat("*", miRNA.getEndPos()-miRNA.getStartPos()+1) 
                        + StringUtils.repeat(" ", 
                                getEndPos()-this.getStartPos()+1 
                                        - (miRNA.getEndPos()-this.getStartPos()+1)                                         
                                        + 10) 
                        + miRNA.getReadCount() + "\n";
                plus = plus.concat(thisMiRLine);
            }else{
                String thisMiRLine = StringUtils.repeat(" ", miRNA.getStartPos()-this.getStartPos()) 
                        + StringUtils.repeat("*", miRNA.getEndPos()-miRNA.getStartPos()+1) 
                        + StringUtils.repeat(" ", 
                                getEndPos()-this.getStartPos()+1 
                                        - (miRNA.getEndPos()-this.getStartPos()+1)                                         
                                        + 10) 
                        + miRNA.getReadCount() + "\n";
                minus = minus.concat(thisMiRLine);
            }
        }
        minus = minus.concat("+" + StringUtils.repeat("-", getEndPos()-this.getStartPos()+1) + "+\n");
        
        return plus + "+" + StringUtils.repeat("-", getEndPos()-this.getStartPos()+1) + "+\n" + minus;
    }
 
    

    
    /**
     * provides a summary of how the mapped reads are distributed across the
     * pre-miRNA
     * 
     * We want a list of start, stop, read count
     * this is something that we can plot in a pretty format
     * 
     * @return String pretty plot of reads including the sequence
     */
    public String prettyPlotReadsWithSeq(){
        String plus = "+" + StringUtils.repeat("-", getEndPos()-this.getStartPos()+1) + "+\n";
        String minus = "";

        if(this.getSeq() == null) 
            return "no sequence defined for this pre-miRNA";
        
        if(this.getSeq().isEmpty())
            return "no sequence defined for this pre-miRNA";
        
        for(MiRNA miRNA: this.getMiRNAList()){
            /*
            
            */
            String seq ;
            int preStart;
            int preEnd;
            if(miRNA.getStartPos()-this.getStartPos()-1<0)
                preStart = 0;
            else
                preStart = miRNA.getStartPos()-this.getStartPos()-1;
            
            if(miRNA.getEndPos()-this.getStartPos()-1 > this.getSeq().length()-1){
                preEnd = this.getSeq().length()-1;
            }else
                preEnd = miRNA.getEndPos()-this.getStartPos()-1;
                        
            // + this.getSequence(miRNA.getStartPos()-this.getStartPos()-1, miRNA.getEndPos()-this.getStartPos()-1)    
            if(miRNA.getStrand().equals(Strand.PLUS.toString())){
                String thisMiRLine = StringUtils.repeat(" ", miRNA.getStartPos()-this.getStartPos()) 
                        + this.getSequence(preStart, preEnd)
                        + StringUtils.repeat(" ", 
                                getEndPos()-this.getStartPos()+1 
                                        - (miRNA.getEndPos()-this.getStartPos()+1)                                         
                                        + 10) 
                        + miRNA.getReadCount() + " ("+ (miRNA.getEndPos()-miRNA.getStartPos()+1) + ")\n";
                plus = plus.concat(thisMiRLine);
            }else{
                String thisMiRLine = StringUtils.repeat(" ", miRNA.getStartPos()-this.getStartPos()) 
//                        + this.getSequence(miRNA.getStartPos()-this.getStartPos()-1, miRNA.getEndPos()-this.getStartPos()-1)
                        + this.getSequence(preStart, preEnd)
                        + StringUtils.repeat(" ", 
                                getEndPos()-this.getStartPos()+1 
                                        - (miRNA.getEndPos()-this.getStartPos()+1)                                         
                                        + 10) 
                        + miRNA.getReadCount() + " ("+ (miRNA.getEndPos()-miRNA.getStartPos()+1) + ") " + miRNA.getName() + "\n";
                minus = minus.concat(thisMiRLine);
            }
        }
        minus = minus.concat("+" + StringUtils.repeat("-", getEndPos()-this.getStartPos()+1) + "+\n");
        
        return plus + "+" + this.getSeq() + "+\n" + minus;
    }
 
    

    
    /**
     * find the range of start and stop positions of mapped reads
     * in the pre-miRNA
     * 
     * @return String summarizes the distribution of read starts along the pre-miRNA
     */
    public  String summarizeStartPositions(){
        
        long veryBigInt = 1000000000000L;
        FeatureLocationRange featureLocRange = new FeatureLocationRange();
        featureLocRange.setMinusMinStart(veryBigInt);
        featureLocRange.setMinusMaxStart(-1);
        featureLocRange.setMinusMinEnd(veryBigInt);
        featureLocRange.setMinusMaxEnd(-1);
        featureLocRange.setPlusMinStart(veryBigInt);
        featureLocRange.setPlusMaxStart(-1);
        featureLocRange.setPlusMinEnd(veryBigInt);
        featureLocRange.setPlusMaxEnd(-1);
        
        int plusSum = 0;
        int minusSum = 0;
        
        for(MiRNA miRNA: this.getMiRNAList()){
            if(miRNA.getStrand().equals(Strand.PLUS.toString())){
                plusSum += (miRNA.getStartPos()-this.getStartPos())*miRNA.getReadCount();
                if(miRNA.getStartPos()-this.getStartPos()<featureLocRange.getPlusMinStart())
                    featureLocRange.setPlusMinStart(miRNA.getStartPos()-this.getStartPos());
                if(miRNA.getStartPos()-this.getStartPos()>featureLocRange.getPlusMaxStart())
                    featureLocRange.setPlusMaxStart(miRNA.getStartPos()-this.getStartPos());
                if(miRNA.getEndPos()-this.getStartPos()<featureLocRange.getPlusMinEnd())
                    featureLocRange.setPlusMinEnd(miRNA.getEndPos()-this.getStartPos());
                if(miRNA.getEndPos()-getEndPos()>featureLocRange.getPlusMaxEnd())
                    featureLocRange.setPlusMaxEnd(miRNA.getEndPos()-this.getStartPos());
            }else{
                minusSum += (miRNA.getStartPos()-this.getStartPos())*miRNA.getReadCount();
                if(miRNA.getStartPos()-this.getStartPos()<featureLocRange.getMinusMinStart())
                    featureLocRange.setMinusMinStart(miRNA.getStartPos()-this.getStartPos());
                if(miRNA.getStartPos()-this.getStartPos()>featureLocRange.getMinusMaxStart())
                    featureLocRange.setMinusMaxStart(miRNA.getStartPos()-this.getStartPos());
                if(miRNA.getEndPos()-this.getStartPos()<featureLocRange.getMinusMinEnd())
                    featureLocRange.setMinusMinEnd(miRNA.getEndPos()-this.getStartPos());
                if(miRNA.getEndPos()-this.getStartPos()>featureLocRange.getMinusMaxEnd())
                    featureLocRange.setMinusMaxEnd(miRNA.getEndPos()-this.getStartPos());                    
            }
        }
        
        
        if(featureLocRange.getPlusMinStart()==veryBigInt){
            featureLocRange.setPlusMinStart(-1);
            featureLocRange.setPlusMaxStart(-1);
            featureLocRange.setPlusAverageStart(-1);
        }else{
            if(this.countReads()==0)
                featureLocRange.setPlusAverageStart(-1);
            else
                featureLocRange.setPlusAverageStart((double)plusSum/(double)this.countReads());            
        }

        if(featureLocRange.getMinusMinStart()==veryBigInt){
            featureLocRange.setMinusMinStart(-1);
            featureLocRange.setMinusMaxStart(-1);
            featureLocRange.setMinusAverageStart(-1);
        }else{
            if(this.countReads()==0)
                featureLocRange.setMinusAverageStart(-1);
            else
                featureLocRange.setMinusAverageStart((double)minusSum/(double)this.countReads());            
        }
            
        DecimalFormat df = new DecimalFormat("#.###");
        df.setRoundingMode(RoundingMode.CEILING);
        
        return featureLocRange.getPlusMinStart() + "\t"
                + featureLocRange.getPlusMaxStart() + "\t"
                + df.format(featureLocRange.getPlusAverageStart()) + "\t"
                + featureLocRange.getMinusMinStart() + "\t"
                + featureLocRange.getMinusMaxStart() + "\t"
                + df.format(featureLocRange.getMinusAverageStart()) 
                ;
    }

    

    
    /**
     * return a breakdown of the read length distribution
     * 
     * @return String summarizing the read length distribution
     * @param minLen - shortest length to report
     * @param maxLen - longest length to report
     */
    public Frequency detailStartRange(int minLen, int maxLen){
        
        Frequency startPosFreqs = new Frequency();
        for(int l=minLen; l<maxLen; l++){
            startPosFreqs.addValue(l);
            startPosFreqs.incrementValue(l, -1);
        }
        for(MiRNA miRNA: this.getMiRNAList()){
            if(miRNA.getStartPos()-this.getStartPos()>=minLen 
                    && miRNA.getStartPos()-this.getStartPos()<=maxLen){
                startPosFreqs.incrementValue(miRNA.getStartPos()-this.getStartPos(), miRNA.getReadCount());
            }
        }
        return startPosFreqs;
        
    }
    
    
    
    
    /**
     * return a breakdown of the read length distribution on the 5' side
     * 
     * @return String summary of read length distribution on the 5' side
     * @param minLen - shortest length to report
     * @param maxLen - longest length to report
     */
    public Frequency detailStartRangeFiveArm(int minLen, int maxLen){
        
        Frequency upperStartPosFreqs = new Frequency();
        for(int l=minLen; l<maxLen; l++){
            upperStartPosFreqs.addValue(l);
            upperStartPosFreqs.incrementValue(l, -1);
        }
        for(MiRNA miRNA: this.getMiRNAList()){
            if(((miRNA.getStartPos()-this.getStartPos()) < this.getLength()/2) && miRNA.getStartPos()-this.getStartPos()>=minLen 
                    && miRNA.getStartPos()-this.getStartPos()<=maxLen){
                upperStartPosFreqs.incrementValue(miRNA.getStartPos()-this.getStartPos(), miRNA.getReadCount());
            }
        }
        return upperStartPosFreqs;
        
    }
    
    
    
    
    /**
     * return a breakdown of the read length distribution
     * 
     * @return Frequency read length distribution
     * @param minLen - shortest length to report
     * @param maxLen - longest length to report
     */
    public Frequency detailStartRangeThreeArm(int minLen, int maxLen){
        
        Frequency lowerStartPosFreqs = new Frequency();
        for(int l=minLen; l<maxLen; l++){
            lowerStartPosFreqs.addValue(l);
            lowerStartPosFreqs.incrementValue(l, -1);
        }
        for(MiRNA miRNA: this.getMiRNAList()){
            if(((miRNA.getStartPos()-this.getStartPos()) >= this.getLength()/2) && miRNA.getStartPos()-this.getStartPos()>=minLen 
                    && miRNA.getStartPos()-this.getStartPos()<=maxLen){
                lowerStartPosFreqs.incrementValue(miRNA.getStartPos()-this.getStartPos(), miRNA.getReadCount());
            }
        }
        return lowerStartPosFreqs;
        
    }
    
    
    
    
    /**
     * find the length range of reads mapped to this pre-miRNA
     * 
     * @return String summarizing length range of reads mapped to this pre-miRNA
     */
    public String summarizeLengthRange(){
        int veryBigInt = 1000000;
        int plusMinReadLength = veryBigInt;
        int plusMaxReadLength = -1;
        int plusReadSum = 0;
        int minusMinReadLength = veryBigInt;
        int minusMaxReadLength = -1;
        int minusReadSum = 0;
        for(MiRNA miRNA: this.getMiRNAList()){
            if(miRNA.getStrand().equals(Strand.PLUS.toString())){
                plusReadSum += (miRNA.getEndPos()-miRNA.getStartPos()+1)*miRNA.getReadCount();
                if(miRNA.getEndPos()-miRNA.getStartPos()+1<plusMinReadLength)
                    plusMinReadLength = miRNA.getEndPos()-miRNA.getStartPos()+1;
                if(miRNA.getEndPos()-miRNA.getStartPos()+1>plusMaxReadLength)
                    plusMaxReadLength = miRNA.getEndPos()-miRNA.getStartPos()+1;
            }else{
                minusReadSum += (miRNA.getEndPos()-miRNA.getStartPos()+1)*miRNA.getReadCount();
                if(miRNA.getEndPos()-miRNA.getStartPos()+1<minusMinReadLength)
                    minusMinReadLength = miRNA.getEndPos()-miRNA.getStartPos()+1;
                if(miRNA.getEndPos()-miRNA.getStartPos()+1>minusMaxReadLength)
                    minusMaxReadLength = miRNA.getEndPos()-miRNA.getStartPos()+1;
            }
        }
        double averagePlusReadLength = -1;
        if(plusMinReadLength==veryBigInt){
            plusMinReadLength = -1;            
        }else  
            if(this.countReads()==0)
                averagePlusReadLength = -1;
            else            
                averagePlusReadLength = (double)plusReadSum/(double)this.countReads();
        
        
        double averageMinusReadLength = -1.0;
        if(minusMinReadLength==veryBigInt){
            minusMinReadLength = -1;            
        }else            
            if(this.countReads()==0)
                averageMinusReadLength = -1.0;
            else            
                averageMinusReadLength = (double)minusReadSum/(double)this.countReads();
        
        DecimalFormat df = new DecimalFormat("#.###");
        df.setRoundingMode(RoundingMode.CEILING);
        return plusMinReadLength + "\t" + plusMaxReadLength + "\t" + df.format(averagePlusReadLength) + "\t"
                + minusMinReadLength + "\t" + minusMaxReadLength + "\t" + df.format(averageMinusReadLength)
                ;
    }
    
    
    
    /**
     * return a breakdown of the read length distribution
     * 
     * @return String summarizing read length distribution for this pre-miRNA
     * @param minLen - shortest length to report
     * @param maxLen - longest length to report
     */
    public Frequency detailLengthRange(int minLen, int maxLen){
        
        Frequency lenFreqs = new Frequency();
        for(int l=minLen; l<maxLen; l++){
            lenFreqs.addValue(l);
            lenFreqs.incrementValue(l, -1);
        }
        for(MiRNA miRNA: this.getMiRNAList()){
            if((miRNA.getEndPos()-miRNA.getStartPos()+1)>=minLen 
                    && (miRNA.getEndPos()-miRNA.getStartPos()+1)<=maxLen){
                lenFreqs.incrementValue((miRNA.getEndPos()-miRNA.getStartPos()+1), miRNA.getReadCount());
            }
        }
        return lenFreqs;
    }
    
    
    
    
    /**
     * return a breakdown of the read length distribution along the upper (5') arm
     * of this pre-miRNA
     * 
     * @return String detailing read length distribution along the upper (5') arm
     *               of this pre-miRNA
     * @param minLen - shortest length to report
     * @param maxLen - longest length to report
     */
    public Frequency detailLengthRangeFiveArm(int minLen, int maxLen){
        
        Frequency fiveLenFreqs = new Frequency();
        for(int l=minLen; l<maxLen; l++){
            fiveLenFreqs.addValue(l);
            fiveLenFreqs.incrementValue(l, -1);
        }
        for(MiRNA miRNA: this.getMiRNAList()){
            if(((miRNA.getStartPos()-this.getStartPos())< this.getLength()/2) && (miRNA.getEndPos()-miRNA.getStartPos()+1)>=minLen 
                    && (miRNA.getEndPos()-miRNA.getStartPos()+1)<=maxLen){
                fiveLenFreqs.incrementValue((miRNA.getEndPos()-miRNA.getStartPos()+1), miRNA.getReadCount());
            }
        }
        return fiveLenFreqs;
    }
    
    
    
    
    /**
     * return a breakdown of the read length distribution of the lower (3') arm 
     * of this pre-miRNA
     * 
     * @return String detailing he read length distribution of the lower (3') arm 
     * of this pre-miRNA
     * @param minLen - shortest length to report
     * @param maxLen - longest length to report
     */
    public Frequency detailLengthRangeThreeArm(int minLen, int maxLen){
        
        Frequency threeLenFreqs = new Frequency();
        for(int l=minLen; l<maxLen; l++){
            threeLenFreqs.addValue(l);
            threeLenFreqs.incrementValue(l, -1);
        }
        for(MiRNA miRNA: this.getMiRNAList()){
            if(((miRNA.getStartPos()-this.getStartPos())>= this.getLength()/2)&& (miRNA.getEndPos()-miRNA.getStartPos()+1)>=minLen 
                    && (miRNA.getEndPos()-miRNA.getStartPos()+1)<=maxLen){
                threeLenFreqs.incrementValue((miRNA.getEndPos()-miRNA.getStartPos()+1), miRNA.getReadCount());
            }
        }
        return threeLenFreqs;
    }
    
    
    
    
    /**
     * add a child miRNA that is associated with this pre-miRNA
     * 
     * @param miRNA MiRNA
     */
    public void addProduct(MiRNA miRNA){
        getMiRNAList().add(miRNA);
        index+=1;
    }

    
    
    
    /**
     * append the supplied string to the note field
     * 
     * @param deltaNote String
     * @return String the updated note
     */
    public String appendNote(String deltaNote){
        this.setNote(this.getNote().concat(deltaNote));
        return this.getNote();
    }
    
    
    
    
    /**
     * add pubmed reference to this pre-miRNA
     * 
     * @param newRef ShortPubMedEntry
     */
    public void addPubmedRef(ShortPubMedEntry newRef){
        getPubmedRefList().add(newRef);
    }

    /**
     * return number of miRNAs associated with this pre-miRNA
     * 
     * @return int number of miRNAs associated with this pre-miRNA
     */
    public int SizeOfProduct(){
        return getMiRNAList().size();
    }

    
    /**
     * because the list of features we might calculate can change, we also store
     * in a hash map to make it simpler to collectively pass information
     */
    public void buildFeatureSet(){

        
        getFeatureSet().put("preRNA_sequence",               this.getSeq());
        getFeatureSet().put("preRNA_structure",              this.getStructureStr());
        getFeatureSet().put("preRNA_energy",                 this.getEnergy());
        getFeatureSet().put("preRNA_size",                   this.getLength());
        getFeatureSet().put("preRNA_GC_content",             super.GCfraction() );
        getFeatureSet().put("preRNA_A_content",              this.Afraction());
        getFeatureSet().put("preRNA_U_content",              this.Ufraction());
        getFeatureSet().put("preRNA_G_content",              this.Gfraction());
        getFeatureSet().put("preRNA_C_content",              this.Cfraction());
    }

    public HashMap getFeatureSet(){
        return featureSet;
    }


    public void setUpperStart(int upperStart) {
        this.upperStart = upperStart;
    }

    /**
     * @return the upperStart
     */
    public int getUpperStart() {
        return upperStart;
    }

    /**
     * @return the upperEnd
     */
    public int getUpperEnd() {
        return upperEnd;
    }

    /**
     * @param upperEnd the upperEnd to set
     */
    public void setUpperEnd(int upperEnd) {
        this.upperEnd = upperEnd;
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

    /**
     * @return the host3code
     */
    public String getHost3code() {
        return host3code;
    }

    /**
     * @param host3code the host3code to set
     */
    public void setHost3Lettercode(String host3code) {
        this.host3code = host3code;
    }

    /**
     * @return the notes
     */
    public String getNote() {
        return note;
    }

    /**
     * @param notes the notes to set
     */
    public void setNote(String notes) {
        this.note = notes;
    }

    /**
     * @return the dbxrefs
     */
    public String getDbxrefs() {
        return dbxrefs;
    }

    /**
     * @param dbxrefs the dbxrefs to set
     */
    public void setDbxrefs(String dbxrefs) {
        this.dbxrefs = dbxrefs;
    }

    /**
     * @return the mimatID
     */
    public String getMiID() {
        return miID;
    }

    /**
     * @param mimatID the mimatID to set
     */
    public void setMimatID(String mimatID) {
        this.miID = mimatID;
    }

    /**
     * @return the parent
     */
    public String getParent() {
        return parent;
    }

    /**
     * @param parent the parent to set
     */
    public void setParent(String parent) {
        this.parent = parent;
    }

    /**
     * @return the miRNAList
     */
    public ArrayList<MiRNA> getMiRNAList() {
        return miRNAList;
    }

    /**
     * @return the readCount
     */
    public long getReadCount() {
        return readCount;
    }

    /**
     * @param readCount the readCount to set
     */
    public void setReadCount(long readCount) {
        this.readCount = readCount;
    }


    /**
     * add to the total read count
     * 
     * @param  dReadCount long
     */
    public void addToReadCount(long dReadCount){
        this.setReadCount(this.getReadCount() + dReadCount);
    }

    /**
     * @return the pubmedRefList
     */
    public ArrayList<ShortPubMedEntry> getPubmedRefList() {
        return pubmedRefList;
    }

    /**
     * @return the mirGeneDBGeneName
     */
    public String getMirGeneDBGeneName() {
        return mirGeneDBGeneName;
    }

    /**
     * @param mirGeneDBGeneName the mirGeneDBGeneName to set
     */
    public void setMirGeneDBGeneName(String mirGeneDBGeneName) {
        this.mirGeneDBGeneName = mirGeneDBGeneName;
    }

    /**
     * @return the mirGeneDBFamily
     */
    public String getMirGeneDBFamily() {
        return mirGeneDBFamily;
    }

    /**
     * @param mirGeneDBFamily the mirGeneDBFamily to set
     */
    public void setMirGeneDBFamily(String mirGeneDBFamily) {
        this.mirGeneDBFamily = mirGeneDBFamily;
    }

    /**
     * @return the mirGeneDBNodeOfOriginFamily
     */
    public String getMirGeneDBNodeOfOriginFamily() {
        return mirGeneDBNodeOfOriginFamily;
    }

    /**
     * @param mirGeneDBNodeOfOriginFamily the mirGeneDBNodeOfOriginFamily to set
     */
    public void setMirGeneDBNodeOfOriginFamily(String mirGeneDBNodeOfOriginFamily) {
        this.mirGeneDBNodeOfOriginFamily = mirGeneDBNodeOfOriginFamily;
    }

    /**
     * @return the mirGeneDBNodeOfOriginGene
     */
    public String getMirGeneDBNodeOfOriginGene() {
        return mirGeneDBNodeOfOriginGene;
    }

    /**
     * @param mirGeneDBNodeOfOriginGene the mirGeneDBNodeOfOriginGene to set
     */
    public void setMirGeneDBNodeOfOriginGene(String mirGeneDBNodeOfOriginGene) {
        this.mirGeneDBNodeOfOriginGene = mirGeneDBNodeOfOriginGene;
    }

    /**
     * @return the terminalLoopStart
     */
    public int getTerminalLoopStart() {
        return terminalLoopStart;
    }

    /**
     * @param terminalLoopStart the terminalLoopStart to set
     */
    public void setTerminalLoopStart(int terminalLoopStart) {
        this.terminalLoopStart = terminalLoopStart;
    }

    /**
     * @return the terminalLoopEnd
     */
    public int getTerminalLoopEnd() {
        return terminalLoopEnd;
    }

    /**
     * @param terminalLoopEnd the terminalLoopEnd to set
     */
    public void setTerminalLoopEnd(int terminalLoopEnd) {
        this.terminalLoopEnd = terminalLoopEnd;
    }


}

