/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mirna;

import java.util.ArrayList;
import java.util.HashMap;
import no.uio.medisin.bag.core.sequence.ChromosomeString;
import no.uio.medisin.bag.core.sequence.SimpleRNASequence;
import no.uio.medisin.bag.core.sequence.SimpleSeq;
import no.uio.medisin.bag.core.sequence.StrandString;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * stores an miRNA entry, while this may seem a cumbersome implementation, the most
 * miRNAs found in a mammalian genome is less than 3000, so we shouldn't be too worried 
 * about overheads
 * 
 * This needs to be made comparable so the entries stored in a list can be sorted
 * to speed things up
 * 
 * Added the ability to store the predicted secondary structure for the miRNA
 * in the duplex
 * 
 * other parameters pertaining to the properties of the secondary structure are
 * stored in a hashTable as the desired properties may change according to need
 * 
 * 
 * @author sr
 */
public class MiRNA extends SimpleRNASequence{

    static Logger               logger                      = LogManager.getLogger();
    
    private static final int        QUALIFIER_START_COL     =   21;
    private static final String     QUALIFIER_ACCESSION     =   "/accession";
    private static final String     QUALIFIER_PRODUCT       =   "/product";
    private static final String     QUALIFIER_EVIDENCE      =   "/evidence";
    private static final String     QUALIFIER_EXPERIMENT    =   "/experiment";
    /**
     * parses out an miRNA entry line from MirGeneDB (http://mirgenedb.org/)
     * The entries don't appear to be exportable, so need to view a specific species
     * and then select and copy. Based on this, an entry would be tab delimited 
     * with the following fields
     * 
     * Gene name	
     * MiRBase ID
     * Family	
     * 5p accession	
     * 3p accession	
     * Chromosome	
     * Start	
     * End	
     * Strand	
     * Node of origin (gene)	
     * Node of origin (family)

     * 
     * @param mirLine 
     */
    
    

    private String              mimatID;
    private String              mimatAlias;
    private String              note;
    //private String              name;
    private String              parent;

    //private String              chromosome;
    //private int                 start;
    //private int                 end;
    //private String              strand;
    private String              isomiRString;
    private HashMap             featureSet = new HashMap();     // stores additional characteristics of this miRNA
    
    private String              evidence="";
    private String              references="";                     // string of ";" delimited pubmed IDs
    
    /*
    private String              structLine1;
    private String              structLine2;
    private String              structLine3;
    private String              structLine4;
    */
    private String              host;
    
    private int                 largestInternalBulge;
    private int                 numOfBulges;
    private int                 numOfUnpairedBases;
    private double              fractOfUnpairedBases;

    private char                firstBase;
    private double              stability;
    private char                dangleBaseOne;
    private char                dangleBaseTwo;

    private int                 miStart;
    private int                 miEnd;

    private long                readCount;
    
    // for storing isomiR information
    //(name + ";" + start + ";" + cigar + ";" + md + ";" + seq + "\t")
    
    private static final int    NAMECOL    = 0;
    private static final int    STARTCOL   = 1;
    private static final int    CIGARCOL   = 2;
    private static final int    MDCOL      = 3;
    private static final int    SEQCOL     = 4;
    
    

    public MiRNA(){
        super();
    }
    
    /**
     * Constructor specifying name, location and sequence of an miRNA
     * (ideal for predicted miRNAs)
     * 
     * @param cName     String      name
     * @param cChrom    String      chromosome
     * @param cStart    int         start position
     * @param cEnd      int         end position
     * @param cStrand   String      Strand
     * 
     */
    public MiRNA(String cName, String cChrom, int cStart, int cEnd, String cStrand){
        this(cName, cChrom, cStart, cEnd, cStrand, "", "", "");
    }
    
    
    
    /**
     * Constructor specifying name, location and sequence of an miRNA
     * (ideal for miRBase entry miRNAs)
     * 
     * @param cName     String      name
     * @param cChrom    String      chromosome
     * @param cStart    int         start position
     * @param cEnd      int         end position
     * @param cStrand   String      Strand
     * @param cMIMAT    String      MIMAT (miRBase) ID
     * @param cParentID String      MI (miRBase) Parent ID
     * @param cSeq      String      sequence
     * 
     */
    public MiRNA(String cName, String cChrom, int cStart, int cEnd, String cStrand, String cMIMAT, String cParentID, String cSeq){
        super.setName(cName);
        super.setChromosome(cChrom);
        super.setStartPos(cStart);
        super.setEndPos(cEnd);
        super.setStrand(cStrand);
        super.setSeq(cSeq);
        parent = cParentID;
        mimatID =cMIMAT;
        isomiRString = "";
        
    }
    
    
    
    /**
     * copy constructor
     * 
     * @param mrnaf MiRNA
     */
    public MiRNA(MiRNA mrnaf){
        
        super(mrnaf);
        //name            = mrnaf.name;
        super.setChromosome(mrnaf.getChromosome());
        super.setStartPos(mrnaf.getStartPos());
        super.setEndPos(mrnaf.getEndPos());
        super.setStrand(mrnaf.getStrand());
        super.setSeq(mrnaf.getSeq());
        parent          = mrnaf.parent;
        mimatID         = mrnaf.mimatID;
        isomiRString    = mrnaf.isomiRString;
        
    }
    
    
    
    /**
     * parse out an miRBase entry for both the pre-miRNA and the child miRNAs
     * 
     * @param line String - a single line in a miRBase GFF file
     * 
     */
    public void parseMiRBaseGFFEntry(String line){

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
            super.setChromosome(ChromosomeString.getChromosomeAsString(ChromosomeString.GetChromNumber(line.split("\t")[PreMiRNA.COL_CHR].trim())));
            super.setStartPos(Integer.parseInt(line.split("\t")[PreMiRNA.COL_STARTPOS].trim()));
            super.setEndPos(Integer.parseInt(line.split("\t")[PreMiRNA.COL_ENDPOS].trim()));
            this.setStrand(StrandString.guessStrand(line.split("\t")[PreMiRNA.COL_STRAND].trim()).toString());

            super.setName("");
            String attribs[] = line.split("\t")[PreMiRNA.COL_ATTRIBUTES].split(";");

            for (String attribStr: attribs){
                String attribType = attribStr.split("=")[0].trim();
                String attribValue = attribStr.split("=")[1].trim();
                switch (attribType){
                    case "ID":
                        mimatID = attribValue;
                        break;

                    case "Alias":
                    case "accession_number":
                        mimatAlias = attribValue;
                        break;

                    case "Name":
                        super.setName(attribValue);
                        break;

                    case "Derives_from":
                    case "derives_from":
                        parent = attribValue;
                        break;

                    default:
                        logger.warn("unknown attribute in parsing miRNA entry \n" + line);
                        logger.info("unknown attribute in parsing miRNA entry \n" + line);
                        break;
                }       
            }
            
        }
        catch(RuntimeException ex){
            logger.error("error parsing miRBase gff reference file miRNA entry \n" + line);
            logger.info("error parsing miRBase gff reference file miRNA entry \n" + line);
            logger.info(ex.toString());
            throw new RuntimeException("error parsing miRBase gff reference file miRNA entry \n" + line);
        }
        
    }
    
    
    

    
    
    
    
    /**
     * Add feature to set
     * 
     * @param key String
     * @param value String
     * @return int FeatureSet size
     */
    public int addFeature(String key, String value){
        featureSet.put(key, value);
        return featureSet.size();
    }
    
    
    
    
    
    /**
     * checks whether Chromosome strings are the same, while attempting
     * to allow for the presence or absence of a variation on the 'Chr' 
     * prefix
     * 
     * @param queryChr String
     * @return Boolean T/F match
     */
    public Boolean chromosomeMatch(String queryChr){
        
        return MiRNA.removeChromosomePrefix(super.getChromosome()).equals(MiRNA.removeChromosomePrefix(queryChr));
        
    }
    
    
    
    
    /**
     * attempt to remove any prefix of the form 'Chr' from the chromosome string
     * 
     * @param chrString String
     * @return String the chromosome string with the prefix (hopefully) removed
     */
    public static String removeChromosomePrefix(String chrString){
        
        if(chrString.contains("chr")){
            chrString = chrString.replace("chr", "");
        }
        else{
            if(chrString.contains("CHR")){
                chrString = chrString.replace("CHR", "");
            }
            else{
                if(chrString.contains("Chr")){
                    chrString = chrString.replace("Chr", "");
                }
            }
            
        }
        return chrString;
        
    }
        
    
    
    
    
    
    
    
    
    
    
    
    
    /**
     * add information to define an isomiR for this entry
     * 
     * @param name String
     * @param start int
     * @param cigar String
     * @param md String
     * @param seq String
     * 
     */
    public void addIsomiR(String name, int start, String cigar, String md, String seq){
        md = md.replace("\t", "|").trim();
        if (isomiRString == null){
            isomiRString = name + ";" + start + ";" + cigar + ";" + md + ";" + seq + "\t";
        }
        else{
             isomiRString = isomiRString.concat(name + ";" + start + ";" + cigar + ";" + md + ";" + seq + "\t");
        }
        
        //this concat() could dhandle null string
//        isomiRString = isomiRString.concat(name + ";" + start + ";" + cigar + ";" + md + ";" + seq + "\t");
        
    }
    
    /**
     * delete the isomiR string
     * 
     */
    public void removeIsomiRs(){
        isomiRString = "";
    }

    
    /**
     * write out isomiRs.
     * 
     * only report reads that are above a baseline, defined in terms of the fraction 
     * of the total number of reads for the miRNA. e.g. if there are 100 reads, and 
     * baseline is 5, then only isomiRs with more than 5 reads will be reported
     * 
     * @param baselinePercent : int     only report isomiRs with reads that
     * @param minCounts       : int     total counts for isomiR must be greater
     *                                  than this value
     * @return String isomiR as printable string
     */
    public String reportIsomiRs(int baselinePercent, int minCounts){
        String reportStr = this.getName() + ":\t[" + this.getChromosome() + "]\t" +this.getStrand()+"\t"
                + this.getStartPos() + "\t" + this.getEndPos() + "\t"+ this.getSeq() + "\n";
        String [] isomiRs = isomiRString.split("\t");
        
        int totalCounts = this.getTotalCounts();
        if (totalCounts < minCounts) return "";
        reportStr = reportStr.concat("Total Counts = " + totalCounts + "\n");
        
        String isomiRStr = "";
        for(String isomiR: isomiRs){
            String[] values = isomiR.split(";");
            if(Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){
                isomiRStr = isomiRStr.concat("name: " + values[0] + "\t"
                            + "start: " + values[STARTCOL] + "\t"
                            + "cigar: " + values[CIGARCOL] + "\t"
                            + "MD: " + values[MDCOL] + "\t" 
                            + "SQ: " + values[SEQCOL] + "\n"
                );
            }
        }
        
        if (isomiRStr.equals("")) return "";
        return reportStr.concat(isomiRStr);
    }
    
    
    
    
    /**
     * report the isomiR in a manner that is visually appealing
     * 
     * @param baselinePercent remove isomiRs with a contribution below this cut off
     * @param minCounts don't report any isomiRs if total counts are below this cut off
     * @return String pretty summary of isomiRs
     */
    public String prettyReportIsomiRs(int baselinePercent, int minCounts){
        
        String reportStr = this.getName() + "|" + this.getMimatID() + " :\tchr" + this.getChromosome() + "\t" +this.getStrand()+"\t"
                + this.getStartPos() + " --> " + this.getEndPos() + " (" + this.getStrand() + ") : " + this.getSeq() + "\n";
        
        int totalCounts = this.getTotalCounts();
        reportStr = reportStr.concat("Total Counts = " + totalCounts + "\n");
        logger.debug(reportStr);
        if (totalCounts < minCounts) return "";
        
    
        String [] isomiRs = isomiRString.split("\t");
        
        int longestName = this.getName().length();
        int longestSeq = this.getSeq().length();
        int longestCounts = 0;
        int longestMD = 0;
        
        for(String isomiR: isomiRs){
            
            if(Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts >= (double)baselinePercent/100.0){
                if(isomiR.split(";")[NAMECOL].length() > longestName)
                    longestName = isomiR.split(";")[NAMECOL].length();

                if(isomiR.split(";")[NAMECOL].split("-")[1].length() > longestCounts)
                    longestCounts = isomiR.split(";")[NAMECOL].split("-")[1].length();

                if(isomiR.split(";")[SEQCOL].length() > longestSeq)
                    longestSeq = isomiR.split(";")[SEQCOL].length();
                if(this.getStrand().equals("+")){
                    if(isomiR.split(";")[SEQCOL].equals(this.getSeq().replace("U", "T")))
                        longestName++;
                }
                else{
                    if(isomiR.split(";")[SEQCOL].equals(SimpleSeq.complement(this.getSeq()).replace("U", "T")))
                        longestName++;                    
                }
                
                    

                if(isomiR.split(";")[MDCOL].length() > longestMD)
                    longestMD = isomiR.split(";")[MDCOL].length();
                
                
            }          

        }    
        
        int leftMargin = 10;
        int ColMargin = 5;
        
        String isoString = "";
        for(String isomiR: isomiRs){
            
            if(Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){

                isoString           = isoString.concat(StringUtils.repeat(" ", leftMargin));
                String isoName      = isomiR.split(";")[NAMECOL];
                if(this.getStrand().equals("+")){
                    if(isomiR.split(";")[SEQCOL].equals(this.getSeq().replace("U", "T")))
                        isoName = "*".concat(isoName);
                }
                else{
                    if(isomiR.split(";")[SEQCOL].equals(SimpleSeq.complement(this.getSeq()).replace("U", "T")))
                        isoName = "*".concat(isoName);
                }
                isoString           = isoString.concat(StringUtils.repeat(" ", longestName - isoName.length()) + isoName);
                
                int isoStart        = Integer.parseInt(isomiR.split(";")[STARTCOL]);
                String isoSeq       = isomiR.split(";")[SEQCOL];
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin + (isoStart - this.getIsomiRMinStart())));
                isoString           = isoString.concat(isoSeq);
                isoString           = isoString.concat(StringUtils.repeat(" ", (this.getIsomiRMaxStop()-(isoStart + isoSeq.length())) + ColMargin));
                
                String isoMD        = isomiR.split(";")[MDCOL];
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin) + isoMD + StringUtils.repeat(" ", (longestMD - isoMD.length()) + ColMargin));
                
                String isoCounts    = isomiR.split(";")[NAMECOL].split("-")[1];
                int isoCountsINT    = Integer.parseInt(isoCounts); 
                float isoPercentage = isoCountsINT * 100.0f / totalCounts;
                isoString           = isoString.concat(StringUtils.repeat(" ", ColMargin + (longestCounts - isoCounts.length())) + isoCounts + StringUtils.repeat(" ", ColMargin)+isoPercentage+ StringUtils.repeat(" ", ColMargin) + "\n");
            
            }
            
        }
        
        if (isoString.equals("")) return "";
        
        return reportStr.concat(isoString + "\n\n\n");

    }
    
    
    
    
     /**
     * an isomiR can be classified as 5´ modification , 3´ modification or polymorphic
     * 
     * @param baselinePercent don't characterize an isomiR that contributes less than 
     *                          this percent of the total counts for the isomiR
     * 
     * @return ArrayList    : list of isomiR points
     */
    public ArrayList characterizeIsomiRs(int baselinePercent){
        
        ArrayList isomiRPts = new ArrayList<>();
        // we can identify 5´ modification from start position
        int totalCounts = this.getTotalCounts();
        String [] isomiRStrings = isomiRString.split("\t");
        String isoString = "";
        logger.debug("this mirRNA name --- "+super.getName());

        for(String isomiR: isomiRStrings){

            
            if(Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts > (double)baselinePercent/100.0){
                logger.debug("isomiR name and strand" +isomiR.split(";")[NAMECOL] + "\t" +this.getStrand() +  "\n" + " : percentage - " + Double.parseDouble(isomiR.split(";")[NAMECOL].split("-")[1]) / (double) totalCounts);


                logger.debug(isomiR);
                int isoStart        = Integer.parseInt(isomiR.split(";")[STARTCOL]);
                String isoSeq       = isomiR.split(";")[SEQCOL];
                int isoEnd          = isoStart + isoSeq.length() - 1;

                int noOf5pSteps = 0;
                int noOf3pSteps = 0;
                int noOfPolySteps = 0;                

                int wStart = 0;
                int wEnd = 0;
                int iStart = 0;
                int iEnd = 0;

                
                if(this.getStrand().equals("+")){


                    logger.debug("isomiR + seq -- "+this.getSeq() + "\n" + isomiR.split(";")[SEQCOL]);

                    if(isoStart != this.getStartPos()){                   
                        noOf5pSteps = isoStart - this.getStartPos();                   
                    }

                    if(isoEnd != this.getEndPos()){
                        noOf3pSteps = isoEnd - this.getEndPos();
                    }
                    logger.debug("noOf5pSteps + -- "+noOf5pSteps + ", " + noOf3pSteps);
                    logger.debug("start + end position -- "+isoStart + ", " + isoEnd + ", " + (this.getStartPos()- isoStart) + ", " + this.getSeq().length());

                    // since we are using wStart, iStart, wEnd and iEnd as an index to get subString, 
                    // it should be -1 before function as index to get subString
                    if(this.getStartPos() >= isoStart ){
                        wStart = 0;
                        iStart = this.getStartPos()- isoStart;
                    }
                    else{
                        wStart = isoStart - this.getStartPos();
                        iStart = 0;
                    }

                    if(this.getEndPos() < isoEnd){
                        wEnd = this.getSeq().length(); 
                        iEnd = isoSeq.length() - (isoEnd - this.getEndPos());
                    }
                    else{
                        wEnd = this.getSeq().length() - (this.getEndPos() - isoEnd);
                        iEnd = isoSeq.length();
                    }
                    
                    if (wStart !=0) wStart= wStart-1;
                    if (iStart !=0) iStart= iStart-1;
                    if (wEnd !=0) wEnd= wEnd-1;
                    if (iEnd !=0) iEnd= iEnd-1;
                   

                    if(this.getSeq().substring(wStart, wEnd).replace("U", "T").equals(isoSeq.substring(iStart, iEnd))==false){
                        for(int b=0; b<wEnd-wStart; b++){
                            if(this.getSeq().charAt(wStart + b) != isoSeq.charAt(iStart + b)) noOfPolySteps++;
                        }
                    }
                }
                else{

                    logger.debug("this.seq - position -- " +this.getSeq() + "\t" + super.getStartPos() + "\t" + super.getEndPos() + "\n");
                    logger.debug(SimpleSeq.complement(isomiR.split(";")[SEQCOL]) + "\t" + isoStart + "\t" + isoEnd + "\n");
                    
                    if(isoStart != this.getStartPos()){                   
                        noOf5pSteps = isoStart - this.getStartPos();                   
                    }

                    if(isoEnd != this.getEndPos()){
                        noOf3pSteps = isoEnd - this.getEndPos();
                    }
                    logger.debug("isomiR- noOf5pSteps -- " +noOf5pSteps + ", " + noOf3pSteps);
                    logger.debug(isoStart + ", " + isoEnd + ", " + (this.getStartPos()- isoStart) + ", " + this.getSeq().length());

                    if(this.getStartPos() >= isoStart ){
                        wStart = 0;
                        iStart = this.getStartPos()- isoStart;
                    }
                    else{
                        wStart = isoStart - this.getStartPos();
                        iStart = 0;
                    }

                    if(this.getEndPos() < isoEnd){
                        wEnd = this.getSeq().length(); 
                        iEnd = isoSeq.length() - (isoEnd - this.getEndPos());
                    }
                    else{
                        wEnd = this.getSeq().length() - (this.getEndPos() - isoEnd);
                        iEnd = isoSeq.length();
                    }

                    if (wStart !=0) wStart= wStart-1;
                    if (iStart !=0) iStart= iStart-1;
                    if (wEnd !=0) wEnd= wEnd-1;
                    if (iEnd !=0) iEnd= iEnd-1;
                    
                    //here it will get an error: java.lang.StringIndexOutOfBoundsException: String index out of range
                    if(this.getSeq().substring(wStart, wEnd).replace("U", "T").equals(SimpleSeq.complement(isoSeq.substring(iStart, iEnd)))==false){
                        String complementWTSeq = SimpleSeq.complement(this.getSeq());



                        for(int b=0; b<wEnd-wStart; b++){
                            if(complementWTSeq.charAt(wStart + b) != isoSeq.charAt(iStart + b)) noOfPolySteps++;
                        }
                    }  
                    
                }


                String isoCounts    = isomiR.split(";")[NAMECOL].split("-")[1];
                
                HashMap isomiRPt = new HashMap();
                isomiRPt.put("5p", noOf5pSteps);
                isomiRPt.put("3p", noOf3pSteps);
                isomiRPt.put("poly", noOfPolySteps);
                isomiRPt.put("fraction", Double.parseDouble(isoCounts)/(double)totalCounts);
                
                isomiRPts.add(isomiRPt);

            }

        }
        
        return isomiRPts;

    }
    

    
    
    
    /**
     * sum counts for this isomiR
     * 
     * @return int total counts
     */
    public int getTotalCounts(){
        
        int totalCounts = 0;
        if(isomiRString.isEmpty()) return 0;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){            
            totalCounts += Integer.parseInt(isomiR.split(";")[NAMECOL].split("-")[1]);
        }
        return totalCounts;
    }
    
    
    /**
     * find the min start position within the isomiRs
     * 
     * @return int min start position 
     */
    public int getIsomiRMinStart(){
        
        int isomiRMinStart = 1000000000;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){         
            if(Integer.parseInt(isomiR.split(";")[STARTCOL]) < isomiRMinStart) 
                isomiRMinStart = Integer.parseInt(isomiR.split(";")[STARTCOL]);
        }
        return isomiRMinStart;
        
    }
    
    
    
    
    /**
     * find the max stop position within the isomiRs
     * 
     * @return int Max stop position
     */
    public int getIsomiRMaxStop(){
        
        int isomiRMaxStop = -1;
        String [] isomiRs = isomiRString.split("\t");

        for(String isomiR: isomiRs){         
            if(Integer.parseInt(isomiR.split(";")[STARTCOL]) + isomiR.split(";")[SEQCOL].length() > isomiRMaxStop) 
                isomiRMaxStop = Integer.parseInt(isomiR.split(";")[STARTCOL]) + isomiR.split(";")[SEQCOL].length() - 1;
        }
        return isomiRMaxStop;
        
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
    	this.setName(headerLine.split(headerLine)[0].trim());
    	this.setMimatID(headerLine.split(headerLine)[2].trim());
    }
    
    /**
     * parse out miRBase name from header line of FastA entry and return
     * @param headerLine
     * @return
     */
    public static String parseNameFromMiRBaseFAHeaderLine(String headerLine) {
    	return headerLine.split(headerLine)[0].trim();
    }
    
    
    /**
     * parse out MIMATID from header line of FastA entry and return
     * @param headerLine
     * @return
     */
    public static String parseIDFromMiRBaseFAHeaderLine(String headerLine) {
    	return headerLine.split(headerLine)[2].trim();
    }
    
    
    /**
     * this will only parse miRNA entries defined in EMBL format in the miRNA.dat
     * file released by miRBase. i.e., only a subset of qualifies are recognized
     * e.g. 
     * 
     *  FT   miRNA           17..38
     *  FT                   /accession="MIMAT0000001"
     *  FT                   /product="cel-let-7-5p"
     *  FT                   /evidence=experimental
     *  FT                   /experiment="cloned [1-3], Northern [1], PCR [4], 454 [5],
     *  FT                   Illumina [6], CLIPseq [7]"
     * 
     * @param emblLines ArrayList of String entries
     */
    public void parseMiRBaseMiRNAEMBLlines(ArrayList<String> emblLines){
        String lastQualifier = "";
        for(String emblLine:emblLines){
            // header line
            if(emblLine.split("\\s+")[1].trim().equals("miRNA")){
                this.miStart= Integer.parseInt(emblLine.substring(QUALIFIER_START_COL, emblLine.indexOf("..")));
                this.miEnd=Integer.parseInt(emblLine.substring(emblLine.indexOf("..")+2).trim());
                continue;
            }
            
            // key / value entries
            String qualifier;
            if(emblLine.contains("="))
                qualifier = emblLine.substring(QUALIFIER_START_COL, emblLine.indexOf("="));
            else
                qualifier="";
            
            switch(qualifier){
                case QUALIFIER_ACCESSION:
                    this.mimatID = emblLine.substring(emblLine.indexOf("=")+1).replaceAll("\"", "");
                    lastQualifier = QUALIFIER_ACCESSION;
                    break;
                    
                case QUALIFIER_PRODUCT:
                    super.setName(emblLine.substring(emblLine.indexOf("=")+1).replaceAll("\"", ""));
                    lastQualifier = QUALIFIER_PRODUCT;
                    break;
                    
                case QUALIFIER_EVIDENCE:
                    this.evidence = emblLine.substring(emblLine.indexOf("=")+1).replaceAll("\"", "");
                    lastQualifier = QUALIFIER_EVIDENCE;
                    break;
                    
                case QUALIFIER_EXPERIMENT:
                    this.references = emblLine.substring(emblLine.indexOf("=")+1).replaceAll("\"", "");
                    lastQualifier = QUALIFIER_EXPERIMENT;
                    break;
                    
                default:
                    concatKeyValue(emblLine, lastQualifier);
            }
        }
    }
    
    
    /**
     * this handles continuation lines in the key/value entry. It should only be 
     * called from @see parseMiRBaseMiRNAlines(ArrayList<String> emblLines). it 
     * is only included in a separate method for readability
     * 
     * @param emblLine
     * @param lastQualifier
     * @return String
     */
    private String concatKeyValue(String emblLine, String lastQualifier){
        
        switch(lastQualifier){
            case QUALIFIER_ACCESSION:
                this.mimatID = mimatID.concat(emblLine.substring(QUALIFIER_START_COL).replaceAll("\"", "") + " ");
                lastQualifier = QUALIFIER_ACCESSION;
                break;

            case QUALIFIER_PRODUCT:
                super.setName(super.getName().concat(emblLine.substring(QUALIFIER_START_COL).replaceAll("\"", "") + " "));
                lastQualifier = QUALIFIER_ACCESSION;
                break;

            case QUALIFIER_EVIDENCE:
                this.evidence = evidence.concat(emblLine.substring(QUALIFIER_START_COL).replaceAll("\"", "") + " ");
                lastQualifier = QUALIFIER_ACCESSION;
                break;

            case QUALIFIER_EXPERIMENT:
                this.references = references.concat(emblLine.substring(QUALIFIER_START_COL).replaceAll("\"", "") + " ");
                lastQualifier = QUALIFIER_ACCESSION;
                break;

            default:
      
        }                    
      
        return lastQualifier;
    }
    
    
    
    
    /**
     * Compare this miRNA with query
     * @param miRFeat miRNA
     * @return int match outcome
     */
    @Deprecated
    public int compareTo(MiRNA miRFeat) {

        int thisChr = -1;
        int queryChr = -1;
        if(this.getChromosome().contains("chr")){
            thisChr = Integer.parseInt(this.getChromosome().replace("chr", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("chr", ""));
        }
        
        if(this.getChromosome().contains("CHR")){
            thisChr = Integer.parseInt(this.getChromosome().replace("CHR", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("CHR", ""));
        }
        
        if(this.getChromosome().contains("Chr")){
            thisChr = Integer.parseInt(this.getChromosome().replace("Chr", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("Chr", ""));
        }
        
     return thisChr - queryChr;

     }
   
    
    
    
    /**
     * base equality on the mimatID, which should be unique
     * 
     * @param qObject Object
     * @return Boolean match T/F
     */
    @Override
    public boolean equals(Object qObject){
        if (qObject != null && qObject instanceof MiRNA)
        {
            return (this.mimatID.equals(((MiRNA) qObject).mimatID));
        }
        return false;

    }
    
    /**
     * generate hashCode for this miRNA for comparison
     * @return int HashCode
     */
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((mimatID == null) ? 0 : mimatID.hashCode());
        return result;
    }    
    
    
    
    
    /**
     * 
     * count the number of GU-pairs by parsing out pretty plot.
     * This works based on the assumption there can only be G-C and G-U pairing,
     * which is (it seems) an assumption of RNAFold. 
     * There are 4 lines in the pretty plot, but the pairing can only occur in the
     * middle two.
     * 
     * The differences in the C and G counts gives us the number of GU pairs
     * 
     * e.g.
     * 
     *     CCGCCAU
     *     |||||||
     *     GGUGGUG
     * 
     * @return int : the number of GU-wobble pairs
     * 
     */
    public int CountGUPairsFromPrettyPlot(){
        int n=upperStructLine2.length()+lowerStructLine4.length();
        int nc=n-(SimpleSeq.NTcount(upperStructLine2,'C')
                +SimpleSeq.NTcount(upperStructLine2,'c')
                +SimpleSeq.NTcount(lowerStructLine4,'C')
                +SimpleSeq.NTcount(lowerStructLine4,'c'));
        int ng=n-(SimpleSeq.NTcount(upperStructLine2,'G')
                +SimpleSeq.NTcount(upperStructLine2,'g')
                +SimpleSeq.NTcount(lowerStructLine4,'G')
                +SimpleSeq.NTcount(lowerStructLine4,'g'));
        return nc-ng;
    }

    @Override
    public void setID(String id){
        super.setID(id);
    }
    /**
     * @return String mimatID
     */
    public String getMimatID() {
        return mimatID;
    }

    /**
     * @param mimatID String mimatID value
     */
    public void setMimatID(String mimatID) {
        this.mimatID = mimatID;
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
     * @return the isomiRString
     */
    public String getIsomiRString() {
        return isomiRString;
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

    /**
     * @return the miStart
     */
    public int getMiStart() {
        return miStart;
    }

    /**
     * @param miStart the miStart to set
     */
    public void setMiStart(int miStart) {
        this.miStart = miStart;
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
    public void setMiEnd(int miEnd) {
        this.miEnd = miEnd;
    }

    /**
     * @return the featureSet
     */
    public HashMap getFeatureSet() {
        return featureSet;
    }

    /**
     * @return the evidence
     */
    public String getEvidence() {
        return evidence;
    }

    /**
     * @param evidence the evidence to set
     */
    public void setEvidence(String evidence) {
        this.evidence = evidence;
    }

    /**
     * @return the references
     */
    public String getReferences() {
        return references;
    }

    /**
     * @param references the references to set
     */
    public void setReferences(String references) {
        this.references = references;
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
     * @return the mimatAlias
     */
    public String getMimatAlias() {
        return mimatAlias;
    }
    
    /**
     * @param newMimatAlias String 
     * change to new mimatAlias when find duplicate record in miRBase gtf file
     */
    public void setMimatAlias(String newMimatAlias) {
        this.mimatAlias = newMimatAlias;
    }


    public static void main (String args[]){
        String miRBaseLine = "chr1\t.\tmiRNA\t3580220\t3580240\t.\t+\t.\tID=MIMAT0011792;"
                        + "Alias=MIMAT0011792;Name=bta-miR-2284i;Derives_from=MI0011294";
        MiRNA mirseq = new MiRNA();
        mirseq.parseMiRBaseGFFEntry(miRBaseLine);
        
        mirseq.addIsomiR("isomiR1-1000", 3580220, "18M", "XA:i:2\t MD:Z:1C9C6\t NM:i:2", "GAGCGCACCGGACGGGCC");
        mirseq.addIsomiR("isomiR2-100",  3580220, "16M", "XA:i:2\t MD:Z:1C9C6\t NM:i:2", "GTGCGCACCGGACGGGCC");
        mirseq.addIsomiR("isomiR3-60",   3580221, "18M", "XA:i:2\t MD:Z:1C9C6\t NM:i:2",  "AGCGCACCGGACGGGCCT");
        mirseq.addIsomiR("isomiR4-25",   3580221, "19M", "XA:i:2\t MD:Z:1C9C6\t NM:i:2", "AAGCGCACCGGACGGGCC");
        
        logger.error(mirseq.prettyReportIsomiRs(5, 200));

        logger.error(mirseq.getIsomiRMinStart());
        logger.error(mirseq.getIsomiRMaxStop());
        
        logger.error(mirseq.prettyReportIsomiRs(5, 200));
        
    }

}
