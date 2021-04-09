/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.reference;

import no.uio.medisin.bag.core.sequence.ChromosomeString;
import no.uio.medisin.bag.core.sequence.StrandString;
import no.uio.medisin.bag.core.sequence.Strand;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * stores a single entry from a GFF line
 * Undefined fields are replaced with the "." character, as described in the original GFF spec.
 *
 * Column 1: "seqid" (e.g. chromosome)
 * Column 2: "source" (e.g. program name)
 * Column 3: "type"  (feature type name, e.g. Gene, Variation, Similarity)
 * Columns 4 & 5: "start" and "end"
 * Column 6: "score"
 * Column 7: "strand"
 * Column 8: "phase"
 * Column 9: "attributes"
 * A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons. 
 * These tags have predefined meanings:
 * 
 *   ID            - Indicates the ID of the feature. must be unique within the scope of the GFF file. 
 *   Name          - Display name for the feature. not necessarily unique
 *   Alias         - A secondary name for the feature. e.g. locus names and accession numbers. not necessarily unique
 *   Parent        - Indicates the parent of the feature. 
 *   Target        - Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. 
 *   Gap           - The alignment of the feature to the target if the two are not collinear 
 *   Derives_from  - Used to disambiguate the relationship between one feature and another 
 *   Note          - A free text note.
 *   Dbxref        - A database cross reference. 
 *   Ontology_term - A cross reference to an ontology term. 
 *   Is_circular   - A flag to indicate whether a feature is circular. 
 * 
 * @author simonray
 */
public class GFFEntry implements Comparable<GFFEntry>{
    
    static                  Logger logger    = LogManager.getLogger();
    
    public static final    int GFF_SEQID  = 0;
    public static final    int GFF_SRC    = 1;
    public static final    int GFF_TYPE   = 2;
    public static final    int GFF_START  = 3;
    public static final    int GFF_STOP   = 4;
    public static final    int GFF_SCORE  = 5;
    public static final    int GFF_STRAND = 6;
    public static final    int GFF_PHASE  = 7;
    public static final    int GFF_ATTR   = 8;
        
    static String   CR = System.getProperty("line.separator");
    
    
    private                 String  refSeqID;
    private                 String  featureID;
    private                 String  src;
    private                 String  type;
    private                 int     start;
    private                 int     stop;
    private                 float   score;
    private                 Strand  strand;
    private                 int     phase;
    private                 String  attrString;
    
    

    public GFFEntry(String line) throws RuntimeException{
        
        try{
            refSeqID    = line.split("\t")[GFF_SEQID];
            src         = line.split("\t")[GFF_SRC];
            type        = line.split("\t")[GFF_TYPE];
        }
        catch(Exception ex){
            logger.error("exception while parsing seqID/src/type values in GFF entry: " + line);
            logger.error(ex);
            throw new RuntimeException("exception while parsing seqID/src/type values in GFF entry" + line);
        }
        try{
            start   = Integer.parseInt(line.split("\t")[GFF_START]);
            stop    = Integer.parseInt(line.split("\t")[GFF_STOP]);
        }
        catch(Exception ex){
            logger.error("exception while parsing start/stop values in GFF entry" + line);
            logger.error(ex);      
            start = -1;
            stop  = -1;
            throw new RuntimeException("exception while parsing seqID/src/type values in GFF entry" + line);
        }
        try{
            score   = Float.parseFloat(line.split("\t")[GFF_SCORE]);
        }
        catch(NumberFormatException exNF){
            score = 0;
        }
        
        strand  = StrandString.guessStrand(line.split("\t")[GFF_STRAND]);
        
        try{
            phase   = Integer.parseInt(line.split("\t")[GFF_PHASE]);
        }
        catch(NumberFormatException exNF){
            phase = 0;
        }
        
        if(line.split("\t")[GFF_ATTR].isEmpty() == false){
            attrString    = line.split("\t")[GFF_ATTR];
            String attribs[] = attrString.split(";");
            for(String a:attribs){
                if(a.toUpperCase().contains("ID=")){
                    featureID=a.split("=")[1];
                }
            }
        }

        
        
        
    }
    
    
    
    
    
    
    /**
     * create a new GFF Entry.
     * We require all fields to be specified because there are too many possible
     * subsets of parameters that might be specified
     * 
     * @param nRefSeqID String
     * @param nSource String
     * @param nType String
     * @param nStart int
     * @param nStop int
     * @param nScore String
     * @param nStrand String
     * @param nPhase String
     * @param nAttrString  String
     * 
     */
    public GFFEntry(String nRefSeqID, String nSource, String nType, int nStart, int nStop, String nScore, String nStrand, String nPhase, String nAttrString){
        
        refSeqID = nRefSeqID;
        src = nSource;
        type = nType;
        start = nStart;
        stop = nStop;
        try{
            score   = Float.parseFloat(nScore);
        }
        catch(NumberFormatException exNF){
            score = 0;
        }
        strand = StrandString.guessStrand(nStrand);
        try{
            phase   = Integer.parseInt(nPhase);
        }
        catch(NumberFormatException exNF){
            phase = 0;
        }
        attrString = nAttrString;
        if(nAttrString.isEmpty() == false){
            String attribs[] = attrString.split(";");
            for(String a:attribs){
                if(a.toUpperCase().contains("ID=")){
                    featureID=a.split("=")[1];
                }
            }
        }
    }

 
    
    
    
    /**
     * check whether the supplied line is a comment line (begins with '#')
     * 
     * @param line String
     * @return Boolean T/F is comment line
     */
    public static Boolean isCommentLine(String line){
        return line.startsWith("#");
    }
     
    
    
    
    
    /**
     * compare by refSeq, strand and then position
     * 
     * @param queryGFFEntry GFFEntry
     * @return int outcome of comparison
     */
    @Override
    public int compareTo(GFFEntry queryGFFEntry){
        
        if(queryGFFEntry.getStrand() != getStrand()){
            if(getStrand() == Strand.PLUS && queryGFFEntry.getStrand() == Strand.MINUS){
                return -1;
            }else{
                return +1;
            }            
        }

        if (queryGFFEntry.getRefSeqID().equals(getRefSeqID()))
        {
            return getStart() - queryGFFEntry.getStart() ;
        }
        
            
        //int l = Math.min(queryGFFEntry.refSeqID.length(), refSeqID.length());
        
        /*
        int i=0;
        while(queryGFFEntry.refSeqID.charAt(i)==(refSeqID.charAt(i)) && i<l){
            i++;
        }
        return  refSeqID.charAt(i) - queryGFFEntry.refSeqID.charAt(i) ;
        */
        
        return getStart() - queryGFFEntry.getStart();
    }
    
    


    /**
     * does the specified region overlap this gffEntry and featureType ?
     * 
     * @param qStart int
     * @param qStop int
     * @param qStrand Strand
     * @param qSeqID String
     * @param featureName String
     * @param bleed int
     * @return Boolean T/F region overlaps
     */
    public Boolean doesRegionOverlap(int qStart, int qStop, Strand qStrand, String qSeqID, String featureName, int bleed){
        
        if(!this.getType().toUpperCase().equals(featureName.toUpperCase())) return false;
        if(this.getStrand() != qStrand) return false;
        if(ChromosomeString.GetChromNumber(this.getSeqID())==ChromosomeString.GetChromNumber(qSeqID) == false) return false;
        return (qStart -this.getStart() > 0 && this.getStop() - qStop > 0);
       

    }
    
    

    
    
    /**
     * does the specified region overlap this gffEntry and featureType ?
     * 
     * @param qStart int
     * @param qStop int
     * @param qStrand Strand
     * @param qChr String
     * @param featureName String
     * @return Boolean T/F feature overlaps region
     */
    public Boolean doesFeatureContainRegion(int qStart, int qStop, Strand qStrand, String qChr, String featureName){
        return 
          this.getType().toUpperCase().equals(featureName.toUpperCase())  
            && this.getStrand() == qStrand
            && this.getSeqID().equals(qChr)
            && this.getStart() <= qStart 
            && this.getStop() >= qStop;
    }




    
    /**
     * does the specified region overlap this gffEntry and featureType ?
     * 
     * @param queryGFFEntry GFFEntry
     * @return  Boolean T/F feature overlaps region
     */
    public Boolean doesFeatureContainRegion(GFFEntry queryGFFEntry){
        return 
           this.getStrand() == queryGFFEntry.getStrand()
            && this.getSeqID().equals(queryGFFEntry.getSeqID())
            && this.getStart() <= queryGFFEntry.getStart()
            && this.getStop() >= queryGFFEntry.getStop();
    }




    
    /**
     * does the specified region overlap this gffEntry regardless of featureType?
     * 
     * @param qStart int
     * @param qStop int 
     * @param qStrand Strand
     * @param qChr String
     * @param bleed int
     * @return  Boolean T/F region overlaps
     */
    public Boolean doesRegionOverlap(int qStart, int qStop, Strand qStrand, String qChr, int bleed){
        return 
            this.getStrand() == qStrand
            && this.getSrc().equals(qChr)
            && Math.abs(this.getStart()- qStart) < bleed 
            && Math.abs(this.getStop() - qStop) < bleed;
    }
    
    
    /**
     * 
     * @return String GFF entry as String
     */
    public String toGFF3String(){
        /*
            chr1            chromosome
            source          n/a here
            miRNA           feature type (n/a)
            start pos
            end pos
            score           n/a here               
            strand          (+/-)
            frame           n/a here
            attributes      e.g. ID=MIMAT0027619;Alias=MIMAT0027619;Name=hsa-miR-6859-3p;Derives_from=MI0022705

        */
        String gff3String = this.getRefSeqID() + "\t"
                + this.getSrc() + "\t"
                + this.getType() + "\t"
                + getStart() + "\t"
                + getStop() + "\t"
                + "." + "\t"
                + getStrand() + "\t"
                + "." + "\t"
//                + "ID=" + featureID + ";"  + "Seq=" + this.getSequence() + "\t"
                + this.getAttr();
        return gff3String;
    }
        
    

    /**
     * add attribute to attribute string
     * 
     * @param attrKey String
     * @param attrVal String
     */
    public void addAttr(String attrKey, String attrVal){
        if(attrString.isEmpty())
            attrString = attrKey + "=" + attrVal;            
        else
            attrString = attrString.concat(";" + attrKey + "=" + attrVal );            
    }   
    
    
    /**
     * remove specified attribute
     * 
     * @param  attrKey String 
     * @return String updated attribute String
     */
    public String removeAttr(String attrKey){
        String[] attrbs = attrString.split(";");
        String newAttr = "";
        for(String attr: attrbs){
            if(attr.split("=")[0].trim().equals(attrKey)==false){
                newAttr = newAttr.concat(attr);
            }
        }
        attrString = newAttr;
        
        return attrString;
    }
    
    
    
    
    /**
     * 
     * clear the attribute String
     * @return String (should be empty)
     * 
     */
    public String clearAttrString(){
        attrString = "";
        return attrString;
    }
    
    
    
    
    /**
     * return the specified Attribute value
     * 
     * @param attrKey String
     * @return String the requested Attribute
     */
    public String getAttrValue(String attrKey){
        String attrs[] = attrString.split(";");
        for (String attr: attrs){
            if(attr.contains(attrKey)){
                return attr.split("=")[1].trim();
            }
        }
        return null;
    }
    
    
    /**
     * @return int seqID
     */
    public String getSeqID() {
        return getRefSeqID();
    }

    /**
     * @param seqID String new value
     */
    public void setSeqID(String seqID) {
        this.setRefSeqID(seqID);
    }

    /**
     * @return String src
     */
    public String getSrc() {
        return src;
    }

    /**
     * @return int type
     */
    public String getType() {
        return type;
    }

    /**
     * @return int start
     */
    public int getStart() {
        return start;
    }

    /**
     * @return int stop
     */
    public int getStop() {
        return stop;
    }

    /**
     * @return float score
     */
    public float getScore() {
        return score;
    }

    /**
     * @return Strand strand
     */
    public Strand getStrand() {
        return strand;
    }

    /**
     * @return int phase
     */
    public int getPhase() {
        return phase;
    }

    /**
     * @return String attr
     */
    public String getAttr() {
        return attrString;
    }
    
    
    
    /**
     * add a sequence by appending to the attribute string
     * 
     * @param seq String
     */
    public void setSequence(String seq){
        attrString = attrString.concat(";seq=" + seq);            
    }
    
    
    
    /**
     * extract sequence from the attribute string
     * 
     * @return String sequence
     */
    public String getSequence(){
        if(attrString.contains("seq=")){
            int startPos = attrString.indexOf("seq=")+4;
            int stopPos = attrString.indexOf(";", startPos);
            if(stopPos==-1)
                stopPos = attrString.length() - 1;
            return attrString.substring(startPos, stopPos);
        }
        return "";
    }
    
    
    
    
    /**
     * write the entry as in FASTA format as used by MiRBase
     * for consistency with MiRBAse the header line must have the format
     * >cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
     * 
     * @return String miRNA as FASTA
     */
    public String toMirbaseFastAString(){
        if(attrString.contains("seq=")){            
            int startPos = attrString.indexOf("seq=")+4;
            int stopPos = attrString.indexOf(";", startPos);
            if (stopPos==-1)
                stopPos=attrString.length()-1;
            return ">" + getFeatureID() + "|" + this.getRefSeqID()  + ":" + getStart() + "-" + getStop() + "(" + getStrand() + ")" + CR + attrString.substring(startPos, stopPos);
        }
        return "";
    }

    /**
     * @return String featureID
     */
    public String getFeatureID() {
        return featureID;
    }

    /**
     * @param featureID String the new featureIDt
     */
    public void setFeatureID(String featureID) {
        this.featureID = featureID;
    }

    /**
     * @return String refSeqID
     */
    public String getRefSeqID() {
        return refSeqID;
    }

    /**
     * @param refSeqID String the new refSeqID
     */
    public void setRefSeqID(String refSeqID) {
        this.refSeqID = refSeqID;
    }
}
