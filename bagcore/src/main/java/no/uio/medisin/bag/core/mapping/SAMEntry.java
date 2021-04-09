/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mapping;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.core.sequence.Strand;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * parse a SAM line
 * 
 * the line can be two types. A read entry or a data line (starts with a @)
 * 
 *
    1   QNAME	   Query template NAME
    2   FLAG	   bitwise FLAG
    3   RNAME	   Reference sequence NAME
    4   POS	   1-based leftmost mapping POSition
    5   MAPQ	   MAPping Quality
    6   CIGAR	   CIGAR string
    7   RNEXT	   Ref. name of the mate/next read
    8   PNEXT	   Position of the mate/next read
    9   TLEN	   observed Template LENgth
    10  SEQ	   segment SEQuence
    11  QUAL	   ASCII of Phred-scaled base QUALity+33
 *                       
 *
 * we don't implement Comparable as there is too much information to consider
 * {@link no.uio.medisin.bag.core.mapping.MappedRead} is a lightweight version that only contains position information
 * and counts that is comparable
 * 
 * @author simonray
 */
public final class SAMEntry {
    
    static Logger                       logger                      = LogManager.getLogger();
    
    
    public static final int             COL_QNAME      = 0;
    public static final int             COL_FLAG       = 1;
    public static final int             COL_RNAME      = 2;
    public static final int             COL_POS        = 3;
    public static final int             COL_MAPQ       = 4;
    public static final int             COL_CIGAR      = 5;
    public static final int             COL_RNEXT      = 6;
    public static final int             COL_PNEXT      = 7;
    public static final int             COL_TLEN       = 8;
    public static final int             COL_SEQ        = 9;
    public static final int             COL_QUAL       = 10;
    public static final int             COL_TAGSTR     = 11;
    
    public static final short           FLAG_MULTI_SEGS         = 0x01;     // template having multiple segments in sequencing
    public static final short           FLAG_PROPER_ALN         = 0x02;     // each segment properly aligned according to the aligner
    public static final short           FLAG_UNMAPPED           = 0x04;     // segment unmapped
    public static final short           FLAG_NEXTSEG_UNMAP      = 0x08;     // next segment in the template unmapped
    public static final short           FLAG_REVCOMP            = 0x10;     // SEQ being reverse complemented
    public static final short           FLAG_NEXTSEG_REVCOMP    = 0x20;     // SEQ of the next segment in the template being reverse complemented
    public static final short           FLAG_FIRSTSEG           = 0x40;     // the first segment in the template
    public static final short           FLAG_LASTSEG            = 0x80;     // the last segment in the template
    public static final short           FLAG_SECONDARY_ALN      = 0x100;    // secondary alignment
    public static final short           FLAG_FILTER_FAIL        = 0x200;    // not passing filters, such as platform/vendor quality controls
    public static final short           FLAG_DUPLICATE          = 0x400;    // PCR or optical duplicate
    public static final short           FLAG_SUPP_ALN           = 0x800;    // supplementary alignment
    

    /*
    Bowtie generates the following tags in column 12
    (https://www.biostars.org/p/257629/)
    Tag Meaning
      NM     Edit distance
      MD     Mismatching positions/bases
      AS     Alignment score
      BC     Barcode sequence
      X0     Number of best hits
      X1     Number of suboptimal hits found by BWA
      XN     Number of ambiguous bases in the reference
      XM     Number of mismatches in the alignment
      XO     Number of gap opens
      XG     Number of gap extentions
      XT     Type: Unique/Repeat/N/Mate-sw
      XA     Alternative hits; format: (chr,pos,CIGAR,NM;)*
      XS     Suboptimal alignment score
      XF     Support from forward/reverse alignment
      XE     Number of supporting seeds

    */
    
    /*
    
    the SAM entry may also contain MD tags
      https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files    
    
      The MD tag is for SNP/indel calling without looking at the reference. 
      It does this by carrying information about the reference that 
      the read does not carry, for a particular alignment. A SNP's alternate 
      base is carried in the read, but without the MD tag or use of the 
      alignment reference, it's impossible to know what the reference base was. 
      Thus, this information is carried in the MD tag. A SNP looks like:

        10A3T0T10    
    
      Here, there are three SNPs:

        11 bases in from the aligned portion of the read (the reference has an A) (excluding softclips).
        15 bases in there is s a T in the reference.
        16 bases in there is a T in the reference. 
    
      Note that 0s are use used to indicate positions of neighboring SNPs.

      Likewise, a reference would be necessary to know the bases deleted from 
      the reference in an alignment. The MD tag stores this information too:  
    
        85^A16
    
      Here, there are 85 matches, 1 deletion from the reference 
      (the reference has an A and the read doesn't), followed by 16 matches.
    */
    
    private String                      samLine;
    private String                      qName;
    private short                       flags;
    private String                      rName;
    private int                         startPos;
    private int                         endPos;
    private Strand                      strand;
    private int                         mapQ;
    private String                      cigar;
    private String                      rNext;
    private int                         pNext;
    private int                         tLen;
    private String                      seq;
    private String                      qual;
    private String                      bowtieTagString;
    private String                      mdString;
    private Boolean                     header;
    private Boolean                     mapped;
    
    private Boolean                     fullInfo = true;
    
    
    public SAMEntry(String sLine){
      samLine = sLine;
      this.parseLine();
    }
    
    /**
     * is line a header line?
     * 
     * @param line samLine to test
     * @return Boolean whether this is a HeaderLine
     */
    public static final Boolean isHeaderLine(String line){
        return line.startsWith("@");
    }
    
    
    /**
     * extract information from input line
     * at the moment, this only extracts the information I am interested in
     * 
     * @param samLine String to be parsed
     * 
     */
    public void parseLine(){

        try{
            if(samLine.startsWith("@")) {
                header = true;
                mapped = false;
                return;
               
            }

            qName = samLine.split("\t")[COL_QNAME];
 
            header = false;
            flags = Short.parseShort(samLine.split("\t")[COL_FLAG]);
            if ((flags & FLAG_UNMAPPED) == FLAG_UNMAPPED) {
                mapped = false;
                return;                
            }
            
            mapped = true;


            
            if(isThisFlagSet(FLAG_REVCOMP)){
                strand = Strand.MINUS;
            }
            else{
                strand = Strand.PLUS;
            }
            if(mapped==true){

                mapQ = Integer.parseInt(samLine.split("\t")[COL_MAPQ]);
                setStartPos(Integer.parseInt(samLine.split("\t")[COL_POS]));
                String cigarStr = samLine.split("\t")[COL_CIGAR].replace("M", "").trim();
                setEndPos(getStartPos() + samLine.split("\t")[COL_SEQ].length() - 1);
                //endPos = getStartPos() + Integer.parseInt(cigarStr) - 1;
                rName= samLine.split("\t")[COL_RNAME].trim();
                //this.extractTags();
                setTagString(samLine.split("\t", COL_TAGSTR+1)[COL_TAGSTR]);
                cigar = samLine.split("\t")[COL_CIGAR];
                rNext = samLine.split("\t")[COL_RNEXT];
                pNext = Integer.parseInt(samLine.split("\t")[COL_PNEXT]);
                tLen = Integer.parseInt(samLine.split("\t")[COL_TLEN]);
                if(fullInfo){
                    seq = samLine.split("\t")[COL_SEQ];
                    qual = samLine.split("\t")[COL_QUAL];
                }
                
            }
        }
        catch(Exception ex){
            logger.error("error parsing samLine " + samLine);
            logger.error(ex);
        }
        
    }
    

    /**
     * Find the String associated with the specified tag.
     * Attributes are in the format 
     * 
     *   TWO LETTER TAG:TAG_TYPE:VALUE
     *   where 
     *      TAG_TYPE is Z, i, or ?
     * 
     * @param tagID the tagID to be extracted
     * @return String the key/value pair for the tagID
     */
    public String getTag(String tagID){
        String tags[] = getBowtieTagString().split("\\s+");
        for(String tag:tags){
            if(tag.contains(tagID)){
                return tag;
            }
        }
        return null;
    }
    
    
    
    /**
     * Check whether the specified flag has been set.
     * 
     * @param qFlag the query flag
     * 
     * @return Boolean true if the flag is set
     * 
     */
    public Boolean isThisFlagSet(int qFlag){
        
        switch (qFlag){
            
            case FLAG_PROPER_ALN:
                return  (this.flags & FLAG_PROPER_ALN) == FLAG_PROPER_ALN ; 
                
            case FLAG_UNMAPPED:
                return  (this.flags & FLAG_UNMAPPED) == FLAG_UNMAPPED ; 
                
            case FLAG_NEXTSEG_UNMAP:
                return  (this.flags & FLAG_NEXTSEG_UNMAP) == FLAG_NEXTSEG_UNMAP ; 
                
            case FLAG_REVCOMP:
                return  (this.flags & FLAG_REVCOMP) == FLAG_REVCOMP ; 
                
            case FLAG_NEXTSEG_REVCOMP:
                return  (this.flags & FLAG_NEXTSEG_REVCOMP) == FLAG_NEXTSEG_REVCOMP ; 
                
            case FLAG_FIRSTSEG:
                return  (this.flags & FLAG_FIRSTSEG) == FLAG_FIRSTSEG ; 
                
            case FLAG_LASTSEG:
                return  (this.flags & FLAG_LASTSEG) == FLAG_LASTSEG ; 
                
            case FLAG_SECONDARY_ALN:
                return  (this.flags & FLAG_SECONDARY_ALN) == FLAG_SECONDARY_ALN ; 
                
            case FLAG_FILTER_FAIL:
                return  (this.flags & FLAG_FILTER_FAIL) == FLAG_FILTER_FAIL ; 
                
            case FLAG_DUPLICATE:
                return  (this.flags & FLAG_DUPLICATE) == FLAG_DUPLICATE ; 
                
            case FLAG_SUPP_ALN:
                return  (this.flags & FLAG_SUPP_ALN) == FLAG_SUPP_ALN ; 
                
            default:
                return null;
                
        }
    }




    
    /**
     * Find the String associated with the specified tag
     * Attributes are in the format 
     * 
     *   TWO LETTER TAG:TAG_TYPE:VALUE
     *   where 
     *      TAG_TYPE is Z, i, or ?
     * 
     * @param tagID the tag to search for
     * @return String the value for the tagID key
     */
    public String getTagValue(String tagID){
        String tags[] = getBowtieTagString().split("\\s+");
        for(String tag:tags){
            if(tag.contains(tagID)){
                return tag.split(":")[2];
            }
        }
        return null;
    }

    
    
    /**
     * break the CIGAR String into individual components
     * @return 
     */
    public ArrayList<String>  splitCigarString(){
        //([0-9]+)([MIDNSHPX=])
        ArrayList<String> cigarElements = new ArrayList();
        Pattern pStart = Pattern.compile("[\\d]+[a-zA-Z|=]");
        Matcher m = pStart.matcher(this.cigar);
        while (m.find()) {
            cigarElements.add(m.group());
        }
        return cigarElements;
    }

    /**
     * @return the qName
     */
    public String getqName() {
        return qName;
    }

    /**
     * @param qName the qName to set
     */
    public void setqName(String qName) {
        this.qName = qName;
    }

    /**
     * @return the flags
     */
    public short getFlags() {
        return flags;
    }

    /**
     * @return the rName
     */
    public String getrName() {
        return rName;
    }

    /**
     * @return the startPos
     */
    public int getStartPos() {
        return startPos;
    }

    /**
     * @return the endPos
     */
    public int getEndPos() {
        return endPos;
    }

    /**
     * @return the strand
     */
    public Strand getStrand() {
        return strand;
    }

    /**
     * @return the mapQ
     */
    public int getMapQ() {
        return mapQ;
    }

    /**
     * @return the rNext
     */
    public String getrNext() {
        return rNext;
    }

    /**
     * @return the tLen
     */
    public int gettLen() {
        return tLen;
    }

    /**
     * @return the seq
     */
    public String getSeq() {
        return seq;
    }

    /**
     * @return the qual
     */
    public String getQual() {
        return qual;
    }

    /**
     * @return the bowtie tag string that may be in col 11
     */
    public String getBowtieTagString() {
        return bowtieTagString;
    }

    /**
     * @return the cigar
     */
    public String getCigar() {
        return cigar;
    }

    /**
     * @return the pNext
     */
    public int getpNext() {
        return pNext;
    }

    /**
     * @return the header
     */
    public Boolean isHeaderLine() {
        return header;
    }

    /**
     * @return the mapped
     */
    public Boolean isMappedRead() {
        return mapped;
    }

    /**
     * @param startPos the startPos to set
     */
    public void setStartPos(int startPos) {
        this.startPos = startPos;
    }

    /**
     * @param endPos the endPos to set
     */
    public void setEndPos(int endPos) {
        this.endPos = endPos;
    }

  /**
   * @return the mdString
   */
  public String getMDString() {
    return mdString;
  }

  /**
   * @param mdString the mdString to set
   */
  public void setMDString(String mdString) {
    this.mdString = mdString;
  }

  /**
   * @param bowtieTagString the bowtieTagString to set
   */
  public void setTagString(String mdStr) {
    this.mdString = mdStr;
  }


  public static void main(String args[]){
      String samString = "SRR1295425.91775966\t337\tchr2\t70252903\t0\t200H27M1D23M"
              + "\tchr19\t49169091\t0\t"
              + "AGCTAATTTTTGTATTCTTTGTAGAGAGGAGGTCTCCCTATGTTGCCCAG\t"
              + "B<<BBBBBBBB<BBBBBBBB<BB<BBB<<B<<BBBBBB<BB<BBBB<<<'\tXA:Z:chr20";

      SAMEntry samEntry = new SAMEntry(samString);
      samEntry.getCigar();
      samEntry.splitCigarString();

  }



}

