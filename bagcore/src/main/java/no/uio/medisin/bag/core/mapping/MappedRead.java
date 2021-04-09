/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mapping;

import java.util.Objects;
import no.uio.medisin.bag.core.sequence.StrandString;
import no.uio.medisin.bag.core.sequence.ChromosomeString;
import no.uio.medisin.bag.core.sequence.Strand;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * A MappedRead instance contains the basic information required to specify the mapped location
 * i.e. 
 *   a referenceID (e.g.Chromosome)
 *   a start position
 *   a stop position
 *   a strand
 * 
 * and optionally, 
 *   a read count
 * 
 * It is a subset of a @see SAMEntry but only contains location information
 * and implements Comparable.
 * The chromosome is stored as an int, it is assumed that the call which created 
 * the instance handled the mapping from String to int. 
 * 
 * 
 * @author simonray
 */
public class MappedRead implements Comparable<MappedRead>{
    static Logger logger = LogManager.getLogger();

    private int startPos;
    private int endPos;
    private int chr;
    private Strand strand;
    private int count;

    
    
    public MappedRead(int s, int e, String c, String t, int n){
        startPos = s;
        endPos = e;
        chr = ChromosomeString.GetChromNumber(c);
        strand = StrandString.guessStrand(t);
        count = n;
    }
    
    
    
    @Override
    public boolean equals(Object o){
     if (this == o)
        return true;
    // null check
    if (o == null)
        return false;
    // type check and cast
    if (getClass() != o.getClass())
        return false;
    MappedRead m = (MappedRead) o;
    // field comparison
    return Objects.equals(this.chr, m.chr)
            && Objects.equals(this.startPos, m.startPos)
            && Objects.equals(this.endPos, m.endPos)
            && Objects.equals(this.strand, m.strand);       
    }
    
    
    /**
     * we don't include count in the hashCode since we only consider location
     * when deciding if a mapped read is equivalent
     * 
     * @return int hashcode
     */
    @Override
    public int hashCode() {
        int prime = 31;
        int result = startPos*prime + endPos*prime*prime + chr*prime*prime*prime;

        result = prime * result +  strand.hashCode();
        result = prime * result + ((strand == null) ? 0 : strand.hashCode());
        return result;
    }
    
    
    /**
     * sort by chromosome, strand and then position
     * 
     * @param mappedRead the query MappedRead
     * @return int outcome of comparison
     */
    @Override
    public int compareTo(MappedRead mappedRead){
        
        int i=0;
        //int l = Math.min(mappedRead.chr.length(), chr.length());
        try{
            if (mappedRead.chr==chr && mappedRead.strand == strand)
            {
                return startPos - mappedRead.startPos ;
            }

            if(mappedRead.strand != strand){
                if(strand == Strand.PLUS && mappedRead.strand == Strand.MINUS){
                    return -1;
                }else{
                    return +1;
                }            
            }
            
            /*
            trying to compare chromosomes is a little more tricky because they
            may be stored in different ways 
            If the string contains 'chr', then remove this and see if we are left
            with a number. Then we can do a simple compare - NumberFormatException
            */
            /*
            int mapI;
            int thisI;
            
            if(mappedRead.chr.toUpperCase().contains("CHR")){
                int mStart = mappedRead.chr.toUpperCase().indexOf("CHR");
                String mapChr = mappedRead.chr.toUpperCase().substring(mStart+3, mappedRead.chr.toUpperCase().length());
                int iStart = mappedRead.chr.toUpperCase().indexOf("CHR");
                String thisChr = this.chr.toUpperCase().substring(iStart+3, this.chr.toUpperCase().length());
                try{
                    mapI = Integer.parseInt(mapChr);
                    thisI = Integer.parseInt(thisChr);
                    return (thisI - mapI);
                }
                catch(RuntimeException exTR){
                    logger.info(mappedRead.toString() + "|" + this.toString() + "\t(" + l + ")");
                    logger.info(exTR);
                }
            }
            else{
                try{
                    mapI = Integer.parseInt(mappedRead.chr);
                    thisI = Integer.parseInt(this.chr);
                    return (thisI - mapI);
                }
                catch(RuntimeException exTR){
                    logger.info(mappedRead.toString() + "|" + this.toString() + "\t(" + l + ")");
                    logger.info(exTR);
                }

            } 
            while(i<l && mappedRead.chr.charAt(i)==(chr.charAt(i)) ){
                i++;
            }
            if(i==l) 
                i--;
                    */
        }
        catch(RuntimeException exAr){
            logger.info("exception comparing: ");
            logger.info(mappedRead.toString());
            logger.info(" &");
            logger.info(this.toString());
            //logger.info("i= " + i  + "\t" + this.chr  + " <-> " + mappedRead.chr + "\tmin length = " + l );
            logger.info(exAr);
        }
        //return  chr.charAt(i) - mappedRead.chr.charAt(i) ;
        return mappedRead.chr - chr;
    }
    
    
    
    
    /**
     * print the mapped read information
     * 
     * @return MappedRead as a String
     */
    @Override
    public String toString(){
        return    "chr "    + chr + "\t" 
                + "start "  + startPos + "\t"
                + "end "    + endPos + "\t"
                + "strand " + strand + "\t"
                + "count "  + getCount();
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
     * @return the chr
     */
    public String getChr() {
        return ChromosomeString.getChromosomeAsString(chr);
    }

    /**
     * @return the strand
     */
    public Strand getStrand() {
        return strand;
    }

    /**
     * @return the count
     */
    public int getCount() {
        return count;
    }
    
}
