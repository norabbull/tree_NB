/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.sequence;

import java.util.ArrayList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * a simple sequence object that contains some basic functionality
 * and stores basic sequence features
 * 
 * @author weibo
 */
public class SimpleSeq {
    
    static Logger logger = LogManager.getRootLogger();
    
    private String          id="";
    private String          seq="";
    private int             length=0;
    private int             start=0;
    private int             end=0;
    private String          chromosome="";
    private String          strand;
    private String          name="";
    private String          accessionNumber;


    /**
     * Empty Class Constructor
     */
    public SimpleSeq(){

    }
    
    /**
     * copy constructor
     * 
     * @param nSimpleSeq 
     */
    public SimpleSeq(SimpleSeq nSimpleSeq){
        this.chromosome         = nSimpleSeq.chromosome;
        this.accessionNumber    = nSimpleSeq.accessionNumber;
        this.strand             = nSimpleSeq.strand;
        this.start              = nSimpleSeq.start;
        this.end                = nSimpleSeq.end;
        this.length             = nSimpleSeq.length;
        this.name               = nSimpleSeq.name;
        this.id                 = nSimpleSeq.id;
        this.seq                = nSimpleSeq.seq;
    }
    

    
    /**
     * Constructor 
     * 
     * @param id
     * @param seq 
     */
    public SimpleSeq(String id,String seq){
        this.name   = id;
        this.id     = id;
        this.seq    = seq;
        this.length = seq.length();
        
    }
    
    
    /*
    return sequence in FASTA format
    */
    public String toFastA(){
        return ">" + this.getName() + "|" + this.getId() + System.lineSeparator()
                + this.getSeq() + System.lineSeparator();
    }
    
    /**
     * summarize the sequence properties
     * 
     * @return 
     */
    @Override
    public String toString(){
        String str = "ID:" + this.getId() + "\t"
                    + "Name:" + this.getName() + "\t"
                    + "Start:" + this.getStartPos() + "\t"
                    + "End:" + this.getEndPos() + "\t"
                    + "Len:" + this.getLength() + "\t"
          + "\n";
        return str;
    }


    /**
     * return subsequence at specified position
     * 
     * @param start
     * @param stop
     * @return 
     */
    public String getSequence(int start, int stop){
        if(start>=0 && stop < seq.length())
            return seq.substring(start, stop+1);

        return null;
    }

    
    /**
     * convert RNA sequence to DNA
     * 
     * @param seq
     * @return 
     */
    public static String rna2dna(String rnaseq){
        return rnaseq.replace("u", "t").replace("U", "T");
    }
    
    
    /**
     * convert DNA sequence to RNA
     * 
     * @param dnaseq
     * @return 
     */
    public static String dna2rna(String dnaseq){
        return dnaseq.replace("t", "u").replace("T", "U");
    }
    
    
    /**
     * Reverse complement RNA/DNA sequence
     * if sequence contains 'U' or 'u' it is assumed the sequence is RNA
     * 
     * @param seqIn
     * @return 
     */
    public static String complement(String seqIn){
        
        StringBuilder Complement = new StringBuilder();
        char [] strReversed = new StringBuilder(seqIn).reverse().toString().toCharArray();

        for (char nt: strReversed) {
            switch (nt){
                case 'a':
                    if(seqIn.toLowerCase().contains("u"))
                        Complement.append("u");
                    else
                        Complement.append("t");
                    break;
                    
                case 'A':
                    if(seqIn.toLowerCase().contains("u"))
                        Complement.append("U");
                    else
                        Complement.append("T");
                    break;
                    
                case 'c':
                    Complement.append("g");
                    break;
                    
                case 'C':
                    Complement.append("G");
                    break;
                    
                case 'g':
                    Complement.append("c");
                    break;
                    
                case 'G':
                    Complement.append("C");
                    break;
                    
                case 't':
                    Complement.append("a");
                    break;
                    
                case 'T':
                    Complement.append("A");
                    break;
                    
                case 'u':
                    Complement.append("a");
                    break;
                    
                case 'U':
                    Complement.append("A");
                    break;
                    
                default:
                    Complement.append("N");
                    break;
                    
            }

        }
        return Complement.toString();
    }


    /**
     * count occurrence of a specific nucleotide within the specified sequence
     * 
     * @param qSeq
     * @param nt
     * @return 
     */
    public static int NTcount(String qSeq, char nt){
        char[] cs=qSeq.toCharArray();
        int n=0;
        for(char c:cs){
            if(c==nt){
                n++;
            }
        }
        return n;
    }
    
    
    
    /**
     * calculate the Shannon entropy of the specified sequence
     * 
     * @param col
     * @return 
     */
    public static double shannonEntropy(String col){
        double aCount = SimpleSeq.NTcount(col.toUpperCase(), 'A');
        double cCount = SimpleSeq.NTcount(col.toUpperCase(), 'C');
        double gCount = SimpleSeq.NTcount(col.toUpperCase(), 'G');
        double tCount = SimpleSeq.NTcount(col.toUpperCase(), 'T') 
                + SimpleSeq.NTcount(col.toUpperCase(), 'U');
        
        double e = 0.0;
        if(aCount > 0)
            e += aCount*(Math.log(aCount)/Math.log(2.0));
        if(cCount > 0)
            e += cCount*(Math.log(cCount)/Math.log(2.0));
        if(gCount > 0)
            e += gCount*(Math.log(gCount)/Math.log(2.0));
        if(tCount > 0)
            e += tCount*(Math.log(tCount)/Math.log(2.0));

        return e;
    }
    
    
    /**
     * split the sequence into fragments according to specified step and window size
     * 
     * @param step
     * @param window
     * @return 
     */
    public ArrayList<SimpleSeq> splitSequence(int step, int window){
        
        ArrayList<SimpleSeq> fragmentList = new ArrayList<>();
        length=seq.length();
        int n=(length-window)/step+1;
        if(n<1) n=1;

        int fragStart=0;
        int fragEnd; 
        for(int i=0;i<n;i++){
            
            if(fragStart>=length) break;
            fragEnd = fragStart + window;
            if(fragEnd>length) fragEnd=length;
            String fragID = name + "_" + (fragStart+1) + "-"+fragEnd;
            String subseq = seq.substring(fragStart,fragEnd);

            SimpleSeq frag = new SimpleSeq(fragID,subseq);
            frag.setAbsStartInQuerySeq(fragStart+1);// count from 1
            frag.setAbsEndInQuerySeq(fragEnd); //count from 1
            frag.setName(fragID);
            fragmentList.add(frag);
            fragStart+= step;
            
        }
        return fragmentList;
    }
    
    
    /**
     * calculate A fraction in the query sequence 
     * @param qSeq
     * @return 
     */
    public static double Afraction(String qSeq){
        return (double)(NTcount(qSeq, 'A') + NTcount(qSeq, 'a'))/(double)qSeq.length();
    }
    
    
    /**
     * calculate A fraction in the sequence 
     * @return 
     */
    public double Afraction(){
        return (double)(NTcount('A') + NTcount('a'))/(double)this.getLength();
    }
    
    
    
    
    /**
     * calculate C fraction in the query sequence 
     * @param qSeq
     * @return 
     */
    public static double Cfraction(String qSeq){
        return (double)(NTcount(qSeq, 'C') + NTcount(qSeq, 'c'))/(double)qSeq.length();
    }
    
    
    /**
     * calculate C fraction in the sequence 
     * @return 
     */
    public double Cfraction(){
        return (double)(NTcount('C') + NTcount('c'))/(double)this.getLength();
    }
    
    
    
    
    /**
     * calculate G fraction in the query sequence 
     * @param qSeq
     * @return 
     */
    public static double Gfraction(String qSeq){
        return (double)(NTcount(qSeq, 'G') + NTcount(qSeq, 'g'))/(double)qSeq.length();
    }
    
    
    /**
     * calculate G fraction in the sequence 
     * @return 
     */
    public double Gfraction(){
        return (double)(NTcount('G') + NTcount('g'))/(double)this.getLength();
    }
    
    
    /**
     * calculate GC fraction
     * @return 
     */
    public double GCfraction(){
        return Cfraction() + Gfraction();
    }
    
    
    /**
     * calculate GC fraction
     * @param qSeq
     * @return 
     */
    public static double GCfraction(String qSeq){
        return (float) (SimpleSeq.NTcount(qSeq, 'g') + SimpleSeq.NTcount(qSeq, 'G') 
                + SimpleSeq.NTcount(qSeq, 'c') + SimpleSeq.NTcount(qSeq, 'C'))
                /(float) qSeq.length();
    }
    
    
    /**
     * calculate T fraction in the query sequence 
     * @param qSeq
     * @return 
     */
    public static double Tfraction(String qSeq){
        return (double)(NTcount(qSeq, 'T') + NTcount(qSeq, 't'))/(double)qSeq.length();
    }
    
    
    /**
     * calculate T fraction in the sequence 
     * @return 
     */
    public double Tfraction(){
        return (double)(NTcount('T') + NTcount('t'))/(double)this.getLength();
    }
    
    
    
   
    /**
     * calculate U fraction in the sequence 
     * @return 
     */
    public double Ufraction(){
        return (double)(NTcount('U') + NTcount('u'))/(double)this.getLength();
    }
    
    
    
   
    /**
     * count occurrence of a specific nucleotide within the defined sequence
     * 
     * @param nt
     * @return 
     */
    public int NTcount(char nt){
        return NTcount(seq, nt);
    }

    /**
     * @return the id
     */
    public String getId() {
        return id;
    }

    /**
     * @param id 
     */
    public void setID(String id) {
        this.id = id;
    }

    /**
     * @return the seq
     */
    public String getSeq() {
        return seq;
    }

    /**
     * @param seq 
     */
    public void setSeq(String seq) {
        this.seq = seq;
        this.setLength(seq.length());
    }

    /**
     * @return length
     */
    public int getLength() {
        return length;
    }

    /**
     * @param length 
     */
    public void setLength(int length) {
        this.length = length;
    }

    /**
     * @return start
     */
    public int getStartPos() {
        return start;
    }

    /**
     * @param start
     */
    public void setAbsStartInQuerySeq(int start) {
        this.setStartPos(start);
    }

    /**
     * @return end
     */
    public int getEndPos() {
        return end;
    }

    /**
     * @param end
     */
    public void setAbsEndInQuerySeq(int end) {
        this.setEndPos(end);
    }

    /**
     * @return name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name 
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the accessionNumber
     */
    public String getAccessionNumber() {
        return accessionNumber;
    }

    /**
     * @param accessionNumber the accessionNumber to set
     */
    public void setAccessionNumber(String accessionNumber) {
        this.accessionNumber = accessionNumber;
    }

    /**
     * @return the chromosome
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * @param chromosome the chromosome to set
     */
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * @return the strand
     */
    public String getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(String strand) {
        this.strand = strand;
    }

    
    public static void main(String args[]){
        
        SimpleSeq simpleSeq = new SimpleSeq();
        simpleSeq.setSeq("AAAACGCGCGCAAAA");
        logger.info((simpleSeq.Cfraction()+simpleSeq.Gfraction())/(float)simpleSeq.length);
         
        
    }

    /**
     * @param start the start to set
     */
    public void setStartPos(int start) {
        this.start = start;
        if(this.getLength()>0)
            this.setEndPos(start+this.getLength()-1);
    }

    /**
     * @param end the end to set
     */
    public void setEndPos(int end) {
        this.end = end;

    }
    


}