/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.sequence;

import java.util.ArrayList;

/**
 * 
 * A class to hold a small alignment and extract slices in either direction
 * This isn't designed for handling large sequence sets or long alignments
 * 
 * @author simonray
 */
public class SimpleSmallAlignment {
    
    int alignmentLen;
    SimpleSequenceSet simpleSeqSet;
    
    
    /**
     * add a set of SimpleSequences to the existing alignment
     * 
     * @param newSeqs
     * @return 
     */
    public int addSequences(SimpleSequenceSet newSeqs){
        newSeqs.getSeqs().stream().forEach((sSeq) -> {
           simpleSeqSet.addSequence(sSeq);
        });
       return simpleSeqSet.getSeqs().size();
    } 

    
    
    /**
     * add a set of SimpleSequences to the existing alignment
     * 
     * @param newSeq
     * @return 
     */
    public int addSequence(SimpleSeq newSeq){
           simpleSeqSet.addSequence(newSeq);
        
       return simpleSeqSet.getSeqs().size();
    } 

    
    /**
     * return sequence in specified row
     * @param row
     * @return 
     */
    public SimpleSeq getSequence(int row){
        if(row<simpleSeqSet.getSeqs().size())
            return simpleSeqSet.getSeqs().get(row);
        return null;
    }
    
    
    /**
     * we assume that only A,C,G, T/U are present
     * If the base is absent it isn't counted
     * if there is no consensus sequence we mark the site with '?'
     * 
     * @return consensusSequence
     * 
     */
    public String findConsensusSequence(){
        String consensus = "";
        
        for(int nt=0; nt<this.getAlignmentLength(); nt++){
            consensus = consensus.concat(getConsensusBase(this.getCol(nt)));
        }
        return consensus;
    }
    
    
    /**
     * pull out specified column from the alignment
     * 
     * @param col
     * @return 
     */
    private SimpleSeq getCol(int col){
        String colString = "";
        for(int s=0; s<this.getNumberOfSequences(); s++){
            colString = colString.concat(this.simpleSeqSet.getSeqs().get(s).getSequence(col, col+1));
        }
        return new SimpleSeq(Integer.toString(col), colString);
    }
    
    
    
    /**
     * get the consensus base for the supplied column
     * This is stored as a SimpleSeq, but it is representing a column within the alignment
     * When there is a tie, return the IUPAC representation.
     * @param seq
     * @return consensus nucleotide
     */
    private String getConsensusBase(SimpleSeq seq){
        /*
    
            IUPAC nucleotide code	Base
            R	A or G
            Y	C or T
            S	G or C
            W	A or T
            K	G or T
            M	A or C
            B	C or G or T
            D	A or G or T
            H	A or C or T
            V	A or C or G
            N	any base
            . or -	gap    
        */
        
        
        double aFrac = seq.Afraction();
        double cFrac = seq.Cfraction();
        double gFrac = seq.Gfraction();
        double tFrac = seq.Tfraction();
        
        String maxBase = "";
        double maxFrac = 0.0;
        if(aFrac > maxFrac){
            maxFrac = aFrac;
            maxBase = "A";
        }

        if(cFrac >= maxFrac){
            maxFrac = cFrac;
            if(cFrac == maxFrac){
                maxBase = maxBase.concat("C");
            }
        }

        if(gFrac >= maxFrac){
            maxFrac = gFrac;
            if(gFrac == maxFrac){
                maxBase = maxBase.concat("G");
           }
         }

        if(tFrac >= maxFrac){
            maxFrac = tFrac;
            if(tFrac == maxFrac){
                maxBase = maxBase.concat("T");                
            }
        }

        
        if(maxBase.length()==1)
            return maxBase;
        
        if(maxBase.length()==2)
        {
            switch(maxBase){
                case "AC":
                    maxBase = "M";
                    break;
                
                case "AG":
                    maxBase = "R";
                    break;
                
                case "AT":
                    maxBase = "W";
                    break;
                
                case "CG":
                    maxBase = "S";
                    break;
                
                case "CT":
                    maxBase = "Y";
                    break;
                
                case "GT":
                    maxBase = "K";
                    break;
                
                default:
                    maxBase = "?";
                    break;
            }
            return maxBase;
        }
        
        if(maxBase.length()==3){
            switch(maxBase){
                case "ACG":
                    maxBase = "V";
                    break;
                    
                case "ACT":
                    maxBase = "H";
                    break;
                    
                case "AGT":
                    maxBase = "D";
                    break;
                    
                case "CGT":
                    maxBase = "B";
                    break;
                    
            }
            return maxBase;
        }
        
        if(maxBase.length()==4)
            maxBase="N";
        
        return maxBase;
        
    }
    
    
    /**
     * Return Alignment in FASTA format
     * 
     * @return 
     */
    public String asFastA(){
        String fastAString = "";
        for(SimpleSeq seq:simpleSeqSet.getSeqs()){
            fastAString = fastAString.concat(seq.toFastA());
        }
        
        return fastAString;
    }
    
    
    
    /**
     * return sequences
     * 
     * @return 
     */
    public ArrayList<SimpleSeq> getSequences(){
        return simpleSeqSet.getSeqs();
    }
    
    
    /**
     * return the number of sequences in the alignment
     * 
     * @return 
     */
    public int getNumberOfSequences(){
        return simpleSeqSet.getSeqs().size();
    }
    
    
    /**
     * 
     * get the alignment length
     * 
     * @return 
     */
    public int getAlignmentLength(){
        if(this.getNumberOfSequences()>0)
            return simpleSeqSet.getSeqs().get(0).getLength();
        return -1;
    }
}




