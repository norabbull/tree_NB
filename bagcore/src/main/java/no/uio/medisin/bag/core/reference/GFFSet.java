/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.reference;

import no.uio.medisin.bag.core.mapping.MappedRead;
import no.uio.medisin.bag.core.sequence.Strand;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Iterator;
import org.apache.commons.math3.stat.Frequency;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Stores an array of GFF Entries and allows searching the set for reads that match
 * 
 * Note:
 * This isn't a very efficient implementation. It would be better to store
 * this data in a tree to allow faster searching
 * 
 * @author simonray
 */
public class GFFSet {
    static Logger logger = LogManager.getLogger();
    private final ArrayList<GFFEntry> GFFEntries;
    
    public GFFSet(){
        GFFEntries = new ArrayList<>();
    }
    
    
    /**
     * add GFFEntry to the set
     * 
     * @param gffEntry GFFEntry
     */
    public void addEntry(GFFEntry gffEntry){
        getGFFEntries().add(gffEntry);
    }
    
    
    /**
     * read specified GFF file
     * 
     * @param filename String
     * @return int number of lines read
     * 
     * @throws IOException if error encountered reading the GFF file
     */
    public int readGFF(String filename) throws IOException, RuntimeException{
        String line = null;
        int lineCount = 0;
        try{
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            while ((line = br.readLine()) != null) {
                if(GFFEntry.isCommentLine(line)==false) 
                    this.addEntry(new GFFEntry(line));
                lineCount ++;
            }
        }
        catch(IOException exIO){
            logger.error("error parsing GFF file " + filename);
            logger.error("exception thrown on line " + lineCount);
            logger.error(line);
            logger.error(exIO);
            throw new IOException("error parsing GFF file " + filename + "\n" 
            + "exception was thrown on line " + lineCount);  
        }
        return lineCount;
    }
    
    
    

    /**
     * read specified GFF file
     * 
     * @param filename String
     * @param featureList ArrayList of Strings
     * 
     * @return int number of lines read
     * 
     * @throws IOException if error encountered reading the GFF file
     */
    public int readGFF(String filename, ArrayList<String> featureList) throws IOException, RuntimeException{
        String line = null;
        int lineCount = 0;
        try{
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            while ((line = br.readLine()) != null) {
                if(GFFEntry.isCommentLine(line)==false && featureList.contains(line.split("\t")[GFFEntry.GFF_TYPE].trim()))
                    this.addEntry(new GFFEntry(line));
                lineCount ++;
            }
        }
        catch(IOException exIO){
            logger.error("error parsing GFF file " + filename);
            logger.error("exception thrown on line " + lineCount);
            logger.error(line);
            logger.error(exIO);
            throw new IOException("error parsing GFF file " + filename + "\n" 
            + "exception was thrown on line " + lineCount);  
        }
        return lineCount;
    }
    
    
    

    /**
     * 
     * return entry within the specified region
     * 
     * @param start int
     * @param stop int
     * @param strand String
     * @param chr String
     * @param bleed int
     * @return GFFEntry hit to GFFEntry
     */
    public GFFEntry findMatch(int start, int stop, String strand, String chr, int bleed){
        Iterator itGF = getGFFEntries().iterator();
        while(itGF.hasNext()){
            GFFEntry gffEntry = (GFFEntry)itGF.next();
            if(gffEntry.getStrand().equals(strand)
                && gffEntry.getAttr().equals(chr)
                && Math.abs(gffEntry.getStart()-start) < bleed 
                && Math.abs(gffEntry.getStop() - stop) < bleed
                    ){
                return gffEntry;
            }
        }
        return null;
    }
    
    
    
    
    /**
     * seek matching entry based on SeqID
     * 
     * @param sID String 
     * @return GFFEntry hit to GFFEntry
     */
    public GFFEntry findEntryByID(String sID){
        for(GFFEntry gffEntry: getGFFEntries()){
            if(gffEntry.getFeatureID().equals(sID)) return gffEntry;
        }
        return null;
    }
    
    
    
    
    
    /**
     * Does the specified region overlap any feature?
     * 
     * @param start int
     * @param stop int 
     * @param strand Strand
     * @param chr String
     * @param bleed int
     * @param featureType String
     * @return Boolean T/F region contains feature
     */
    public Boolean doesRegionContainFeature(int start, int stop, Strand strand, String chr, String featureType, int bleed){
        
        for(GFFEntry gffEntry: getGFFEntries()){
            if(gffEntry.doesRegionOverlap(start, stop, strand, chr, featureType, bleed)){
                return true;
            }
        }
        return false;
        
    }
    
    
    
    /**
     * 
     * @param queryGFFEntry GFFEntry query
     * @return Boolean T/F region contains feature
     */
    public Boolean doesGFFSetContainFeature(GFFEntry queryGFFEntry){
        for(GFFEntry thisEntry:this.getGFFEntries()){
            if(thisEntry.doesFeatureContainRegion(queryGFFEntry))
                return true;
        }
        return false;
    }
    
    
    /**
     * Does the specified region overlap a feature regardless of its type?
     * 
     * @param start int
     * @param stop int 
     * @param strand Strand
     * @param chr String
     * @param bleed int
     * @return Boolean T/F region contains feature
     */
    public Boolean doesRegionContainFeature(int start, int stop, Strand strand, String chr, int bleed){
        
        for(GFFEntry gffEntry: getGFFEntries()){
            if(gffEntry.doesRegionOverlap(start, stop, strand, chr, bleed)){
                return true;
            }
        }
        return false;
        
    }
    
    
    
    
    
    /**
     * find if specified region overlaps a gffEntry in this set
     * 
     * @param start int
     * @param stop int
     * @param strand Strand
     * @param chr String
     * @param bleed int 
     * @param featureType String
     * 
     * @return GFFEntry GFFEntry hit
     */
    public GFFEntry findOverlappingFeature(int start, int stop, Strand strand, String chr, String featureType, int bleed){
        
        for(GFFEntry gffEntry: getGFFEntries()){
            if(gffEntry.doesRegionOverlap(start, stop, strand, chr, featureType, bleed)){
                return gffEntry;
            }
        }
        
        return null;
    }
    
    
    
    
    
    /**
     * find if specified mappedRead overlaps a gffEntry in this set
     * 
     * @param queryRead MappedRead
     * @param bleed String
     * @param featureType String
     * @return GFFEntry GFFEntry hit
     */
    public GFFEntry findOverlappingFeature(MappedRead queryRead, String featureType, int bleed){
        
        for(GFFEntry gffEntry: getGFFEntries()){
            if(gffEntry.doesRegionOverlap(queryRead.getStartPos(), queryRead.getEndPos(), queryRead.getStrand(), queryRead.getChr(), featureType, bleed)){
                return gffEntry;
            }
        }
        
        return null;       
    }
    
    
    
    
    
    /**
     * find if specified mappedRead overlaps a gffEntry in this set
     * 
     * @param queryRead MappedRead
     * @param featureType String
     * @return GFFEntry GFFEntry hit
     */
    public GFFEntry doesFeatureContainRegion(MappedRead queryRead, String featureType){
        
        for(GFFEntry gffEntry: getGFFEntries()){
            if(gffEntry.doesFeatureContainRegion(queryRead.getStartPos(), queryRead.getEndPos(), queryRead.getStrand(), queryRead.getChr(), featureType)){
                return gffEntry;
            }
        }
        
        return null;       
    }
    
    
    
    
    
    /**
     * 
     * @param bwLD BufferedWriter
     * @param start int
     * @param stop int 
     * @throws IOException throws exception if error encountered writing file
     */
    public void writeLengthDistribution(BufferedWriter bwLD, int start, int stop) throws IOException{
        Frequency freqDist = new Frequency();
        Iterator itGF = getGFFEntries().iterator();
        while(itGF.hasNext()){
            GFFEntry gffEntry = (GFFEntry)itGF.next();            
            freqDist.addValue(gffEntry.getStop() - gffEntry.getStart() + 1);
        }

        for(int l=start; l<=stop; l++){
            bwLD.write(l + "\t" + freqDist.getCount(l) + "\n");
        }
        
    }
    
    
    /**
     * write out the features in GFF format
     * 
     * @param bwFT BufferedWriter
     * @param ID String
     * @throws IOException throws exception if error encountered while writing GFF file
     */
    public void writeFeaturesAsGFF3(BufferedWriter bwFT, String ID) throws IOException{
            bwFT.write("# created " + new Timestamp((new java.util.Date()).getTime()));
            bwFT.write("# from GFFSet");
            bwFT.write("# " + ID);
            bwFT.write("# ");
            bwFT.write("# \n");
        
        Iterator itGF = getGFFEntries().iterator();
        while(itGF.hasNext()){
            GFFEntry gffEntry = (GFFEntry)itGF.next();            
            bwFT.write(gffEntry.toGFF3String() + "\n");
        }
    }
    
    
    
    
    /**
     * write features in FASTA format, one / line
     * 
     * @param bwFA BufferedWriter 
     * @throws IOException throws exception if error encountered while writing FASTA file
     */
    public void writeFeaturesAsFastA(BufferedWriter bwFA) throws IOException{
        Iterator itGF = getGFFEntries().iterator();
        while(itGF.hasNext()){
            GFFEntry gffEntry = (GFFEntry)itGF.next();            
            logger.info(gffEntry.getFeatureID());
            bwFA.write(gffEntry.toMirbaseFastAString() + "\n");
        }
    }
    
    
    /**
     * return the number of entries in the Array
     * 
     * @return int NumberofEntries
     */
    public int getNoOfEntries(){
        return this.getGFFEntries().size();
    }

    /**
     * @return ArrayList  of GFFEntries
     */
    public ArrayList<GFFEntry> getGFFEntries() {
        return GFFEntries;
    }
}
