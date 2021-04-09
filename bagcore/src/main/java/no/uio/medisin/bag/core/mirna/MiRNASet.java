/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.mirna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * Stores a complete miRBase release of miRNA data. It doesn't load any information
 * about the parent pre-miRNA
 * 
 * @author sr
 */
public class MiRNASet {
    
    static  Logger               logger                              = LogManager.getLogger();
    
    private List<MiRNA>          miRBaseMiRNAList                    = new ArrayList<>();

    
    public void loadMiRNAsFromFA(String fastaFilePath)throws IOException{
            
        try{
            logger.info("reading miRBase fasta reference file <" +  fastaFilePath + ">");
            BufferedReader brFA = new BufferedReader(new FileReader(new File(fastaFilePath)));
            String lineFA = null;
            while ((lineFA = brFA.readLine())!=null){

                String seq = brFA.readLine().trim();
                MiRNA miRNA = new MiRNA();
                miRNA.setName(lineFA);
                miRNA.setSeq(seq);
                this.miRBaseMiRNAList.add(miRNA);

            }
            brFA.close();
            logger.info("read " + miRBaseMiRNAList.size() + " entries");
            logger.info("done");

        }
        catch(IOException ex){
            logger.error("error reading miRBase fasta reference file <" +  fastaFilePath + ">\n" + ex.toString());
            throw new IOException("error reading miRBase fasta reference file <" +  fastaFilePath + ">");
        }
    }
    /**
     * 1. 
     * load miRNA specs (name and Chromosome position) from GFF file
     * downloaded from miRBase
     * Because different releases of miRBase use different releases of the
     * reference genome, we have to track both miRBase and genome reference IDs
     * 2.
     * Load Sequence from Mature.fa file
     * 
     * The pre-miRNA data is not loaded
     * 
     * @param host
     * @param gffFilePath           
     * @param fastaFilePath         
     * @throws IOException
     * 
     */
    public void loadMiRBaseData(String host, String gffFilePath, String fastaFilePath) throws IOException{

        HashMap <String, String> miRBaseSeq = null;
        if(!host.isEmpty()){
          miRBaseSeq = new HashMap();
          try{
              logger.info("reading miRBase fasta reference file <" +  fastaFilePath + ">");
              BufferedReader brFA = new BufferedReader(new FileReader(new File(fastaFilePath)));
              String lineFA = null;
              while ((lineFA = brFA.readLine())!=null){

                  String seq = brFA.readLine().trim();
                  String entryHost = lineFA.split(" ")[0].substring(1).split("-")[0].trim();
                  if(entryHost.equals(host)){
                      String mimatID = lineFA.split(" ")[1].trim();
                      miRBaseSeq.put(mimatID, seq);
                  }

              }
              brFA.close();
              logger.info("read " + miRBaseSeq.size() + " entries");
              logger.info("done");

          }
          catch(IOException ex){
              logger.error("error reading miRBase fasta reference file <" +  fastaFilePath + ">\n" + ex);
              throw new IOException("error reading miRBase fasta reference file <" +  fastaFilePath + ">");
          }
          
        }else{
          logger.info("No fasta reference file specified ");          
        }
        
        
        
        try{
            String line = null;
            logger.info("reading miRBase gff reference file <" +  gffFilePath + ">");
            BufferedReader brMiR = new BufferedReader(new FileReader(new File(gffFilePath)));
                while((line = brMiR.readLine())!= null){

                    if(line.startsWith("#")) continue;
                    if(line.contains("miRNA_primary_transcript")) continue;
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
                    //logger.debug(line);
                    MiRNA miRNA = new MiRNA();
                    miRNA.parseMiRBaseGFFEntry(line);
                    // this may need revising. In some cases we dont care if there is no sequence, we still want the feature
                    String alias =miRNA.getMimatAlias();

                    if(miRBaseSeq != null){
                      String seq = miRBaseSeq.get(alias);
                      miRNA.setSeq(seq);
                      if(seq == null) {
                          logger.warn("no sequence specified for entry <" + alias + ">. Skipping");
                      }
                    }else{
                      
                    }
                    String mirName = miRNA.getName();
                    //setSeq() is not definied in MiRNA class! so miRNA.getSeuence() will return null
                    //for merging, change back to setSeq()

                    boolean isDuplicate = checkDuplicate(mirName, this.miRBaseMiRNAList);
                    if (isDuplicate){
                        miRNA.setName(mirName+".2");
                    }
                    this.miRBaseMiRNAList.add(miRNA);

                }
            brMiR.close();
            logger.info("read " + miRBaseMiRNAList.size() + " miRNA entries");
        }
        catch(IOException ex){
            logger.error("error reading miRBase gff reference file <" +  gffFilePath + ">\n" + ex.toString());
            throw new IOException("error reading miRBase gff reference file <" +  gffFilePath + ">");
        }
        
    }
    

    
    
    
    /**
     * Does the read sufficiently overlap a defined miRNA entry?
     * 
     * @param start
     * @param stop
     * @param chr
     * @param strand
     * @param bleed     int : specifies how much a read can 'miss' an entry
     *                        and still be counted
     * 
     * @return MiRNAFeature
     */
    public MiRNA doesReadOverlapKnownMiRNA(int start, int stop, String chr, String strand, int bleed){
        
        for(MiRNA miRBaseEntry: this.getMiRBaseMiRNAList()){
            if (miRBaseEntry.chromosomeMatch(chr)){
                if(strand.equals(miRBaseEntry.getStrand())){
                    
                    if( java.lang.Math.abs(start - miRBaseEntry.getStartPos()) <= bleed){
                        
                        if( java.lang.Math.abs(stop - miRBaseEntry.getEndPos()) <= bleed){
                            return miRBaseEntry;
                        }

                    }

                }
            }
        }
        
        return null;
        
    }
    
    
    
    
    /**
     * search the miRNA list by miRNA name
     * 
     * @param querymiRNA
     * @return 
     */
    public MiRNA findMiRNAbyName(String querymiRNA){
        for(MiRNA miRBaseEntry: this.getMiRBaseMiRNAList()){
            if(miRBaseEntry.getName().equals(querymiRNA))
                return miRBaseEntry;
        }        
        return null;
    }
    
    
    
    /**
     * get the number of entries in this release
     * 
     * @return 
     */
    public int getNumberOfEntries(){
        return getMiRBaseMiRNAList().size();
    }

    /**
     * @return the miRBaseMiRNAList
     */
    public List<MiRNA> getMiRBaseMiRNAList() {
        return miRBaseMiRNAList;
    }

    private boolean checkDuplicate(String mirnaName, List<MiRNA> miRBaseMiRNAList) {
        boolean isDuplicate = false;
        for(MiRNA miRBaseEntry: this.getMiRBaseMiRNAList()){
            if (mirnaName.equals(miRBaseEntry.getName())) {
                isDuplicate = true;
            }
        }
        return isDuplicate;
    }
    
}
