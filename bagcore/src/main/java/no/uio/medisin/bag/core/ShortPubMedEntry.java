/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core;

import java.io.IOException;
import java.util.ArrayList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author simonray
 */
public class ShortPubMedEntry {
    static Logger                   logger                  = LogManager.getLogger();


    private static final String     REFERENCE_AUTHOR        =   "RA";
    private static final String     REFERENCE_COMMENT       =   "RC";
    private static final String     REFERENCE_LOCATION      =   "RL";
    private static final String     REFERENCE_NUMBER        =   "RN";
    private static final String     REFERENCE_TITLE         =   "RT";
    private static final String     REFERENCE_XF            =   "RX";
    
    private static final String     PUBMED_MARKER           =   "PUBMED;";
    
    private int         pubmedID;
    private String      authors;
    private String      title;
    private String      journalRef;
    private String      refNum;




    /**
     * initialize the entry by parsing out the array of lines
     * 
     * @param lines ArrayList of Strings
     * @throws IOException when error is thrown trying to parse EMBL entries
     * 
     */
    public ShortPubMedEntry(ArrayList<String> lines) throws Exception{
        this();
        parseEMBLLines(lines);
    }
    
    
    
    /**
     * empty constructor
     * 
     */
    public ShortPubMedEntry() {
        pubmedID    = -1;
        authors     = "";
        title       = "";
        journalRef  = "";
        refNum      = "";

    }


    /**
     * input will be a series of lines with EMBL two letter field codes corresponding to
     * reference entries.  
     * 
     * For example:
     * 
     *      RN   [7]
     *      RX   PUBMED; 20062054.
     *      RA   Zisoulis DG, Lovci MT, Wilbert ML, Hutt KR, Liang TY, Pasquinelli AE, Yeo
     *      RA   GW;
     *      RT   "Comprehensive discovery of endogenous Argonaute binding sites in
     *      RT   Caenorhabditis elegans";
     *      RL   Nat Struct Mol Biol. 17:173-179(2010).
     * 
     * 
     * @param lines 
     */
    private void parseEMBLLines(ArrayList<String> lines) throws Exception{
    
        for(String line: lines){

            switch(line.trim().substring(0, 2)){
                
                case REFERENCE_AUTHOR:
                    this.setAuthors(this.getAuthors().concat(" " + line.substring(2).trim()));
                    break;
                    
                case REFERENCE_COMMENT: // for now, we dont store this.
                    break;
                    
                case REFERENCE_LOCATION:
                    this.setJournalRef(this.getJournalRef().concat(" " + line.substring(2).trim()));
                    break;
                    
                case REFERENCE_NUMBER: 
                    this.setRefNum(line.substring(3).trim());
                    break;
                    
                case REFERENCE_TITLE:
                    this.setTitle(this.getTitle().concat(" " + line.substring(2).trim()));
                    break;
                    
                case REFERENCE_XF:                    
                    try{
                        String pubmedStr = line.substring(line.trim().indexOf(PUBMED_MARKER) + PUBMED_MARKER.length() + 1, line.length()-1);
                        this.setPubmedID(Integer.parseInt(line.substring(line.trim().indexOf(PUBMED_MARKER) + PUBMED_MARKER.length() + 1, line.length()-1).trim()));
                    }
                    catch(Exception exEx){
                        logger.info("couldn't parse pubmedID <" + line + "> to integer");
                        logger.error("couldn't parse pubmedID <" + line + "> to integer");
                        throw new NumberFormatException("\"couldn't parse pubmedID <\" + line + \"> to integer\"");
                    }
                    break;
                    
                default:
                    logger.info("couldn't parse line <" + line + "> \n  - unrecognized two letter code '" 
                            + line.trim().substring(1, 2) + ">" );
                    logger.error("couldn't parse line <" + line + "> \n  - unrecognized two letter code '" 
                            + line.trim().substring(1, 2) + ">" );
                    throw new Exception("couldn't parse line <" + line + "> \n  - unrecognized two letter code '" 
                            + line.trim().substring(1, 2) + ">");
                    
                    
            }
        }
        this.setTitle(this.getTitle().substring(2, this.getTitle().length()-2));
        
        
    }
    
    
    
    /**
     * @return the authors
     */
    public String getAuthors() {
        return authors;
    }

    /**
     * @param authors the authors to set
     */
    public void setAuthors(String authors) {
        this.authors = authors;
    }

    /**
     * @return the title
     */
    public String getTitle() {
        return title;
    }

    /**
     * @param title the title to set
     */
    public void setTitle(String title) {
        this.title = title;
    }

    /**
     * @return the journal
     */
    public String getJournalRef() {
        return journalRef;
    }

    /**
     * @param journal the journal to set
     */
    public void setJournalRef(String journal) {
        this.journalRef = journal;
    }

    /**
     * @return the pubmedID
     */
    public int getPubmedID() {
        return pubmedID;
    }

    /**
     * @param pubmedID the pubmedID to set
     */
    public void setPubmedID(int pubmedID) {
        this.pubmedID = pubmedID;
    }

    /**
     * @return the refNum
     */
    public String getRefNum() {
        return refNum;
    }

    /**
     * @param refNum the refNum to set
     */
    public void setRefNum(String refNum) {
        this.refNum = refNum;
    }
}
