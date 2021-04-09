/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.sequence;

import no.uio.medisin.bag.core.mirna.PreMiRNASet;

/**
 *
 * @author simonray
 */
public class GFF2SeqGrabber {
  private String                      organismCode;
  private String                      genomeFAFilepath;
  private String                      gffFilepath;
  private Boolean                     singleSeqFile = false;

  private PreMiRNASet                 preMiRSet;
    
    
  public GFF2SeqGrabber(){
    preMiRSet = new PreMiRNASet();
  }
  
  
  public void gff2seq(){
    preMiRSet.setOrganismCode(organismCode);
    preMiRSet.setGenomeFAFilepath(genomeFAFilepath);
    preMiRSet.setGffFilepath(gffFilepath);
    
    
  }
  
  
  
  /**
   * @return the genomeFAFilepath
   */
  public String getGenomeFAFilepath() {
    return genomeFAFilepath;
  }

  /**
   * @param genomeFAFilepath the genomeFAFilepath to set
   */
  public void setGenomeFAFilepath(String genomeFAFilepath) {
    this.genomeFAFilepath = genomeFAFilepath;
  }

  /**
   * @return the gffFilepath
   */
  public String getGffFilepath() {
    return gffFilepath;
  }

  /**
   * @param gffFilepath the gffFilepath to set
   */
  public void setGffFilepath(String gffFilepath) {
    this.gffFilepath = gffFilepath;
  }

  /**
   * @return the singleSeqFile
   */
  public Boolean getSingleSeqFile() {
    return singleSeqFile;
  }

  /**
   * @param singleSeqFile the singleSeqFile to set
   */
  public void setSingleSeqFile(Boolean singleSeqFile) {
    this.singleSeqFile = singleSeqFile;
  }

  /**
   * @return the organismCode
   */
  public String getOrganismCode() {
    return organismCode;
  }

  /**
   * @param organismCode the organismCode to set
   */
  public void setOrganismCode(String organismCode) {
    this.organismCode = organismCode;
  }

  /**
   * @return the preMiRSet
   */
  public PreMiRNASet getPreMiRSet() {
    return preMiRSet;
  }

  /**
   * @param preMiRSet the preMiRSet to set
   */
  public void setPreMiRSet(PreMiRNASet preMiRSet) {
    this.preMiRSet = preMiRSet;
  }
}
