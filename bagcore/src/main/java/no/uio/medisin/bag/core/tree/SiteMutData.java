/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

/**
 *
 * @author sr
 */
public class SiteMutData {
  private int ntPos;
  private int nt;
  private int codonPos;
  private int codonStart;
  private int consNT;
  private String codon;
  private String aAcid;
  private String consCodon;
  private String consAAcid;
  private double conservation;
  private SubstitutionType subType;


  public String printSummary()
  {
    String s = ntPos + "\t" + nt + "\t" + codonPos + "\t"
            + aAcid + "\t" + consAAcid + "\t" + codon + "\t" + consCodon;
    if(subType == no.uio.medisin.bag.core.tree.SubstitutionType.SYNONYMOUS)
      s = s.concat("S");
    else
      s = s.concat("N");
    return s;
  }

  /**
   * @return the ntPos
   */
  public int getNtPos() {
    return ntPos;
  }

  /**
   * @param ntPos the ntPos to set
   */
  public void setNtPos(int ntPos) {
    this.ntPos = ntPos;
  }

  /**
   * @return the nt
   */
  public int getNt() {
    return nt;
  }

  /**
   * @param nt the nt to set
   */
  public void setNt(int nt) {
    this.nt = nt;
  }

  /**
   * @return the codonPos
   */
  public int getCodonPos() {
    return codonPos;
  }

  /**
   * @param codonPos the codonPos to set
   */
  public void setCodonPos(int codonPos) {
    this.codonPos = codonPos;
  }

  /**
   * @return the conservation
   */
  public double getConservation() {
    return conservation;
  }

  /**
   * @param conservation the conservation to set
   */
  public void setConservation(double conservation) {
    this.conservation = conservation;
  }

  /**
   * @return the subType
   */
  public SubstitutionType getSubType() {
    return subType;
  }

  /**
   * @param subType the subType to set
   */
  public void setSubType(SubstitutionType subType) {
    this.subType = subType;
  }

  /**
   * @return the codon
   */
  public String getCodon() {
    return codon;
  }

  /**
   * @param codon the codon to set
   */
  public void setCodon(String codon) {
    this.codon = codon;
  }

  /**
   * @return the aAcid
   */
  public String getaAcid() {
    return aAcid;
  }

  /**
   * @param aAcid the aAcid to set
   */
  public void setaAcid(String aAcid) {
    this.aAcid = aAcid;
  }

  /**
   * @return the consCodon
   */
  public String getConsCodon() {
    return consCodon;
  }

  /**
   * @param consCodon the consCodon to set
   */
  public void setConsCodon(String consCodon) {
    this.consCodon = consCodon;
  }

  /**
   * @return the consAAcid
   */
  public String getConsAAcid() {
    return consAAcid;
  }

  /**
   * @param consAAcid the consAAcid to set
   */
  public void setConsAAcid(String consAAcid) {
    this.consAAcid = consAAcid;
  }

  /**
   * @return the codonStart
   */
  public int getCodonStart() {
    return codonStart;
  }

  /**
   * @param codonStart the codonStart to set
   */
  public void setCodonStart(int codonStart) {
    this.codonStart = codonStart;
  }

  /**
   * @return the consNT
   */
  public int getConsNT() {
    return consNT;
  }

  /**
   * @param consNT the consNT to set
   */
  public void setConsNT(int consNT) {
    this.consNT = consNT;
  }
}
