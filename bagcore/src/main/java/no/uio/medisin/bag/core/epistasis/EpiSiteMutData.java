/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.epistasis;

import no.uio.medisin.bag.core.tree.NodeType;
import no.uio.medisin.bag.core.tree.SiteMutData;


/**
 *
 * @author sr
 */
public class EpiSiteMutData {
  private String seqName;
  private SiteMutData seqSiteAMutData;
  private SiteMutData seqSiteBMutData;
  private NodeType nodeLocation;
  private double alignDist;

  public EpiSiteMutData(SiteMutData epiA, SiteMutData epiB)
  {
    seqSiteAMutData = epiA;
    seqSiteBMutData = epiB;
  }

  public EpiSiteMutData()
  {
    seqSiteAMutData = new SiteMutData();
    seqSiteBMutData = new SiteMutData();
  }


  public String printMutationSummary()
  {
    String s = seqName + "\t" + alignDist + "\t"
            + seqSiteAMutData.printSummary() + "\t"
            + seqSiteBMutData.printSummary();

    return s;
  }


  /**
   * returns true if the sequence is epistatic
   *
   * @param s
   * @return
   */
  public Boolean isSeqEpistatic()
  {
     if(seqSiteAMutData.getNt() == seqSiteAMutData.getConsNT() &&
      seqSiteBMutData.getNt() == seqSiteBMutData.getConsNT())
       return false;    // same as the consensus sequence
     else
       if(seqSiteAMutData.getNt() != seqSiteAMutData.getConsNT() &&
        seqSiteBMutData.getNt() == seqSiteBMutData.getConsNT())
         return false;  // only first site has mutated
       else
       if(seqSiteAMutData.getNt() == seqSiteAMutData.getConsNT() &&
        seqSiteBMutData.getNt() != seqSiteBMutData.getConsNT())
         return false;  // only second site has mutated
       else
         return true;   // both sites have mutated
  }

  /**
   * @return the seqAEpiData
   */
  public SiteMutData getSeqAEpiData() {
    return seqSiteAMutData;
  }

  /**
   * @return the seqBEpiData
   */
  public SiteMutData getSeqBEpiData() {
    return seqSiteBMutData;
  }


  /**
   * @return the alignDist
   */
  public double getAlignDist() {
    return alignDist;
  }

  /**
   * @param alignDist the alignDist to set
   */
  public void setAlignDist(double alignDist) {
    this.alignDist = alignDist;
  }

  /**
   * @param seqAEpiData the seqAEpiData to set
   */
  public void setSeqAEpiData(SiteMutData seqAEpiData) {
    this.seqSiteAMutData = seqAEpiData;
  }

  /**
   * @param seqBEpiData the seqBEpiData to set
   */
  public void setSeqBEpiData(SiteMutData seqBEpiData) {
    this.seqSiteBMutData = seqBEpiData;
  }

  /**
   * @return the seqName
   */
  public String getSeqName() {
    return seqName;
  }

  /**
   * @param seqName the seqName to set
   */
  public void setSeqName(String seqName) {
    this.seqName = seqName;
  }

  /**
   * @return the nodeLocation
   */
  public NodeType getNodeLocation() {
    return nodeLocation;
  }

  /**
   * @param nodeLocation the nodeLocation to set
   */
  public void setNodeLocation(NodeType nodeLocation) {
    this.nodeLocation = nodeLocation;
  }

}
