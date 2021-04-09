/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.epistasis;

import no.uio.medisin.bag.core.aminoacid.AminoPairList;

/**
 *
 * @author sr
 */
public class EpiMutationStats {

  private AminoPairList aaPairs;

  private int siteAcodon1count;
  private int siteAcodon2count;
  private int siteAcodon3count;

  private int nonSynACount;
  private int synACount;

  private int siteBcodon1count;
  private int siteBcodon2count;
  private int siteBcodon3count;

  private int nonSynBCount;
  private int synBCount;


  public EpiMutationStats()
  {
    aaPairs = new AminoPairList();
  }

  /**
   * add new amino acid pair to the pair list
   * 
   * @param newAApair
   * @return
   */
  public int addPair(EpistaticAAPair newAApair)
  {

    aaPairs.add(newAApair);
    return ((EpistaticAAPair)
            aaPairs.getAApair(aaPairs.getAaPairs().size()-1)).getCount();
  }

  /**
   * increment the SiteA count for Codon1 mutations
   *
   * @return
   */
  public int incSiteACodon1Count()
  {
    siteAcodon1count++;
    return siteAcodon1count;
  }

  /**
   * increment the SiteA count for Codon2 mutations
   *
   * @return
   */
  public int incSiteACodon2Count()
  {
    siteAcodon2count++;
    return siteAcodon2count;
  }


  /**
   * increment the SiteA count for Codon3 mutations
   *
   * @return
   */
  public int incSiteACodon3Count()
  {
    siteAcodon3count++;
    return siteAcodon3count;
  }


  /**
   * increment the SiteB count for Codon1 mutations
   *
   * @return
   */
  public int incSiteBCodon1Count()
  {
    siteBcodon1count++;
    return siteBcodon1count;
  }


  /**
   * increment the SiteB count for Codon2 mutations
   *
   * @return
   */
  public int incSiteBCodon2Count()
  {
    siteBcodon2count++;
    return siteBcodon2count;
  }


  /**
   * increment the SiteB count for Codon3 mutations
   *
   * @return
   */
  public int incSiteBCodon3Count()
  {
    siteBcodon3count++;
    return siteBcodon3count;
  }


  /**
   * increment the Site A count for non synonymous substitutions
   *
   * @return
   */
  public int incSiteANonSynCount()
  {
    this.nonSynACount++;
    return nonSynACount;
  }

  /**
   * increment the Site A count for synonymous substitutions
   *
   * @return
   */
  public int incSiteASynCount()
  {
    this.synACount++;
    return synACount;
  }

  /**
   * increment the Site B count for non synonymous substitutions
   *
   * @return
   */
  public int incSiteBNonSynCount()
  {
    this.nonSynBCount++;
    return nonSynBCount;
  }

  /**
   * increment the Site B count for synonymous substitutions
   *
   * @return
   */
  public int incSiteBSynCount()
  {
    this.synBCount++;
    return synBCount;
  }


  /**
   * @return the siteAcodon1count
   */
  public int getSiteAcodon1count() {
    return siteAcodon1count;
  }

  /**
   * @param siteAcodon1count the siteAcodon1count to set
   */
  public void setSiteAcodon1count(int siteAcodon1count) {
    this.siteAcodon1count = siteAcodon1count;
  }

  /**
   * @return the siteAcodon2count
   */
  public int getSiteAcodon2count() {
    return siteAcodon2count;
  }

  /**
   * @param siteAcodon2count the siteAcodon2count to set
   */
  public void setSiteAcodon2count(int siteAcodon2count) {
    this.siteAcodon2count = siteAcodon2count;
  }

  /**
   * @return the siteAcodon3count
   */
  public int getSiteAcodon3count() {
    return siteAcodon3count;
  }

  /**
   * @param siteAcodon3count the siteAcodon3count to set
   */
  public void setSiteAcodon3count(int siteAcodon3count) {
    this.siteAcodon3count = siteAcodon3count;
  }

  /**
   * @return the siteBcodon1count
   */
  public int getSiteBcodon1count() {
    return siteBcodon1count;
  }

  /**
   * @param siteBcodon1count the siteBcodon1count to set
   */
  public void setSiteBcodon1count(int siteBcodon1count) {
    this.siteBcodon1count = siteBcodon1count;
  }

  /**
   * @return the siteBcodon2count
   */
  public int getSiteBcodon2count() {
    return siteBcodon2count;
  }

  /**
   * @param siteBcodon2count the siteBcodon2count to set
   */
  public void setSiteBcodon2count(int siteBcodon2count) {
    this.siteBcodon2count = siteBcodon2count;
  }

  /**
   * @return the siteBcodon3count
   */
  public int getSiteBcodon3count() {
    return siteBcodon3count;
  }

  /**
   * @param siteBcodon3count the siteBcodon3count to set
   */
  public void setSiteBcodon3count(int siteBcodon3count) {
    this.siteBcodon3count = siteBcodon3count;
  }

  /**
   * @return the aaPairs
   */
  public AminoPairList getAaPairs() {
    return aaPairs;
  }

  /**
   * @return the nonSynACount
   */
  public int getNonSynACount() {
    return nonSynACount;
  }

  /**
   * @param nonSynACount the nonSynACount to set
   */
  public void setNonSynACount(int nonSynACount) {
    this.nonSynACount = nonSynACount;
  }

  /**
   * @return the synACount
   */
  public int getSynACount() {
    return synACount;
  }

  /**
   * @param synACount the synACount to set
   */
  public void setSynACount(int synACount) {
    this.synACount = synACount;
  }

  /**
   * @return the nonSynBCount
   */
  public int getNonSynBCount() {
    return nonSynBCount;
  }

  /**
   * @param nonSynBCount the nonSynBCount to set
   */
  public void setNonSynBCount(int nonSynBCount) {
    this.nonSynBCount = nonSynBCount;
  }

  /**
   * @return the synBCount
   */
  public int getSynBCount() {
    return synBCount;
  }

  /**
   * @param synBCount the synBCount to set
   */
  public void setSynBCount(int synBCount) {
    this.synBCount = synBCount;
  }
}
