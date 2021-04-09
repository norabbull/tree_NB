/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

/**
 *
 * @author sr
 */
public class MutationPair {
  private int site1;
  private int site2;
  private String base1;
  private String base2;

  public MutationPair(int s1, int s2, String b1, String b2)
  {
    site1 = s1;
    site2= s2;
    base1 = b1;
    base2 = b2;
  }
  /**
   * @return the site1
   */
  public int getSite1() {
    return site1;
  }

  /**
   * @param site1 the site1 to set
   */
  public void setSite1(int site1) {
    this.site1 = site1;
  }

  /**
   * @return the site2
   */
  public int getSite2() {
    return site2;
  }

  /**
   * @param site2 the site2 to set
   */
  public void setSite2(int site2) {
    this.site2 = site2;
  }

  /**
   * @return the base1
   */
  public String getBase1() {
    return base1;
  }

  /**
   * @param base1 the base1 to set
   */
  public void setBase1(String base1) {
    this.base1 = base1;
  }

  /**
   * @return the base2
   */
  public String getBase2() {
    return base2;
  }

  /**
   * @param base2 the base2 to set
   */
  public void setBase2(String base2) {
    this.base2 = base2;
  }

}
