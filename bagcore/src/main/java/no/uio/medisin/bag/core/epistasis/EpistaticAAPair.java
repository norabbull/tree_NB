/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.epistasis;

/**
 *
 * @author sr
 */
public class EpistaticAAPair {
  private String aa1;
  private String aa2;
  private int count;


  public EpistaticAAPair(String acid1, String acid2)
  {
    aa1 = acid1;
    aa2 = acid2;
    count = 1;
  }

  /**
   * @return the aa1
   */
  public String getAa1() {
    return aa1;
  }

  /**
   * @param aa1 the aa1 to set
   */
  public void setAa1(String aa1) {
    this.aa1 = aa1;
  }

  /**
   * @return the aa2
   */
  public String getAa2() {
    return aa2;
  }

  /**
   * @param aa2 the aa2 to set
   */
  public void setAa2(String aa2) {
    this.aa2 = aa2;
  }

  /**
   * @return the count
   */
  public int getCount() {
    return count;
  }

  /**
   * @param count the count to set
   */
  public void setCount(int count) {
    this.count = count;
  }



}
