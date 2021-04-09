/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;

/**
 *
 * @author sr
 */
public class AAProperties {

  private String codon;
  private String shortName;
  private String threeLetterName;
  private String longName;

  public AAProperties(String c, String m, String s, String l)
  {
    codon = c;
    shortName = s;
    threeLetterName = m;
    longName = l;
  }
  /**
   * @return the codon
   */
  public String getCodon() {
    return codon;
  }

  /**
   * @return the shortName
   */
  public String getShortName() {
    return shortName;
  }

  /**
   * @return the threeLetterName
   */
  public String getThreeLetterName() {
    return threeLetterName;
  }

  /**
   * @return the longName
   */
  public String getLongName() {
    return longName;
  }


}
