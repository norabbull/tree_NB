/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;

/**
 *
 * @author sr
 */
public class entropy implements java.io.Serializable {

  double h;

  public void setEntropy(double e)
  {
    h = e;
  }

  public double getEntropy()
  {
    return h;
  }

  public void incEntropy(double i)
  {
    h += i;
  }
}
