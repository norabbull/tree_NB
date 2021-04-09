/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.epistasis;

import no.uio.medisin.bag.core.tree.MutationPair;


/**
 *
 * @author sr
 */
public class EpistasisSet {
  private java.util.ArrayList<MutationPair> epistasisSet;

  public EpistasisSet()
  {
    epistasisSet = new java.util.ArrayList<MutationPair>();
  }

  /**
   * @return the epistasisSet
   */
  public java.util.ArrayList<MutationPair> getEpistasisSet() {
    return epistasisSet;
  }


}
