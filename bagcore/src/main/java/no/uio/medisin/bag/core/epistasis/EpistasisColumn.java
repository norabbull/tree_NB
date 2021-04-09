package no.uio.medisin.bag.core.epistasis;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author sr
 */
public class EpistasisColumn {
  private java.util.ArrayList<EpistasisSet> epistasisCol;

  public EpistasisColumn()
  {
    epistasisCol = new java.util.ArrayList<EpistasisSet>();
  }

  /**
   * @return the epistasisCol
   */
  public java.util.ArrayList<EpistasisSet> getEpistasisCol() {
    return epistasisCol;
  }

}
