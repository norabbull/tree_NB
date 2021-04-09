/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

/**
 *
 * @author sr
 */
public class BranchMutualPair {
  private double branchLength;
  private double mutualInfo;

  public BranchMutualPair(double bl, double mi)
  {
    branchLength = bl;
    mutualInfo = mi;
  }

  /**
   * @return the branchLength
   */
  public double getBranchLength() {
    return branchLength;
  }

  /**
   * @param branchLength the branchLength to set
   */
  public void setBranchLength(double branchLength) {
    this.branchLength = branchLength;
  }

  /**
   * @return the mutualInfo
   */
  public double getMutualInfo() {
    return mutualInfo;
  }

  /**
   * @param mutualInfo the mutualInfo to set
   */
  public void setMutualInfo(double mutualInfo) {
    this.mutualInfo = mutualInfo;
  }

}
