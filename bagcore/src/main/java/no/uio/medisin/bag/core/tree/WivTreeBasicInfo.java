/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

import no.uio.medisin.bag.core.tree.TreeType;

/**
 *
 * @author sr
 */
public class WivTreeBasicInfo {
  private TreeType treeType;
  private String treeFilename;
  private String treeLine;
  private double treeLength;
  static final String eol = System.getProperty("line.separator");

  /**
   * @return the treeType
   */
  public TreeType getTreeType() {
    return treeType;
  }

  /**
   * @param treeType the treeType to set
   */
  public void setTreeType(TreeType treeType) {
    this.treeType = treeType;
  }

  /**
   * @return the treeFilename
   */
  public String getTreeFilename() {
    return treeFilename;
  }

  /**
   * @param treeFilename the treeFilename to set
   */
  public void setTreeFilename(String treeFilename) {
    this.treeFilename = treeFilename;
  }

  /**
   * @return the treeLine
   */
  public String getTreeLine() {
    return treeLine;
  }

  /**
   * @param treeLine the treeLine to set
   */
  public void setTreeLine(String treeLine) {
    this.treeLine = treeLine;
  }

  /**
   * @return the treeLength
   */
  public double getTreeLength() {
    return treeLength;
  }

  /**
   * @param treeLength the treeLength to set
   */
  public void setTreeLength(double treeLength) {
    this.treeLength = treeLength;
  }


}
