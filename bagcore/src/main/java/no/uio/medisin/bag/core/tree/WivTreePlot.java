/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

import java.awt.Color;
import no.uio.medisin.bag.core.tree.SuperTreeNode;
import no.uio.medisin.bag.core.tree.TreeNode;

/**
 * This is an extended version of the tree, it contains additional information
 * about how to plot the tree.  It's the same basic structure as WivTree, just the nodes
 * are different.
 * 
 * This is no longer used as all the plotting information is stored in a HashMap
 * in the WivTree class.
 *
 * @author sr
 */
@Deprecated 
public class WivTreePlot{

  private WivTreeBasicInfo treeBasicInfo;
  private java.util.ArrayList<SuperTreeNode> sNodes;

  
  private int defaultLineWidth = 1;
  private int selectLineWidth = 4;
  
  private Color defaultLineColour = Color.blue;
  private Color highlightLineColour = Color.red;
  

  public WivTreePlot()
  {
    treeBasicInfo  = new WivTreeBasicInfo();
    sNodes = new java.util.ArrayList<SuperTreeNode>();
  }

  public WivTreePlot(WivTree wTree)
  {
    treeBasicInfo  = new WivTreeBasicInfo();
    treeBasicInfo.setTreeFilename(wTree.getTreeFilename());
    treeBasicInfo.setTreeLength(wTree.getTreeLength());
    treeBasicInfo.setTreeType(wTree.getTreeType());

    // copy over the nodes
    java.util.ListIterator itN = wTree.getNodes().listIterator();
    while(itN.hasNext())
    {
      TreeNode n = (TreeNode)itN.next();
      sNodes.add(new SuperTreeNode(n.getNodeID(), treeBasicInfo.getTreeType()));
    }
  }


  /**
   * copy over the nodes from a standard WivTree
   * @param wTree
   */
  public void copyTree(WivTree wTree)
  {
    treeBasicInfo  = new WivTreeBasicInfo();
    treeBasicInfo.setTreeFilename(wTree.getTreeFilename());
    treeBasicInfo.setTreeLength(wTree.getTreeLength());
    treeBasicInfo.setTreeType(wTree.getTreeType());

    // copy over the nodes
    java.util.ListIterator itN = wTree.getNodes().listIterator();
    while(itN.hasNext())
    {
      TreeNode n = (TreeNode)itN.next();
      sNodes.add(new SuperTreeNode(n, treeBasicInfo.getTreeType()));
    }
  }

  /**
   * return the parentID for this node
   *
   * @param thisNode
   * @return
   */
  public int getParentID(TreeNode thisNode)
  {
    return thisNode.getParentNode();
  }


  /**
   * return the node at index int
   * @param ID
   * @return
   */
  public SuperTreeNode getNode(int ID)
  {
    return this.sNodes.get(ID);
  }


  /**
   * @return the treeBasicInfo
   */
  public WivTreeBasicInfo getTreeBasicInfo() {
    return treeBasicInfo;
  }

  /**
   * @return the sNodes
   */
  public java.util.ArrayList<SuperTreeNode> getNodes() {
    return sNodes;
  }

  /**
   * @return the defaultLineWidth
   */
  public int getDefaultLineWidth() {
    return defaultLineWidth;
  }

  /**
   * @param defaultLineWidth the defaultLineWidth to set
   */
  public void setDefaultLineWidth(int defaultLineWidth) {
    this.defaultLineWidth = defaultLineWidth;
  }

  /**
   * @return the selectLineWidth
   */
  public int getSelectLineWidth() {
    return selectLineWidth;
  }

  /**
   * @param selectLineWidth the selectLineWidth to set
   */
  public void setSelectLineWidth(int selectLineWidth) {
    this.selectLineWidth = selectLineWidth;
  }

  /**
   * @return the defaultLineColour
   */
  public Color getDefaultLineColour() {
    return defaultLineColour;
  }

  /**
   * @param defaultLineColour the defaultLineColour to set
   */
  public void setDefaultLineColour(Color defaultLineColour) {
    this.defaultLineColour = defaultLineColour;
  }

  /**
   * @return the highlightLineColour
   */
  public Color getHighlightLineColour() {
    return highlightLineColour;
  }

  /**
   * @param highlightLineColour the highlightLineColour to set
   */
  public void setHighlightLineColour(Color highlightLineColour) {
    this.highlightLineColour = highlightLineColour;
  }

  
}
