/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

import java.util.HashMap;
import java.util.Map;

/**
 * adds additional functionality for selecting and plotting trees.
 * This is no longer used, as the TreeNode now includes a HashMap for storing
 * parameters related to plotting
 * @author sr
 */
@Deprecated 
public class SuperTreeNode extends TreeNode{


  private int lineWidth;
  private int lineColour;
  private int fontColour;
  private String fontName;
  private int fontStyle;

  private HashMap<String, String>  otherParams  = new HashMap<>();;


  public SuperTreeNode(int n, TreeType tType)
  {
    super(n, tType);

  }

  public SuperTreeNode(TreeNode node, TreeType tType)
  {
    super(node.getNodeID(), TreeType.NEXUSPLUSPLUS);

    // need to copy everything over
   setParentNode(node.getParentNode());
   setName(node.getName());
   setBranchLength(node.getBranchLength());
   this.setBootStrap(node.getBootStrap());
   this.setX(node.getX());
   this.setY(node.getY());
   this.setxBounds(node.getxBounds());
   this.setyBounds(node.getyBounds());
   this.setLeftChild(node.getLeftChild());
   this.setRightChild(node.getRightChild());
   for(int child: node.getChildren())
   {
     this.getChildren().add(child);
   }
   java.util.Iterator itH = this.getOtherParams().entrySet().iterator();
   while(itH.hasNext())
   {
     Map.Entry p = (Map.Entry) itH.next();
     this.getOtherParams().put(p.getKey(), p.getValue());
   }
  }

  /**
   * @return the lineColour
   */
  public int getLineColour() {
    return lineColour;
  }

  /**
   * @param lineColour the lineColour to set
   */
  public void setLineColour(int lineColour) {
    this.lineColour = lineColour;
  }

  /**
   * @return the fontColour
   */
  public int getFontColour() {
    return fontColour;
  }

  /**
   * @param fontColour the fontColour to set
   */
  public void setFontColour(int fontColour) {
    this.fontColour = fontColour;
  }

  /**
   * @return the fontName
   */
  public String getFontName() {
    return fontName;
  }

  /**
   * @param fontName the fontName to set
   */
  public void setFontName(String fontName) {
    this.fontName = fontName;
  }

  /**
   * @return the fontStyle
   */
  public int getFontStyle() {
    return fontStyle;
  }

  /**
   * @param fontStyle the fontStyle to set
   */
  public void setFontStyle(int fontStyle) {
    this.fontStyle = fontStyle;
  }

  /**
   * @return the lineWidth
   */
  public int getLineWidth() {
    return lineWidth;
  }

  /**
   * @param lineWidth the lineWidth to set
   */
  public void setLineWidth(int lineWidth) {
    this.lineWidth = lineWidth;
  }

  /**
   * @return the otherParams
   */
  public HashMap getOtherParams() {
    return otherParams;
  }

}
