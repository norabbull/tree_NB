/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.tree;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * contains a group of related nodes. 
 * For example, a set of nodes that form a sub-tree, or a set of nodes that
 * share a particular property such as a shared geographical origin or host.
 * 
 * @author simonrayner
 */
public class TreeNodeGroup {
  private int groupID;
  private String groupName;
  private HashMap<String, String> otherParams = new HashMap<>();
  private ArrayList<Integer> leafNodesInGroup = new ArrayList<>();      // Added by Nora, nodeIDs

  
  public TreeNodeGroup(int id, String name){
    groupID = id;
    groupName = name;
  }

  /**
   * @return the groupID
   */
  public int getGroupID() {
    return groupID;
  }

  /**
   * @param groupID the groupID to set
   */
  public void setGroupID(int groupID) {
    this.groupID = groupID;
  }

  /**
   * @return the groupName
   */
  public String getGroupName() {
    return groupName;
  }

  /**
   * @param groupName the groupName to set
   */
  public void setGroupName(String groupName) {
    this.groupName = groupName;
  }

  /**
   * @return the otherParams
   */
  public HashMap<String, String> getOtherParams() {
    return otherParams;
  }
  
  public String getParameterValue(String key){
    return this.otherParams.get(key);
  }  
  
  public void setParameter(String key, String value){
    otherParams.put(key, value);
  }

  // Added by Nora
  /**
   *
   * @param ID ID of child node that is added to group
   */
  public void addLeafNodeToGroup(Integer ID) {
    this.leafNodesInGroup.add(ID);
  }

  // Added by Nora
  public ArrayList<Integer> getLeafNodesInGroup() { return leafNodesInGroup; }
  
}
