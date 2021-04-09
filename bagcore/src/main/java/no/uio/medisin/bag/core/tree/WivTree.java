/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

import java.awt.Color;
import java.io.*;
import java.lang.reflect.Array;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * a tree parsing class
 * 
 * The class is supposed to be able to handle all forms of the Newick format
 * although this seems to a little varied
 * the details are provided here (http://evolution.genetics.washington.edu/phylip/newick_doc.html)
 * and BEAST trees appear to follow this format
 * 
 * attention should be paid to the following
 * 
 * Notes:
 *  Unquoted labels may not contain blanks, parentheses, square brackets,
 *       single_quotes, colons, semicolons, or commas.
 *  Underscore characters in unquoted labels are converted to blanks.
 *  Single quote characters in a quoted label are represented by two single
 *       quotes.
 *  Blanks or tabs may appear anywhere except within unquoted labels or
 *       branch_lengths.
 *  Newlines may appear anywhere except within labels or branch_lengths.
 *  Comments are enclosed in square brackets and may appear anywhere
 * 
 * I could use regex to search for matching parentheses, but it seems that this can be
 * unreliable, so i perform my own (slower) searching as my method is more trustworthy
 * 
 * 
 * @author sr
 */
public class WivTree {


  public static final String DEFAULT_LINECOLOUR_KEYSTRING = "default_line_colour";
  public static final String DEFAULT_LINEWIDTH_KEYSTRING = "default_line_width";

  static Logger logger = LogManager.getRootLogger();

  private TreeType treeType;
  private String treeFilename;
  private String sampleListFilename;
  private java.util.ArrayList<TreeNode> nodes;
  private java.util.ArrayList<String> taxaNames;
  private java.util.ArrayList<String> translateNames;
  private java.util.ArrayList<TreeNodeGroup> treeGroups;
  private SampleClassificationList sampleClassificationList = new SampleClassificationList();

  private String treeLine;
  private double treeLength;

  private HashMap<String, String> pigTreeParams = new HashMap<>(); // store PigTree parameters related to display properties
  private HashMap<String, String> figTreeParams = new HashMap<>(); // store figTree parameters related to display properties
  // Noras adding OBS! ELIMINATE THESE
  public double SDV;              // Any point in making these public or not?


  public WivTree() {
    nodes = new java.util.ArrayList<TreeNode>();
    taxaNames = new java.util.ArrayList<String>();
    translateNames = new java.util.ArrayList<String>();
    treeGroups = new ArrayList<TreeNodeGroup>();

    this.setPigtreeParameter(WivTree.DEFAULT_LINECOLOUR_KEYSTRING, Color.BLACK.toString());
    this.setPigtreeParameter(WivTree.DEFAULT_LINEWIDTH_KEYSTRING, Integer.toString(2));
  }

  /**
   * read and parse a tree file
   *
   * @param
   * @throws java.io.IOException
   */
  public void readTree(String fileName) throws java.io.IOException {
    /*
     * check tree type by reading the first line
     */
    setTreeLine("");
    BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
    String line = br.readLine();
    br.close();
    if (line.equals("#NEXUS")) {
      this.readNexusTree(fileName);
    } else {
      readNewickTree(fileName);
    }


  }

  /**
   * Nexus Trees contain a series of 'Blocks' defined by 'begin' and 'end' lines
   * <p>
   * I'm not sure of the complete set of allowed Blocks. For now, we recognize:
   * <p>
   * begin taxa;     end;
   * begin trees;    end;
   * begin figtree;  end;  - display information for FigTree
   * begin pigtree;  end;  - display information for PigTree
   *
   * @param nexusFileName
   * @throws java.io.IOException
   */
  public void readNexusTree(String nexusFileName) throws java.io.IOException {

    String line = "";

    BufferedReader brNx = new BufferedReader(new FileReader(new File(nexusFileName)));
    while ((line = brNx.readLine()) != null) {
      if (line.toLowerCase().contains("begin taxa;")) {
        brNx = this.readTaxaBlock(brNx);
      }

      if (line.toLowerCase().contains("begin trees;")) {
        brNx = this.readTreeBlock(brNx);
      }

      if (line.toLowerCase().contains("begin figtree;")) {
        brNx = this.readFigtreeBlock(brNx);
      }

      if (line.toLowerCase().contains("begin pigtree;")) {
        brNx = this.readPigtreeBlock(brNx);
      }

    }

    brNx.close();

  }

  /**
   * parse the 'taxa' block. This is delimited by 'begin' / 'end' lines
   *
   * @param brNx
   * @return
   * @throws java.io.IOException
   */
  private BufferedReader readTaxaBlock(BufferedReader brNx) throws java.io.IOException {

    String taxaLine = "";
    int numTaxa = 0;
    while ((taxaLine = brNx.readLine()) != null) {
      if (taxaLine.toLowerCase().contains("taxlabels")) {
        continue;
      }
      if (taxaLine.trim().contentEquals(";")) {
        continue;
      }
      if (taxaLine.toLowerCase().contains("dimensions")) {
        taxaLine = taxaLine.split(";")[0].trim(); // may be a ";" at the end of the line
        numTaxa = Integer.valueOf(taxaLine.split("=")[1].trim());
        continue;
      }
      if (taxaLine.toLowerCase().contains("end;")) {
        if (numTaxa != taxaNames.size()) {
          logger.warn("--parsing Taxlabels in Taxa block");
          logger.warn("--expected <" + numTaxa + "> Taxlabels, but found <" + taxaNames.size() + ">");
          int t = 1;
          for (String label : taxaNames) {
            logger.warn("----" + "[" + t + "]" + label);
            t++;
          }
        }
        return brNx;
      }
      taxaNames.add(taxaLine.trim());
    }

    throw new IOException("parsing error reading the taxa block. Didn't find the end of block");
  }


  /**
   * parse the 'tree' block. This is delimited by 'begin' / 'end' lines
   *
   * @param brNx
   * @return
   * @throws java.io.IOException
   */
  private BufferedReader readTreeBlock(BufferedReader brNx) throws java.io.IOException {

    String nxTreeLine = "";
    int treeCount = 0;
    while ((nxTreeLine = brNx.readLine()) != null) {
      if (nxTreeLine.toLowerCase().contains("end;")) {
        /*
         * if this is a BEAST tree, there will be extra characters at the front of the
         * tree string that need to be removed.
         * I'm still not entirely sure what these extra characters can include, but
         * can be something like
         *  tree TREE1 = [&R]
         *  where [&R] may be marking the Root node and [&U] indicates unrooted
         */
        String nxTreeName = this.treeLine.substring(0, this.treeLine.indexOf("="));
        nxTreeName = nxTreeName.split("tree ")[1].trim();
        this.treeLine = this.treeLine.substring(this.treeLine.indexOf("("));
        //this.setTreeLine(nxTreeLine);
        if (treeCount > 1) {
          logger.warn("there was more than one tree found in the trees block.\n"
                  + "PigTree doesn't currently handle multiple trees.\n"
                  + "Only the last tree in the block <" + nxTreeName + "> was saved");
        }
        return brNx;
      } else {
        if (nxTreeLine.contains("tree")) {
          this.setTreeLine(nxTreeLine);
          treeCount++;
        } else {
          if (nxTreeLine.toLowerCase().contains("translate")) {
            this.parseTranslateBlock(brNx);
          }
        }
      }
    }

    throw new IOException("parsing error reading the trees block. Didn't find the end of block");
  }


  /**
   * parse out TaxLabels in taxa block within Nexus tree
   *
   * @param brNx
   * @return
   * @throws java.io.IOException
   */
  private BufferedReader parseTranslateBlock(BufferedReader brNx) throws java.io.IOException {
    String nxTreeLine = "";

    while ((nxTreeLine = brNx.readLine()) != null) {
      if (nxTreeLine.trim().contains(";")) {
        logger.info("--found <" + translateNames.size() + "> translateLines");
        return brNx;
      } else {
        translateNames.add(nxTreeLine.split(",")[0].trim());
      }
    }
    return brNx;
  }


  /**
   * parse the 'figtree' block. This is delimited by 'begin' / 'end' lines
   *
   * @param brNx
   * @return
   * @throws java.io.IOException
   */
  private BufferedReader readFigtreeBlock(BufferedReader brNx) throws java.io.IOException {

    logger.warn("Found a FigTree block.\n"
            + "PigTree doesn't currently process this information.\n");
    String nxTreeLine = "";
    int lineCount = 1;
    while ((nxTreeLine = brNx.readLine()) != null) {
      if (nxTreeLine.contains("nd;")) {
        logger.info("found <" + (lineCount - 1) + "> lines in this block");
        return brNx;
      }
      lineCount++;
    }

    throw new IOException("parsing error reading the FigTree block. Didn't find the end of block");
  }


  /**
   * parse the 'figtree' block. This is delimited by 'begin' / 'end' lines
   *
   * @param brNx
   * @return
   * @throws java.io.IOException
   */
  private BufferedReader readPigtreeBlock(BufferedReader brNx) throws java.io.IOException {

    logger.warn("Found a FigTree block.\n"
            + "PigTree doesn't currently process this information.\n");
    String nxTreeLine = "";
    int lineCount = 1;
    while ((nxTreeLine = brNx.readLine()) != null) {
      if (nxTreeLine.contains("nd;")) {
        logger.info("found <" + (lineCount - 1) + "> lines in this block");
        return brNx;
      }
      lineCount++;
    }

    throw new IOException("parsing error reading the FigTree block. Didn't find the end of block");
  }


  /**
   * (I think) Newick trees are a single line containing the tree
   *
   * @param newickFileName
   * @throws java.io.IOException
   */
  public void readNewickTree(String newickFileName) throws java.io.IOException {
    BufferedReader brNw = new BufferedReader(new FileReader(new File(newickFileName)));
    this.treeLine = brNw.readLine();
    brNw.close();

  }


  /**
   * iterate through all class types (e.g. Super and Sub populations)
   * and class instances (e.g. AFR, EUR... and YOR, GIH ....)
   * and calculate the pairwise distance between each
   */
  public void getPairDistancesForClassificationTypes(String filename) {
    ArrayList<String> classificationTypes = this.getSampleClassificationList().getClassificationTypes();
    for (String classificationType : classificationTypes) {

    }
  }


  /**
   * iterate through all class types (e.g. Super and Sub populations)
   * and class instances (e.g. AFR, EUR... and YOR, GIH ....)
   * and calculate the pairwise distance between each
   */
  public void getPairDistancesForClassificationNames(String filename) throws IOException {
    // open file for writing
    BufferedWriter bwDP = new BufferedWriter(new FileWriter(new File(filename)));
    bwDP.write("#" + this.calcTreeLength() + System.lineSeparator());
    ArrayList<String> classificationTypes = this.getSampleClassificationList().getClassificationTypes();
    for (String classificationType : classificationTypes) {
      ArrayList<String> classificationNames = this.getSampleClassificationList().getAllEntriesForClassificationType(classificationType);

      for (String classificationName : classificationNames) {
        ArrayList<Integer> nodeIDSet = this.getMatchingLeafNodes(classificationName);
        ArrayList<Double> pairDistances = this.getPairDistancesForAllPairsInSet(nodeIDSet);

        for (Double pairDist : pairDistances) {
          bwDP.write(pairDist + "\t" + classificationType + "\t" + classificationName + System.lineSeparator());
        }

      }
    }
    bwDP.close();
  }

  /**
   * return the subset of leaf nodes that contain the sampleID
   * (this is not an ideal implementation as it searches for the presence of the
   * specified sanpleID in the sample name. So there may be multiple hits)
   *
   * @param sampleID query ID
   * @return list of node IDs that matched the query in their sampleName
   */
  public ArrayList<Integer> getMatchingLeafNodes(String sampleID) {
    ArrayList<Integer> nodeIDset = new ArrayList<>();
    for (TreeNode node : this.getNodes()) {
      if (node.getChildren().size() == 0 && node.getName().contains(sampleID)) {
        nodeIDset.add(node.getNodeID());
      }
    }
    return nodeIDset;
  }


  /**
   * calculate the distance between each pairwise combination of the nodes
   *
   * @param sampleIDs the list of leaf nodes to be processed
   * @return
   */
  public ArrayList<Double> getPairDistancesForAllPairsInSet(ArrayList<Integer> sampleIDs) {
    ArrayList<Double> pairDistances = new ArrayList<>();

    for (int i = 0; i < sampleIDs.size(); i++) {
      for (int j = i + 1; j < sampleIDs.size(); j++) {
        double pairDistance = this.findPairDistance(this.getNode(sampleIDs.get(i)).getName(), this.getNode(sampleIDs.get(j)).getName());
        pairDistances.add(pairDistance);
      }
    }

    return pairDistances;
  }

  /**
   * find the distance between two nodes by adding up the branch lengths.
   * We do this by tracing the two modes up to the top node and then comparing
   * the paths to find where they intersect.
   * e.g. if Node1 has the path 8->6->3->1
   * and Node2 has the path 10->9->5->3->1
   * <p>
   * then the path betwix the nodes is 8->6->3->5->9->10
   *
   * @param node1
   * @param node2
   * @return
   */
  public double findPairDistance(String node1, String node2) {
    if (node1.equals(node2))
      return 0;

    // first of all need to find the nodes
    int node1no = this.getNodeID(node1);
    if (node1no < 0)
      return -1.0;
    int node2no = this.getNodeID(node2);
    if (node2no < 0)
      return -1.0;

    // trace Node1 up to the top node
    java.util.ArrayList<Integer> path1 = new java.util.ArrayList<>();
    TreeNode currentTN = this.nodes.get(node1no);
    while (1 == 1) {
      path1.add(currentTN.getNodeID());
      currentTN = nodes.get(currentTN.getParentNode());
      if (currentTN.getNodeID() == 0) {
        path1.add(currentTN.getNodeID());
        break;
      }
    }

    // Trace Node2 up to the point where it intersects Node1
    java.util.ArrayList<Integer> path2 = new java.util.ArrayList<>();
    currentTN = this.nodes.get(node2no);
    int matchingNode = 0;
    while (1 == 1) {
      path2.add(currentTN.getNodeID());
      if (path1.contains(currentTN.getNodeID())) {
        matchingNode = path1.indexOf(currentTN.getNodeID());
        break;
      }
      currentTN = nodes.get(currentTN.getParentNode());
    }

    // Now join up the two paths and calculate the branch length
    // for path1 start from where we matched the nodes in the two lists
    // and work down
    double branchLength = 0.0;
    int n1 = 0;
    String path = "";
    while (n1 != matchingNode) {
      path = path.concat(path1.get(n1) + "->");
      branchLength +=
              ((TreeNode) nodes.get((Integer) path1.get(n1))).getBranchLength();
      n1++;
    }

    // For path2 start from the bottom of the list and go to the top.
    int n2 = path2.size() - 2;
    while (n2 > 0) {
      path = path.concat(path2.get(n2) + "->");
      branchLength += ((TreeNode) nodes.get((Integer) path2.get(n2))).getBranchLength();
      n2--;
    }
    path = path.concat(path2.get(n2).toString());
    branchLength += ((TreeNode) nodes.get((Integer) path2.get(n2))).getBranchLength();

    return branchLength;
  }

  /**
   * summarize the nodes in the tree
   */
  public void listNodes() {
    ListIterator itN = nodes.listIterator();
    while (itN.hasNext()) {
      TreeNode node = (TreeNode) itN.next();
      logger.debug(node.getParentNode()
              + "<--" + node.getName()
              + "(" + node.getNodeID() + ")"
              + "[" + node.getBranchLength() + "]" + "\t"
              + "\t-->[0]" + node.getLeftChild()
              + "\t-->[1]" + node.getRightChild()
      );
    }
  }

  /**
   * calculate longest distance from root to node.
   *
   * @return
   */
  public double calcTreeLength() {
    treeLength = 0.0;
    ListIterator itN = nodes.listIterator();
    while (itN.hasNext()) {
      TreeNode node = (TreeNode) itN.next();
      if (node.getNodeID() == 0) continue;
//      logger.debug("current node is " + node.getName());
      double thisBranch = node.getBranchLength();
      while (true) {
        node = nodes.get(node.getParentNode());
//        logger.debug(node.getName() + "," + node.getBranchLength());
        thisBranch += node.getBranchLength();
        if (node.getNodeID() == 0) break;
      }
      if (thisBranch > getTreeLength()) {
//        logger.debug("longest branch is now " + thisBranch);
        treeLength = thisBranch;
      }
    }
    return getTreeLength();
  }

  /**
   * find the distance from the specified node to the root of the tree
   * (even if the tree is unrooted, there is a parent node with branch length 0
   * which is always node 0)
   *
   * @param node
   * @return
   */
  public double calcTipToRootDist(TreeNode node) {
    double dist = node.getBranchLength();
    if (node.getNodeID() == 0) return 0.0; // this is the root node
    TreeNode currNode = node;
    while (true) {
      currNode = nodes.get(currNode.getParentNode());
      dist += currNode.getBranchLength();
      if (currNode.getNodeID() == 0) break;
    }
    return dist;
  }


  /**
   * Find all the nodes underneath this one.
   *
   * @param queryNode
   * @return
   */
  public ArrayList<Integer> findChildren(TreeNode queryNode) {

    ArrayList<Integer> childList = new ArrayList<>();
    for (int child : queryNode.getChildren()) {
      if (this.getNode(child).getChildren().isEmpty()) {
        childList.add(child);
      } else {
        childList.addAll(findChildren(this.getNode(child)));
      }

    }
    return childList;
  }


  /**
   * calculates the positions of the nodes within the tree assuming the
   * dimensions are 1 x 1.
   * <p>
   * this is only really useful for plotting the tree.
   */
  public void calcNodePositions() {
    double dY = 2.0 / (nodes.size() + 2.0);
    double y = dY;
    double xScale = 1.0 / this.calcTreeLength();
    /*
     * first run through the nodes and assign X and Y positions to the leaf nodes.
     * This is easy because the leaf nodes are uniformly distributed in Y
     * and the X positions are the distances from tip to root
     *
     */
    java.util.ListIterator itN = nodes.listIterator();
    while (itN.hasNext()) {
      TreeNode currNode = (TreeNode) itN.next();
      if (currNode.getChildren().isEmpty()) {
        // we have a leaf node
        currNode.setX(this.calcTipToRootDist(currNode) * xScale);
        currNode.setY(y);
        y += dY;
      }
    }
    /*
     * now we have placed the leaf nodes, we can work up the tree and
     * place the internal nodes
     */
    while (nodes.get(0).getX() < 0.0) {
      /*
       * to assign X and Y positions to a node, both children need (X,Y) assigned
       */
      itN = nodes.listIterator();
      while (itN.hasNext()) {
        TreeNode currNode = (TreeNode) itN.next();
//        logger.debug(currNode.getNodeID() + " | " + currNode.getName());
        if (currNode.getX() < 0.0) {
         /*   logger.debug("--Internal Node: children are "
                    + "NODE "+ currNode.getChildren().get(0)
                    + " & "
                    + "NODE "+ currNode.getChildren().get(1));*/
          // we have an internal node
          if (nodes.get(currNode.getChildren().get(0)).getY() >= 0
                  &&
                  nodes.get(currNode.getChildren().get(1)).getY() >= 0
          ) {
            //logger.debug("--set X & Y");
            currNode.setX(this.calcTipToRootDist(currNode) * xScale);
            currNode.setY(
                    (nodes.get(currNode.getChildren().get(0)).getY()
                            +
                            nodes.get(currNode.getChildren().get(1)).getY())
                            / 2.0);
            // logger.debug("-->"+ currNode.getX() + "\t" + currNode.getY());
          }
        } else {
          //logger.debug("-- positions set");
        }
      }

    }
  }

  /**
   * list the x, y positions of all nodes
   */
  public void reportNodePositions() {
    java.util.ListIterator itN = nodes.listIterator();
    while (itN.hasNext()) {
      TreeNode currNode = (TreeNode) itN.next();
      logger.info(currNode.getNodeID() + " | " + currNode.getName()
              + "[" + currNode.getX() + ", " + currNode.getY() + "]");
    }
  }


  /**
   * find Node by name
   *
   * @param nodeName
   * @return
   */
  public int getNodeID(String nodeName) {
    java.util.ListIterator itTr = nodes.listIterator();
    int n = 0;
    while (itTr.hasNext()) {
      TreeNode tn = (TreeNode) itTr.next();
      if (tn.getName().trim().equals(nodeName.trim()))
        return n;
      n++;
    }
    return -1;
  }


  /**
   * return the node at index int
   *
   * @param ID
   * @return
   */
  public TreeNode getNode(int ID) {
    return this.nodes.get(ID);
  }

  /**
   * adds a new group to the Tree
   * We need to ensure
   * (1) we generate a unique ID for the new group.
   * (2) the GroupName doesn't already exist in the list
   *
   * @return : the ID of the new TreeNodeGroup
   */
  public int addNewNodeGroup(String groupName) {
    int lastGroupID = -1;
    int g = 0;
    while (g < treeGroups.size()) {
      TreeNodeGroup treeNodeGroup = (TreeNodeGroup) treeGroups.get(g);
      if (treeNodeGroup.getGroupID() > lastGroupID)
        lastGroupID = treeNodeGroup.getGroupID();

      if (groupName.equals(treeNodeGroup.getGroupName())) {
        logger.warn("WARNING:group name <" + groupName + "> already exists in the group list");
        groupName = groupName.concat("_1");
        logger.warn("changing new group name to <" + groupName + ">");
        g = -1;
      }

      g++;
    }

    treeGroups.add(new TreeNodeGroup(lastGroupID + 1, groupName));
    return lastGroupID + 1;

  }

  /**
   * check whether the node is part of a existing group
   *
   * @param groupName
   * @return
   */
  public int doesNodeGroupExist(String groupName) {

    int g = 0;
    while (g < treeGroups.size()) {
      TreeNodeGroup treeNodeGroup = (TreeNodeGroup) treeGroups.get(g);
      if (treeNodeGroup.getGroupName().equals(groupName))
        return g;
    }
    return -1;
  }


  /**
   * remove the group from the TreeNodeGroupList
   * 1.  the groupID from all relevant nodes
   * 2.  delete the Group from the TreeNodeGroupList
   *
   * @param
   * @return
   */
  public int removeNodeGroup(String groupID) {
    java.util.ListIterator itN = this.getNodes().listIterator();
    int nodeCount = 0;
    while (itN.hasNext()) {
      TreeNode currNode = (TreeNode) itN.next();
      if (currNode.getGroupName() == null)
        continue;
      if (currNode.getGroupName().equals(groupID)) {
        currNode.setGroupName(null);
      }
      nodeCount++;
    }

    int g = 0;
    java.util.ListIterator itG = this.treeGroups.listIterator();
    while (itG.hasNext()) {
      TreeNodeGroup treeGroup = (TreeNodeGroup) itG.next();
      if (treeGroup.getGroupName().equals(groupID)) {
        break;
      }
      g++;
    }

    this.treeGroups.remove(g);
    return nodeCount;

  }


  /**
   * get a list of all the nodes in the specified group
   *
   * @param groupID
   * @return
   */
  public ArrayList<Integer> getNodesInGroup(String groupID) {
    ArrayList<Integer> nodeIDList = new ArrayList();

    java.util.ListIterator itN = this.getNodes().listIterator();
    int nodeCount = 0;
    while (itN.hasNext()) {
      TreeNode currNode = (TreeNode) itN.next();
      if (currNode.getGroupName() == null) {
        nodeCount++;
        continue;
      }
      if (currNode.getGroupName().equals(groupID)) {
        nodeIDList.add(nodeCount);
      }
      nodeCount++;
    }
    return nodeIDList;

  }


  /**
   * look up TreeNodeGroup entry by groupID
   *
   * @param queryGroupID
   * @return TreeNodeGroup entry or null
   */
  public TreeNodeGroup findGroupNodeByID(int queryGroupID) {

    for (TreeNodeGroup tng : treeGroups) {
      if (tng.getGroupID() == queryGroupID)
        return tng;
    }
    return null;
  }


  /**
   * look up TreeNodeGroup entry by groupID
   *
   * @param queryGroupName
   * @return TreeNodeGroup entry or null
   */
  public TreeNodeGroup findGroupNodeByName(String queryGroupName) {

    for (TreeNodeGroup tng : treeGroups) {
      if (tng.getGroupName().equals(queryGroupName))
        return tng;
    }
    return null;
  }


  /**
   * count how many leafNodes are under this node
   *
   * @param queryNode
   * @return
   */
  public int countLeafNodes(TreeNode queryNode) {

    int leafCount = 0;
    for (int child : queryNode.getChildren()) {
      if (this.getNode(child).isLeafNode()) {
        leafCount++;
      } else {
        leafCount += countLeafNodes(this.getNode(child));
      }
    }
    return leafCount;

  }


  /**
   * set the specified parameters for all child nodes (internal and leaf) of
   * the given TreeNode
   *
   * @param currentNode
   * @return the number of nodes that were affected
   */
  public int setParameterInChildren(TreeNode currentNode, String key, String value) {
    int nodeCount = 0;
    for (int child : currentNode.getChildren()) {
      this.getNode(child).setParameter(key, value);
      logger.debug("set highlight parameters for node <" + this.getNode(child).getNodeID()
              + "|" + this.getNode(child).getName() + "> to (" + key + "," + value + ")");

      nodeCount++;
      if (this.getNode(child).isInternalNode()) {
        nodeCount += setParameterInChildren(this.getNode(child), key, value);
      }
    }

    return nodeCount;
  }


  /**
   * set child nodes to specified GroupID
   *
   * @param currentNode
   * @param groupID
   * @return
   */
  public int setGroupNameInChildren(TreeNode currentNode, String groupID) {
    int nodeCount = 0;
    for (int child : currentNode.getChildren()) {
      this.getNode(child).setGroupName(groupID);
      logger.debug("set groupID for node <" + this.getNode(child).getNodeID()
              + "|" + this.getNode(child).getName() + "> to " + groupID);

      nodeCount++;
      if (this.getNode(child).isInternalNode()) {
        nodeCount += setGroupNameInChildren(this.getNode(child), groupID);
      }
    }

    return nodeCount;
  }


  /**
   * Tree parsing methods
   */


  public void parseNewickTree() {
    this.removeWhiteSpaceFromTreeString();
    try {
      String remainingTreeLine = this.findNewickParentNode();
      logger.debug(remainingTreeLine);
      logger.debug(this.getNode(0).printParameters());
      remainingTreeLine = remainingTreeLine.substring(1, remainingTreeLine.length() - 1);
      this.parseOutNewickChildren(remainingTreeLine, 0, 0);
      logger.debug("\n" + this.printNodes());
      logger.info("\n" + this.printNodes());
      logger.info("finished");
    } catch (IOException exIO) {
      logger.info(exIO.toString());
    }

  }

  /**
   * Find parent node for the tree. the search method depends on whether we have NEXUS, NEWICK or BEAST tree
   *
   * @return
   * @throws IOException
   */
  private String findParentNode() throws IOException {
    String remainingTreeString = "";
    if (this.getTreeType() == TreeType.BEAST) {
      remainingTreeString = findBeastParentNode();
    } else {
      if (this.getTreeType() == TreeType.NEWICK) {
        remainingTreeString = findNewickParentNode();
      }
    }
    return remainingTreeString;
  }


  /**
   * There are four options for the root node:
   * 1. ()[];
   * 2. ''[], ();
   * 3. ()[], ''[]
   * 4. throw an exception
   * <p>
   * moreover, some trees place a 0.0 branch length after the last []
   *
   * @return remainder of the treeString that needs to be parsed
   * @throws IOException
   */
  private String findBeastParentNode() throws IOException {
    logger.debug("--<findBeastParentNode>");
    if (getTreeLine().endsWith(";"))
      setTreeLine(StringUtils.chop(getTreeLine()));

    TreeNode rootNode = new TreeNode(0, getTreeType());
    rootNode.setParentNode(-1);
    rootNode.setParameter(TreeNode.DEFAULT_LINECOLOUR_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINECOLOUR_KEYSTRING));
    rootNode.setParameter(TreeNode.DEFAULT_LINEWIDTH_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINEWIDTH_KEYSTRING));

    /*
     * this is the same basic logic as findNewickParentNode, except we have to take into account the []
     */
    if (getTreeLine().startsWith("(") && getTreeLine().endsWith("]") || getTreeLine().startsWith("(") && getTreeLine().endsWith("0")) {
      // not sure if this logic is correct. can you have ()[]? (i.e. parameters for the unspecified root node)
      // yes you can, and some trees place a zero branch length afterwards too, i.e. []:0.0;
      // Case 1: ()
      int lastLSquareParenPos = getTreeLine().lastIndexOf("[");
      String nodeString = getTreeLine().substring(lastLSquareParenPos, getTreeLine().length());
      rootNode.parseNewickNodeString(nodeString);
      rootNode.setName("NODE:0");
      this.getNodes().add(rootNode);
      return getTreeLine();
    } else {
      if (getTreeLine().endsWith(")")) {
        // Case 2: N[], () , so need to find the first occurrence of ']'
        String nodeString = getTreeLine().substring(0, getTreeLine().indexOf("]"));
        rootNode.parseNewickNodeString(nodeString);
        rootNode.printParameters();
        // NAME[]:dist
        this.getNodes().add(rootNode);
        return getTreeLine().substring(getTreeLine().indexOf("("));
      } else {
        if (getTreeLine().startsWith("(")) {
          // Case 3: (), N[], so need to find the last occurrence of '[' and then step back to the previous ','
          int lastLSquareParenPos = getTreeLine().lastIndexOf("[");
          int lastNodeStartPos = getTreeLine().lastIndexOf(",", lastLSquareParenPos);
          String nodeString = getTreeLine().substring(lastNodeStartPos, getTreeLine().length() - 1);
          rootNode.parseNewickNodeString(nodeString);
          rootNode.printParameters();
          this.getNodes().add(rootNode);
          return getTreeLine().substring(0, getTreeLine().indexOf(")"));
        } else {
          throw new IOException("unrecognized tree structure while searching for parent node");
        }
      }

    }


  }


  /**
   * There are four options for the root node:
   * 1. (); or ()[];
   * 2. NODE, ();
   * 3. (), NODE;
   * 4. throw an exception
   *
   * @return remainder of the treeString that needs to be parsed
   * @throws IOException
   */
  private String findNewickParentNode() throws IOException {
    logger.debug("--<findNewickParentNode>");
    if (getTreeLine().endsWith(";"))
      setTreeLine(StringUtils.chop(getTreeLine()));

    TreeNode rootNode = new TreeNode(0, getTreeType());
    rootNode.setParentNode(-1);
    rootNode.setParameter(TreeNode.DEFAULT_LINECOLOUR_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINECOLOUR_KEYSTRING));
    rootNode.setParameter(TreeNode.DEFAULT_LINEWIDTH_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINEWIDTH_KEYSTRING));

    if (getTreeLine().startsWith("(") && getTreeLine().endsWith(")")) {
      // Case 1: ()
      rootNode.setName("NODE:0");
      this.getNodes().add(rootNode);
      return getTreeLine();
    } else {
      if (getTreeLine().endsWith(")")) {
        // Case 2: N, ()
        String nodeString = getTreeLine().substring(0, getTreeLine().indexOf("("));
        rootNode.parseNewickNodeString(nodeString);
        rootNode.printParameters();
        // NAME[]:dist
        this.getNodes().add(rootNode);
        return getTreeLine().substring(getTreeLine().indexOf("("));
      } else {
        if (getTreeLine().startsWith("(")) {
          // Case 3: (), N
          String nodeString = getTreeLine().substring(getTreeLine().indexOf(")") + 1, getTreeLine().length() - 1);
          rootNode.parseNewickNodeString(nodeString);
          rootNode.printParameters();
          this.getNodes().add(rootNode);
          return getTreeLine().substring(0, getTreeLine().indexOf(")"));
        } else {
          throw new IOException("unrecognized tree structure while searching for parent node");
        }
      }

    }


  }

  /**
   * try to figure out the TreeType from the tree string. If we can't work this out, then we can't parse
   * out the string.
   * This is a bit of a lazy search method. For example, if a name was '[sample]', the tree would be identified incorrectly as BEAST
   *
   * @return TreeType
   */
  public static final TreeType findTreeType(String queryTreeString) {
    logger.info("--<findTreeType>");
    TreeType queryTreeType = TreeType.UNKNOWN;
    if (queryTreeString.contains("[")) {
      queryTreeType = TreeType.BEAST;
      return queryTreeType;
    }
    // if we have ')number:number' then we have a (bootstrap:branch length)/NEWICKBOOTSTRAP tree 
    Pattern pattern = Pattern.compile("\\)\\d+:\\d+", Pattern.CASE_INSENSITIVE);
    Matcher matcher = pattern.matcher(queryTreeString);
    if (matcher.find())
      return TreeType.NEWICKBOOTSTRAP;
    else
      return TreeType.NEWICK;
    
    /*
	  int startPos = queryTreeString.indexOf(")");
	  int endPos = queryTreeString.indexOf(",", startPos);
	  String testString = queryTreeString.substring(startPos, endPos);
	  if(testString.split(":").length == 2){
		  queryTreeType = TreeType.NEWICKBOOTSTRAP;
	  }else {
		  if(testString.split(":").length == 1) {
			  queryTreeType = TreeType.NEWICK;
		  }else {
			  queryTreeType = TreeType.NEWICK;
		  }
	  }
    */
    //return queryTreeType;

  }


  /**
   * remove all whitespace from the tree (leaving it in confuses the parser)
   *
   * @return number of whitespace characters that were removed
   */
  public int removeWhiteSpaceFromTreeString() {
    int oldLength = this.getTreeLine().length();
    this.setTreeLine(this.getTreeLine().replaceAll("\\s+", ""));
    return oldLength - this.getTreeLine().length();
  }

  /**
   * possibilities are Leaf Node or Internal Node
   * Internal Nodes are delimited by () and can contain additional combinations of L and I
   * Combinations are
   * (L,L....)
   * ((I),L...)
   * (L,(I)....)
   * ((I),(I)...)
   *
   * @param treeStringFrag
   * @param parentID
   * @return
   */
  private String parseOutNewickChildren(String treeStringFrag, int parentID, int level) {
    logger.debug("--<parseOutChildren2:[" + level + "]>");
    logger.info("--<parseOutChildren2:[" + level + "]>");
    while (treeStringFrag.length() > 0) {
      if (this.isNextNodeLeaf(treeStringFrag)) {
        String leafNodeString = this.getNextNewickLeafNodeString(treeStringFrag, 0);

        // create a new leaf node from this string
        TreeNode leafNode = new TreeNode(nodes.size(), getTreeType());
        leafNode.setParentNode(parentID);
        leafNode.setParameter(TreeNode.DEFAULT_LINECOLOUR_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINECOLOUR_KEYSTRING));
        leafNode.setParameter(TreeNode.DEFAULT_LINEWIDTH_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINEWIDTH_KEYSTRING));
        leafNode.parseNewickNodeString(leafNodeString);
        nodes.get(parentID).getChildren().add(nodes.size()); // add this as a child of the parent node
        nodes.add(leafNode);
        logger.info("-- Added Leaf Node <" + leafNode.getNodeID() + ">|" + leafNode.getName() + "-->" + parentID);
        logger.debug("-- Added Leaf Node <" + leafNode.getNodeID() + ">|" + leafNode.getName() + "-->" + parentID);

        // remove this leafNodeString from the treeString and pass on what is left.
        if (leafNodeString.equals(treeStringFrag)) {
          treeStringFrag = "";
        } else {
          treeStringFrag = treeStringFrag.substring(leafNodeString.length() + 1, treeStringFrag.length()).trim();
        }
      } else {
        // create a new internal node
        TreeNode intNode = new TreeNode(nodes.size(), getTreeType());
        intNode.setParentNode(parentID);
        nodes.get(parentID).getChildren().add(nodes.size()); // add this as a child of the parent node
        intNode.setName("INTNODE:" + nodes.size());
        intNode.setParameter(TreeNode.DEFAULT_LINECOLOUR_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINECOLOUR_KEYSTRING));
        intNode.setParameter(TreeNode.DEFAULT_LINEWIDTH_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINEWIDTH_KEYSTRING));
        nodes.add(intNode);
        logger.info("-- Added Internal Node <" + intNode.getNodeID() + ">|" + intNode.getName() + "-->" + parentID);
        logger.debug("-- Added Internal Node <" + intNode.getNodeID() + ">|" + intNode.getName() + "-->" + parentID);
        String intNodeStringAndParams = this.getNextInternalNodeString(treeStringFrag, 0);
        logger.debug("----intNodeStringAndParams <" + intNodeStringAndParams + ">");

        // This internal node string needs to be parsed out
        // first need to parse the internal node data such as bootstrapValue and branchLength
        String intNodeStringNoParams = this.parseOutInternalNodeParametersAndTrim(intNode.getNodeID(), intNodeStringAndParams);
        logger.debug("----intNodeStringNoParams <" + intNodeStringNoParams + ">");

        // now parse out the contents inside the parenthesis
        // The nodeID of this internal node will be the parentID of any children
        level = level + 1;
        this.parseOutNewickChildren(intNodeStringNoParams, intNode.getNodeID(), level);

        /* remove this internal node from the string and pass on what is left.
         * If we have additional nodes within the string that need parsing then
         *   because we removed the left parenthesis but kept the right one, we need to add +1 or +2 to the start pos
         *   +1 for the right parenthesis and +1 for the "," that follows the end of the internal node
         * otherwise, if it's the last node, then we just empty the string
         *
         * what happens here depends on what the next node is
         * If it's a leaf node, remove 2, if it's another closing ')', then remove 1
         */
        if ((treeStringFrag.length() - intNodeStringAndParams.length()) > 1) {
//				  logger.debug("--<treeStringSplit[" + treeStringFrag.substring(intNodeStringAndParams.length()-3, 
//						  intNodeStringAndParams.length()+ 15) +  "]>");
//				  logger.debug("--<treeStringSplit[--|---------|-----]>");

          treeStringFrag = treeStringFrag.substring(intNodeStringAndParams.length() + 2, treeStringFrag.length()).trim();
          logger.debug(treeStringFrag);
          logger.debug(treeStringFrag);
        } else {
          treeStringFrag = ""; // last
        }
      }

    }


    logger.debug("--<finished[" + level + "]>");
    return "";
  }


  /**
   * possibilities are Leaf Node or Internal Node
   * Internal Nodes are delimited by () and can contain additional combinations of L and I
   * Combinations are
   * (L,L....)
   * ((I),L...)
   * (L,(I)....)
   * ((I),(I)...)
   *
   * @param treeStringFrag
   * @param parentID
   * @return
   */
  private String parseOutBeastChildren(String treeStringFrag, int parentID, int level) {
    logger.debug("--<parseOutBeastChildren:[" + level + "]>");
    while (treeStringFrag.length() > 0) {
      if (this.isNextNodeLeaf(treeStringFrag)) {
        String leafNodeString = this.getNextBeastLeafNodeString(treeStringFrag, 0);

        // create a new leaf node from this string
        TreeNode leafNode = new TreeNode(nodes.size(), getTreeType());
        leafNode.setParentNode(parentID);
        leafNode.setParameter(TreeNode.DEFAULT_LINECOLOUR_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINECOLOUR_KEYSTRING));
        leafNode.setParameter(TreeNode.DEFAULT_LINEWIDTH_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINEWIDTH_KEYSTRING));
        leafNode.parseBeastParameters(leafNodeString);
        nodes.get(parentID).getChildren().add(nodes.size()); // add this as a child of the parent node
        nodes.add(leafNode);
        logger.info("-- Added Leaf Node <" + leafNode.getNodeID() + ">|" + leafNode.getName() + "-->" + parentID);
        logger.debug("-- Added Leaf Node <" + leafNode.getNodeID() + ">|" + leafNode.getName() + "-->" + parentID);

        // remove this leafNodeString from the treeString and pass on what is left.
        if (leafNodeString.equals(treeStringFrag)) {
          treeStringFrag = "";
        } else {
          treeStringFrag = treeStringFrag.substring(leafNodeString.length() + 1, treeStringFrag.length()).trim();
        }
      } else {
        // create a new internal node
        TreeNode intNode = new TreeNode(nodes.size(), getTreeType());
        intNode.setParentNode(parentID);
        nodes.get(parentID).getChildren().add(nodes.size()); // add this as a child of the parent node
        intNode.setName("INTNODE:" + nodes.size());
        intNode.setParameter(TreeNode.DEFAULT_LINECOLOUR_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINECOLOUR_KEYSTRING));
        intNode.setParameter(TreeNode.DEFAULT_LINEWIDTH_KEYSTRING, this.getPigtreeParameterValue(WivTree.DEFAULT_LINEWIDTH_KEYSTRING));
        nodes.add(intNode);
        logger.info("-- Added Internal Node <" + intNode.getNodeID() + ">|" + intNode.getName() + "-->" + parentID);
        logger.debug("-- Added Internal Node <" + intNode.getNodeID() + ">|" + intNode.getName() + "-->" + parentID);
        String intNodeStringAndParams = this.getNextBeastInternalNodeString(treeStringFrag, 0);
        logger.debug("----intNodeStringAndParams <" + intNodeStringAndParams + ">");

        // This internal node string needs to be parsed out
        // first need to parse the internal node data such as bootstrapValue and branchLength
        String intNodeStringNoParams = this.parseOutInternalNodeParametersAndTrim(intNode.getNodeID(), intNodeStringAndParams);
        logger.debug("----intNodeStringNoParams <" + intNodeStringNoParams + ">");

        // now parse out the contents inside the parenthesis
        // The nodeID of this internal node will be the parentID of any children
        level = level + 1;
        this.parseOutBeastChildren(intNodeStringNoParams, intNode.getNodeID(), level);

        /* remove this internal node from the string and pass on what is left.
         * If we have additional nodes within the string that need parsing then
         *   because we removed the left parenthesis but kept the right one, we need to add +1 or +2 to the start pos
         *   +1 for the right parenthesis and +1 for the "," that follows the end of the internal node
         * otherwise, if it's the last node, then we just empty the string
         *
         * what happens here depends on what the next node is
         * If it's a leaf node, remove 2, if it's another closing ')', then remove 1
         */
        if ((treeStringFrag.length() - intNodeStringAndParams.length()) > 1) {
//				  logger.debug("--<treeStringSplit[" + treeStringFrag.substring(intNodeStringAndParams.length()-3, 
//						  intNodeStringAndParams.length()+ 15) +  "]>");
//				  logger.debug("--<treeStringSplit[--|---------|-----]>");

          treeStringFrag = treeStringFrag.substring(intNodeStringAndParams.length() + 2, treeStringFrag.length()).trim();
          logger.debug(treeStringFrag);
          logger.debug(treeStringFrag);
        } else {
          treeStringFrag = ""; // last
        }
      }

    }


    logger.debug("--<finished[" + level + "]>");
    return "";
  }


  /**
   * The next leaf node in the query string is delimited by a ','
   * This is only placed in a separate method for readability
   *
   * @param queryString
   * @param startPos
   * @return
   */
  private String getNextNewickLeafNodeString(String queryString, int startPos) {
    return queryString.split(",")[0].trim();
  }

  /**
   * The next leaf node in the query string is delimited by a ','
   * This is only placed in a separate method for readability
   *
   * @param queryString
   * @param startPos
   * @return
   */
  private String getNextBeastLeafNodeString(String queryString, int startPos) {
    int sqRightParenPos = queryString.indexOf("]");
    int nextCommaPos = queryString.indexOf(",", sqRightParenPos);
    if (nextCommaPos > 0)
      return queryString.substring(0, nextCommaPos);
    else
      return queryString;
  }

  /**
   * return position of next '('.  allowing for the possibility that the sample name may also contain parentheses
   * At the moment, this hasn't been implemented as there isn't much need at present to be able to handle this much
   * flexibility in the sample name
   *
   * @return
   */
  private int getNextLeftParenthesisPos() {
    return -1;
  }

  /**
   * parse out internal node by searching for matching parenthesis '()'
   * also need to check for optional branch length and bootstrap values after the right parenthesis
   *
   * @param queryString
   * @param startPos
   * @return
   */
  private String getNextInternalNodeString(String queryString, int startPos) {
    logger.debug("--<getNextNewickInternalNodeString>");
    /*
     * start with a simple check. If there is
     *   (i) only one '(' and ')'
     *   and we don't have something like
     *   (...)branchLength
     * then we have the node string already.
     *
     */
    if (StringUtils.countMatches(queryString, ")") == 1 && StringUtils.countMatches(queryString, ")") == 1) {
      int nextrightParenLoc = queryString.indexOf(")", startPos + 1);
      int nextCommaLoc = queryString.indexOf(",", nextrightParenLoc);
      if (nextCommaLoc == -1) {
        logger.debug("--<return simple match>");
        return queryString.substring(startPos + 1, queryString.length());
      }
    }

    int depth = 0;

    //if we find another "(" before we find a ")" then we have to go deeper
    //( but we also need to check for the presence of ' ')  <--
    int start = startPos;
    int nextleftParenLoc = queryString.indexOf("(", start + 1);
    int nextrightParenLoc = queryString.indexOf(")", start + 1);
    Boolean match = false;
    while (match == false) {
      logger.debug("Depth<" + depth + ">");
      if (nextleftParenLoc >= 0)
        logger.debug(" LEFT:[" + nextleftParenLoc + "]");
      else
        logger.debug(" LEFT:[" + nextleftParenLoc + "\t" + "NO MORE MATCHES");

      if (nextrightParenLoc >= 0)
        logger.debug("RIGHT:[" + nextrightParenLoc + "]");
      else
        logger.debug("RIGHT:[" + nextrightParenLoc + "\t" + "NO MORE MATCHES");

      if (nextleftParenLoc < nextrightParenLoc && nextleftParenLoc >= 0) {
        depth = depth + 1;
        nextleftParenLoc++;
        nextleftParenLoc = queryString.indexOf("(", nextleftParenLoc);
        logger.debug("+Depth[" + depth + "] left|right" + nextleftParenLoc + "|" + nextrightParenLoc);
      } else {
        if (depth == 0) {
          match = true;
        } else {
          depth--;
          nextrightParenLoc++;
          nextrightParenLoc = queryString.indexOf(")", nextrightParenLoc);
          logger.debug("-Depth[" + depth + "] left|right" + nextleftParenLoc + "|" + nextrightParenLoc);
        }
      }

    }
    if (nextleftParenLoc == -1)
      nextleftParenLoc = 0;
    if (nextrightParenLoc == -1)
      nextrightParenLoc = queryString.length();
    /*
     * lastly, we need to check whether we have branch length and bootstrap information
     * following on from the right ')'
     *   if next char is ')' or ',' we don't have any following information
     *
     *   For the following:
     *   if next char is ":" we have ):branchLength
     *   if next char is a digit, we have )bootstrapValue:branchLength
     *   if next char is "'", we have a MEGA style label for the node
     *     )'label':branchLength or  )'label':bootstrapValue:branchLength
     *
     *  for all these cases, we just look for the next ')' or ',', whichever comes first
     *  we don't need to know what we have, we just want to return the string
     *
     *
     */
    if (queryString.substring(nextrightParenLoc + 1, nextrightParenLoc + 2).equalsIgnoreCase(")")
            || queryString.substring(nextrightParenLoc + 1, nextrightParenLoc + 2).equalsIgnoreCase(",")) {
      logger.debug("--<return string>");
      return queryString.substring(startPos + 1, nextrightParenLoc);
    }

    /**
     * check for BEAST format '[]'
     * need to check this case first as there may be ',' within the '[]'
     */

    /*
     * :branchLength or bootstrap:branchLength
     * Options are:
     * number terminated by ')' or ',', or by the end of the string
     * 1. )number)
     * 2. )number,
     * 3. )number
     */
    int nextNextrightParenLoc = queryString.indexOf(")", nextrightParenLoc + 3);
    int nextNextrightComma = queryString.indexOf(",", nextrightParenLoc + 3);
    int endOfNumber = -1;
    if (nextNextrightComma > 0 && nextNextrightParenLoc > 0) {
      if (nextNextrightComma < nextNextrightParenLoc)
        endOfNumber = nextNextrightComma;
      else
        endOfNumber = nextNextrightParenLoc;
    } else {
      if (nextNextrightComma > 0) {
        endOfNumber = nextNextrightComma;
      } else {
        if (nextNextrightParenLoc > 0) {
          endOfNumber = nextNextrightParenLoc;
        } else {
          // if we get here something odd happened and we should throw an error
          endOfNumber = queryString.length() - 1;
          //logger.error("couldn't make sense of the tree string <" + queryString + ">");

        }
      }
    }

    logger.debug("--<return branch length>");
    return queryString.substring(startPos + 1, endOfNumber);


  }


  /**
   * parse out internal node by searching for matching parenthesis '()'
   * also need to check for optional branch length and bootstrap values after the right parenthesis
   *
   * @param queryString
   * @param startPos
   * @return
   */
  private String getNextBeastInternalNodeString(String queryString, int startPos) {
    logger.debug("--<getNextBeastInternalNodeString>");
    /*
     * start with a simple check. If there is
     *   (i) only one '(' and ')'
     *   and we don't have something like
     *   (...)branchLength
     * then we have the node string already.
     *
     */
    if (StringUtils.countMatches(queryString, "(") == 1 && StringUtils.countMatches(queryString, ")") == 1) {
      int nextrightParenLoc = queryString.indexOf(")", startPos + 1);
      int nextCommaLoc = queryString.indexOf(",", nextrightParenLoc);
      if (nextCommaLoc == -1) {
        logger.debug("--<return simple match>");
        return queryString.substring(startPos + 1, queryString.length());


      }
    }

    int depth = 0;

    //if we find another "(" before we find a ")" then we have to go deeper
    //( but we also need to check for the presence of ' ')  <--
    int start = startPos;
    int nextleftParenLoc = queryString.indexOf("(", start + 1);
    int nextrightParenLoc = queryString.indexOf(")", start + 1);
    Boolean match = false;
    while (match == false) {
      logger.debug("Depth<" + depth + ">");
      if (nextleftParenLoc >= 0)
        logger.debug("Next LEFT:[" + nextleftParenLoc + "]");
      else
        logger.debug("Next LEFT:[" + nextleftParenLoc + "\t" + "NO MORE MATCHES");

      if (nextrightParenLoc >= 0)
        logger.debug("Next RIGHT:[" + nextrightParenLoc + "]");
      else
        logger.debug("Next RIGHT:[" + nextrightParenLoc + "\t" + "NO MORE MATCHES");

      if (nextleftParenLoc < nextrightParenLoc && nextleftParenLoc >= 0) {
        depth = depth + 1;
        nextleftParenLoc++;
        nextleftParenLoc = queryString.indexOf("(", nextleftParenLoc);
        logger.debug("+Depth[" + depth + "] (left|right)" + nextleftParenLoc + "|" + nextrightParenLoc);
      } else {
        if (depth == 0) {
          match = true;
        } else {
          depth--;
          nextrightParenLoc++;
          nextrightParenLoc = queryString.indexOf(")", nextrightParenLoc);
          logger.debug("-Depth[" + depth + "] (left|right)" + nextleftParenLoc + "|" + nextrightParenLoc);
        }
      }

    }
    if (nextleftParenLoc == -1)
      nextleftParenLoc = 0;
    if (nextrightParenLoc == -1)
      nextrightParenLoc = queryString.length();


    /**
     * 1. check for BEAST format '[]'
     * need to check this case first as there may be ',' within the '[]'
     *
     * (I think) the next character must be a '[', but there may or may not be a ':number' afterwards
     *
     */
    if (queryString.substring(nextrightParenLoc + 1, nextrightParenLoc + 2).equals("[")) {
      // get position of ']'
      int nextrightSquareParenLoc = queryString.indexOf("]", nextrightParenLoc + 3);
      if (queryString.substring(nextrightSquareParenLoc + 1, nextrightSquareParenLoc + 2).equals(":")) {
        // get position of ']'
        /*
         * :branchLength (but i don't think we get bootstrap:branchLength)
         * Options are:
         * number terminated by ')' or ',', or by the end of the string
         * 1. )[]number)
         * 2. )[]number,
         * 3. )[]number
         */
        int nextNextrightParenLoc = queryString.indexOf(")", nextrightSquareParenLoc + 1);
        int nextNextrightComma = queryString.indexOf(",", nextrightSquareParenLoc + 1);
        int endOfNumber = -1;
        if (nextNextrightComma > 0 && nextNextrightParenLoc > 0) {
          if (nextNextrightComma < nextNextrightParenLoc)
            endOfNumber = nextNextrightComma;
          else
            endOfNumber = nextNextrightParenLoc;
        } else {
          if (nextNextrightComma > 0) {
            endOfNumber = nextNextrightComma;
          } else {
            if (nextNextrightParenLoc > 0) {
              endOfNumber = nextNextrightParenLoc;
            } else {
              // if we get here something odd happened and we should throw an error
              endOfNumber = queryString.length() - 1;
              logger.error("couldn't make sense of the tree string <" + queryString + ">");

            }
          }
        }
        logger.debug("--<return branch length>");
        logger.debug("--" + queryString.substring(startPos + 1, endOfNumber));
        logger.debug("\n\n");
        return queryString.substring(startPos + 1, endOfNumber);
      }

    }


    /*
     * 2. check whether we have branch length and bootstrap information
     * following on from the right ')'
     *   if next char is ')' or ',' we don't have any following information
     *
     *   For the following:
     *   if next char is ":" we have ):branchLength
     *   if next char is a digit, we have )bootstrapValue:branchLength
     *   if next char is "'", we have a MEGA style label for the node
     *     )'label':branchLength or  )'label':bootstrapValue:branchLength
     *
     *  for all these cases, we just look for the next ')' or ',', whichever comes first
     *  we don't need to know what we have, we just want to return the string
     *
     *
     */
    if (queryString.substring(nextrightParenLoc + 1, nextrightParenLoc + 2).equalsIgnoreCase(")")
            || queryString.substring(nextrightParenLoc + 1, nextrightParenLoc + 2).equalsIgnoreCase(",")) {
      logger.debug("--<return string>");
      return queryString.substring(startPos + 1, nextrightParenLoc);
    }

    logger.debug(queryString.substring(startPos + 1, nextrightParenLoc));

    /*
     * :branchLength or bootstrap:branchLength
     * Options are:
     * number terminated by ')' or ',', or by the end of the string
     * 1. )number)
     * 2. )number,
     * 3. )number
     */
    int nextNextrightParenLoc = queryString.indexOf(")", nextrightParenLoc + 3);
    int nextNextrightComma = queryString.indexOf(",", nextrightParenLoc + 3);
    int endOfNumber = -1;
    if (nextNextrightComma > 0 && nextNextrightParenLoc > 0) {
      if (nextNextrightComma < nextNextrightParenLoc)
        endOfNumber = nextNextrightComma;
      else
        endOfNumber = nextNextrightParenLoc;
    } else {
      if (nextNextrightComma > 0) {
        endOfNumber = nextNextrightComma;
      } else {
        if (nextNextrightParenLoc > 0) {
          endOfNumber = nextNextrightParenLoc;
        } else {
          // if we get here something odd happened and we should throw an error
          endOfNumber = queryString.length() - 1;
          logger.error("couldn't make sense of the tree string <" + queryString + ">");

        }
      }
    }

    logger.debug("--<return branch length>");
    return queryString.substring(startPos + 1, endOfNumber);


  }


  /**
   * Need to check for the following options:
   * if next char is "[" we have BEAST )[]:branchLength
   * if next char is ":" we have ):branchLength
   * if next char is a digit, we have )bootstrapValue:branchLength
   * if next char is "'", we have a MEGA style label for the node
   * )'label':branchLength or  )'label':bootstrapValue:branchLength
   *
   * @param nodeID
   * @param internalNodeString
   */
  private String parseOutInternalNodeParametersAndTrim(int nodeID, String internalNodeString) {
    logger.debug("--<parseOutInternalNodeParametersAndTrim>");
    int lastRightParenPos = internalNodeString.lastIndexOf(")");
    String paramString = internalNodeString.substring(lastRightParenPos + 1).trim();

    String firstChar = paramString.substring(0, 1);
    if (firstChar.equals("[")) {
      this.getNode(nodeID).parseBeastParameters(paramString);
    }


    // branchLength only
    if (firstChar.equals(":")) {
      Double branchLength = Double.parseDouble(paramString.substring(1, paramString.length()));
      this.getNode(nodeID).setBranchLength(branchLength);
    }

    if (firstChar.equals("'")) {
      int nextQuotePos = paramString.indexOf("'", 1);
      String nodeName = paramString.substring(0, nextQuotePos);
      this.getNode(nodeID).setName(nodeName);
      /*
       * next, we need to figure out whether we have bootstrapValue and branchLength values
       */

      int firstColonPos = paramString.indexOf(":", nextQuotePos + 1);
      int secondColonPos = paramString.indexOf(":", firstColonPos + 1);

      if (secondColonPos > 0) {
        // we have both bootstrapValue and branchLength
        Double bootstrapValue = Double.parseDouble(paramString.substring(firstColonPos + 1, secondColonPos));
        Double branchLength = Double.parseDouble(paramString.substring(secondColonPos + 1));
        this.getNode(nodeID).setBootStrap(bootstrapValue);
        this.getNode(nodeID).setBranchLength(branchLength);

      } else {
        Double branchLength = Double.parseDouble(paramString.substring(firstColonPos + 1));
        this.getNode(nodeID).setBranchLength(branchLength);
      }
    }

    if (StringUtils.isNumeric(firstChar)) {
      int firstColonPos = paramString.indexOf(":", 1);
      Double bootstrapValue = Double.parseDouble(paramString.substring(0, firstColonPos));
      Double branchLength = Double.parseDouble(paramString.substring(firstColonPos + 1));
      this.getNode(nodeID).setBootStrap(bootstrapValue);
      this.getNode(nodeID).setBranchLength(branchLength);

    }

    logger.debug("--<finished>");
    return internalNodeString.substring(0, lastRightParenPos);
  }


  /**
   * added for readability
   * If the string starts with "(" it's an internal node
   * for example: ():1,L,():1,L,():1
   *
   * @param treeStringFrag
   * @return
   */
  private Boolean isNextNodeLeaf(String treeStringFrag) {
    return !treeStringFrag.startsWith("(");
  }


  /**
   * print out all nodes in the tree in short form (no parameter information)
   *
   * @return
   */
  public String printNodes() {
    String nodeString = "";

    for (TreeNode node : this.getNodes()) {
      nodeString = nodeString.concat(node.printNodeShort() + "+" + StringUtils.repeat("-", 40) + "+\n");
    }
    return nodeString;
  }


  /**
   * remove the nodes and reset parameters
   */
  public void dropTree() {
    nodes.clear();
    setTreeLine("");
    treeLength = 0;
    taxaNames.clear();
  }


  /**
   * return the parentID for this node
   *
   * @param thisNode
   * @return
   */
  public int getParentID(TreeNode thisNode) {
    return thisNode.getParentNode();
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
   * @return the nodes
   */
  public java.util.ArrayList<TreeNode> getNodes() {
    return nodes;
  }

  /**
   * @return the treeLength
   */
  public double getTreeLength() {
    return treeLength;
  }

  /**
   * @return the treeType
   */
  public TreeType getTreeType() {
    return treeType;
  }


  public void setPigtreeParameter(String key, String value) {
    pigTreeParams.put(key, value);
  }

  public String getPigtreeParameterValue(String key) {
    return this.pigTreeParams.get(key);
  }

  public void setFigtreeParameter(String key, String value) {
    figTreeParams.put(key, value);
  }

  public String getFigtreeParameterValue(String key) {
    return this.figTreeParams.get(key);
  }

  private void setTreeType(TreeType treeType) {
    this.treeType = treeType;
  }

  /**
   * @return the sampleClassificationList
   */
  public SampleClassificationList getSampleClassificationList() {
    return sampleClassificationList;
  }

  /**
   * @param sampleClassificationList the sampleClassificationList to set
   */
  public void setSampleClassificationList(SampleClassificationList sampleClassificationList) {
    this.sampleClassificationList = sampleClassificationList;
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
   * @return the sampleListFilename
   */
  public String getSampleListFilename() {
    return sampleListFilename;
  }

  /**
   * @param sampleListFilename the sampleListFilename to set
   */
  public void setSampleListFilename(String sampleListFilename) {
    this.sampleListFilename = sampleListFilename;
  }


  // Noras adding

  /**
   * Initiate TreeNodeGroups.
   * All nodes are added to a group (not only leafs).
   * Both superpops (types) and subpops (names) are added as individual groups.
   */
  public void createTreeNodeGroups() {

    // This was tested thoroughly 23.03.
    ArrayList<String> classificationTypes = this.getSampleClassificationList().getClassificationTypes();
    for (String classificationType : classificationTypes) {
      ArrayList<String> classificationNames = this.getSampleClassificationList().getAllEntriesForClassificationType(classificationType);
      for (String classificationName : classificationNames) {
        int nodeGroupID = this.addNewNodeGroup(classificationName);   // Add group and retrieve ID
        TreeNodeGroup currentGroupNode = findGroupNodeByID(nodeGroupID);  // fetch object to set parameter
        currentGroupNode.setParameter("ClassificationType", classificationType);   // Set parameter. May be redundant.
      }
    }
    // Add nodes to groups
    for (TreeNode node : this.getNodes()) {
      if (node.isLeafNode()) {

        String classificationSuper = node.getName().split("___")[0];   // Eg. EUR
        String classificationSub = node.getName().split("___")[1];   // Eg. GBR

        // Obtain matching group objects for node
        TreeNodeGroup groupNodeSuper = findGroupNodeByName(classificationSuper);   // Lager du ny eller henter du gammel?
        TreeNodeGroup groupNodeSub = findGroupNodeByName(classificationSub);

        // Update list of node IDs for group objects
        groupNodeSuper.addLeafNodeToGroup(node.getNodeID());
        groupNodeSub.addLeafNodeToGroup(node.getNodeID());
      }
    }
  }

  /**
   * Calculates SDR for all groups in tree (both super and sub).
   * SDR = mean within-subtype pairwise distance / mean between-subtype pairwise
   * Idea: Add option to calculate for only super pops or all.
   */

  public void calculateSDRs () {

    this.createTreeNodeGroups();

    // Calculate mean-within pairwise distances
    HashMap<String, Double> meanWithinDistances = new HashMap<>();

    for (TreeNodeGroup treeGroup : this.treeGroups) {
      ArrayList<Integer> IDs = treeGroup.getLeafNodesInGroup();
      ArrayList<Double> pairDistances = getPairDistancesForAllPairsInSet(IDs);

      // Test print out
      if (treeGroup.getGroupName().equals("ESN")) {
        System.out.println("ESN ----------------");
        System.out.println("IDs: " + IDs.toString());
        System.out.println("pairDistances: " + pairDistances.toString());
      }

      double sumDist = 0;
      if (!pairDistances.isEmpty()) {
        for (Double d : pairDistances) sumDist += d;
        meanWithinDistances.putIfAbsent(treeGroup.getGroupName(), sumDist / pairDistances.size());
      }
      /*  Test print
      System.out.println("--------------------mean within distances-------------------------");
      meanWithinDistances.forEach((key, value) -> System.out.println(key + ":" + value));
      System.out.println("----------------------------------------------------------");

      */
    }

      // Make the between distances claculations
      Map<String, Double> fullDistanceBetween = new HashMap<>();
      Map<String, Integer> countsBetweenComparisons = new HashMap<>();

      // Fetch list of child nodes for comparison
      ArrayList<TreeNode> childNodes = new ArrayList<>();
      for (TreeNode node : this.getNodes()) {
        if (node.isLeafNode()) {
          childNodes.add(node);
        }
      }

      for (int i = 0; i <= (childNodes.size()); i++) {
        for (int j = i + 1; j <= (childNodes.size() - 1); j++) {

          TreeNode node1 = childNodes.get(i);
          TreeNode node2 = childNodes.get(j);

          String superpop1 = node1.getName().split("___")[0];
          String superpop2 = node2.getName().split("___")[0];

          String subpop1 = node1.getName().split("___")[1];
          String subpop2 = node2.getName().split("___")[1];

          double pairDistance = this.findPairDistance(node1.getName(), node2.getName());
          System.out.println("Pair distance for nodes: " + node1.getName() + " and " + node2.getName());
          System.out.println("Pair distance: " + pairDistance);

          if (!subpop1.equals(subpop2)) {      // Superpop will also be equal
            if (superpop1.equals(superpop2)) {    // Superpops equal but not subpops

              for (String pop : new String[]{subpop1, subpop2}) {     // Must add to both subpops

                fullDistanceBetween.computeIfPresent(pop, (key, val) -> val + pairDistance);  // sjekk
                countsBetweenComparisons.computeIfPresent(pop, (key, val) -> val + 1);
                fullDistanceBetween.putIfAbsent(pop, pairDistance);
                countsBetweenComparisons.putIfAbsent(pop, 1);

              }

            } else {      // All different, must add to between-map for all

              for (String pop : new String[]{subpop1, subpop2, superpop1, superpop2}) {
                fullDistanceBetween.computeIfPresent(pop, (key, val) -> val + pairDistance);
                countsBetweenComparisons.computeIfPresent(pop, (key, val) -> val + 1);
                fullDistanceBetween.putIfAbsent(pop, pairDistance);
                countsBetweenComparisons.putIfAbsent(pop, 1);
              }
            }
          }
        }
      }

      // Calculate SDRs if possible:
      if (countsBetweenComparisons.isEmpty()) {
        logger.debug("Cannot calculate SDR for a single cluster");
        System.out.println("Cannot calculate SDR for a single cluster.");   // Remove at some point
      } else {

        HashMap<String, Double> meanBetweenDistances = new HashMap<>();

        // Tested thoroughly
        for (Map.Entry<String, Double> entry : fullDistanceBetween.entrySet()) {
          String popKey = entry.getKey();
          meanBetweenDistances.putIfAbsent(popKey, entry.getValue() / countsBetweenComparisons.get(popKey));
        }

        // Test print total distances
        System.out.println("Print within and between lists check -------------------------");
        System.out.println("Between: " + meanBetweenDistances.toString());
        System.out.println("Within: " + meanWithinDistances.toString());
        System.out.println("--------------------------------------------------------------");


        // Add SDR to group objects
        for (Map.Entry<String, Double> entry : meanWithinDistances.entrySet()) {
          String popKey = entry.getKey();
          if (meanBetweenDistances.containsKey(popKey)) {
            double SDR = entry.getValue() / meanBetweenDistances.get(popKey);
            TreeNodeGroup group = findGroupNodeByName(popKey);
            group.setParameter("SDR", Double.toString(SDR));
          }
        }
      }

  }

  public static void main(String args[]){
      WivTree testTree = new WivTree();
      // All possible combinations of other parameters
      // testTree.setTreeLine("F[&lr={0.19,0.82},a=8.88-16,b=3.14,&cat={0.31,0.42}]:0.9,((B[]:0.2,(C[]:0.3,D[]:0.4),E[]:0.5),A[]:0.1,(C2[]:0.3,D2[]:0.4));");

      String beastTestString = "F[a=8]:0.9,((B[b=1]:0.2,(C[c=2]:0.3,D[d=4]:0.4),E[e=5]:0.5),A[a2=6]:0.1,(C2[c2=7]:0.3,D2[d2=9]:0.4));";
      String bootSimpleString  = "((((chimp:0.0608,bonobo:0.0556):0.103,homo_sapiens:0.204):0.0864,gorilla:0.226)'demoLabel2':0.245,(orangutan:0.119,sumatran:0.099):0.432,gibbon:0.903);";
      String bootTestString = "((raccoon:19.19,bear:6.80)50:0.84,((sea_lion:11.99, seal:12.00)100:7.52,((monkey:100.85,cat:47.14)80:20.59, weasel:18.87)75:2.09)50:3.87,dog:25.46);";
      String recurseString = "((('L1:shortname,long_name':0.01, L2:0.02)0.11:0.1,L3:0.03, (L4:0.04, L5:0.05)0.22:0.2)0.33:0.3,L6:0.06);";

      // String beastTreeFile = "/Users/simonray/Dropbox/DResearch/phylogeny/sample_trees/Primates.MCC.tree";
      // String bayesTreeFile = "/Users/simonray/Dropbox/data/1KGenomes/3utr/pathway/R_HSA_70326/PFKL/ENST00000349048/MAF_0.1/aligned_PFKL_ENST00000349048_MAF_0.1_consensus.nex.con.tre";
      // String ninjaTreeFile = "/Users/simonray/Dropbox/dropData/ninja_trees/ENSG00000137154___RS6_HUMAN_aln.tre";
      String ninjaTreeFile = "C:\\Users\\norab\\MasterDisaster\\javacode\\data\\NINJA3rounded.tre";
      // String popClassificationFile = "/Users/simonray/Dropbox/dropData/1KGenomes/population_reference/population_classes.tsv";

      // Tested with nwk tree.
      // String ninjaTreeFile = "C:\\Users\\norab\\MasterDisaster\\javacode\\data\\gene_tree.nwk";
      // String ninjaTreeFile = "/Users/simonray/Dropbox/DResearch/phylogeny/sample_trees/NINJA.tre";
      String popClassificationFile = "C:\\Users\\norab\\MasterDisaster\\javacode\\data\\phydist_population_classes.tsv";
      testTree.setSampleListFilename(popClassificationFile);

      try{
        //logger.debug("debug");
        //logger.info("info");
        testTree.getSampleClassificationList().readSampleGroupFile(popClassificationFile);
        testTree.readTree(ninjaTreeFile);
        //testTree.setTreeLine(recurseString);
        //testTree.removeWhiteSpaceFromTreeString();
        TreeType currentTreeType = WivTree.findTreeType(testTree.getTreeLine());
        testTree.setTreeType(currentTreeType);

        if(testTree.getTreeType()==TreeType.NEWICK || testTree.getTreeType()==TreeType.NEWICKBOOTSTRAP) {
          String remainingTreeLine=testTree.findNewickParentNode();
          //logger.debug(remainingTreeLine);
          //logger.debug(testTree.getNode(0).printParameters());
          remainingTreeLine = remainingTreeLine.substring(1,  remainingTreeLine.length()-1);
          testTree.parseOutNewickChildren(remainingTreeLine, 0, 0);
          String pairDistancesFile = Paths.get(FilenameUtils.getFullPath(popClassificationFile),
                  FilenameUtils.getBaseName(popClassificationFile) + "_pairdists.tsv").toString();
          testTree.getPairDistancesForClassificationNames(pairDistancesFile);
        } else {
          if(testTree.getTreeType()==TreeType.BEAST) {
            String remainingTreeLine=testTree.findBeastParentNode();
            //logger.debug(remainingTreeLine);
            //logger.debug(testTree.getNode(0).printParameters());
            remainingTreeLine = remainingTreeLine.substring(1,  remainingTreeLine.length()-1);
            testTree.parseOutBeastChildren(remainingTreeLine, 0, 0);
          }
        }
        //logger.debug("\n"+ testTree.printNodes());
        testTree.calculateSDRs();
        //testTree.getSDR();
        System.out.println("finished");
      }catch(IOException exIO){
        System.out.println(exIO.toString());
        System.out.println("This did not work out...");
      }


    /*
      WivTree testTree = new WivTree();
      // All possible combinations of other parameters
      // testTree.setTreeLine("F[&lr={0.19,0.82},a=8.88-16,b=3.14,&cat={0.31,0.42}]:0.9,((B[]:0.2,(C[]:0.3,D[]:0.4),E[]:0.5),A[]:0.1,(C2[]:0.3,D2[]:0.4));");

      String beastTestString = "F[a=8]:0.9,((B[b=1]:0.2,(C[c=2]:0.3,D[d=4]:0.4),E[e=5]:0.5),A[a2=6]:0.1,(C2[c2=7]:0.3,D2[d2=9]:0.4));";
      String bootSimpleString  = "((((chimp:0.0608,bonobo:0.0556):0.103,homo_sapiens:0.204):0.0864,gorilla:0.226)'demoLabel2':0.245,(orangutan:0.119,sumatran:0.099):0.432,gibbon:0.903);";
      String bootTestString = "((raccoon:19.19,bear:6.80)50:0.84,((sea_lion:11.99, seal:12.00)100:7.52,((monkey:100.85,cat:47.14)80:20.59, weasel:18.87)75:2.09)50:3.87,dog:25.46);";
      String recurseString = "((('L1:shortname,long_name':0.01, L2:0.02)0.11:0.1,L3:0.03, (L4:0.04, L5:0.05)0.22:0.2)0.33:0.3,L6:0.06);";

      String beastTreeFile = "/Users/simonray/Dropbox/DResearch/phylogeny/sample_trees/Primates.MCC.tree";
      String bayesTreeFile = "/Users/simonray/Dropbox/data/1KGenomes/3utr/pathway/R_HSA_70326/PFKL/ENST00000349048/MAF_0.1/aligned_PFKL_ENST00000349048_MAF_0.1_consensus.nex.con.tre";
      //String ninjaTreeFile = "/Users/simonray/Dropbox/dropData/ninja_trees/ENSG00000137154___RS6_HUMAN_aln.tre";
      //String popClassificationFile = "/Users/simonray/Dropbox/dropData/1KGenomes/population_reference/population_classes.tsv";
      String ninjaTreeFile = "/Users/simonray/Dropbox/DResearch/phylogeny/sample_trees/mega_bootstrap2.nwk";
      //String ninjaTreeFile = "/Users/simonray/Dropbox/DResearch/phylogeny/sample_trees/NINJA.tre";
      String popClassificationFile = "/Users/simonray/Dropbox/dropData/1KGenomes/population_reference/population_classes.tsv";
      testTree.setSampleListFilename(popClassificationFile);

      try{
        logger.debug("debug");
        logger.info("info");
        testTree.getSampleClassificationList().readSampleGroupFile(popClassificationFile);
        testTree.readTree(ninjaTreeFile);
        //testTree.setTreeLine(recurseString);
        //testTree.removeWhiteSpaceFromTreeString();
        TreeType currentTreeType = WivTree.findTreeType(testTree.getTreeLine());
        testTree.setTreeType(currentTreeType);

	      if(testTree.getTreeType()==TreeType.NEWICK || testTree.getTreeType()==TreeType.NEWICKBOOTSTRAP) {
          String remainingTreeLine=testTree.findNewickParentNode();
  	      logger.debug(remainingTreeLine);
  	      logger.debug(testTree.getNode(0).printParameters());
  	      remainingTreeLine = remainingTreeLine.substring(1,  remainingTreeLine.length()-1);
		      testTree.parseOutNewickChildren(remainingTreeLine, 0, 0);	      	
          String pairDistancesFile = Paths.get(FilenameUtils.getFullPath(popClassificationFile),  
                  FilenameUtils.getBaseName(popClassificationFile) + "_pairdists.tsv").toString();
          testTree.getPairDistancesForClassificationNames(pairDistancesFile);
	      } else {
		      if(testTree.getTreeType()==TreeType.BEAST) {
	          String remainingTreeLine=testTree.findBeastParentNode();
	  	      logger.debug(remainingTreeLine);
	  	      logger.debug(testTree.getNode(0).printParameters());
	  	      remainingTreeLine = remainingTreeLine.substring(1,  remainingTreeLine.length()-1);
			      testTree.parseOutBeastChildren(remainingTreeLine, 0, 0);	      	
		      }	      	
	      }
	      logger.debug("\n"+ testTree.printNodes());
	      System.out.println("finished");
      }catch(IOException exIO){
          System.out.println(exIO.toString());
      }
      */
  }
 }
