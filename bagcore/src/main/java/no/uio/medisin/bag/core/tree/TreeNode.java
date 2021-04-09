package no.uio.medisin.bag.core.tree;

import java.text.ParseException;
import java.util.HashMap;

import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *
 * @author sr
 */
public class TreeNode {

  static Logger				logger = LogManager.getLogger();
  public static final String DEFAULT_LINECOLOUR_KEYSTRING = "default_node_line_colour";
  public static final String DEFAULT_LINEWIDTH_KEYSTRING = "default_node_line_width";
  
  
  private int nodeID;
  private int parentNode;
  private String name;
  private double branchLength;
  private double bootStrap;
  private HashMap<String, String>  otherParams  = new HashMap<>();
  private double x;
  private double y;
  private double xBounds;
  private double yBounds;
  private int leftChild;
  private int rightChild;
  private java.util.ArrayList<Integer> children;
  private String groupName;    //  To be removed


  public TreeNode()
  {
  }

  public TreeNode(int n, TreeType tType)
  {
    nodeID = n;
    leftChild = -1;
    rightChild = -1;
    x = -1.0;
    y = -1.0;
    children = new java.util.ArrayList<Integer> ();
    groupName = null;   //  To be removed

  }

  public TreeNode(int id, String treeString, TreeType tType)
  {
	  this(id, tType);
	  parseTreeString(treeString, tType);
  }


  
  /**
   * parse out the supplied string into TreeNodes
   * If the string contains parenthesis, then there are child nodes to parse
   * In this case, the node needs to create subnodes and they need to parse out the relevant substring
   * 
   * @param treeString
   * @return
   */
  public String parseTreeString(String treeString, TreeType tType) {
	  String returnString = "";
	  
	  if(treeString.startsWith("(") && treeString.endsWith(")")) {
		  treeString = treeString.substring(1,  treeString.length()-1);
		  while(treeString.length()>0) {
			  if(this.isNextNodeLeafNode(treeString)) {
				  String nextNodeString = grabNextLeafNode(treeString);
				  this.parseNewickNodeString(nextNodeString);
				  treeString = treeString.substring(nextNodeString.length(), treeString.length()-1);
			  }else {
				  // internal node - create new node
				  
			  }
		  }
		  
	  }else {
		  logger.error("unmatched parenthesis");
	  }
	  
	  return returnString;
  }
  
  /**
   * grab the first leaf node from the front queryString.
   * This assumes the first node is a leafNode,
   * if there is an internal node, the behaviour will be unpredictable
   * 
   * @param queryTreeString
   * @return
   */
  private String grabNextLeafNode(String queryTreeString) {
	  int nextCommaLoc = queryTreeString.indexOf(",");
	  String currentNodeString = queryTreeString.substring(0, nextCommaLoc);
	  logger.debug("-- nextLeadNodeString = <" + currentNodeString + ">");
	  return currentNodeString;
  }
  
  /**
   * grab the next internal node from the front of queryString.
   * This assumes the next node is an internal node
   * if there is a leaf node first, the behaviour will be unpredictable
   * 
   * internal node is delimited by (), but may contain further nested ()
   * may be branch length and bootstrap values in the format
   *   ()bootstrap:branchlength
   *   
   * @param queryTreeString
   * @return
   */
  private String grabNextInternalNode(String queryTreeString) throws ParseException {
	  /*
	   * As we assume the first node in the string is an internal node,
	   * we expect the string to begin with a "("
	   */
	  if(!queryTreeString.startsWith("(")) {
		  String exceptionString = "Exception parsing queryTreeString\n"
		  		+ "Expected string to start with a '(', but found <\n"
		  		+ queryTreeString.substring(0,1) + ">\n"
		  				+ "Full string is \n" + queryTreeString;
		  throw new ParseException("Exception parsing queryTreeString", 0);
	  }
	  
	  return "";
	  
  }
  
  
  private static String getNextInternalNodeString(String queryString, int startPos) {
	  int depth = 0;
	  
	  //if we find another "(" before we find a ")" then we have to go deeper
	  int start = startPos;
	  int nextleftParenLoc = queryString.indexOf("(", start+1);
	  int nextrightParenLoc = queryString.indexOf(")", start+1);
	  Boolean match = false;
	  while(match==false) {
		  logger.debug("Depth<" +depth + ">");
		  if(nextleftParenLoc>=0)
			  logger.debug(" LEFT:[" + nextleftParenLoc + "\t" + queryString.substring(nextleftParenLoc));
		  else
			  logger.debug(" LEFT:[" + nextleftParenLoc + "\t" + "NO MORE MATCHES");
			  
		  if(nextrightParenLoc>=0)
			  logger.debug("RIGHT:[" + nextrightParenLoc + "\t" + queryString.substring(nextrightParenLoc));
		  else
			  logger.debug("RIGHT:[" + nextrightParenLoc + "\t" + "NO MORE MATCHES");
			  
		  if(nextleftParenLoc<nextrightParenLoc && nextleftParenLoc>=0) {
			  depth = depth +1;
			  nextleftParenLoc ++;
			  nextleftParenLoc = queryString.indexOf("(", nextleftParenLoc);
			  logger.debug(nextleftParenLoc + "|" + nextrightParenLoc);
		  }else {
			  if(depth==0) {
				  match = true;
			  }else {
				  depth--;
				  nextrightParenLoc++;
				  nextrightParenLoc = queryString.indexOf(")", nextrightParenLoc);
			  }
		  }

	  }
	  if(nextleftParenLoc==-1) 
		  nextleftParenLoc=0;
	  if(nextrightParenLoc==-1)
		  nextrightParenLoc=queryString.length();
	  return queryString.substring(startPos+1, nextrightParenLoc);
	  
	  
  }
  
  /**
   * there are only two options, 
   *   1. internal node delimited by () + optional BranchLength and BootstrapValue
   *   2. leaf node name + optional BranchLength
   *   
   * if next "," comes before the next "(" next node is leaf node
   * 
   * ---> For now, we only consider NEWICK format, i.e., not BEAST
   * 
   * @param treeStringFrag
   * @return
   */
  private Boolean isNextNodeLeafNode(String treeStringFrag) {
	  int nextCommaLoc = treeStringFrag.indexOf(",");
	  int nextRightParLoc = treeStringFrag.indexOf(")");
	  if(nextCommaLoc<nextRightParLoc) {
		  return true;
	  }
	  return false;	  
  }
  
  
  //private parseLeafNode()
  
  
  
  /**
   * parse out the information contained in the supplied string
   * @param nodeString
   * 
   * Depending on the tree type, different information may be supplied
   * we have :
   * 
   * 	Case 1: NodeName:[additional parameters]:distance
   * 	Case 2: NodeName:param1,param2  (usually Branch Length, Bootstrap) but not necessarily in that order
   * 	Case 3: Something simple (need to figure out what)
   *    Case 4: throw Exception
   * 
   */
  public void parseNewickNodeString(String nodeString) {
	  logger.debug("--<parseNodeString>");
	  
	  /* Case 1:
	     it is possible someone uses this combination in a Sample name, but we don't attempt to handle this
	   */
	  if(nodeString.contains("[")){		  
		  parseBeastParameters(nodeString);
	  }else {
		  parseNewickString(nodeString);
	  }
	  logger.debug("done");
  }
  
  /**
   * parse out parameters from nodeString in Newick Format
   * the information contained in the string will depend on whether it is a
   * leaf node or internal node
   * An internal node can possibly have two numerical parameters '(...)num1:num2' 
   * A leaf node can only have one numerical parameter 'leafnode_name:num1'
   * 
   * @param nodeString
   */
  private void parseNewickString(String nodeString){
    
    // do we have a leaf node?
    try {
      Double.parseDouble(nodeString.split(":")[0].trim());
      if(nodeString.substring(0, 1).equals(":")){
        this.setBranchLength(Double.parseDouble(nodeString.split(":")[1].trim()));		  
      }else{
        try {
          double d = Double.parseDouble(nodeString.split(":")[0].trim());
            this.setBootStrap(d);
        } catch (NumberFormatException nfe) {
            this.setName(nodeString.split(":")[0].trim());
        }
        this.setBranchLength(Double.parseDouble(nodeString.split(":")[1].trim()));
      }
    } catch (NumberFormatException nfe) {
      logger.debug("leaf node");
      this.setName(nodeString.split(":")[0].trim());
      this.setBranchLength(Double.parseDouble(nodeString.split(":")[1].trim()));
    }
    

	  logger.debug(nodeString);
	  logger.debug(this.getName() + "|" + this.getBranchLength() + "|" + this.getBootStrap());
	  logger.debug("--<finished parsing NodeString>");
  }
  
  
  
  /**
   * parse out the extra parameters associated with the node (used by BEAST)
   * These will be in the format []:number
   * and within the [], the parameters can be in the format
   * 'param1=number' or  'param2_range={rangeLow, rangeHigh}'
   * 
   * 
   * @param extraParamString
   */
  public void parseBeastParameters(String extraParamString) {

	  logger.debug("--<parseBeastParameters>");
	  logger.debug(extraParamString);

	  /*
	   * check for node name
	   */
	  if(!extraParamString.startsWith("[")) {
	  	this.setName(extraParamString.split("\\[")[0]);
	  	extraParamString = extraParamString.split("\\[")[1];
	  }
	  
	  
	  /**
	   * first step, is there a number after the [] ?
	   */

	  if(extraParamString.contains(":")) {
	  	String numberString = extraParamString.split(":")[1].trim();
	  	this.setBranchLength(Double.parseDouble(numberString));
	  }
	  extraParamString = extraParamString.split(":")[0].trim();
	  
	  if(extraParamString.startsWith("[")) {
	  	extraParamString=extraParamString.substring(1);
	  }
	  if(extraParamString.endsWith("]")) {
	  	extraParamString=StringUtils.chop(extraParamString);
	  }
	  
	  String bits[] = extraParamString.split(",");
	  
	  int b=0;
	  while(b<bits.length) {
	  	String pair = "";
	  	String paramString = "";
	  	String paramValue = "";
	  	if(bits[b].contains("{")) {
	  		pair = bits[b].concat(bits[b+1]);
	  		b+=2;	  		
	  	}else {
	  		pair = bits[b].trim();
	  		b++;
	  	}
  		paramString=pair.split("=")[0].trim();
  		paramValue=pair.split("=")[1].trim();
  		logger.debug("--|" + Integer.toString(b) + "|<" + paramString + "," + paramValue + ">");
  		otherParams.put(paramString, paramValue);
	  }
	  logger.debug("--<finished>");

	  
  }
  
  /**
   * print out details of the node, but not the parameter information
   * @return
   */
  public String printNodeShort() {
	  String nodeString = "NODE <" + this.getNodeID()+ ">\n"
			  + "\tPARENT_ID:\t" + this.getParentNode() + "\n"
			  + "\tNAME:\t" + this.getName() + "\n"
			  + "\tBRANCHLEN:\t" + this.getBranchLength() + "\n"
			  + "\tCHILDREN:\t" + "\n";
	  for(Integer nodeID: this.getChildren()) {
		  nodeString = nodeString.concat("\t--\t<" + nodeID + ">\n");
	  }
	  
	  return nodeString;
			  
			  
	  
  }
  
  
  /**
   * print other parameters in pretty format
   * @return other parameters as a single string
   */
  public String printParameters() {
	  String  paramString  = "";
		for (HashMap.Entry<String, String> entry : otherParams.entrySet()) {
			String k = entry.getKey();
			String v = entry.getValue();
			paramString = paramString.concat("[" + k + "," + v + "]\n");
		}
	  return paramString;
  }
  
  
  /**
   * is this a Leaf Node
   * @return 
   */
  public Boolean isLeafNode(){
    return getChildren().isEmpty();
  }
  
  /**
   * is this an internal Node
   * @return 
   */
  public Boolean isInternalNode(){
    return !getChildren().isEmpty();
  }
  
  
  public void setParameter(String key, String value){
    otherParams.put(key, value);
  }
  
  public String getParameterValue(String key){
    return this.otherParams.get(key);
  }

  
  
  /**
   * @return the nodeID
   */
  public int getNodeID() {
    return nodeID;
  }

  /**
   * @param nodeID the nodeID to set
   */
  public void setNodeID(int nodeID) {
    this.nodeID = nodeID;
  }

  /**
   * @return the name
   */
  public String getName() {
    return name;
  }

  /**
   * @param name the name to set
   */
  public void setName(String name) {
    this.name = name;
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
   * @return the parentNode
   */
  public int getParentNode() {
    return parentNode;
  }

  /**
   * @param parentNode the parentNode to set
   */
  public void setParentNode(int parentNode) {
    this.parentNode = parentNode;
  }

  /**
   * @return the leftChild
   */
  public int getLeftChild() {
    return leftChild;
  }

  /**
   * @param leftChild the leftChild to set
   */
  public void setLeftChild(int leftChild) {
    this.leftChild = leftChild;
  }

  /**
   * @return the rightChild
   */
  public int getRightChild() {
    return rightChild;
  }

  /**
   * @param rightChild the rightChild to set
   */
  public void setRightChild(int rightChild) {
    this.rightChild = rightChild;
  }

  /**
   * @return the children
   */
  public java.util.ArrayList<Integer> getChildren() {
    return children;
  }

  /**
   * @return the x
   */
  public double getX() {
    return x;
  }

  /**
   * @return the y
   */
  public double getY() {
    return y;
  }

  /**
   * @param x the x to set
   */
  public void setX(double x) {
    this.x = x;
  }

  /**
   * @param y the y to set
   */
  public void setY(double y) {
    this.y = y;
  }

  /**
   * @return the xBounds
   */
  public double getxBounds() {
    return xBounds;
  }

  /**
   * @param xBounds the xBounds to set
   */
  public void setxBounds(double xBounds) {
    this.xBounds = xBounds;
  }

  /**
   * @return the yBounds
   */
  public double getyBounds() {
    return yBounds;
  }

  /**
   * @param yBounds the yBounds to set
   */
  public void setyBounds(double yBounds) {
    this.yBounds = yBounds;
  }

  /**
   * @return the bootStrap
   */
  public double getBootStrap() {
    return bootStrap;
  }

  /**
   * @param bootStrap the bootStrap to set
   */
  public void setBootStrap(double bootStrap) {
    this.bootStrap = bootStrap;
  }

  /**
   * @return the groupName
   */
  public String getGroupName() { return groupName; }   // To be removed

  /**
   * @param groupName the groupName to set
   */
  public void setGroupName(String groupName) { this.groupName = groupName; }  // To be removed
  
  public static void main(String args[]) {
	  String treeString = "(((L1, L2),L3, (L4, L5)),L6);";
	  String intNodeString = "";
	  logger.debug(treeString);
	  logger.debug("start at <0>");
	  intNodeString = TreeNode.getNextInternalNodeString(treeString, 0);
	  logger.debug("<" + intNodeString + ">\n\n");
	  logger.debug("start at <1>");
	  intNodeString = TreeNode.getNextInternalNodeString(treeString, 1);
	  logger.debug("<" + intNodeString + ">");
	  logger.debug("start at <2>");
	  intNodeString = TreeNode.getNextInternalNodeString(treeString, 2);
	  logger.debug("<" + intNodeString + ">");
	  logger.debug("start at <15>");
	  intNodeString = TreeNode.getNextInternalNodeString(treeString, 15);
	  logger.debug("<" + intNodeString + ">");
	  logger.debug("finished");
  }
  
}
