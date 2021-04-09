package no.uio.medisin.bag.core.blast;

import java.io.IOException;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import no.uio.medisin.bag.core.mirna.PreMiRNASet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.w3c.dom.traversal.DocumentTraversal;
import org.w3c.dom.traversal.TreeWalker;
import org.xml.sax.SAXException;

import no.uio.medisin.bag.exceptions.BLASTException;

/**
 * 
 * @author simon rayner
 * A lightweight parser to extract elements of interest from a BLAST XML result file
 * This for small BLAST result files, which is why I use DOM rather than SAX
 * 
 * The basic structure of the XML appears to be:
 * 
 * <BlastOutput>
 *   <BlastOutput_iterations>
 *     <Iteration>
 *       <Iteration_iter-num>1</Iteration_iter-num>
 *       <Iteration_query-len>47</Iteration_query-len>
 *       <Iteration_hits>
 *         <Hit>
 *           <Hit_num>1</Hit_num>
 *           <Hit_hsps>
 *             <Hsp>
 *             </Hsp>
 *             <Hsp>
 *             </Hsp>
 *           </Hit_hsps>
 *        </Hit>
 *      </Iteration_hits>
 *      <Iteration_stat>
 *        <Statistics>
 *          <Statistics_entropy>0.85</Statistics_entropy>
 *        </Statistics>
 *      </Iteration_stat>
 *    </Iteration>
 *   </BlastOutput_iterations>
 * </BlastOutput>
 *
 * 
 *
 */
public class BLASTXMLParser {
	
	static Logger logger = LogManager.getLogger();
	DocumentBuilderFactory dbFactory;
	DocumentBuilder dBuilder;
	Document docBLAST;
	DocumentTraversal traversal;
  private String  xmlFile;
	
	private BLASTRunParameters blastParams;
	private BLASTQuery 			blastQuery;
  private ArrayList<BLASTHit> blastHitList;
		
	
	public BLASTXMLParser() {
		blastParams 	= new BLASTRunParameters();
		blastQuery 		= new BLASTQuery();
    blastHitList  = new ArrayList<BLASTHit>();
	}
  
	/**
	 * initialize the DOM classes
   * @throws ParserConfigurationException
	 */
	public void initialize() throws ParserConfigurationException{
		try {
			dbFactory = DocumentBuilderFactory.newInstance();
			dBuilder = dbFactory.newDocumentBuilder();

		}catch(ParserConfigurationException exPs) {
      logger.error("error initializing DocumentBuilderFactory");
      throw new ParserConfigurationException("error initializing DocumentBuilderFactory");
		}
	}
  
  public void clearData(){
    blastHitList.clear();
  }
	
	public static void main( String args[]) {

    String xmlFile = "/Users/simonray/Dropbox/data/testdata/seqExtract/MI0022143_vs_ssc.XML";
    BLASTXMLParser blastXMLParser = new BLASTXMLParser();
    try {
    	blastXMLParser.initialize();
    	blastXMLParser.loadXML(xmlFile);
    	blastXMLParser.parseOutBLASTParameters();
    	blastXMLParser.parseAllHits();
    	
    	
    }catch(Exception exEx) {
    	exEx.printStackTrace();
    }

	}	
	
	/**
	 * load the XML
	 * @param xmlFile
	 * @throws ParserConfigurationException (if parameters aren't found)
	 * 
	 */
	public void loadXML(String xmlFile) throws ParserConfigurationException, SAXException {
		try {
		  docBLAST = dBuilder.parse(xmlFile);
		  docBLAST.getDocumentElement().normalize();			
      DocumentTraversal traversal = (DocumentTraversal) docBLAST;
		}catch(IOException exIO) {
      logger.error("IO error trying to parse XML file <" + xmlFile + ">");
      throw new ParserConfigurationException("IO error trying to parse XML file <" + xmlFile + ">");
		}catch(SAXException exSx) {
      logger.error("SAXerror trying to parse XML file <" + xmlFile + ">");
      throw new SAXException("SAXerror trying to parse XML file <" + xmlFile + ">");
		}
		
	}

	/**
	 * parse out BLAST parameter values stored within the <BlastOutput_param></BlastOutput_param> keys
	 * 
	 * Not sure what the full set is, but this is what we get in default XML output
	 * 	    <Parameters>
	 * 	      <Parameters_expect>10</Parameters_expect>
	 * 	      <Parameters_sc-match>1</Parameters_sc-match>
	 * 	      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
	 * 	      <Parameters_gap-open>0</Parameters_gap-open>
	 * 	      <Parameters_gap-extend>0</Parameters_gap-extend>
	 * 	      <Parameters_filter>L;m;</Parameters_filter>
	 * 	    </Parameters>
	 */
	public void parseOutBLASTParameters() throws BLASTException{
		NodeList nParamList = docBLAST.getElementsByTagName("Parameters");		
		Node pNode = nParamList.item(0);
    	
    System.out.println("\nCurrent Element :" + pNode.getNodeName());

  	if (pNode.getNodeType() == Node.ELEMENT_NODE) {
  		Element eElement = (Element) pNode;
  		getBlastParams().setExpectationValue(Double.parseDouble(eElement.getElementsByTagName("Parameters_expect").item(0).getTextContent()));
  		getBlastParams().setScoreMatch(Integer.parseInt(eElement.getElementsByTagName("Parameters_sc-match").item(0).getTextContent()));
  		getBlastParams().setScoreMismatch(Integer.parseInt(eElement.getElementsByTagName("Parameters_sc-mismatch").item(0).getTextContent()));
  		getBlastParams().setGapOpen(Integer.parseInt(eElement.getElementsByTagName("Parameters_gap-open").item(0).getTextContent()));
  		getBlastParams().setGapExtend(Integer.parseInt(eElement.getElementsByTagName("Parameters_gap-extend").item(0).getTextContent()));
  		getBlastParams().setFilterString(eElement.getElementsByTagName("Parameters_filter").item(0).getTextContent());
    }else {
    	logger.error("didn't find any BLAST parameters in the XML file");
    	throw new BLASTException("didn't find any BLAST parameters in the XML file");
    }
	}
	
	
	/**
	 * parse out the hits that are associated with a single query
	 * This will generally involve more than one reference sequence (i.e., multiple chromosomes or viral strains)
	 * 
	 * @return
	 */
	public int parseAllHits() {
		
		Node nQueryInfo = docBLAST.getElementsByTagName("Iteration").item(0);
  	if (nQueryInfo.getNodeType() == Node.ELEMENT_NODE) {
  		logger.debug("--" + nQueryInfo.getNodeName());
  		Element eQueryInfo = (Element) nQueryInfo;
  		logger.debug(eQueryInfo.getElementsByTagName("Iteration_iter-num").item(0).getTextContent());
  		logger.debug(eQueryInfo.getElementsByTagName("Iteration_query-ID").item(0).getTextContent());
  		logger.debug(eQueryInfo.getElementsByTagName("Iteration_query-def").item(0).getTextContent());
  		logger.debug(eQueryInfo.getElementsByTagName("Iteration_query-len").item(0).getTextContent());
  		getBlastQuery().setQueryNum(Integer.parseInt(eQueryInfo.getElementsByTagName("Iteration_iter-num").item(0).getTextContent()));
  		getBlastQuery().setQueryID(eQueryInfo.getElementsByTagName("Iteration_query-ID").item(0).getTextContent());
  		getBlastQuery().setQueryDef(eQueryInfo.getElementsByTagName("Iteration_query-def").item(0).getTextContent());
  		getBlastQuery().setQueryLength(Integer.parseInt(eQueryInfo.getElementsByTagName("Iteration_query-len").item(0).getTextContent()));
  	}
  	/*
  	 * then cycle through the hits
  	 */
		NodeList itHitsList = docBLAST.getElementsByTagName("Iteration_hits");	
		for(int it = 0; it < itHitsList.getLength(); it++) { 
			NodeList hitList = ((Element)itHitsList.item(it)).getElementsByTagName("Hit");
      if(hitList.item(0)==null)
        return -1;
			NodeList hitListChildren = hitList.item(0).getChildNodes();
			/**
			 * we divide by 6 because there are 6 child nodes/hit
			 */
			for(int h = 0; h < hitListChildren.getLength()/6; h++) { 
				Node hit = hitList.item(h);
        if(hit == null)
          break;
  			if( hit.getNodeType() == Node.ELEMENT_NODE) {
          BLASTHit				blastHit = new BLASTHit();
  				/* Sample Entry
  		    <Hit_num>1</Hit_num>
  			  <Hit_id>gnl|BL_ORD_ID|1</Hit_id>
  			  <Hit_def>gi|155573622|ref|NC_006273.2| Human herpesvirus 5 strain Merlin, complete genome</Hit_def>
  			  <Hit_accession>1</Hit_accession>
  			  <Hit_len>235646</Hit_len>
  				*/
		  		Element eHitInfo = (Element) hit;
		  		logger.debug("\n-- Hit_num --> " + eHitInfo.getElementsByTagName("Hit_num").item(0).getTextContent());  		
		  		logger.debug("-- Hit_id  --> " + eHitInfo.getElementsByTagName("Hit_id").item(0).getTextContent());  		
		  		logger.debug(eHitInfo.getElementsByTagName("Hit_def").item(0).getTextContent());  		
		  		logger.debug(eHitInfo.getElementsByTagName("Hit_accession").item(0).getTextContent());  		
		  		logger.debug(eHitInfo.getElementsByTagName("Hit_len").item(0).getTextContent());  		
		  		blastHit.setHitNum(Integer.parseInt(eHitInfo.getElementsByTagName("Hit_num").item(0).getTextContent()));
		  		blastHit.setHitID(eHitInfo.getElementsByTagName("Hit_id").item(0).getTextContent());
		  		blastHit.setHitDef(eHitInfo.getElementsByTagName("Hit_def").item(0).getTextContent());
		  		blastHit.setAccNum(eHitInfo.getElementsByTagName("Hit_accession").item(0).getTextContent());
		  		blastHit.setHitLength(Integer.parseInt(eHitInfo.getElementsByTagName("Hit_len").item(0).getTextContent()));
		  		

          NodeList hit_hspsList = ((Element)hit).getElementsByTagName("Hit_hsps");
          NodeList hit_hspsListChildren = hit_hspsList.item(0).getChildNodes();
          for(int iHsp = 0; iHsp < hit_hspsListChildren.getLength(); iHsp++) { 
            Node current = hit_hspsListChildren.item(iHsp);
            if(current.getNodeType() == Node.ELEMENT_NODE) {
              BLASTHighScoringPair bHSP = new BLASTHighScoringPair();
              Element eCurrent = (Element) current;
              logger.debug("\n---- Hsp-num ----> " + eCurrent.getElementsByTagName("Hsp_num").item(0).getTextContent());
              logger.debug(eCurrent.getElementsByTagName("Hsp_bit-score").item(0).getTextContent());
              logger.debug(eCurrent.getElementsByTagName("Hsp_hit-from").item(0).getTextContent());
              logger.debug(eCurrent.getElementsByTagName("Hsp_hit-to").item(0).getTextContent());
              bHSP.setID(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_num").item(0).getTextContent()));
              bHSP.setBitScore(Double.parseDouble(eCurrent.getElementsByTagName("Hsp_bit-score").item(0).getTextContent()));
              bHSP.setScore(Double.parseDouble(eCurrent.getElementsByTagName("Hsp_score").item(0).getTextContent()));
              bHSP.seteValue(Double.parseDouble(eCurrent.getElementsByTagName("Hsp_evalue").item(0).getTextContent()));
              bHSP.setQueryFrom(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_query-from").item(0).getTextContent()));
              bHSP.setQueryTo(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_query-to").item(0).getTextContent()));
              bHSP.setQueryFrame(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_query-frame").item(0).getTextContent()));
              bHSP.setHitFrom(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_hit-from").item(0).getTextContent()));
              bHSP.setHitTo(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_hit-to").item(0).getTextContent()));
              bHSP.setHitFrame(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_hit-frame").item(0).getTextContent()));
              bHSP.setIdentity(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_identity").item(0).getTextContent()));
              bHSP.setPositive(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_positive").item(0).getTextContent()));
              bHSP.setGaps(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_gaps").item(0).getTextContent()));
              bHSP.setAlignLength(Integer.parseInt(eCurrent.getElementsByTagName("Hsp_align-len").item(0).getTextContent()));
              bHSP.setQueryString(eCurrent.getElementsByTagName("Hsp_qseq").item(0).getTextContent());
              bHSP.setHitString(eCurrent.getElementsByTagName("Hsp_hseq").item(0).getTextContent());
              bHSP.setMidLine(eCurrent.getElementsByTagName("Hsp_midline").item(0).getTextContent());
              blastHit.addHighScoringPair(bHSP);

            }
          
          }
          getBlastHitList().add(blastHit);
		  		/*
			      <Hsp_num>1</Hsp_num>
			      <Hsp_bit-score>40.14</Hsp_bit-score>
			      <Hsp_score>20</Hsp_score>
			      <Hsp_evalue>3.45415e-06</Hsp_evalue>
			      <Hsp_query-from>1</Hsp_query-from>
			      <Hsp_query-to>20</Hsp_query-to>
			      <Hsp_hit-from>27992</Hsp_hit-from>
			      <Hsp_hit-to>28011</Hsp_hit-to>
			      <Hsp_query-frame>1</Hsp_query-frame>
			      <Hsp_hit-frame>1</Hsp_hit-frame>
			      <Hsp_identity>20</Hsp_identity>
			      <Hsp_positive>20</Hsp_positive>
			      <Hsp_gaps>0</Hsp_gaps>
			      <Hsp_align-len>20</Hsp_align-len>
			      <Hsp_qseq>TAACTAGCCTTCCCGTGAGA</Hsp_qseq>
			      <Hsp_hseq>TAACTAGCCTTCCCGTGAGA</Hsp_hseq>
			      <Hsp_midline>||||||||||||||||||||</Hsp_midline>
	  			*/
	  			
	  		}
        logger.debug("finished parsing HSPs");
			}
		}				
		return 0;
	}
  
  
	
	private static void traverseLevel(TreeWalker walker, String indent, int level) {

		Node node = walker.getCurrentNode();
    logger.info("Node[" + level + "] " + indent + node.getNodeName());
		

		if (node.getNodeType() == Node.ELEMENT_NODE) {
      logger.info("Element Node[" + level + "] " + indent + node.getNodeName());
		}

		if (node.getNodeType() == Node.TEXT_NODE) {
      String content_trimmed = node.getTextContent().trim();
      if (content_trimmed.length() > 0) {
      	logger.info("Content Text" + indent + "%s%n", content_trimmed);
      }
		}
		level++;
		for (Node n = walker.firstChild(); n != null; n = walker.nextSibling()) {
      traverseLevel(walker, indent + "  ", level);
		}

		walker.setCurrentNode(node);
	}

  /**
   * @return the xmlFile
   */
  public String getXmlFile() {
    return xmlFile;
  }

  /**
   * @param xmlFile the xmlFile to set
   */
  public void setXmlFile(String xmlFile) {
    this.xmlFile = xmlFile;
  }

  /**
   * @return the blastHitList
   */
  public ArrayList<BLASTHit> getBlastHitList() {
    return blastHitList;
  }

  /**
   * @param blastHitList the blastHitList to set
   */
  public void setBlastHitList(ArrayList<BLASTHit> blastHitList) {
    this.blastHitList = blastHitList;
  }

  /**
   * @return the blastParams
   */
  public BLASTRunParameters getBlastParams() {
    return blastParams;
  }

  /**
   * @param blastParams the blastParams to set
   */
  public void setBlastParams(BLASTRunParameters blastParams) {
    this.blastParams = blastParams;
  }

  /**
   * @return the blastQuery
   */
  public BLASTQuery getBlastQuery() {
    return blastQuery;
  }

  /**
   * @param blastQuery the blastQuery to set
   */
  public void setBlastQuery(BLASTQuery blastQuery) {
    this.blastQuery = blastQuery;
  }


  
}
