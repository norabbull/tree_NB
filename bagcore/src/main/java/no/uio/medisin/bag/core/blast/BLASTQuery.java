package no.uio.medisin.bag.core.blast;

/**
 * Stores information for a single query (as specified within the <Iteration></Iteration> tags in the output XML
 * 
 * Not sure what the full set of identifies are, but found the following in a sample XML file:
 * 
 *   <Iteration_iter-num>1</Iteration_iter-num>
 *   <Iteration_query-ID>Query_1</Iteration_query-ID>
 *   <Iteration_query-def>hcmv-miR-UL22A-5p MIMAT0001574 Human cytomegalovirus miR-UL22A-5p</Iteration_query-def>
 *   <Iteration_query-len>20</Iteration_query-len>
 *
 * 
 * @author simonray
 *
 */
public class BLASTQuery {
	/*
 */
	private int			queryNum;
	private String	queryID;
	private String	queryDef;
	private int			queryLength;
	
	
	
	
	public int getQueryNum() {
		return queryNum;
	}
	public void setQueryNum(int queryNum) {
		this.queryNum = queryNum;
	}
	public String getQueryID() {
		return queryID;
	}
	public void setQueryID(String queryID) {
		this.queryID = queryID;
	}
	public String getQueryDef() {
		return queryDef;
	}
	public void setQueryDef(String queryDef) {
		this.queryDef = queryDef;
	}
	public int getQueryLength() {
		return queryLength;
	}
	public void setQueryLength(int queryLength) {
		this.queryLength = queryLength;
	}
}
