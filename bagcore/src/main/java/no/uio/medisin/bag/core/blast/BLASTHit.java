package no.uio.medisin.bag.core.blast;

import java.util.ArrayList;

/**
 * 
 * Stores information for a single query (as specified within the <Hit></Hit> tags in the output XML
 * 
 * Not sure what the full set of identifies are, but found the following in a sample XML file:
 * 
 * <Hit>
 *   <Hit_num>1</Hit_num>
 *   <Hit_id>gnl|BL_ORD_ID|1</Hit_id>
 *   <Hit_def>gi|155573622|ref|NC_006273.2| Human herpesvirus 5 strain Merlin, complete genome</Hit_def>
 *   <Hit_accession>1</Hit_accession>
 *   <Hit_len>235646</Hit_len>
 *   <Hit_hsps>
 *     <Hsp></Hsp>
 *     <Hsp></Hsp>
 *   </Hit_hsps>
 * </Hit>
 * 
 * @author simonray
 *
 */
public class BLASTHit {
	private int 														hitNum;
	private String 													hitID;
	private String 													hitDef;
	private String 													accNum;
	private int 														hitLength;
	private ArrayList<BLASTHighScoringPair> hspList;
	
	
	public BLASTHit() {
		hspList = new ArrayList<>();
	}
	
	
	/**
	 * add the specified high scoring pair to the hit list
	 * 
	 * @param newHSP
	 * @return
	 */
	public int addHighScoringPair(BLASTHighScoringPair  newHSP)
	{
		hspList.add(newHSP);
		return hspList.size();
	}
	
	
	public int getHitNum() {
		return hitNum;
	}
	public void setHitNum(int hitNum) {
		this.hitNum = hitNum;
	}
	public String getHitDef() {
		return hitDef;
	}
	public void setHitDef(String hitDef) {
		this.hitDef = hitDef;
	}
	public String getAccNum() {
		return accNum;
	}
	public void setAccNum(String accNum) {
		this.accNum = accNum;
	}
	public int getHitLength() {
		return hitLength;
	}
	public void setHitLength(int hitLength) {
		this.hitLength = hitLength;
	}
	public ArrayList<BLASTHighScoringPair> getHspList() {
		return hspList;
	}


	public String getHitID() {
		return hitID;
	}


	public void setHitID(String hitID) {
		this.hitID = hitID;
	}
	
}
