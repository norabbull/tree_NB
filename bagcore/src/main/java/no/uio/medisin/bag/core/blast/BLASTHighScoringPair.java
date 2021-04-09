package no.uio.medisin.bag.core.blast;

/**
 * contains all the information about a single BLAST hit 
 * @author simon rayner
 *
 *		  		
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
public class BLASTHighScoringPair {
	private int 		ID;
	private double	bitScore;
	private double	score;
	private double	eValue;
	private int			queryFrom;
	private int			queryTo;
	private int			hitFrom;
	private int			hitTo;
	private int			queryFrame;
	private int			hitFrame;
	private int			identity;
	private int			positive;
	private int			gaps;
	private int			alignLength;
	private String	queryString;
	private	String	hitString;
	private String	midLine;
	
	
	
	public int getID() {
		return ID;
	}
	public void setID(int iD) {
		ID = iD;
	}
	public double getBitScore() {
		return bitScore;
	}
	public void setBitScore(double bitScore) {
		this.bitScore = bitScore;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public double geteValue() {
		return eValue;
	}
	public void seteValue(double eValue) {
		this.eValue = eValue;
	}
	public int getQueryFrom() {
		return queryFrom;
	}
	public void setQueryFrom(int queryFrom) {
		this.queryFrom = queryFrom;
	}
	public int getQueryTo() {
		return queryTo;
	}
	public void setQueryTo(int queryTo) {
		this.queryTo = queryTo;
	}
	public int getHitFrom() {
		return hitFrom;
	}
	public void setHitFrom(int hitFrom) {
		this.hitFrom = hitFrom;
	}
	public int getHitTo() {
		return hitTo;
	}
	public void setHitTo(int hitTo) {
		this.hitTo = hitTo;
	}
	public int getQueryFrame() {
		return queryFrame;
	}
	public void setQueryFrame(int queryFrame) {
		this.queryFrame = queryFrame;
	}
	public int getHitFrame() {
		return hitFrame;
	}
	public void setHitFrame(int hitFrame) {
		this.hitFrame = hitFrame;
	}
	public int getIdentity() {
		return identity;
	}
	public void setIdentity(int identity) {
		this.identity = identity;
	}
	public int getPositive() {
		return positive;
	}
	public void setPositive(int positive) {
		this.positive = positive;
	}
	public int getGaps() {
		return gaps;
	}
	public void setGaps(int gaps) {
		this.gaps = gaps;
	}
	public int getAlignLength() {
		return alignLength;
	}
	public void setAlignLength(int alignLength) {
		this.alignLength = alignLength;
	}
	public String getQueryString() {
		return queryString;
	}
	public void setQueryString(String queryString) {
		this.queryString = queryString;
	}
	public String getHitString() {
		return hitString;
	}
	public void setHitString(String hitString) {
		this.hitString = hitString;
	}
	public String getMidLine() {
		return midLine;
	}
	public void setMidLine(String midLine) {
		this.midLine = midLine;
	}
	
	
}
