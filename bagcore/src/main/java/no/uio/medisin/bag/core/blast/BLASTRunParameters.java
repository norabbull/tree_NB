package no.uio.medisin.bag.core.blast;

/**
 * the run files that were used in a completed run
 * @author simonray
 *
 */
public class BLASTRunParameters {
	private String 	blastProgram;
	private String	blastVersion; 
	private double	expectationValue 	;
	private int			scoreMatch				;
	private	int			scoreMismatch 		;
	private int			gapOpen						;
	private int			gapExtend					;
	private String	filterString			;
	
	
	
	
	
	
	
	public String getBlastProgram() {
		return blastProgram;
	}
	public void setBlastProgram(String blastProgram) {
		this.blastProgram = blastProgram;
	}
	public String getBlastVersion() {
		return blastVersion;
	}
	public void setBlastVersion(String blastVersion) {
		this.blastVersion = blastVersion;
	}
	public double getExpectationValue() {
		return expectationValue;
	}
	public void setExpectationValue(double expectationValue) {
		this.expectationValue = expectationValue;
	}
	public int getScoreMatch() {
		return scoreMatch;
	}
	public void setScoreMatch(int scoreMatch) {
		this.scoreMatch = scoreMatch;
	}
	public int getScoreMismatch() {
		return scoreMismatch;
	}
	public void setScoreMismatch(int scoreMismatch) {
		this.scoreMismatch = scoreMismatch;
	}
	public int getGapOpen() {
		return gapOpen;
	}
	public void setGapOpen(int gapOpen) {
		this.gapOpen = gapOpen;
	}
	public int getGapExtend() {
		return gapExtend;
	}
	public void setGapExtend(int gapExtend) {
		this.gapExtend = gapExtend;
	}
	public String getFilterString() {
		return filterString;
	}
	public void setFilterString(String filterString) {
		this.filterString = filterString;
	}
}
