package no.uio.medisin.bag.exceptions;


/**
 * for handling errors related to parsing BLAST input/output
 * @author simonray
 *
 */
public class BLASTException extends Exception {

	private static final long serialVersionUID = -2511080310126748855L;
		public BLASTException() { super(); }
	  public BLASTException(String message) { super(message); }
	  public BLASTException(String message, Throwable cause) { super(message, cause); }
	  public BLASTException(Throwable cause) { super(cause); }
	
}
