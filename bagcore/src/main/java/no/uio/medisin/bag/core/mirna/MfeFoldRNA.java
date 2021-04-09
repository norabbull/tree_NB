
package no.uio.medisin.bag.core.mirna;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Performs RNA sequence folding using the RNAFold library
 * @author weibo
 */
public class MfeFoldRNA {

    static Logger                       logger = LogManager.getLogger();
    
    private String pathToLibrary;
    private static native float fold(String seq, float t);
    private static native void initIDs();
    
    
    static{
        logger.info(System.getProperty("java.library.path"));
        if(System.getProperty("os.name").toUpperCase().contains("MAC")){
            System.loadLibrary("RNAFold");            
        }else{
            if(System.getProperty("os.name").toUpperCase().contains("LINUX")){
                System.loadLibrary("RNAFold");                        
            }else
            {
                logger.info("unsupported OS: <" + System.getProperty("os.name") + ">");
                logger.error("unsupported OS: <" + System.getProperty("os.name") + ">");
                throw new RuntimeException("unsupported OS: <" + System.getProperty("os.name") + ">");
            }
        }
        
	initIDs();
        
    }

    private static final float  temperature=37;
    private static String       structure="";



    public static float foldSequence(String sequence){
        structure = "";
        return MfeFoldRNA.fold(sequence, temperature);
    }
    
    
    
    public static float foldSequence(String sequence, String str){
        structure = str;
        return MfeFoldRNA.fold(sequence, temperature);
    }

    
    
    public static String getStructure(){
        return structure;
    }

}
