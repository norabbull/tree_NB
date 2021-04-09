/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

/**
 *this class use two different threads to handle standard input stream reader, 
 * get stdIn and stdError message to avoid deadlock.
 * having the TRUE statement in while loop is important, without it, pipeline crashed again.
 * maybe pipeline hang problem is coming from here, other than Runtime.exec()  or thread
 * need look more into this.
 * @author joey
 */
public class ReadStdIOStream implements Runnable {

    String msgType;
    InputStream inputStream;
    Thread thread;
    static Logger logger = LogManager.getLogger();
    String stepID;
    private static final String STD_ERROR = "stderr";
    private static final String STD_IN = "stdin";
    
    private ArrayList<String> streamData;


    public ReadStdIOStream(String msgType, InputStream inputStream, String stepIDString) {
        this.msgType = msgType;
        this.inputStream = inputStream;
        this.stepID = stepIDString;
    }

    
    
    
    public void start() {
        thread = new Thread(this);
        thread.start();
    }

    
    
    /**
    * having the TRUE statement in the while loop is important, without it, 
    * pipeline crashes maybe pipeline hang problem is coming from here, 
    * other than Runtime.exec() or thread - need look more into this.
     */
    @Override
    public void run() {
        try {
            InputStreamReader isr = new InputStreamReader(inputStream);
            BufferedReader br = new BufferedReader(isr);
            streamData = new ArrayList<>();
            
 
            if (msgType.toLowerCase().equals(STD_IN)) {
                String line = "";
                while (true && (line = br.readLine()) != null) {
                    logger.info("OUTPUT:" + line);
                    streamData.add("OUTPUT:" + line);
                }
            }


            if (msgType.toLowerCase().equals(STD_ERROR)) {
                int skipCount = 0;

                String line = "";
                while (true && (line = br.readLine()) != null) {
                    if (line.contains("Warning: Skipping") && line.contains("less than")) {
                        skipCount++;
                    } else {
                        logger.info("ERROR:" + line);
                        streamData.add("ERROR:" + line);
                    }
                }
                logger.info(skipCount + " lines were skipped because the read was too short");
                inputStream.close();
            }
        } catch (IOException ex) {
            logger.error("Exception executing command from within " + this.stepID + "class");
            logger.error(ex);
        }
    }

    
    
    
    /**
    * return the data collected for the stream
    * @return streamData 
    */
    public ArrayList<String> getStreamData() {
        return streamData;
    }
}
