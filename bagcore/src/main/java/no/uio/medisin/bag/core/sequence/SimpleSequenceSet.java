/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList; 
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * Class to Serialize, Search & Process a set of SimpleSeqs 
 * 
 * @author weibo / Simon Rayner
 */
public class SimpleSequenceSet {

    static Logger                       logger = LogManager.getLogger();

    private File                        fastaFilename;
    private BufferedReader              brFA;
    private ArrayList<SimpleSeq>        seqs;

    
    
    /**
     * read query sequences in FASTA format
     * 
     * @param fname : query sequence filename

     */
    public SimpleSequenceSet(File fname){
      this.fastaFilename=fname;
    }
    
    public SimpleSequenceSet(){
    	seqs  = new ArrayList<>();
    }
    
    /**
     * Read in a list of sequences in FASTA format.  Stores results in a SimpleSeq List
     * @param filename
     * @throws java.io.IOException
     */
    public void readFromFastaFile(String filename) throws java.io.IOException{
    	this.fastaFilename = new File(filename);
    	readFromFastaFile();
    }
    
    
    /**
    * Read in a list of sequences in FASTA format.  Stores results in a SimpleSeq List
    *
    * @throws java.io.IOException
    */
    public void readFromFastaFile() throws java.io.IOException{
        logger.info("Start loading data " + fastaFilename.getName()+" ...");
        long totalSize=fastaFilename.length();
        long readedSize=0;
        
        seqs=new ArrayList<>();

        String line = "";
        String title = "";
        try{
            brFA = new BufferedReader(new FileReader(fastaFilename));

            StringBuilder seq = new StringBuilder();

            title = brFA.readLine();
            readedSize+=title.length()+1;
            if(title.substring(0, 1).equals(">") == false){
                throw new java.io.IOException("bad FASTA format on first line\n");
            }
            title = title.substring(1, title.length());

            while((line = brFA.readLine()) != null){
                readedSize+=line.length()+1;
                if(line.startsWith(">")==true){
                    line = line.substring(1, line.length());
                    addSequence(title, seq.toString());
                    title = line;
                    seq = new StringBuilder();
                }
                else{
                    line=line.replaceAll("\\s+", "").replace('T', 'u').toUpperCase();
                    seq.append(line);
                }
            }
            addSequence(title,seq.toString());
            brFA.close();
        }
        catch(IOException exIO){
            logger.info("error reading sequence <"+ title + "> from file <" + fastaFilename.getName() + ">");
            logger.error("error reading sequence <"+ title + "> from file <" + fastaFilename.getName() + ">");
            throw new IOException("error reading sequence <"+ title + "> from file <" + fastaFilename.getName() + ">");
        }
        logger.info("completed successfully" + fastaFilename.getName()+" ...");
    }
    
    
    
    /**
     * add this sequence to the sequence set
     * 
     * @param title
     * @param seq 
     */
    private void addSequence(String title, String seq) {

        getSeqs().add(new SimpleSeq(title,seq));        

    }

    
    
    /**
     * add this sequence to the sequence set
     * 
     * 
     * @param newSeq 
     */
    public void addSequence(SimpleSeq newSeq) {

        getSeqs().add(newSeq);        

    }

    
    
    /**
     * return the sequence list
     * 
     * @return ArrayList<SimSeq>: the seqs
     */
    public ArrayList<SimpleSeq> getSeqs() {
        return seqs;
    }


  }
