/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.targetscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

/**
 *
 * @author simonray
 */
public class TargetScanPredictedTargetList {
    
    static Logger                       logger                          = LogManager.getLogger();



    
    private ArrayList<TargetScanPredictedTarget> targetScanPredictedTargets;
    
    
    /**
     * load target entries from Target Scan 
     * @param predictedTargetInfoFilename String
     * @throws IOException thrown if error encountered during file loading
     */
    public void loadPredictedTargetInfo(String predictedTargetInfoFilename) throws IOException{
        targetScanPredictedTargets = new ArrayList<>();
        
        String targetLine ="";
        try(BufferedReader brCF = new BufferedReader(new FileReader(new File(predictedTargetInfoFilename)))){
            brCF.readLine();
            while((targetLine=brCF.readLine())!=null){
                targetScanPredictedTargets.add(new TargetScanPredictedTarget(targetLine));
            }
            brCF.close(); 
        }
        catch(IOException exIO){
            
            logger.error("error reading <Predicted_Targets_Info< file < " + predictedTargetInfoFilename);
            logger.error("error occurred on this line");
            logger.error(targetLine);
            
            throw new IOException("error reading " + "<Predicted_Targets_Info< file < " + predictedTargetInfoFilename + ">\n"
                    + "error occurred on this line\n" + targetLine);
            
        }
        
    }
    
    
    
    /**
     * match the input list of miR Family ids to the PredictedTarget entries
     * 
     * @param qMiRfamilyList ArrayList of Strings
     * @return ArrayList of Strings containing hits to supplied miRNA list
     */
    public ArrayList<String> getGeneTargetHits(ArrayList<String> qMiRfamilyList){
        ArrayList<String> geneTargetHits = new ArrayList<>();
        
        for( TargetScanPredictedTarget predictedTarget: targetScanPredictedTargets){
            if(qMiRfamilyList.contains(predictedTarget.getMiRFamily())){
                String targetString = predictedTarget.getMiRFamily()  + "|" + predictedTarget.getGeneID() 
                        + predictedTarget.getGeneSymbol() + "|" + predictedTarget.getTranscriptID();
                if(geneTargetHits.contains(targetString)==false){
                    geneTargetHits.add(targetString);
                }
            }
        }
        
        return geneTargetHits;
    }
            
    /**
     * return the number of entries in the dataset
     * 
     * @return int number of entries in predicted target list
     */
    public int getNumberOfEntries(){
        return targetScanPredictedTargets.size();
    }
    
}
