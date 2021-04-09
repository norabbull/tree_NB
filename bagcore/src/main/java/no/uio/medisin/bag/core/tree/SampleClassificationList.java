/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.tree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author simonray
 */
public class SampleClassificationList {
  static Logger				logger = LogManager.getLogger();  
  private ArrayList<SampleClassification> sampleClassificationList;  
  public SampleClassificationList(){
    sampleClassificationList = new ArrayList<>();
  }
  
  /**
   * return a list of unique ClassificationTypes in the list
   * For example, all specified Super and Sub populations
   * @return 
   */
  public ArrayList<String> getClassificationTypes(){
    ArrayList<String> uniqClassifications = new ArrayList<>();
    for (SampleClassification sampleClassification: sampleClassificationList){
      if(uniqClassifications.contains(sampleClassification.getClassificationType())==false){
        uniqClassifications.add(sampleClassification.getClassificationType());
      }
    }
    return uniqClassifications;
  }
  
  /**
   * get all entries in the list that match the specified ClassificationType
   * For example, all entries of AFR super-population 
   * @param classificationType
   * @return 
   */
  public ArrayList<String> getAllEntriesForClassificationType(String classificationType){
    ArrayList<String> matchEntries = new ArrayList<>();
    for (SampleClassification sampleClassification: sampleClassificationList){
      if(sampleClassification.getClassificationType().equals(classificationType)){
        matchEntries.add(sampleClassification.getClassificationName());
      }
    }
    return matchEntries;
  }
  
  
  /**
   * load list of sample classes from tab delimited file
   * 
   * @param sampleGroupFile
   * @throws java.io.IOException 
   */
  public void readSampleGroupFile(String sampleGroupFile) throws java.io.IOException{
    logger.info("--reading Sample Group File <" + sampleGroupFile + ">...");
    String sampleLine = "";
    BufferedReader brSG = new BufferedReader(new FileReader(new File(sampleGroupFile)));
    // Skip Header
    sampleLine = brSG.readLine();
      while((sampleLine = brSG.readLine()) !=null){
        this.sampleClassificationList.add(new SampleClassification(
                sampleLine.split("\t")[0].trim(),
                sampleLine.split("\t")[1].trim(),
                sampleLine.split("\t")[2].trim()));
      }
      logger.info("--read <" + this.sampleClassificationList.size() + "> entries");
    brSG.close();
  }  
}
