/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.core.tree;

import java.util.HashMap;

/**
 * This stores an information about a type of sample classification
 * For example, a sample can be part of a SuperPopulation (Africa) and SubPopulation (Yoruba)
 * These would be stored as two SampleClassification instances
 * 
 * For example:
 *   'Super', 'AFR'
 *   'Sub'  , 'YOR'
 * 
 * where these would match the information associated with the samples in the tree
 * 
 * It's a simple solution for now, but makes it possible to iterate through 
 * all 'Super' and 'Sub' types and perform operations on them 
 * (such as counting up the numbers in each class)
 * 
 * @author simonray
 */
public class SampleClassification {
  private String classificationType;        // e.g. 'Super'
  private String classificationName;        // e.g. 'AFR'
  private String classificationDescription; // e.g. 'AFR'
  //  used to store additional information generated about this class
  // for example, average pairwise distance
  private HashMap<String, String>  otherParams  = new HashMap<>();  
  

  public SampleClassification(String cType, String cName, String cDesc ){
    classificationType = cType;
    classificationName = cName;
  }
  /**
   * @return the classificationType
   */
  public String getClassificationType() {
    return classificationType;
  }

  /**
   * @param classificationType the classificationType to set
   */
  public void setClassificationType(String classificationType) {
    this.classificationType = classificationType;
  }

  /**
   * @return the classificationName
   */
  public String getClassificationName() {
    return classificationName;
  }

  /**
   * @param classificationName the classificationName to set
   */
  public void setClassificationName(String classificationName) {
    this.classificationName = classificationName;
  }

  /**
   * @return the classificationDescription
   */
  public String getClassificationDescription() {
    return classificationDescription;
  }

  /**
   * @param classificationDescription the classificationDescription to set
   */
  public void setClassificationDescription(String classificationDescription) {
    this.classificationDescription = classificationDescription;
  }
          
}
