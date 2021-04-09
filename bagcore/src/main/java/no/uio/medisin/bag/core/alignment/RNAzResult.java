/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.alignment;

import java.io.*;


/**
 * This stores the results from running RNAz.
 * This is a program that predicts conserved RNA structures from a sequence
 * alignment and gives an estimate that the predicted structure is meaningful.
 * Result file has the following header
 *
 * ############################  RNAz 2.0  ##############################
 *
 *  Sequences: 6
 *  Columns: 121
 *  Reading direction: forward
 *  Mean pairwise identity:  96.45
 *  Shannon entropy: 0.06800
 *  G+C content: 0.39389
 *  Mean single sequence MFE: -30.95
 *  Consensus MFE: -28.89
 *  Energy contribution: -29.25
 *  Covariance contribution:   0.36
 *  Combinations/Pair:   1.03
 *  Mean z-score:  -1.74
 *  Structure conservation index:   0.93
 *  Background model: mononucleotide
 *  Decision model: sequence based alignment quality
 *  SVM decision value:   2.30
 *  SVM RNA-class probability: 0.991995
 *  Prediction: RNA
 *
 * ######################################################################
 *
 * @author sr
 */
public class RNAzResult {

  // READ DIRECTION
  public final static int FORWARD = 0;
  public final static int REVERSE = 1;
  // PREDICTION
  public final static int RNA = 0;
  public final static int NOTRNA = 1;

  // DECISION MODEL
  public final static int SEQUENCE_BASED_ALIGNMENT_QUALITY = 0;
  public final static int UNKNOWN_DECISION = 1;
  // BACKGROUND MODEL
  public final static int MONONUCLEOTIDE = 0;
  public final static int UNKNOWN_BACKGROUND = 1;
  
  private int noOfSequences;
  private int alignmentLength;
  private int readingDirection;
  private float pairwiseIdentity;
  private float shannonEntropy;
  private float GCcontent;
  private float meanSingleSequenceMFE;
  private float consensusMFE;
  private float energyContribution;
  private float covarianceContribution;
  private float combinationsPerPair;
  private float meanZScore;
  private float structureConservationIndex;
  private int backgroundModel;
  private int decisionModel;
  private float svmDecisionValue;
  private float svmRNAclassProbability;
  private int prediction;

  public RNAzResult()
  {

  }

  /**
   * Reads in the result file from RNAz prediction
   * @param br
   * @throws java.io.IOException
   */
  public void readRNAzResultHeader(BufferedReader br) throws java.io.IOException
  {
    int l=0;
    // One blank lines, one Delimiter line, One blank line.
    String line = br.readLine();
    line = br.readLine();
    line = br.readLine();

    noOfSequences = Integer.parseInt(br.readLine().split("\\:")[1].trim());
    alignmentLength = Integer.parseInt(br.readLine().split("\\:")[1].trim());
    String directionString = br.readLine().split("\\:")[1].trim();
    if(directionString.equals("forward"))
      readingDirection = FORWARD;
    else
      readingDirection = REVERSE;
      
    pairwiseIdentity  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    shannonEntropy  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    GCcontent = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    meanSingleSequenceMFE  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    consensusMFE  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    energyContribution  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    covarianceContribution  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    combinationsPerPair  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    meanZScore  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    structureConservationIndex = Float.parseFloat(br.readLine().split("\\:")[1].trim());

    String modelString = br.readLine().split("\\:")[1].trim();
    if(modelString.equals("mononucleotide"))
      backgroundModel = MONONUCLEOTIDE;
    else
      backgroundModel = UNKNOWN_BACKGROUND;

    String decisionString = br.readLine().split("\\:")[1].trim();
    if(decisionString.equals("sequence based alignment quality"))
      decisionModel = SEQUENCE_BASED_ALIGNMENT_QUALITY;
    else
      decisionModel = UNKNOWN_DECISION;

    svmDecisionValue  = Float.parseFloat(br.readLine().split("\\:")[1].trim());
    svmRNAclassProbability  = Float.parseFloat(br.readLine().split("\\:")[1].trim());


    String predictString = br.readLine().split("\\:")[1].trim();
    if(predictString.equals("RNA"))
      prediction = RNA;
    else
      prediction = NOTRNA;
    
  }

  public void writeHeaderForSummary(BufferedWriter bw) throws java.io.IOException
  {
    bw.write("position\tmeanZScore\tsvmRNAclassProbability\tstructureConservationIndex\t"
            + "GCcontent\tmeanSingleSequenceMFE\tconsensusMFE\n"
            );
  }

  public void writeSummary(BufferedWriter bw) throws java.io.IOException
  {
    bw.write(meanZScore + "\t"
            + svmRNAclassProbability + "\t"
            + structureConservationIndex + "\t"
            + GCcontent + "\t"
            + meanSingleSequenceMFE + "\t"
            + consensusMFE + "\n"
            );
  }
}
