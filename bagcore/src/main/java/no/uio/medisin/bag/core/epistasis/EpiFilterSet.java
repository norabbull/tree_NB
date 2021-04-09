/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.epistasis;

/**
 *
 * @author sr
 */
public class EpiFilterSet {
  private double minAlignDist;
  private double maxAlignDist;
  private int minEpiSeqCount;
  private int maxEpiSeqCount;
  private int minMutSeqCount;
  private int maxMutSeqCount;
  private double minSMetric;
  private double maxSMetric;
  private int minBLcount;
  private int maxBLcount;
  private double minBLspread;
  private double maxBLspread;
  private double minBLmean;
  private double maxBLmean;
  private double minBLsd;
  private double maxBLsd;
  private double minBLmoment;
  private double maxBLmoment;
  private double minBL;
  private double maxBL;
  private Boolean tripped;

  public EpiFilterSet()
  {
    minAlignDist = -100;
    maxAlignDist = 100000000;
    minEpiSeqCount = -100;
    maxEpiSeqCount = 100000000;
    minMutSeqCount = -100;
    maxMutSeqCount = 100000000;
    minSMetric = -100.0;
    maxSMetric = 100000000.0;
    minBLcount = -100;
    maxBLcount = 100000000;
    minBLspread = -100.0;
    maxBLspread = 100000000.0;
    minBLmean = -100.0;
    maxBLmean = 100000000.0;
    minBLsd = -100.0;
    maxBLsd = 100000000.0;
    minBLmoment = -100.0;
    maxBLmoment = 100000000.0;
    minBL = -100.0;
    maxBL = 100000000.0;
    tripped = false; // this is set to true if any of the set methods are called
  }

  public void resetToZero()
  {
    minAlignDist = 0.0;
    maxAlignDist = 0.0;
    minEpiSeqCount = 0;
    maxEpiSeqCount = 0;
    minMutSeqCount = 0;
    maxMutSeqCount = 0;
    minSMetric = 0.0;
    maxSMetric = 0.0;
    minBLcount = 0;
    maxBLcount = 0;
    maxBLspread = 0.0;
    minBLmean = 0.0;
    minBLmean = 0.0;
    maxBLmean = 0.0;
    minBLsd = 0.0;
    maxBLsd = 0.0;
    minBLmoment = 0.0;
    maxBLmoment = 0.0;
    setMinBL(0.0);
    setMaxBL(0.0);
    tripped = false; // this is set to true if any of the set methods are called
  }
  /**
   * @return the minAlignDist
   */
  public double getMinAlignDist() {
    return minAlignDist;
  }

  /**
   * @param minAlignDist the minAlignDist to set
   */
  public void setMinAlignDist(double minAlignDist) {
    this.minAlignDist = minAlignDist;
    tripped = true;
  }

  /**
   * @return the maxAlignDist
   */
  public double getMaxAlignDist() {
    return maxAlignDist;
  }

  /**
   * @param maxAlignDist the maxAlignDist to set
   */
  public void setMaxAlignDist(double maxAlignDist) {
    this.maxAlignDist = maxAlignDist;
    tripped = true;
  }

  /**
   * @return the minEpiSeqCount
   */
  public int getMinEpiSeqCount() {
    return minEpiSeqCount;
  }

  /**
   * @param minEpiSeqCount the minEpiSeqCount to set
   */
  public void setMinEpiSeqCount(int minEpiSeqCount) {
    this.minEpiSeqCount = minEpiSeqCount;
    tripped = true;
  }

  /**
   * @return the maxEpiSeqCount
   */
  public int getMaxEpiSeqCount() {
    return maxEpiSeqCount;
  }

  /**
   * @param maxEpiSeqCount the maxEpiSeqCount to set
   */
  public void setMaxEpiSeqCount(int maxEpiSeqCount) {
    this.maxEpiSeqCount = maxEpiSeqCount;
    tripped = true;
  }

  /**
   * @return the minMutSeqCount
   */
  public int getMinMutSeqCount() {
    return minMutSeqCount;
  }

  /**
   * @param minMutSeqCount the minMutSeqCount to set
   */
  public void setMinMutSeqCount(int minMutSeqCount) {
    this.minMutSeqCount = minMutSeqCount;
    tripped = true;
  }

  /**
   * @return the maxMutSeqCount
   */
  public int getMaxMutSeqCount() {
    return maxMutSeqCount;
  }

  /**
   * @param maxMutSeqCount the maxMutSeqCount to set
   */
  public void setMaxMutSeqCount(int maxMutSeqCount) {
    this.maxMutSeqCount = maxMutSeqCount;
    tripped = true;
  }

  /**
   * @return the minSMetric
   */
  public double getMinSMetric() {
    return minSMetric;
  }

  /**
   * @param minSMetric the minSMetric to set
   */
  public void setMinSMetric(double minSMetric) {
    this.minSMetric = minSMetric;
    tripped = true;
  }

  /**
   * @return the maxSMetric
   */
  public double getMaxSMetric() {
    return maxSMetric;
  }

  /**
   * @param maxSMetric the maxSMetric to set
   */
  public void setMaxSMetric(double maxSMetric) {
    this.maxSMetric = maxSMetric;
    tripped = true;
  }

  /**
   * @return the minBLcount
   */
  public int getMinBLcount() {
    return minBLcount;
  }

  /**
   * @param minBLcount the minBLcount to set
   */
  public void setMinBLcount(int minBLcount) {
    this.minBLcount = minBLcount;
    tripped = true;
  }

  /**
   * @return the maxBLcount
   */
  public int getMaxBLcount() {
    return maxBLcount;
  }

  /**
   * @param maxBLcount the maxBLcount to set
   */
  public void setMaxBLcount(int maxBLcount) {
    this.maxBLcount = maxBLcount;
    tripped = true;
  }

  /**
   * @return the minBLmean
   */
  public double getMinBLmean() {
    return minBLmean;
  }

  /**
   * @param minBLmean the minBLmean to set
   */
  public void setMinBLmean(double minBLmean) {
    this.minBLmean = minBLmean;
    tripped = true;
  }

  /**
   * @return the maxBLmean
   */
  public double getMaxBLmean() {
    return maxBLmean;
  }

  /**
   * @param maxBLmean the maxBLmean to set
   */
  public void setMaxBLmean(double maxBLmean) {
    this.maxBLmean = maxBLmean;
    tripped = true;
  }

  /**
   * @return the minBLsd
   */
  public double getMinBLsd() {
    return minBLsd;
  }

  /**
   * @param minBLsd the minBLsd to set
   */
  public void setMinBLsd(double minBLsd) {
    this.minBLsd = minBLsd;
    tripped = true;
  }

  /**
   * @return the maxBLsd
   */
  public double getMaxBLsd() {
    return maxBLsd;
  }

  /**
   * @param maxBLsd the maxBLsd to set
   */
  public void setMaxBLsd(double maxBLsd) {
    this.maxBLsd = maxBLsd;
    tripped = true;
  }

  /**
   * @return the minBLmoment
   */
  public double getMinBLmoment() {
    return minBLmoment;
  }

  /**
   * @param minBLmoment the minBLmoment to set
   */
  public void setMinBLmoment(double minBLmoment) {
    this.minBLmoment = minBLmoment;
    tripped = true;
  }

  /**
   * @return the maxBLmoment
   */
  public double getMaxBLmoment() {
    return maxBLmoment;
  }

  /**
   * @param maxBLmoment the maxBLmoment to set
   */
  public void setMaxBLmoment(double maxBLmoment) {
    this.maxBLmoment = maxBLmoment;
    tripped = true;
  }

  /**
   * @return the tripped
   */
  public Boolean getTripped() {
    return tripped;
  }

  /**
   * @return the minBL
   */
  public double getMinBL() {
    return minBL;
  }

  /**
   * @param minBL the minBL to set
   */
  public void setMinBL(double minBL) {
    this.minBL = minBL;
  }

  /**
   * @return the maxBL
   */
  public double getMaxBL() {
    return maxBL;
  }

  /**
   * @param maxBL the maxBL to set
   */
  public void setMaxBL(double maxBL) {
    this.maxBL = maxBL;
  }

  /**
   * @return the minBLspread
   */
  public double getMinBLspread() {
    return minBLspread;
  }

  /**
   * @param minBLspread the minBLspread to set
   */
  public void setMinBLspread(double minBLspread) {
    this.minBLspread = minBLspread;
  }

  /**
   * @return the maxBLspread
   */
  public double getMaxBLspread() {
    return maxBLspread;
  }

  /**
   * @param maxBLspread the maxBLspread to set
   */
  public void setMaxBLspread(double maxBLspread) {
    this.maxBLspread = maxBLspread;
  }

}