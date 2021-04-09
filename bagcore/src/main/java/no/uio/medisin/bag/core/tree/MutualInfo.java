/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.BufferedWriter;

/**
 *
 * @author sr
 */
public class MutualInfo implements java.io.Serializable {

  private double prob[];
  private double entropyIJ[];
  private double mInfo[];

  private double maxMutual;

  MutualInfo(int i)
  {
    prob = new double[i];
    entropyIJ = new double[i];
    mInfo = new double[i];
    maxMutual = 0;
  }

  public double[] getEntropy()
  {
    return entropyIJ;
  }

  public void setEntropyIJ(int i, double value)
  {
    entropyIJ[i] = value;
  }

  public double[] getMInfo()
  {
    return mInfo;
  }

  public void setMInfo(int i, double value)
  {
    mInfo[i] = value;
    if(value>maxMutual)
    {
      maxMutual=value;
    }
  }

  public double findMax()
  {
    int i=0;
    while(i<mInfo.length)
    {
      if(mInfo[i]>maxMutual)
      {
        maxMutual=mInfo[i];
      }
      i++;
    }
    return maxMutual;
  }

  public double getMax()
  {
    return maxMutual;
  }

  public double getMInfo(int i)
  {
    return mInfo[i];
  }
  public double[] getProb()
  {
    return prob;
  }

  public void setProb(int i, double value)
  {
    prob[i] = value;
  }

  public double getProb(int i)
  {
    return prob[i];
  }

  public void writeMutInfoBin(DataOutputStream dos, double threshold) throws java.io.IOException
  {
    int j=0;
    while(j< mInfo.length)
    {
      if(mInfo[j] >= threshold)
      {
        dos.writeDouble(mInfo[j++]);
      }
      else
      {
        dos.writeDouble(0.0);
      }
    }
  }

  public void writeMutInfo(BufferedWriter bw, double threshold) throws java.io.IOException
  {
    int j=0;
    while(j< mInfo.length)
    {
      if(mInfo[j] >= threshold)
      {
        bw.write(Double.toString(mInfo[j]));
      }
      else
      {
        bw.write(Double.toString(0.0));
      }
      j++;
    }
  }

  public void readMutInfoBin(DataInputStream dos) throws java.io.IOException
  {
    int j=0;
    mInfo[j++] = dos.readDouble();
  }

//  public void writeEntropy(BufferedWriter bw, double threshold) throws java.io.IOException
//  {
//    bw.write(getEntropy().toString());
//  }

//  public void writeProbability(BufferedWriter bw, double threshold)
//          throws java.io.IOException
//  {
//    bw.write(getEntropy().toString());
//  }

}
