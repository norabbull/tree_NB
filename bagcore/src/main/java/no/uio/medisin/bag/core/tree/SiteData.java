/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.tree;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.BufferedWriter;
import no.uio.medisin.bag.core.aminoacid.aminoAcidListx;
import no.uio.medisin.bag.core.aminoacid.entropy;
//import java.util.*;
/**
 *
 * @author sr
 */
public class SiteData implements java.io.Serializable {
    private int noOfSequences;
    private aminoAcidListx aal;
    private entropy ent;
    private MutualInfo mutInfo;
    private int position;

    public SiteData(int n, String name)
    {
      aal = new aminoAcidListx(name);
      ent = new entropy();
      mutInfo = new MutualInfo(n);
    }

    public aminoAcidListx getAAL()
    {
      return aal;
    }

    public entropy getEntropy()
    {
      return ent;
    }

    public void setEntropy(double e)
    {
      ent.setEntropy(e);
    }

    public int getNoOfSeqs()
    {
      return noOfSequences;
    }

    public int getSitePosition()
    {
      return position;
    }

    public void setSitePosition(int p)
    {
      position = p;
    }

    public MutualInfo getMutualInfo()
    {
      return mutInfo;
    }

    public void writeMutInfoBin(DataOutputStream dos) throws java.io.IOException
    {
      getMutualInfo().writeMutInfoBin(dos, 0.0);
    }

    public void writeMutInfo(BufferedWriter bw) throws java.io.IOException
    {
      getMutualInfo().writeMutInfo(bw, 0.0);
    }

    public void readMutInfoBin(DataInputStream dos) throws java.io.IOException
    {
      int j=0;
      mutInfo.setMInfo(j++, dos.readDouble());
    }

    public void writeEntropy(BufferedWriter bw) throws java.io.IOException
    {
      bw.write(getEntropy().toString());
    }

    public void writeProbability(BufferedWriter bw, double threshold)
            throws java.io.IOException
    {
      bw.write(getEntropy().toString());
    }

    public void setProb(int i, double p)
    {
      mutInfo.setProb(i, p);
    }

    public double getProb(int i)
    {
      return mutInfo.getProb(i);
    }
}
