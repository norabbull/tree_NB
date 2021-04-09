/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;

/**
 *
 * @author sr
 */
  public class AACoords{
    double x;
    double y;
    double z;

    AACoords(double px, double py, double pz)
    {
      x = px;
      y = py;
      z = pz;
    }

    public double getX()
    {
      return x;
    }

    public double getY()
    {
      return y;
    }

    public double getZ()
    {
      return z;
    }
  }
