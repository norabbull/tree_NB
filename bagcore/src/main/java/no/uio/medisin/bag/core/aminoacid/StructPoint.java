/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package no.uio.medisin.bag.core.aminoacid;

import no.uio.medisin.bag.core.aminoacid.AACoords;

/**
 * Stores coordinate information associated with the amino acids in a
 * protein structure.
 * This was really created for overlaying the solved sections of a
 * protein structure on a sequence so they can be aligned.
 * A structure point consists of a 3d position stored in an AACoord
 * instance and a site positiom indicating the location within a reference
 * sequence.
 * @author Lep Rayner
 */
public class StructPoint {

  AACoords coords;
  int site;

  /**
   * Create a new StructPoint at site s and coordinate (px, py, pz)
   * @param px
   * @param py
   * @param pz
   * @param s
   */
  public StructPoint(double px, double py, double pz, int s)
  {
    coords = new AACoords(px, py, pz);
    site = s;
  }

  /**
   * return site
   * @return site
   */
  public int getSite()
  {
    return site;
  }

  /**
   * return coordinates in a AACoord instance
   * @return coords
   */
  public AACoords returnCoords()
  {
    return coords;
  }

}
