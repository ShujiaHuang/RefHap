/*******************************************************************************
 * SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * Copyright 2011 Jorge Duitama
 *
 * This file is part of SingleIndividualHaplotyper.
 *
 *     SingleIndividualHaplotyper is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     SingleIndividualHaplotyper is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with SingleIndividualHaplotyper.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package mpg.molgen.sih.model;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


public class FragmentsFileHandler {
	public List<Fragment> loadFragments(String filename) throws IOException {
		List<Fragment> fragments = new ArrayList<Fragment>();
		FileInputStream fis = new FileInputStream(filename);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		String line=in.readLine();
		for(int i=0;line!=null;i++) {
			String [] items = line.split(" |\t");
			if(items.length > 2) {
				try {
					StringBuilder calls = new StringBuilder();
					String qualityScores = null;
					int startPos = Integer.parseInt(items[2])-1;
					int pos = startPos;
					for(int j=2;j<items.length;j+=2) {
						//Assumes that start positions are one based
						if(j<items.length-1) {
							int nextStart = Integer.parseInt(items[j])-1;
							while(pos<nextStart) {
								calls.append(Fragment.NODATACHAR);
								pos++;
							}
							calls.append(items[j+1]);
							pos+=items[j+1].length();
						} else {
							qualityScores = items[j];
						}
					}
					int k=0;
					double [] extendedProbabilities = null;
					Fragment f;
					if(qualityScores!=null) {
						extendedProbabilities = new double[calls.length()];
						Arrays.fill(extendedProbabilities, 0);
						for(int j=0;j<calls.length();j++) {
							if(calls.charAt(j)!=Fragment.NODATACHAR) {
								extendedProbabilities [j] = calculateProbability(qualityScores.charAt(k));
								k++;
							}
						}
						
						f = new Fragment(items[1],startPos, calls.toString(),extendedProbabilities);
					} else {
						f = new Fragment(items[1],startPos, calls.toString());
					}
					fragments.add(f);
				} catch(Exception e) {
					throw new IOException("Error reading line: "+ line+ " of file "+filename,e);
				}
			}
			
			line=in.readLine();
		}
		fis.close();
		Collections.sort(fragments,new FragmentsComparator());
		return fragments;
	}
	private double calculateProbability(char qual) {
		double score = (double)(qual -33);
		return 1 - Math.pow(10, -score /10);
	}
	public List<Integer> loadVariantPositions(String filename, int posCol) throws IOException {
		List<Integer> variantPos = new ArrayList<Integer>();
		FileInputStream fis = new FileInputStream(filename);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		String line=in.readLine();
		for(int i=0;line!=null;i++) {
			String [] items = line.split(" |\t");
			try {
				variantPos.add(Integer.parseInt(items[posCol]));
			} catch(NumberFormatException e) {
				System.err.println("WARN: Could not read variant id in line: "+line );
			}
			line=in.readLine();
		}
		fis.close();
		return variantPos;
	}
}
