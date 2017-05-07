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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Block {
	private List<Fragment> fragments;
	private int firstPos;
	private int lastPos;
	private String haplotype;
	private int MEC=0;
	private int phased=0;
	private int calls=0;
	private double hapScores [];
	
	
	public Block(List<Fragment> fragments) {
		super();
		this.fragments = fragments;
		this.firstPos = fragments.get(0).getFirstPos();
		this.lastPos = calculateLastPos(fragments); 
	}
	private int calculateLastPos(List<Fragment> fragments) {
		int answer =0;
		for(Fragment f:fragments) {
			if(answer < f.getLastPos()) {
				answer = f.getLastPos();
			}
		}
		return answer;
	}
	/**
	 * @return the solution
	 */
	public String getHaplotype() {
		return haplotype;
	}
	
	public void setHaplotype(String haplotype) {
		assert haplotype.length() == this.length();
		this.haplotype = haplotype;
		updateStatistics();
	}
	
	private void updateStatistics() {
		phased=0;
		for(int i=0;i<haplotype.length();i++) {
			if(haplotype.charAt(i)!=Fragment.NODATACHAR) {
				phased++;
			}
		}
		MEC = 0;
		calls=0;
		hapScores = new double[haplotype.length()];
		int [] countsGood = new int [haplotype.length()];
		int [] countsTotal = new int [haplotype.length()];
		Arrays.fill(countsGood, 0);
		Arrays.fill(countsTotal, 0);
		for(Fragment f:fragments) {
			int nCalls = f.getnCalls();
			calls+=nCalls;
			int d1 = f.getHammingDistance(haplotype,firstPos);
			int overlapping = f.getOverlappingCount(haplotype, firstPos);
			int d2 = overlapping - d1;
			int dMin = Math.min(d1, d2);
			MEC+= dMin;
			boolean complement = (dMin == d2);
			String fcalls = f.getExtendedCalls();
			for(int i=f.getFirstPos();i<=f.getLastPos();i++) {
				int relPos = i-firstPos;
				char hapCall = haplotype.charAt(relPos);
				char call = fcalls.charAt(i-f.getFirstPos());
				if(hapCall!=Fragment.NODATACHAR && call!=Fragment.NODATACHAR) {
					if(complement == (call != hapCall)) {
						countsGood[relPos]++;
					}
					countsTotal[relPos]++;
				}
			}
		}
		for(int i=0;i<hapScores.length;i++) {
			if(countsTotal[i]==0) {
				hapScores[i] = 0;
			} else {
				hapScores[i] = (double)countsGood[i]/countsTotal[i];
			}
		}
	}
	public int getMEC() {
		return MEC;
	}
	public int getPhased() {
		return phased;
	}
	public int getPhased(int first, int last) {
		int answer = 0;
		for(int i=first;i<=last;i++) {
			if(haplotype.charAt(i-firstPos)!=Fragment.NODATACHAR) {
				answer++;
			}
		}
		return answer;
	}
	public int getCalls() {
		return calls;
	}
	/**
	 * @return the fragments
	 */
	public List<Fragment> getFragments() {
		return fragments;
	}
	/**
	 * @return the firstPos
	 */
	public int getFirstPos() {
		return firstPos;
	}
	/**
	 * @return the lastPos
	 */
	public int getLastPos() {
		return lastPos;
	}
	public int length() {
		return lastPos-firstPos+1;
	}
	public List<Fragment> getFragments(int first,int last) {
		List<Fragment> answer = new ArrayList<Fragment>();
		for(Fragment f:fragments) {
			if(f.getFirstPos()<=last && f.getLastPos()>=first) {
				answer.add(f);
			}
		}
		return answer;
	}
	/**
	 * @return the hapScores
	 */
	public double[] getHapScores() {
		return hapScores;
	}
	
}
