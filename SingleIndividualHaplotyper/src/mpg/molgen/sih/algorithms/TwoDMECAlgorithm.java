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
package mpg.molgen.sih.algorithms;

import java.util.List;

import mpg.molgen.sih.model.Block;
import mpg.molgen.sih.model.Fragment;

public class TwoDMECAlgorithm implements SIHAlgorithm {
	private boolean [] cut;
	private String haplotype="";
	@Override
	public void buildHaplotype(Block b) {
		cut = new boolean [b.getFragments().size()];
		initCut(b);
		for (int i=0;i<100;i++) {
			String newhap=CutHaplotypeTranslator.getHaplotype(b, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
			if(newhap.equals(haplotype)) {
				break;
			}
			haplotype = newhap;
			updateCut( b);
		}
		b.setHaplotype(haplotype);
	}
	private void initCut(Block b) {
		List<Fragment> fragments = b.getFragments();
		char [] hap1 = new char [b.length()];
		char [] hap2 = new char [b.length()];
		 
		for(int i=0;i<hap1.length;i++) {
			hap1[i] = hap2[i] = Fragment.NODATACHAR;
		}
		boolean [] assigned = new boolean [fragments.size()];
		for(int i=0;i<assigned.length;i++) {
			assigned[i] = false;
		}
		//Get the two fragments with larger hamming distance
		int maxD=0;
		int maxI=0;
		int maxJ=0;
		for(int i=0;i<fragments.size();i++) {
			Fragment f1 = fragments.get(i);
			for(int j=i+1;j<fragments.size();j++) {
				Fragment f2 = fragments.get(j);
				if(f2.getFirstPos()> f1.getLastPos()) {
					break;
				}
				int d = f1.getHammingDistance(f2);
				if(d>maxD) {
					maxD=d;
					maxI=i;
					maxJ=j;
				}
			}
		}
		Fragment f1 = fragments.get(maxI);
		Fragment f2 = fragments.get(maxJ);
		assigned[maxI] = assigned [maxJ] = true;
		cut[maxI] = false;
		cut[maxJ] = true;
		updateHaplotype(hap1, b.getFirstPos(), f1);
		updateHaplotype(hap2, b.getFirstPos(), f2);
		//Assign the other fragments
		for(int i=0;i<fragments.size()-2;i++) {
			int maxAbsDiff =0;
			maxJ = 0;
			int maxDiff =0;
			//Find the best next fragment
			for(int j=0;j<fragments.size();j++) {
				if(!assigned[j]) {
					Fragment f = fragments.get(j);
					int d1 = f.getHammingDistance(new String(hap1),b.getFirstPos());
					int d2 = f.getHammingDistance(new String(hap2),b.getFirstPos());
					int diff = d1 - d2;
					int absDiff = Math.abs(diff);
					if(absDiff > maxAbsDiff) {
						maxAbsDiff = absDiff;
						maxJ = j;
						maxDiff = diff;
					}
				}
			}
			Fragment f = fragments.get(maxJ);
			assigned[maxJ] = true;
			cut[maxJ] = (maxDiff==maxAbsDiff);
			if(cut[maxJ]) {
				updateHaplotype(hap2, b.getFirstPos(), f);
			} else {
				updateHaplotype(hap1, b.getFirstPos(), f);
			}
		}
	}
	private void updateHaplotype(char[] hap, int start, Fragment f) {
		
		String calls = f.getExtendedCalls();
		for(int i=0;i<calls.length();i++) {
			int relPos = f.getFirstPos() + i - start;
			if(hap[relPos]==Fragment.NODATACHAR && calls.charAt(i)!=Fragment.NODATACHAR) {
				hap[relPos] = calls.charAt(i);
			}
		}
		
	}
	private void updateCut(Block b) {
		List<Fragment> fragments = b.getFragments();
		for(int i=0;i<fragments.size();i++) {
			Fragment f = fragments.get(i);
			int distanceHap1 = f.getHammingDistance(haplotype,b.getFirstPos());
			int overlapping = f.getOverlappingCount(haplotype, b.getFirstPos());
			int distanceHap2 = overlapping - distanceHap1;
			if(distanceHap1!=distanceHap2) {
				cut[i] = distanceHap2 < distanceHap1;
			} else {
				distanceHap1 = f.getHamming2(haplotype,b.getFirstPos());
				distanceHap2 = -distanceHap1;
				if(distanceHap1 != distanceHap2) {
					cut[i] = distanceHap2 < distanceHap1;
				}
				
			}
		}
	}

}
