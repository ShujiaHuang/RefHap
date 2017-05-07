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

public class DGSAlgorithm implements SIHAlgorithm {
	private boolean [] cut;
	private String haplotype="";
	@Override
	public void buildHaplotype(Block b) {
		cut = new boolean [b.getFragments().size()];
		initCut(b);
		for (int i=0;i<1000;i++) {
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
		char [] hap = new char [b.length()];
		 
		for(int i=0;i<hap.length;i++) {
			hap[i] = Fragment.NODATACHAR;
		}
		boolean [] assigned = new boolean [fragments.size()];
		for(int i=0;i<assigned.length;i++) {
			assigned[i] = false;
		}
		//Get the fragment with larger number of calls
		int maxC=0;
		int maxI=0;
		for(int i=0;i<fragments.size();i++) {
			int calls = fragments.get(i).getnCalls();
			if(calls>maxC) {
				maxC=calls;
				maxI=i;
			}
		}
		Fragment f = fragments.get(maxI);
		assigned[maxI] = true;
		cut[maxI] = false;
		updateHaplotype(hap, b.getFirstPos(), f,false);
		//Assign the other fragments
		for(int i=0;i<fragments.size()-1;i++) {
			int maxAbsScore =0;
			int maxJ = 0;
			int maxScore =0;
			//Find the best next fragment
			for(int j=0;j<fragments.size();j++) {
				if(!assigned[j]) {
					f = fragments.get(j);
					int score = f.getHamming2(new String(hap),b.getFirstPos());
					int absScore = Math.abs(score);
					if(absScore > maxAbsScore) {
						maxAbsScore = absScore;
						maxJ = j;
						maxScore = score;
					}
				}
			}
			f = fragments.get(maxJ);
			assigned[maxJ] = true;
			cut[maxJ] = (maxScore==maxAbsScore);
			updateHaplotype(hap, b.getFirstPos(), f, cut[maxJ]);
		}
	}
	private void updateHaplotype(char[] hap, int start, Fragment f, boolean reverse) {
		String calls = f.getExtendedCalls();
		for(int i=0;i<calls.length();i++) {
			int relPos = f.getFirstPos() + i - start;
			if(hap[relPos]==Fragment.NODATACHAR && calls.charAt(i)!=Fragment.NODATACHAR) {
				if(!reverse) {
					hap[relPos] = calls.charAt(i);
				} else if (calls.charAt(i) == Fragment.ALLELE1CHAR) {
					hap[relPos] = Fragment.ALLELE2CHAR;
				} else {
					hap[relPos] = Fragment.ALLELE1CHAR;
				}
				
			}
		}
		
	}
	private void updateCut(Block b) {
		List<Fragment> fragments = b.getFragments();
		for(int i=0;i<fragments.size();i++) {
			Fragment f = fragments.get(i);
			int score = f.getHamming2(haplotype,b.getFirstPos());
			if(score != 0) {
				cut[i] = score > 0;
			}
			
		}
	}

}
