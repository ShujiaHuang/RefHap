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
import java.util.Random;

import mpg.molgen.sih.model.Block;
import mpg.molgen.sih.model.Fragment;

public class SHRThreeAlgorithm implements SIHAlgorithm {

	private Random random = new Random();
	@Override
	public void buildHaplotype(Block b) {
		int minMEC = -1;
		boolean [] minCut = null;
		for(int iter=0;iter < 10;iter++) {
			boolean [] cut = generateCut(b);
			int mec1 = calculateMec(b,cut,false);
			int mec2 = calculateMec(b,cut,true);
			int mec = Math.max(mec1, mec2);
			if(minMEC == -1 || mec < minMEC) {
				minMEC = mec;
				minCut = cut;
			}
		}
		b.setHaplotype(CutHaplotypeTranslator.getHaplotype(b, minCut, CutHaplotypeTranslator.CONSENSUS_GROUP_1));
	}
	private boolean[] generateCut(Block b) {
		
		List<Fragment> fragments = b.getFragments();
		boolean [] cut = new boolean [fragments.size()];
		if(fragments.size()==1) {
			return cut;
		}
		int index1 = random.nextInt(cut.length);
		int index2 = index1;
		while(index2==index1) {
			index2 = random.nextInt(cut.length);
		}
		Fragment f1 = fragments.get(index1);
		Fragment f2 = fragments.get(index2);
		cut[index1] = false;
		cut[index2] = true;
		for(int i=0;i<fragments.size();i++) {
			Fragment f = fragments.get(i);
			if(i!=index1 && i!= index2) {
				int d1 = f.getHammingDistance(f1);
				int d2 = f.getHammingDistance(f2);
				if(d1 < d2) {
					cut[i] = false;
				} else if (d1 > d2) {
					cut[i] = true;
				} else {
					cut[i] = random.nextBoolean();
				}
			}
		}
		return cut;
	}
	private int calculateMec(Block b, boolean[] cut, boolean subset) {
		int mec = 0; 
		List<Fragment> fragments = b.getFragments();
		for(int i=0;i<fragments.size();i++) {
			if(cut[i] == subset) {
				Fragment f1 = fragments.get(i);
				for(int j=i+1;j<fragments.size();j++) {
					if(cut[j] == subset) {
						Fragment f2 = fragments.get(j);
						if(f2.getFirstPos()> f1.getLastPos()) {
							break;
						}
						mec += f1.getHammingDistance(f2);
					}
				}
			}
		}
		return mec;
	}

}
