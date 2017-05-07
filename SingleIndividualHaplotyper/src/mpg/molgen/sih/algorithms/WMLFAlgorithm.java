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

public class WMLFAlgorithm implements SIHAlgorithm {
	private boolean [] cut;
	private String haplotype="";
	private Random random = new Random();
	@Override
	public void buildHaplotype(Block b) {
		cut = new boolean [b.getFragments().size()];
		int bestMEC = -1;
		String bestHaplotype = "";
		for (int iter=0;iter<100;iter++) {
			initCut(b);
			haplotype="";
			for (int i=0;i<100;i++) {
				//TODO: Make get weighted haplotype
				String newhap=CutHaplotypeTranslator.getHaplotype(b, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
				if(newhap.equals(haplotype)) {
					break;
				}
				haplotype = newhap;
				updateCut( b);
			}
			if(bestMEC == -1 || b.getMEC()< bestMEC) {
				bestMEC = b.getMEC();
				bestHaplotype = haplotype;
			}
			
		}
		b.setHaplotype(bestHaplotype);
	}
	private void updateCut(Block b) {
		List<Fragment> fragments = b.getFragments();
		for(int i=0;i<fragments.size();i++) {
			Fragment f = fragments.get(i);
			//TODO: Weighted distance
			int distanceHap1 = f.getHammingDistance(haplotype,b.getFirstPos());
			int overlapping = f.getOverlappingCount(haplotype, b.getFirstPos());
			int distanceHap2 = overlapping - distanceHap1;
			if(distanceHap2 != distanceHap1) {
				cut[i] = distanceHap2 < distanceHap1;
			}
		}
		
	}
	private void initCut(Block b) {
		for(int i=0;i<cut.length;i++) {
			cut[i] = random.nextBoolean();
		}
	}
}
