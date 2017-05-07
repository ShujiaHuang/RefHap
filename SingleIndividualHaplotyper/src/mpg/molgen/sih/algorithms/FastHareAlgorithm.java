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

import java.util.ArrayList;
import java.util.List;

import mpg.molgen.sih.model.Block;
import mpg.molgen.sih.model.Fragment;

public class FastHareAlgorithm implements SIHAlgorithm {

	@Override
	public void buildHaplotype(Block b) {
		List<Fragment> allFragments = b.getFragments();
		boolean cut [] = new boolean [allFragments.size()];
		List<Fragment> f1= new ArrayList<Fragment>();
		List<Fragment> f2= new ArrayList<Fragment>();
		int start1=0;
		int start2=0;
		for(int i=0;i<allFragments.size();i++) {
			Fragment f = allFragments.get(i);
			if(f1.size()==0) {
				f1.add(f);
				cut[i]=false;
				start1=f.getFirstPos();
			} else {
				//TODO: Do this more efficiently
				String hap1 = CutHaplotypeTranslator.getHaplotype(new Block(f1), cut, CutHaplotypeTranslator.CONSENSUS_ALL);
				int d1 = f.getHamming2(hap1, start1);
				String hap2 = "";
				int d2=0;
				if(f2.size()>0) {
					hap2 = CutHaplotypeTranslator.getHaplotype(new Block(f2), cut, CutHaplotypeTranslator.CONSENSUS_ALL);
					d2 = f.getHamming2(hap2, start2);
				}
				if(d1<d2) {
					f1.add(f);
					cut[i] = false;
				} else {
					f2.add(f);
					cut[i] = true;
					if(f2.size()==1) {
						start2 = f.getFirstPos();
					}
				}
			}
		}
		b.setHaplotype(CutHaplotypeTranslator.getHaplotype(b, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED));
	}
	
}