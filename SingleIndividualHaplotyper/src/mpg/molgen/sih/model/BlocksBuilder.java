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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class BlocksBuilder {
	private int [] assignments;
	private int lastBlock=0;
	private List<Join>joins=new ArrayList<Join>(); 
	public List<Block> buildBlocks(Block bigBlock) {
		fillAssignments(bigBlock);
		joinBlocks();
		return distributeFragments(bigBlock.getFragments());
	}
	private void fillAssignments(Block bigBlock) {
		assignments = new int [bigBlock.getLastPos()-bigBlock.getFirstPos()+1];
		Arrays.fill(assignments, 0);
		List<Fragment> fragments = bigBlock.getFragments();
		int firstPos = bigBlock.getFirstPos();
		for(Fragment f:fragments) {
			int fragBlock = 0;
			String calls = f.getExtendedCalls();
			int relPos = f.getFirstPos()-firstPos;
			for(int k=0;k<calls.length();k++) {
				if(calls.charAt(k)!=Fragment.NODATACHAR) {
					int block = assignments[relPos+k];
					if(block > 0) {
						if(fragBlock == 0) {
							fragBlock = block;
						} else if(block > fragBlock) {
							joins.add(new Join(block,fragBlock));
						} else if(fragBlock > block) {
							joins.add(new Join(fragBlock,block));
							fragBlock = block;
						}
					}
				}
			}
			if(fragBlock==0) {
				lastBlock++;
				fragBlock = lastBlock;
			}
			for(int k=0;k<calls.length();k++) {
				if(calls.charAt(k)!=Fragment.NODATACHAR) {
					assignments[relPos+k] = fragBlock;
				}
			}
		}
		
	}
	private void joinBlocks() {
		Collections.sort(joins,new JoinComparator());
		int lastGroup = 0;
		int newGroup = 0; 
		for(Join join:joins ) {
			if(lastGroup != join.group1) {
				newGroup = join.group1;
			}
			if(newGroup > join.group2) {
				for(int j=0;j<assignments.length;j++) {
					if(assignments[j]==newGroup) {
						assignments[j] = join.group2;
					}
				}
				newGroup = join.group2;
			}
			
			lastGroup = join.group1;
		}
	}
	private List<Block> distributeFragments(List<Fragment> fragments) {
		List<List<Fragment>> distFragments = new ArrayList<List<Fragment>>(lastBlock);
		for(int i=0;i<lastBlock;i++) {
			distFragments.add(new ArrayList<Fragment>());
		}
		int firstPos = fragments.get(0).getFirstPos();
		for(Fragment f:fragments) {
			int pos = assignments[f.getFirstPos()-firstPos]-1;
			distFragments.get(pos).add(f);
		}
		List<Block> blocks = new ArrayList<Block>();
		for(List<Fragment> l:distFragments) {
			if(l.size()>0) {
				blocks.add(new Block(l));
			}
		}
		return blocks;
	}
}
class Join {
	int group1;
	int group2;
	public Join(int group1, int group2) {
		super();
		this.group1 = group1;
		this.group2 = group2;
	}
	public int getGroup1() {
		return group1;
	}
	public int getGroup2() {
		return group2;
	}
	
	
}
class JoinComparator implements Comparator<Join> {

	@Override
	public int compare(Join j1, Join j2) {
		if(j1.group1!=j2.group1) {
			return j2.group1-j1.group1;
		}
		return j2.group2-j1.group2;
	}
	
}
