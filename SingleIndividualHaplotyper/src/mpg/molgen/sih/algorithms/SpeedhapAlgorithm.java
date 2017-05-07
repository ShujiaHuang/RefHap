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
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import mpg.molgen.sih.model.Block;
import mpg.molgen.sih.model.Fragment;

public class SpeedhapAlgorithm implements SIHAlgorithm {
	private DataMatrix dataMatrix;
	private ConflictMatrix conflictMatrix;
	private static boolean test = false;
	@Override
	public void buildHaplotype(Block b) {
		dataMatrix = new DataMatrix(b);
		conflictMatrix = new ConflictMatrix(dataMatrix);
		//System.out.println("Initial rows: "+dataMatrix.nrows());
		//System.out.println("Initial cols: "+dataMatrix.ncols());
		preprocessing();
		if(test)System.out.println("Assigned cols preprocessing: "+countAssignedCols());
		int assignedCols = 0;
		for(int i=0;i<100 && assignedCols<dataMatrix.ncols();i++) {
			if(i>0) {
				conflictMatrix = new ConflictMatrix(dataMatrix);
			}
			phase1(b,i==0);
			if(test)System.out.println("Iteration: "+i+". Assigned cols phase1: "+countAssignedCols());
			phase2(b,i==0,true);
			assignedCols = countAssignedCols();
			if(test)System.out.println("Iteration: "+i+". Assigned cols phase2: "+assignedCols);
			
		}
		while(assignedCols < dataMatrix.ncols()) {
			phase2(b,false,false);
			assignedCols = countAssignedCols();
			if(test)System.out.println("Final assignment. Assigned cols: "+assignedCols);
		}
		String hap = buildHaplotype();
		b.setHaplotype(hap);
	}
	private String buildHaplotype() {
		char [] hap = new char [dataMatrix.ncols()];
		Arrays.fill(hap, Fragment.NODATACHAR);
		for(int j=0;j<dataMatrix.ncols();j++) {
			/*for(int i=0;i<dataMatrix.nrows();i++) {
				if(dataMatrix.superProfile[i][j] == 1) {
					hap[j]=dataMatrix.matrix[i][j];
				} else if (dataMatrix.superProfile[i][j] == 2) {
					if(dataMatrix.matrix[i][j] == Fragment.ALLELE1CHAR) {
						hap[j] = Fragment.ALLELE2CHAR;
					} else {
						hap[j] = Fragment.ALLELE1CHAR;
					}
				}
			}*/
			if(dataMatrix.assignedCols[j]==1) {
				hap[j] = Fragment.ALLELE1CHAR;
			} else if (dataMatrix.assignedCols[j]==2) {
				hap[j] = Fragment.ALLELE2CHAR;
			}
		}
		return new String (hap);
	}
	private void preprocessing() {
		//Fill out homozygous columns
		//fillHomozygous();
		boolean changes = true;
		while(changes) {
			changes = false;
			int indexMaxConflicts = conflictMatrix.getIndexMaxConflicts();
			if(indexMaxConflicts>=0) {
				changes = removeConflict(indexMaxConflicts);
			}
			if(changes) conflictMatrix = new ConflictMatrix(dataMatrix);
			else {
				//Fix columns with conflicts no matter if they are not the most conflicting ones
				for(int i=0;i<conflictMatrix.conflictCols.length && !changes;i++) {
					if(conflictMatrix.conflictCols[i].size()>1) {
						changes = removeConflict(i);
					}
				}
			}
			if(changes) conflictMatrix = new ConflictMatrix(dataMatrix);
		}
		//Remove remaining conflicts but keep conflict matrix
		int nconflicts=0;
		for(int i=0;i<conflictMatrix.matrix.length;i++) {
			for(int j=0;j<i;j++) {
				ConflictCell cell = conflictMatrix.matrix[i][j];
				if(cell !=null) {
					IndexList conflictRows = cell.getConflictRows();
					List<Integer> rowsToRemove = conflictRows.getData();
					for(int k:rowsToRemove) {
						dataMatrix.matrix[k][i] = Fragment.NODATACHAR;
						dataMatrix.matrix[k][j] = Fragment.NODATACHAR;
						nconflicts+=2;
					}
				}
			}
		}
		if(test && nconflicts>0)System.out.println("Final Preprocessing. Removed: "+nconflicts+" conflicting entries");
	}
	private boolean removeConflict(int indexConflict) {
		
		List<Integer> conflicts = conflictMatrix.conflictCols[indexConflict].getData();
		List<Integer> rowsToRemove = new ArrayList<Integer>();
		for(Integer i:conflicts) {
			if(conflictMatrix.conflictCols[i].size()==1) {
				rowsToRemove.addAll(conflictMatrix.matrix[indexConflict][i].getConflictRows().getData());
			}
		}
		boolean changes = rowsToRemove.size()>0;
		if(changes) {
			int nconflictsR=0;
			for(int i:rowsToRemove) {
				dataMatrix.matrix[i][indexConflict] = Fragment.NODATACHAR;
				nconflictsR++;
			}
			if(test) System.out.println("Preprocessing. Removed: "+nconflictsR+" conflicting entries from column: "+indexConflict);
			//
		}
		return changes;
	}
	
	/*public void fillHomozygous() {
		char [][] matrix = dataMatrix.matrix;
		for(int j=0;j<matrix[0].length;j++) {
			int count1=0;
			int count2=0;
			for(int i=0;i<matrix.length;i++) {
				if(matrix[i][j]==Fragment.ALLELE1CHAR) {
					count1++;
				} else if(matrix[i][j]==Fragment.ALLELE2CHAR) {
					count2++;
				} 
			}
			if(count1>0 && count2 == 0) {
				dataMatrix.hap[j] = Fragment.ALLELE1CHAR; 
			} else if(count1==0 && count2 > 0) {
				dataMatrix.hap[j] = Fragment.ALLELE2CHAR;
			}
		}
	}*/
	private void phase1(Block b, boolean firstTime) {
		conflictMatrix.calculateConnectedCols(firstTime);
		int selId = -1;
		short h = 1;
		if(firstTime) {
			int [] groupsMap = buildGroups();
			int [] groupSizes = new int [groupsMap.length];
			Arrays.fill(groupSizes, 0);
			for (int i=0;i< groupsMap.length;i++) {
				int gId = groupsMap[i];
				if(gId>=0) groupSizes[gId]++;
				else assert dataMatrix.assignedCols[i]!=0;
			}
			int maxSize =0;
			selId = 0;
			for(int i=0;i<groupSizes.length;i++) {
				if(groupSizes[i]>maxSize) {
					maxSize=groupSizes[i];
					selId = i;
				}
			}
			
		} else {
			int maxOverlap = -1;
			
			ConflictArray ca = new ConflictArray(dataMatrix);
			for(int i=0;i<ca.array.length;i++) {
				if(dataMatrix.assignedCols[i] == 0) {
					ConflictCell cell = ca.array[i];
					int ov1 = cell.rows11.size() + cell.rows22.size();
					int ov2 = cell.rows12.size() + cell.rows21.size();
					int ovT = Math.max(ov1, ov2);
					if(ovT>maxOverlap) {
						maxOverlap = ovT;
						selId = i;
						if(ovT == ov1) h= 1;
						else h= 2;
						if(test) System.out.println("New max column: "+i+". Overlap1: "+ov1+" Overlap 2: "+ov2+" Haplotype: "+h);
					}
				}
			}
		}
		
		assignColumn(selId,h);
		//Build super profiles
		LinkedList<Integer> stack = new LinkedList<Integer>();
		
		stack.push(selId);
		while (stack.size()>0) {
			int index = stack.pop();
			short a1 = dataMatrix.assignedCols[index];
			for(int i:conflictMatrix.connectedCols[index].data) {
				ConflictCell cell = conflictMatrix.matrix[index][i];
				if(dataMatrix.assignedCols[i]==0) {	
					stack.push(i);
					if(cell.rows11.size()>0 || cell.rows22.size()>0) {
						assignColumn(i,a1);
					} else if(a1==1) {
						assignColumn(i,(short)2);
					} else {
						assignColumn(i,(short)1);
					}
				} /*else {
					short a2 = 1;
					if(a1==a2) a2= 2;
					assignEntries(cell.rows11,i,a1);
					assignEntries(cell.rows22,i,a2);
					assignEntries(cell.rows12,i,a1);
					assignEntries(cell.rows21,i,a2);
				}*/
			}
		}
		//Assign fragments
		assignFragments();
	}
	private void assignFragments() {
		for(int i=0;i<dataMatrix.nrows();i++) {
			if(dataMatrix.rowsAssignment[i] == 0) {
				int score = calculateScore(i);
				if(score<0) {
					dataMatrix.rowsAssignment[i] = 1;
				} else if(score > 0){
					dataMatrix.rowsAssignment[i] = 2;
				}
			}
		}
	}
	private int calculateScore(int i) {
		int score =0;
		for(int j=0;j<dataMatrix.ncols();j++) {
			if (dataMatrix.superProfile[i][j]==1) {
				score--;
			} else if (dataMatrix.superProfile[i][j]==2) {
				score++;
			}
		}
		return score;
	}
	public void assignEntries(IndexList rows, int j, short h) {
		for (int i:rows.data) {
			char call = dataMatrix.matrix[i][j];
			if(call!=Fragment.NODATACHAR) {
				if(dataMatrix.superProfile[i][j]==0)dataMatrix.superProfile[i][j]=h;
				else if(dataMatrix.superProfile[i][j]!=h) {
					dataMatrix.superProfile[i][j] = 3;
					if(test) System.out.println("Conflict found in cell: "+i+","+j);
				}
			}
		}
	}
	private void assignColumn(int j, short h) {
		for (int i=0;i<dataMatrix.nrows();i++) {
			char call = dataMatrix.matrix[i][j];
			if(call!=Fragment.NODATACHAR) {
				short cellHap = h ;
				if(call == Fragment.ALLELE2CHAR) {
					//Switch haplotype membership
					if(cellHap == 2) cellHap = 1;
					else cellHap = 2;
				}
				if(dataMatrix.superProfile[i][j]==0)dataMatrix.superProfile[i][j]=cellHap;
				else if(dataMatrix.superProfile[i][j]!=cellHap) {
					dataMatrix.superProfile[i][j] = 3;
					if(test) System.out.println("Conflict found in cell: "+i+","+j);
				}
			}
		}
		if(dataMatrix.assignedCols[j]==0) dataMatrix.assignedCols[j] = h; 
	}
	
	private int[] buildGroups() {
		int [] groupsMap = new int [conflictMatrix.connectedCols.length];
		Arrays.fill(groupsMap, -1);
		Queue<Integer> queue = new LinkedList<Integer>();
		while (true) {
			
			int nextIndex = -1;
			for(int i=0;i<groupsMap.length;i++) {
				if(dataMatrix.assignedCols[i]==0 && groupsMap[i]==-1) {
					nextIndex = i;
					break;
				}
			}
			if(nextIndex == -1) {
				break;
			}
			groupsMap[nextIndex] = nextIndex;
			queue.clear();
			queue.add(nextIndex);
			while (queue.size()>0) {
				int index = queue.remove();	
				for(Integer i:conflictMatrix.connectedCols[index].data) {
					if(groupsMap[i]==-1) {
						groupsMap[i] = groupsMap[index];
						queue.add(i);
					} else {
						assert groupsMap[i] == groupsMap[index]; 
					}
				}
			}
		}
		return groupsMap;
	}
	private void phase2(Block b, boolean checkFullRank, boolean checkConflict) {
		boolean change = true;
		while(change) {
			change = false;
			ConflictArray ca = new ConflictArray(dataMatrix);
			for(int i=0;i<ca.array.length;i++) {
				if(dataMatrix.assignedCols[i] == 0) {
					ConflictCell cell = ca.array[i];
					if(!checkConflict || !cell.isEmpty()) {
						if(!checkConflict || !cell.hasConflict() ) {
							if(!checkFullRank || cell.isFullRank()) {
								if(cell.rows11.size()+cell.rows22.size()>cell.rows12.size()+cell.rows21.size()) {
									assignColumn(i,(short)1);
								} else {
									assignColumn(i,(short)2);
								}
								change = true;
							}
						}
					}
					
				}
			}
			//Assign fragments 
			if (change) assignFragments();
			//Correct errors
			for(int i=0;i<ca.array.length;i++) {
				if(dataMatrix.assignedCols[i] == 0) {
					ConflictCell cell = ca.array[i];
					List<Integer> conflicts = cell.getConflictRows().getData();
					change = change || conflicts.size()>0;
					for(int k:conflicts) {
						dataMatrix.matrix[k][i] = Fragment.NODATACHAR;
					}
					if(test && conflicts.size()>0) System.out.println("Phase 2. Fixed "+conflicts.size()+" conflicts in column: "+i);	
				}
			}
		}
	}
	private int countAssignedCols() {
		int count=0;
		for(int i=0;i<dataMatrix.assignedCols.length;i++) {
			if (dataMatrix.assignedCols[i]!=0) {
				count++;
			}
		}
		return count;
	}
}
class DataMatrix {
	short [] rowsAssignment; 
	//char [] hap;
	char [][] matrix;
	short [][] superProfile;
	short [] assignedCols;
	public DataMatrix (Block b) {
		List<Fragment> fragments = b.getFragments();
		matrix = new char [fragments.size()][b.length()];
		superProfile = new short [fragments.size()][b.length()];
		rowsAssignment = new short [matrix.length];
		Arrays.fill(rowsAssignment, (short)0);
		assignedCols = new short [b.length()];
		Arrays.fill(assignedCols, (short)-1);
		for(int i=0;i<matrix.length;i++) {
			Fragment f = fragments.get(i);
			String calls = f.getExtendedCalls();
			int relStart = f.getFirstPos() - b.getFirstPos();
			int relEnd = relStart + calls.length() - 1;
			for(int j=0;j<matrix[0].length;j++) {
				if(j>=relStart && j<=relEnd) {
					matrix[i][j] = calls.charAt(j-relStart);
					if(matrix[i][j]!=Fragment.NODATACHAR) {
						assignedCols[j] = 0;
					}
				} else {
					matrix[i][j] = Fragment.NODATACHAR;
				}
				superProfile[i][j] = 0;
			}
		}
	}
	public int nrows() {
		return rowsAssignment.length;
	}
	public int ncols () {
		return assignedCols.length; 
	}
}
class ConflictMatrix {
	ConflictCell matrix [][];
	IndexList [] conflictCols;
	IndexList [] connectedCols;
	DataMatrix dM;
	
	public ConflictMatrix(DataMatrix dM) {
		this.dM = dM;
		matrix = new ConflictCell[dM.ncols()][dM.ncols()];
		conflictCols = new IndexList [dM.ncols()];
		connectedCols = new IndexList [dM.ncols()];
		
		for(int i=0;i<matrix.length;i++) {
			conflictCols[i] = new IndexList();
			connectedCols[i] = new IndexList();
			ConflictCell c = new ConflictCell(dM, i, i);
			matrix[i][i] = c;
			if(dM.assignedCols[i]==0) {
				for(int j=0;j<i;j++) {
					if(dM.assignedCols[j]==0) {
						c = new ConflictCell(dM, i, j);
						matrix[i][j]=matrix[j][i]=c;
						if(c.hasConflict()) {
							conflictCols[i].add(j);
							conflictCols[j].add(i);
						}
					}
				}
			}
		}
		for(int i=0;i<conflictCols.length;i++) {
			conflictCols[i].sort();
		}
	}
	public void calculateConnectedCols (boolean checkFullRank) {
		int average = (int)calculateAverageConflictDegree()+1;
		//System.out.println("Averege conflict degree: "+average);
		for(int i=0;i<matrix.length;i++) {
			if(dM.assignedCols[i]==0 && conflictCols[i].size() <= average) {
				for(int j=0;j<i;j++) {
					if(dM.assignedCols[j]==0 &&conflictCols[j].size() <= average) {
						ConflictCell c = matrix[i][j];
						if(!c.isEmpty() && !c.hasConflict() && (!checkFullRank ||c.isFullRank()) ) {
							connectedCols[i].add(j);
							connectedCols[j].add(i);
						}
						//else if (c.hasConflict())System.out.println("Columns: "+i+" and "+j+" have conflict");
						//else if (!c.isFullRank())System.out.println("Columns: "+i+" and "+j+" are not full rank");
					}
				}
			}
		}
		for(int i=0;i<connectedCols.length;i++) {
			connectedCols[i].sort();
		}
	}
	private double calculateAverageConflictDegree() {
		double sum=0;
		int count =0;
		for(int i=0;i<conflictCols.length;i++) {
			if (conflictCols[i].size()>0) {
				sum += conflictCols[i].size();
				count++;
			}
		}
		if(count>0) return sum/count;
		return 0;
	}
	public int getIndexMaxConflicts() {
		int maxConflicts = 1;
		int maxIndex = -1;
		for(int i=0;i<conflictCols.length;i++) {
			if(conflictCols[i].size()>maxConflicts) {
				maxConflicts = conflictCols[i].size();
				maxIndex = i;
			}
		}
		return maxIndex;
	}
}
class ConflictArray {
	ConflictCell [] array;
	public ConflictArray (DataMatrix dM) {
		array = new ConflictCell [dM.ncols()];
		for(int i=0;i<array.length;i++) {
			if(dM.assignedCols[i]==0) {
				array[i] = new ConflictCell(dM, i);
			}
		}
	}
}
class ConflictCell {
	IndexList rows11=new IndexList();
	IndexList rows12=new IndexList();
	IndexList rows21=new IndexList();
	IndexList rows22=new IndexList();
	public ConflictCell (DataMatrix dm, int j1, int j2) {
		char [][] matrix = dm.matrix;
		for(int i=0;i<dm.nrows();i++) {
			if(matrix[i][j1]!=Fragment.NODATACHAR && matrix[i][j2]!=Fragment.NODATACHAR) {
				if(matrix[i][j1]==Fragment.ALLELE1CHAR && matrix[i][j2]==Fragment.ALLELE1CHAR) {
					rows11.add(i);
				} else if(matrix[i][j1]==Fragment.ALLELE1CHAR && matrix[i][j2]==Fragment.ALLELE2CHAR) {
					rows12.add(i);
				} else if(matrix[i][j1]==Fragment.ALLELE2CHAR && matrix[i][j2]==Fragment.ALLELE1CHAR) {
					rows21.add(i);
				} else if(matrix[i][j1]==Fragment.ALLELE2CHAR && matrix[i][j2]==Fragment.ALLELE2CHAR) {
					rows22.add(i);
				}
			}
		}
	}
	public ConflictCell (DataMatrix dm, int j) {
		char [][] matrix = dm.matrix;
		for(int i=0;i<dm.nrows();i++) {
			if(dm.rowsAssignment[i]!=0 && matrix[i][j]!=Fragment.NODATACHAR) {
				if(dm.rowsAssignment[i]==1 && matrix[i][j]==Fragment.ALLELE1CHAR) {
					rows11.add(i);
				} else if(dm.rowsAssignment[i]==1 && matrix[i][j]==Fragment.ALLELE2CHAR) {
					rows12.add(i);
				} else if(dm.rowsAssignment[i]==2 && matrix[i][j]==Fragment.ALLELE1CHAR) {
					rows21.add(i);
				} else if(dm.rowsAssignment[i]==2 && matrix[i][j]==Fragment.ALLELE2CHAR) {
					rows22.add(i);
				}
			}
		}
	}
	public boolean hasConflict () {
		int count=0;
		if(rows11.size()>0) {
			count++;
		}
		if(rows12.size()>0) {
			count++;
		}
		if(rows21.size()>0) {
			count++;
		}
		if(rows22.size()>0) {
			count++;
		}
		return count > 2;
	}
	public boolean isFullRank () {
		if(rows11.size()>0 && rows22.size()>0) {
			return true;
		}
		if(rows12.size()>0 && rows21.size()>0) {
			return true;
		}
		return false;
	}
	public boolean isEmpty() {
		return (rows11.size()+rows22.size()+rows21.size()+rows12.size()==0);
	}
	public IndexList getConflictRows() {
		if(rows11.size()>0 && rows12.size()>0 && rows21.size()>0 && rows22.size()>0) {
			return new IndexList();
		}
		if(!isFullRank()) {
			return new IndexList();
		}
		if(rows11.size()>0 && rows22.size()>0) {
			if(rows12.size() > 0) {
				return rows12;
			} else {
				return rows21;
			}
		} else {
			if(rows11.size() > 0) {
				return rows11;
			} else {
				return rows22;
			}
		}
	}
}
class IndexList {
	List<Integer> data = new ArrayList<Integer>();
	public void add(int index) {
		data.add(index);
	}
	public void add(IndexList i2) {
		data.addAll(i2.data);
	}
	public List<Integer> getData() {
		return data;
	}
	public int size () {
		return data.size();
	}
	public void sort() {
		Collections.sort(data);
	}
}
