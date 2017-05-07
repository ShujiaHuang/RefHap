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
package mpg.molgen.sih.statistics;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import mpg.molgen.sih.model.Block;

public class BlockStatistics {
	private List<BlockItem> blockStats = new ArrayList<BlockItem>();
	private int totalVariants=0;
	private int totalPhased =0;
	private int totalSpan = 0;
	private double totalAdjustedSpan = 0;
	private int totalFragments = 0;
	private int totalCalls = 0;
	private int totalMEC=0;
	private int totalGSMEC=0;
	private int totalSwitchErrors=0;
	private int totalOverlap = 0;
	private DecimalFormat df = new DecimalFormat("##0.00");
	public void addVariants (List<Integer> variantPositions) {
		totalVariants += variantPositions.size();
	}
	public void addBlock (List<Integer> variantPositions, Block b, int se, int overlap, int gsMEC) {
		BlockItem item = createBlock(variantPositions, b.getFirstPos(), b.getLastPos(), b.getPhased());
		item.numFragments = b.getFragments().size();
		item.numCalls = b.getCalls();
		item.mec = b.getMEC();
		item.gsmec = gsMEC;
		item.switchErrors = se;
		item.overlap = overlap;
		totalFragments += item.numFragments;
		totalCalls += item.numCalls;
		totalMEC += item.mec;
		totalGSMEC += item.gsmec;
		totalSwitchErrors += item.switchErrors;
		totalOverlap += item.overlap;
	}
	public void addPartition(List<Integer> variantPositions, Block b, List<Integer> partition) {
		for(int i=0;i<partition.size();i+=2) {
			int firstPos = partition.get(i);
			int lastPos = partition.get(i+1);
			createBlock(variantPositions, firstPos,lastPos,b.getPhased(firstPos,lastPos));
		}
	}
	public BlockItem createBlock (List<Integer> variantPositions, int firstPos, int lastPos, int phased) {	
		int size = lastPos-firstPos+1;
		int firstVarPos = variantPositions.get(firstPos);
		int lastVarPos = variantPositions.get(lastPos);
		int span = lastVarPos - firstVarPos + 1;
		//System.out.println("Phased: "+phased +" firstpos: "+firstPos+" lastpos: "+lastPos);
		BlockItem item = new BlockItem(size, span, phased);
		blockStats.add(item);
		totalPhased+=phased;
		totalSpan += span;
		totalAdjustedSpan += item.adjustedSpan;
		return item;
	}
	
	/**
	 * @return the totalPhased
	 */
	public int getTotalPhased() {
		return totalPhased;
	}
	/**
	 * @return the totalSpan
	 */
	public int getTotalSpan() {
		return totalSpan;
	}
	/**
	 * @return the totalAdjustedSpan
	 */
	public double getTotalAdjustedSpan() {
		return totalAdjustedSpan;
	}
	public int getNumBlocks() {
		return blockStats.size();
	}
	public int getSizeN50() {
		List<BlockItem> sortedStats = new ArrayList<BlockItem>();
		sortedStats.addAll(blockStats);
		Collections.sort(sortedStats,new BlockItemPhasedComparator());
		int sum = 0;
		for(BlockItem item:sortedStats) {
			sum+=item.phased;
			if(sum > totalVariants/2) {
				return item.phased;
			}
		}
		return 0;
	}
	public int getSpanN50() {
		List<BlockItem> sortedStats = new ArrayList<BlockItem>();
		sortedStats.addAll(blockStats);
		Collections.sort(sortedStats,new BlockItemSpanComparator());
		int phased = 0;
		for(BlockItem item:sortedStats) {
			phased+=item.phased;
			if(phased > totalVariants/2) {
				return item.span;
			}
		}
		return 0;
	}
	public double getAdjustedSpanN50() {
		List<BlockItem> sortedStats = new ArrayList<BlockItem>();
		sortedStats.addAll(blockStats);
		Collections.sort(sortedStats,new BlockItemAdjustedSpanComparator());
		int phased = 0;
		for(BlockItem item:sortedStats) {
			phased+=item.phased;
			if(phased > totalVariants/2) {
				return item.adjustedSpan;
			}
		}
		return 0;
	}
	public void printSummaryStatistics(PrintStream out) {
		if(totalCalls> 0) {
			double avgFragCalls = (double)totalCalls/totalFragments;
			double avgCovVariants = (double)totalCalls/totalVariants;
			out.print (""+totalVariants+" "+totalFragments+" "+totalCalls+" "+df.format(avgFragCalls)+" "+df.format(avgCovVariants)+" "+blockStats.size());
			double mecPercentage = (double)totalMEC*100.0/totalCalls;
			double gsMECPercentage = (double)totalGSMEC*100.0/totalCalls;
			double sePercentage = (double)totalSwitchErrors*100.0/totalOverlap;
			out.print (" "+totalMEC+" "+df.format(mecPercentage)+" "+totalGSMEC+" "+df.format(gsMECPercentage)+" "+totalSwitchErrors+" "+df.format(sePercentage));
		}
		int phs50 = getSizeN50();
		int phn50 = getSpanN50();
		double phan50 = getAdjustedSpanN50();
		out.print (" "+phs50+" "+totalPhased);
		out.print (" "+phn50+" "+totalSpan);
		out.print (" "+df.format(phan50)+" "+df.format(totalAdjustedSpan));
	}
	public void printDetailStatistics(PrintStream out) {
		for(BlockItem item:blockStats) {
			out.print(item.numVariants+" "+item.span+" "+df.format(item.adjustedSpan)+" "+item.phased);
			if(item.numCalls> 0) {
				double avgFragCalls = (double)item.numCalls/item.numFragments;
				double avgCovVariants = (double)item.numCalls/item.numVariants;
				out.print (" "+item.numFragments+" "+item.numCalls+" "+item.mec+" "+item.gsmec+" "+item.switchErrors+" "+item.overlap+" "+df.format(avgFragCalls)+" "+df.format(avgCovVariants));
				double mecPercetage = (double)item.mec*100.0/item.numCalls;
				double gsmecPercetage = (double)item.gsmec*100.0/item.numCalls;
				double sePercentage = 0;
				if(item.overlap>0) sePercentage = (double)item.switchErrors*100.0/item.overlap;
				out.print(" "+df.format(mecPercetage)+" "+df.format(gsmecPercetage)+" "+df.format(sePercentage));
			}
			out.println();
		}
	}
}
class BlockItem {
	int numVariants;
	int span;
	double adjustedSpan;
	int phased;
	int numFragments;
	int numCalls;
	int mec;
	int gsmec;
	int switchErrors;
	int overlap;
	public BlockItem(int numVariants, int span, int phased) {
		super();
		this.numVariants = numVariants;
		this.span = span;
		this.phased = phased;
		this.adjustedSpan = (double)span*phased/numVariants;
	}
}
class BlockItemPhasedComparator implements Comparator <BlockItem>{

	@Override
	public int compare(BlockItem i1, BlockItem i2) {
		return i2.phased -i1.phased;
	}
}
class BlockItemSpanComparator implements Comparator <BlockItem>{

	@Override
	public int compare(BlockItem i1, BlockItem i2) {
		return i2.span -i1.span;
	}
}

class BlockItemAdjustedSpanComparator implements Comparator <BlockItem>{

	@Override
	public int compare(BlockItem i1, BlockItem i2) {
		if(i1.adjustedSpan > i2.adjustedSpan) return -1;
		else if (i1.adjustedSpan < i2.adjustedSpan) return 1;
		return 0;
	}
}


