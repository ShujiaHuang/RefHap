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

import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import mpg.molgen.sih.algorithms.FragmentsCutBuilder;
import mpg.molgen.sih.model.Block;
import mpg.molgen.sih.model.BlocksBuilder;
import mpg.molgen.sih.model.Fragment;
import mpg.molgen.sih.model.FragmentsFileHandler;
import mpg.molgen.sih.model.HaplotypesFileHandler;

public class PhasingStatisticsCalculator {
	private BlockStatistics initialStats = new BlockStatistics();
	private BlockStatistics phasingBlockStats = new BlockStatistics();
	private boolean printSEFragments = false;
	
	/*public String getHaplotype (List<Block> phasedBlocks, int length) {
		StringBuilder builder = new StringBuilder();
		int pos = 0;
		for(Block b:phasedBlocks) {
			for(;pos<b.getFirstPos();pos++) {
				builder.append(Fragment.NODATACHAR);
			}
			//In overlapping fragments this is not already set after the previous cycle
			pos = b.getFirstPos();
			String blockHap = b.getHaplotype();
			for(int i=0;i<blockHap.length();i++) {
				char allele = blockHap.charAt(i);
				if(pos<builder.length()) {
					if(builder.charAt(pos)== Fragment.NODATACHAR) {
						builder.setCharAt(pos, allele);
					}
				} else {
					builder.append(allele);
				} 
				pos++;
			}
			
		}
		while(builder.length()<length) {
			builder.append(Fragment.NODATACHAR);
		}
		return builder.toString();
	}*/
	public int calculateSE(String testHaplotype, String realHaplotype) {
		Fragment dummyFrag = new Fragment("", 0, testHaplotype);
		List<Fragment> dummyList = new ArrayList<Fragment>();
		dummyList.add(dummyFrag);
		Block dummyBlock = new Block(dummyList);
		dummyBlock.setHaplotype(testHaplotype);
		return calculateSE(calculateSEPartition(dummyBlock, realHaplotype));  
	}
	public int calculateSE(List<Integer> partition) {
		return Math.max(0, partition.size()/2-1);
	}
	public int calculateOverlapping(String testHap, String realHap) {
		int answer = 0;
		for(int i=0;i<testHap.length()&&i<realHap.length();i++) {
			char alleleHap = testHap.charAt(i);
			char actReal = realHap.charAt(i);
			if(alleleHap!=Fragment.NODATACHAR && actReal!=Fragment.NODATACHAR) {
				answer++;
			}
		}
		return answer;
	}
	public List<Integer> calculateSEPartition(Block b, String realHap) {
		List<Integer> partition = new ArrayList<Integer>();
		
		String haplotype = b.getHaplotype();
		int lastPos = -1;
		for(int i=b.getFirstPos();i<=b.getLastPos();i++) {
			char actReal = realHap.charAt(i);
			char alleleHap = haplotype.charAt(i-b.getFirstPos());
			if(alleleHap!=Fragment.NODATACHAR && actReal!=Fragment.NODATACHAR) {
				if(lastPos>=0) {
					char lastReal = realHap.charAt(lastPos);
					char lastHap = haplotype.charAt(lastPos-b.getFirstPos());
					if(lastReal == actReal && lastHap != alleleHap) {
						partition.add(lastPos);
						partition.add(i);
					} else if (lastReal != actReal && lastHap == alleleHap) {
						partition.add(lastPos);
						partition.add(i);
					}
				} else {
					//First phased pos
					partition.add(i);
				}
				lastPos = i;
			}
		}
		if(lastPos!=-1) partition.add(lastPos);
		return partition;
	}
	public void calculateStatistics(String testHaplotype, List<Integer> variantPositions, String realHaplotype) {
		Fragment dummyFrag = new Fragment("", 0, testHaplotype);
		List<Fragment> dummyList = new ArrayList<Fragment>();
		dummyList.add(dummyFrag);
		Block dummyBlock = new Block(dummyList);
		dummyBlock.setHaplotype(testHaplotype);
		List<Block> dummyBlockList = new ArrayList<Block>();
		dummyBlockList.add(dummyBlock);
		calculateStatistics(dummyBlockList,variantPositions,realHaplotype);
	}
	public String calculateStatistics(List<Block> blocks, List<Integer> variantPositions, String realHaplotype) {
		StringBuilder outputHaplotype = new StringBuilder(realHaplotype);
		initialStats.addVariants(variantPositions);
		phasingBlockStats.addVariants(variantPositions);
		for(Block b: blocks) {
			List<Integer> partition = calculateSEPartition(b, realHaplotype);
			fillOutputHaplotype(b,partition,outputHaplotype);
			if(printSEFragments) printSEFragments (b,variantPositions, partition,realHaplotype,System.out);
			Block b2 = new Block(b.getFragments());
			String realSubhap = realHaplotype.substring(b2.getFirstPos(), b2.getLastPos()+1);
			b2.setHaplotype(realSubhap);
			int gsMEC =b2.getMEC();
			int overlap = calculateOverlapping(b.getHaplotype(), realSubhap);
			initialStats.addBlock(variantPositions, b, calculateSE(partition),overlap, gsMEC);
			phasingBlockStats.addPartition(variantPositions,b, partition);
		} 
		return outputHaplotype.toString();
	}
	
	private void fillOutputHaplotype(Block b, List<Integer> partition,StringBuilder outputHaplotype) {
		String haplotype = b.getHaplotype();
		int prevLast = b.getFirstPos()-1;
		boolean prevComp = false;
		for(int i=0;i<partition.size();i+=2) {
			int first = partition.get(i);
			int last = partition.get(i+1);
			boolean complement = outputHaplotype.charAt(first)!=haplotype.charAt(first-b.getFirstPos());
			if(i==0) {
				first = b.getFirstPos();
			} else {
				if(prevComp == complement) System.err.println("Inconsistent switch error at: "+(first+1));
				if(prevLast+1<first) {
					first = findPosMinScore(b,prevLast+1,first-1);
					fillHaplotype(b, outputHaplotype, haplotype, prevLast+1, first-1, prevComp);
				}
			}
			if(i==partition.size()-1) {
				last = b.getLastPos();
			}
			fillHaplotype(b, outputHaplotype, haplotype, first, last, complement);
			prevLast = last;
			prevComp = complement;
		}
	}
	private void fillHaplotype(Block b, StringBuilder outputHaplotype,String haplotype, int first, int last, boolean complement) {
		for(int j=first;j<=last;j++) {
			char call = haplotype.charAt(j-b.getFirstPos());
			if(call!=Fragment.NODATACHAR) {
				if(complement) {
					if(call == Fragment.ALLELE1CHAR) call = Fragment.ALLELE2CHAR;
					else call = Fragment.ALLELE1CHAR;
				}
				if(outputHaplotype.charAt(j)==Fragment.NODATACHAR) {
					outputHaplotype.setCharAt(j, call);
				} else if(outputHaplotype.charAt(j) != call){
					System.err.println("Error assembling output haplotypes at position: "+(j+1)+". Real call: "+outputHaplotype.charAt(j)+". Hap call: "+call);
				}
			}
		}
	}
	private int findPosMinScore(Block b, int first, int last) {
		double[] hapScores = b.getHapScores();
		double minScore = hapScores[first - b.getFirstPos()];
		int posMinScore = first;
		for(int i=first+1;i<=last;i++) {
			int relPos = i - b.getFirstPos();
			if(hapScores[relPos]<minScore) {
				minScore = hapScores[relPos];
				posMinScore = i;
			}
		}
		return posMinScore;
	}
	public void printSEFragments(Block b, List<Integer> variantPositions, List<Integer> partition, String realHaplotype, PrintStream out) {
		DecimalFormat df = new DecimalFormat("##0.00");
		out.println("Switch errors for block "+b.getFirstPos()+" "+b.getLastPos());
		boolean [] cut = null;
		if(partition.size()>2) {
			FragmentsCutBuilder builder = new FragmentsCutBuilder(b.getFragments());
			builder.calculateMaxCut();
			cut = builder.getCut();
			/*for(int j=0;j<b.getFragments().size();j++) {
				out.println(""+b.getFragments().get(j).getId()+" "+cut[j]);
			}*/
		}
		for(int i=1;i<partition.size()-1;i+=2) {
			int pos = partition.get(i);
			int pos2 = partition.get(i+1);
			List<Fragment> frags = b.getFragments(pos-2,pos2+2);
			boolean [] subcut = calculateSubcut(b.getFragments(),frags,cut);
			int startAln = Math.max(0, pos-10);
			int endAln = Math.min(realHaplotype.length(), pos2+10);
			out.println("Switch error between "+variantPositions.get(pos)+" and "+variantPositions.get(pos2));
			out.println("Relative start alignment "+(startAln+1)+" pos1 "+(pos+1)+" pos2 "+(pos2+1));
			
			out.println(realHaplotype.substring(startAln,endAln));
			String hap = b.getHaplotype();
			int relStart = startAln - b.getFirstPos();
			int relEnd = endAln - b.getFirstPos();
			int relPos = relStart;
			while(relPos<0) {
				out.print(Fragment.NODATACHAR);
				relPos++;
			}
			int end = Math.min(hap.length(),relEnd);
			out.print(hap.substring(relPos,end));
			relPos = end;
			while(relPos<relEnd) {
				out.print(Fragment.NODATACHAR);
				relPos++;
			}
			out.println();
			int j=0;
			for(Fragment f: frags) {
				String calls = f.getExtendedCalls();
				relStart = startAln - f.getFirstPos();
				relEnd = endAln - f.getFirstPos();
				relPos = relStart;
				while(relPos<0) {
					out.print(Fragment.NODATACHAR);
					relPos++;
				}
				end = Math.min(calls.length(),relEnd);
				out.print(calls.substring(relPos,end));
				relPos = end;
				while(relPos<relEnd) {
					out.print(Fragment.NODATACHAR);
					relPos++;
				}
				out.print("\t"+f.getId()+" "+subcut[j]);
				double [] probs = f.getExtendedProbabilities();
				if(probs!=null) {
					for(int k=0;k<probs.length;k++) {
						out.print(" "+df.format(probs[k]));
					}
				}
				out.println();
				j++;
			}
		}
		
	}
	private boolean[] calculateSubcut(List<Fragment> allFrags, List<Fragment> selectedFrags, boolean[] cut) {
		boolean answer [] = new boolean [selectedFrags.size()];
		int i = 0;
		for(int j=0;j<answer.length;j++) {
			Fragment selF = selectedFrags.get(j);
			while(!allFrags.get(i).getId().equals(selF.getId())) {
				i++;
			}
			answer[j] = cut[i];
		}
		return answer;
	}
	public void calculateStatistics (String seqName, String methodName,boolean singleHap) throws IOException {
		FragmentsFileHandler loader = new FragmentsFileHandler();
		HaplotypesFileHandler loaderHaps = new HaplotypesFileHandler();
		String variantsFile = seqName+".allvars";
		String haplotypeFile = seqName+"."+methodName+".phase";
		String realHaplotypeFile = seqName+".real.phase";
		String outputHapFile = seqName+".real_"+methodName+".phase";
		String matrixFile = seqName+".matrix.SORTED";
		List <Integer> variantPositions= loader.loadVariantPositions(variantsFile,2);
		String realHaplotype = loaderHaps.loadHaplotype(realHaplotypeFile,variantPositions);
		if(!singleHap) {
			List<Fragment> fragments = loader.loadFragments(matrixFile);
			BlocksBuilder blocksBuilder = new BlocksBuilder();
			List<Block> blocks = blocksBuilder.buildBlocks(new Block(fragments));
			loaderHaps.loadPhasedBlocks(blocks, haplotypeFile);
			String outputHaplotype = calculateStatistics(blocks, variantPositions, realHaplotype);
			printHaplotype(variantPositions,outputHaplotype,new PrintStream(outputHapFile));
		} else {
			String testHaplotype = loaderHaps.loadHaplotype(haplotypeFile, variantPositions);
			calculateStatistics(testHaplotype, variantPositions, realHaplotype);
		}
	}
	
	private void printHaplotype(List<Integer> variantPositions,String outputHaplotype, PrintStream out) {
		for(int i=0;i<variantPositions.size();i++) {
			char call = outputHaplotype.charAt(i);
			char call2 = call;
			if(call == Fragment.ALLELE1CHAR) call2 = Fragment.ALLELE2CHAR;
			else if (call == Fragment.ALLELE2CHAR) call2 = Fragment.ALLELE1CHAR;
			out.println(""+variantPositions.get(i)+"\t"+call+"\t"+call2);
		}
		out.flush();
		out.close();
	}
	public static void main(String[] args) throws Exception {
		PhasingStatisticsCalculator calculator = new PhasingStatisticsCalculator();
		
		int i=0;
		boolean singleHap = false;
		String seqName = null;
		while(i<args.length && args[i].startsWith("-")) {
			if ("-h".equals(args[i])) {
				singleHap = true;
			} else if ("-s".equals(args[i])) {
				i++;
				seqName = args[i];
			}
			i++;
		}
		String methodName = args[i++];
		if(seqName==null) {
			for(int c=1;c<=22;c++) {
				calculator.calculateStatistics("chr"+c, methodName, singleHap);
				
			}
			calculator.calculateStatistics("chrX", methodName, singleHap);
		} else {
			calculator.calculateStatistics(seqName, methodName, singleHap);
		}
		
		calculator.initialStats.printSummaryStatistics(System.out);
		calculator.phasingBlockStats.printSummaryStatistics(System.out);
		System.out.println();
		if(!singleHap) {
			PrintStream out = new PrintStream("haplotyper."+methodName+".initial.stats");
			calculator.initialStats.printDetailStatistics(out);
			out.flush();
			out.close();
			out = new PrintStream("haplotyper."+methodName+".quality.stats");
			calculator.phasingBlockStats.printDetailStatistics(out);
			out.flush();
			out.close();
		}
		
		//System.out.println(realHaplotype.length());
		//System.out.println("0"+blocks.get(0).getHaplotype());
	}	
}
