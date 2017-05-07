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
package mpg.molgen.sih.main;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import mpg.molgen.sih.algorithms.SIHAlgorithm;
import mpg.molgen.sih.model.Block;
import mpg.molgen.sih.model.BlocksBuilder;
import mpg.molgen.sih.model.Fragment;
import mpg.molgen.sih.model.FragmentsFileHandler;

public class SIH {
	private List<Integer> variantPositions = new ArrayList<Integer>();
	private String algorithmClassName = "mpg.molgen.sih.algorithms.RefhapAlgorithm";
	private SIHAlgorithm algClass;
	private Random r = new Random();
	/**
	 * @return the algorithmClassName
	 */
	public String getAlgorithmClassName() {
		return algorithmClassName;
	}
	/**
	 * @param algorithmClassName the algorithmClassName to set
	 * @throws ClassNotFoundException 
	 * @throws IllegalAccessException 
	 * @throws Exception 
	 */
	public void setAlgorithmClassName(String algorithmClassName) throws Exception {
		this.algorithmClassName = algorithmClassName;
		algClass = (SIHAlgorithm)Class.forName(algorithmClassName).newInstance();
	}
	public List<Block> buildHaplotypes (List<Fragment> fragments) throws Exception {
		BlocksBuilder blocksBuilder = new BlocksBuilder();
		List<Block> blocks = blocksBuilder.buildBlocks(new Block(fragments));
		for(Block b:blocks) {
			buildHaplotype(b);
		}
		return blocks;
	}
	public void buildHaplotype (Block b) throws Exception {
		if(algClass == null) {
			algClass = (SIHAlgorithm)Class.forName(algorithmClassName).newInstance();
		}
		algClass.buildHaplotype(b);
	}
	private void printBlock(Block b, PrintStream out) {
		out.println("BLOCK: offset: "+(b.getFirstPos()+1)+" len: "+b.getHaplotype().length()+" phased: "+b.getPhased());
		String hap = b.getHaplotype();
		boolean printReverse = r.nextBoolean();
		for(int i=0;i<hap.length();i++) {
			char c = hap.charAt(i);
			char c2;
			if(c==Fragment.NODATACHAR) {
				c2=c;
			} else if(c==Fragment.ALLELE1CHAR) {
				c2=Fragment.ALLELE2CHAR;
			} else {
				c2=Fragment.ALLELE1CHAR;
			}
			int pos = b.getFirstPos()+i;
			if(variantPositions.size()>pos) {
				pos = variantPositions.get(pos);
			} else {
				pos++;
			}
			if(printReverse) {
				out.println(""+pos+"\t"+c2+"\t"+c);
			} else {
				out.println(""+pos+"\t"+c+"\t"+c2);
			}
			
		}
		out.println("********");
		out.flush();
		
	}
	public static void main(String[] args) throws Exception {
		long time = System.currentTimeMillis();
		FragmentsFileHandler loader = new FragmentsFileHandler();
		SIH h = new SIH();
		String snpsFile =null;
		int posColumn = 0;
		String algorithmName = "Refhap";
		int i=0;
		
		while(i<args.length && args[i].startsWith("-")) {
			if("-v".equals(args[i])) {
				i++;
				snpsFile = args[i];
			} else if("-c".equals(args[i])) {
				i++;
				posColumn = Integer.parseInt(args[i])-1;
			} else if("-a".equals(args[i])) {
				i++;
				algorithmName = args[i];
			} else if ("-h".equals(args[i])) {
				printUsage(System.out);
				System.exit(0);
			} else {
				System.err.println("Unrecognized option: "+args[i]);
				printUsage(System.err);
				System.exit(1);
			}
			i++;
		}
		if(snpsFile!=null) h.variantPositions = loader.loadVariantPositions(snpsFile,posColumn);
		if(i>args.length-2) {
			System.err.println("Missing required parameters.");
			printUsage(System.err);
			System.exit(1);
		}
		String inputFile = args[i++];
		String outputFile = args[i++];
		h.setAlgorithmClassName("mpg.molgen.sih.algorithms."+algorithmName+"Algorithm");
		List<Fragment> f = loader.loadFragments(inputFile);
		List<Block> blocks = h.buildHaplotypes(f);
		System.out.println("Number of blocks:" +blocks.size());
		PrintStream out = new PrintStream(outputFile);
		int totalPhased=0;
		int totalCalls=0;
		int totalMEC =0;
		for(Block b:blocks) {
			System.out.println("Offset: "+b.getFirstPos()+" Phased: "+b.getPhased()+ " Calls: "+b.getCalls()+ " MEC: "+b.getMEC());
			totalPhased+=b.getPhased();
			totalCalls+=b.getCalls();
			totalMEC+=b.getMEC();
			h.printBlock(b,out);
		}
		out.close();
		double diff = System.currentTimeMillis() - time;
		diff/=1000;
		System.out.println("Phased: "+totalPhased+ " Calls: "+totalCalls+ " MEC: "+totalMEC+" Time(s): "+diff);
	}
	private static void printUsage(PrintStream out) {
		out.println("USAGE");
		out.println("java -cp SIH.jar mpg.molgen.sih.main.SIH <OPTIONS> <INPUT_FILE> <OUTPUT_FILE>");
		out.println("OPTIONS: ");
		out.println("\t-v FILE\t\t: Text file with genomic coordinates of variants.");
		out.println("\t-c INT\t\t: Column in the variants file where coordinates are located.");
		out.println("\t-a STRING\t: Name of the algorithm to run (Refhap, DGS, FastHare).");
	}
	
}
