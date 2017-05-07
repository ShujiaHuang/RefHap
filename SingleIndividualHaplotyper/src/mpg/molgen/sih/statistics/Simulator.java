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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import mpg.molgen.sih.main.SIH;
import mpg.molgen.sih.model.Block;
import mpg.molgen.sih.model.BlocksBuilder;
import mpg.molgen.sih.model.Fragment;
import mpg.molgen.sih.model.FragmentsComparator;

public class Simulator {
	private int avgFragLength = 6;
	private int hapLength = 100;
	private  int numFrags = 100;
	private double errorRate = 0.05;
	private double gapRate= 0.1;
	//private String [] algorithmNames = {"Refhap","DGS","FastHare","SHRThree","Speedhap","TwoDMEC","WMLF"};
	//private String [] algorithmNames = {"DGS","SHRThree","TwoDMEC","WMLF"};
	private String [] algorithmNames = {"Refhap","DGS","FastHare"};
	//private String [] algorithmNames = {"WMLF"};
	
	private Random r = new Random();
	public String createRandomHap () {
		boolean [] calls = new boolean[hapLength];
		for(int i=0;i<hapLength;i++) {
			calls[i] = r.nextBoolean();
		}
		return format(calls);
	}
	
	private List<Fragment> buildFragments(String realHap) {
		List<Fragment> frags = new ArrayList<Fragment>();
		Random r = new Random();
		for(int i=0;i<numFrags;i++) {
			int fragLength = 1;
			while(fragLength < 2 ) {
				fragLength = (int)Math.round(r.nextGaussian()+avgFragLength);
			}
			int start = r.nextInt(hapLength-fragLength+1);
			char [] seq = realHap.substring(start,start+fragLength).toCharArray();
			double [] probabilities = new double [seq.length];
			Arrays.fill(probabilities, 0);
			if(r.nextBoolean()) {
				flipAll(seq);
			}
			for(int j=0;j<seq.length;j++) {
				if(j>0 && j<seq.length-1 && r.nextDouble() < gapRate) {
					seq[j] = '-';
				} else if (errorRate > 0 && r.nextDouble() < errorRate) {
					if(seq[j] == '1') {
						seq[j] = '0'; 
					} else {
						seq[j] = '1';
					}
					double p = r.nextDouble()/5+0.2;
					probabilities[j] = p;
				} else {
					double p = r.nextDouble()/5;
					probabilities[j] = p;
					
				}
			}
			Fragment f = new Fragment(""+(i+1),start,new String(seq),probabilities);
			//System.out.println("Fragment start:"+f.getFirstPos());
			frags.add(f);
		}
		Collections.sort(frags,new FragmentsComparator());
		return frags;
	}
	private void flipAll(char[] seq) {
		for(int j=0;j<seq.length;j++) {
			if(seq[j] == '1') {
				seq[j] = '0'; 
			} else {
				seq[j] = '1';
			}
		}
		
	}
	private String format(boolean[] calls) {
		StringBuilder builder = new StringBuilder(calls.length);
		for(int i=0;i<calls.length;i++) {
			if(calls[i]) {
				builder.append('1');
			} else {
				builder.append('0');
			}
		}
		return builder.toString();
	}
	public void printHap(String hap, PrintStream out) {
		out.println("BLOCK: offset: 1 len: "+hap.length()+" phased: "+hap.length());
		for(int i=0;i<hap.length();i++) {
			char c1 = hap.charAt(i);
			char c2 = '1';
			if(c1 == c2) {
				c2 = '0';
			}
			out.println(""+(i+1)+" "+c1+ " "+c2);
		}
	}
	private void printVariants(int num, PrintStream out) {
		for(int i=0;i<num;i++) {
			out.println(""+(i+1)+" chr1 "+(i+1)+ " A T AT ");
		}
		
	}
	public void printMatrix(List<Fragment> fragments, PrintStream out) {
		out.println(""+(fragments.size()+1)+" "+hapLength);
		for(Fragment f:fragments) {
			int firstPos = f.getFirstPos()+1;
			String [] subseqs = f.getExtendedCalls().split("-");
			StringBuffer finalSeq = new StringBuffer();
			int nParts = 0;
			
			int nextPos = firstPos;
			for(int i=0;i<subseqs.length;i++) {
				String subseq = subseqs[i];
				if(subseq.length() > 0) {
					finalSeq.append(" "+nextPos+" "+subseq);
					nextPos+=subseq.length()+1;
					nParts++;
				} else {
					nextPos++;
				}
			}
			out.print(""+nParts+" "+f.getId()+finalSeq);
			//TODO: Print properly quality scores
			/*String qualScores = f.getExtendedQualityScores(); 
			if(f!=null) {
				for(int j=0;j<qualScores.length();j++) {
					if (j==0) out.print (" ");
					char score = qualScores.charAt(j);
					if(score!=Fragment.NODATACHAR) out.print(qualScores.charAt(j));
				}
					
			} */
			out.println();
		}
		out.flush();
	}
	public void printMatrixHapler(List<Fragment> fragments, PrintStream out) {
		for(Fragment f:fragments) {
			int i;
			for(i=0;i<f.getFirstPos();i++) {
				out.print(Fragment.NODATACHAR);
			}
			String calls = f.getExtendedCalls(); 
			out.print(calls);
			i+=calls.length();
			for(;i<hapLength;i++) {
				out.print(Fragment.NODATACHAR);
			}
			out.println();
		}
		out.flush();
	}
	public void printScoresHapler(List<Fragment> fragments, PrintStream out) {
		for(Fragment f:fragments) {
			out.println(""+f.getFirstPos()+" "+f.getLastPos());
			for(int i=0;i<f.length();i++) {
				out.print(" "+(1-errorRate)+" ");
			}
			out.println();
		}
		out.flush();
	}
	
	
	public void runSimulation (int times, PrintStream out) throws Exception {
		DecimalFormat df = new DecimalFormat("##0.00");
		PhasingStatisticsCalculator calculator = new PhasingStatisticsCalculator();
		double [] mecs= new double[algorithmNames.length];
		double [] mecsSquared= new double[algorithmNames.length];
		double [] switchErrors= new double[algorithmNames.length];
		double [] switchErrorsSquared= new double[algorithmNames.length];
		double [] runtimes= new double[algorithmNames.length];
		double [] runtimesSquared= new double[algorithmNames.length];
		Arrays.fill(mecs,0);
		Arrays.fill(mecsSquared,0);
		Arrays.fill(switchErrors,0);
		Arrays.fill(switchErrorsSquared,0);
		Arrays.fill(runtimes,0);
		Arrays.fill(runtimesSquared,0);
		double calls = 0; 
		double phasedVars = 0;
		double numBlocks = 0;
		for(int i=0;i<times;i++) {
			String realHap = createRandomHap();
			List<Fragment> fragments = buildFragments(realHap);
			SIH hap = new SIH();
			BlocksBuilder blocksBuilder = new BlocksBuilder();
			List<Block> blocks = blocksBuilder.buildBlocks(new Block(fragments));
			numBlocks+=blocks.size();
			for(int j=0;j<algorithmNames.length;j++) {
				String algClassName = "mpg.molgen.sih.algorithms."+algorithmNames[j]+"Algorithm";
				hap.setAlgorithmClassName(algClassName);
				int mec=0;
				int se=0;
				long start = System.currentTimeMillis();
				for(Block b:blocks) {
					hap.buildHaplotype(b);
					if(j==0) {
						calls+=b.getCalls();
						phasedVars+=b.getPhased();
					}
					mec+=b.getMEC();
					se+=Math.max(0,calculator.calculateSEPartition(b, realHap).size()/2-1);
				}
				long runtime = System.currentTimeMillis()-start;
				mecs[j]+=mec;
				mecsSquared[j]+=mec*mec;
				//String answer = hap.getHaplotype(blocks, realHap.length());
				switchErrors[j]+=se;
				switchErrorsSquared[j]+=(se*se);
				runtimes[j]+=runtime;
				runtimesSquared[j]+=runtime*runtime;
			}
		}
		numBlocks/=times;
		calls/=times;
		phasedVars/=times;
		out.print (""+hapLength+" "+numFrags+" "+avgFragLength+" "+df.format(errorRate)+" "+df.format(gapRate));
		out.print (" "+numBlocks+" "+df.format(phasedVars)+" "+df.format(calls)+" "+df.format(calls/hapLength));
		for(int j=0;j<switchErrors.length;j++) {
			double average = switchErrors[j]/times;
			double variance = (switchErrorsSquared[j]-switchErrors[j]*switchErrors[j]/times)/(times-1);
			out.print (" "+df.format(average)+" "+df.format(variance));
			average = runtimes[j]/times;
			variance = (runtimesSquared[j]-runtimes[j]*runtimes[j]/times)/(times-1);	
			out.print (" "+df.format(average)+" "+df.format(variance));
			average = mecs[j]/times;
			variance = (mecsSquared[j]-mecs[j]*mecs[j]/times)/(times-1);	
			out.print (" "+df.format(average)+" "+df.format(variance));
		}
		out.println();
		out.flush();
	}
	
	public static void main(String[] args) throws Exception {
		Simulator s = new Simulator();
		int i=0;
		boolean print = false;
		boolean one = false;
		int times = 100;
		for(;i<args.length && args[i].charAt(0)=='-';i++) {
			if("-fl".equalsIgnoreCase(args[i])) {
				i++;
				s.avgFragLength = Integer.parseInt(args[i]);
			} else if("-hl".equalsIgnoreCase(args[i])) {
				i++;
				s.hapLength = Integer.parseInt(args[i]);
			} else if("-n".equalsIgnoreCase(args[i])) {
				i++;
				s.numFrags = Integer.parseInt(args[i]);
			} else if("-e".equalsIgnoreCase(args[i])) {
				i++;
				s.errorRate = Double.parseDouble(args[i]);
			} else if("-g".equalsIgnoreCase(args[i])) {
				i++;
				s.gapRate = Double.parseDouble(args[i]);
			} else if("-p".equalsIgnoreCase(args[i])) {
				print = true;
			} else if("-1".equalsIgnoreCase(args[i])) {
				one = true;
			} else if("-t".equalsIgnoreCase(args[i])) {
				i++;
				times = Integer.parseInt(args[i]);
			} else if("-a".equalsIgnoreCase(args[i])) {
				i++;
				s.algorithmNames = args[i].split(",");
			} 
		}
		
		if(print) {
			String realHap = s.createRandomHap();
			List<Fragment> fragments = s.buildFragments(realHap);
			PrintStream out = new PrintStream("sim.matrix.SORTED");
			s.printMatrix(fragments,out);
			out.close();
			out = new PrintStream("sim.variants");
			s.printVariants(realHap.length(), out);
			out.close();
			out = new PrintStream("sim.real.phase");
			s.printHap(realHap,out);
			out.close();
			out = new PrintStream("sim.random.phase");
			s.printHap(s.createRandomHap(),out);
			out.close();
			out = new PrintStream("sim.hapler.matrix");
			s.printMatrixHapler(fragments, out);
			out.close();
			out = new PrintStream("sim.hapler.weights");
			s.printScoresHapler(fragments, out);
			out.close();
		} else if(!one ){
			int fragmentLengths[] = {3,6,10,15};
			int coverages [] = {4,6,8,10,12};
			double errorRates [] = {0,0.01,0.02,0.05,0.1};
			for(int f=0;f<fragmentLengths.length;f++) {
				s.avgFragLength=fragmentLengths[f];
				for(int c=0;c<coverages.length;c++) {
					s.setCoverage(coverages[c]);
					for(int e=0;e<errorRates.length;e++) {
						s.errorRate = errorRates[e];
						s.runSimulation(times, System.out);
					}
				}
			}
		} else {
			s.runSimulation(times, System.out);
		}
	}

	

	private void setCoverage(int coverage) {
		double expectedCalls = (double)(coverage * hapLength);
		if(avgFragLength==2) {
			numFrags = (int) Math.round( expectedCalls / avgFragLength);
		} else if (avgFragLength == 3) {
			numFrags = (int) Math.round( expectedCalls / (avgFragLength*(1-0.5*gapRate)));
		} else {
			numFrags = (int)Math.round(expectedCalls / (avgFragLength*(1-gapRate)));
		}
			
	}
}