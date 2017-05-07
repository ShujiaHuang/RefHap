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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

public class HaplotypesFileHandler {
	public void loadPhasedBlocks(List<Block> fragmentBlocks, String filename) throws IOException {
		FileInputStream fis = new FileInputStream(filename);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		String line=in.readLine();
		StringBuilder lastHap =null;
		int i=0;
		while(line!=null) {
			if(line.startsWith("BLOCK")) {
				if(lastHap !=null) {
					fragmentBlocks.get(i).setHaplotype(lastHap.toString());
					i++;
				}
				lastHap = new StringBuilder();
			} else {
				String [] items = line.split(" |\t");
				if(items.length==3) {
					lastHap.append(items[1].charAt(0));
				}
			}
			line=in.readLine();
		}
		fis.close();
		fragmentBlocks.get(i).setHaplotype(lastHap.toString());
	}
	public String loadHaplotype (String filename,List<Integer> variantPositions) throws IOException {
		StringBuilder haplotype = new StringBuilder();
		FileInputStream fis = new FileInputStream(filename);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		String line=in.readLine();
		int i = 0;
		while(line!=null) {
			String [] items = line.split(" |\t");
			if(items.length==3) {
				try {
					int variantPos = Integer.parseInt(items[0]);
					while(i<variantPositions.size() && variantPositions.get(i)<variantPos) {
						haplotype.append(Fragment.NODATACHAR);
						i++;
					}
					if(i<variantPositions.size() ) {
						if(variantPositions.get(i) == variantPos) {
							char allele = items[1].charAt(0);
							haplotype.append(allele);
							i++;
						}
					}
				} catch (NumberFormatException e) {		
				}
			}
			line=in.readLine();
		}
		fis.close();
		while(i<variantPositions.size()) {
			haplotype.append(Fragment.NODATACHAR);
			i++;
		}
		return haplotype.toString();
	}
}
