package progistar.pXg.train;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Random;

public class AADistribution {
	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/resources/Sequences/Human_UP000005640_202204.fasta");
		double[][][] aaTable = pairCompositionTableFrompFasta(file, 1000, 1000);
		
		System.out.print("Num");
		for(int i=0; i<26; i++) {
			if((char)(i+'A') == 'B' || (char)(i+'A') == 'J' || (char)(i+'A') == 'O' || (char)(i+'A') == 'U' || (char)(i+'A') == 'Z') {
				continue;
			}
			for(int j=0; j<26;j ++) {
				if((char)(j+'A') == 'B' || (char)(j+'A') == 'J' || (char)(j+'A') == 'O' || (char)(j+'A') == 'U' || (char)(j+'A') == 'Z') {
					continue;
				}
				System.out.print("\t"+(char)(i+'A')+""+(char)(j+'A'));
			}
			
		}
		for(int num = 0; num<aaTable.length; num++) {
			System.out.println();
			System.out.print(num+1);
			
			double sum = 0;
			for(int i=0; i<aaTable[num].length; i++) {
				if((char)(i+'A') == 'B' || (char)(i+'A') == 'J' || (char)(i+'A') == 'O' || (char)(i+'A') == 'U' || (char)(i+'A') == 'Z') {
					continue;
				}
				for(int j=0; j<aaTable[num][i].length; j++) {
					if((char)(j+'A') == 'B' || (char)(j+'A') == 'J' || (char)(j+'A') == 'O' || (char)(j+'A') == 'U' || (char)(j+'A') == 'Z') {
						continue;
					}
					sum += aaTable[num][i][j]; 
				}
			}
			
			for(int i=0; i<aaTable[num].length; i++) {
				if((char)(i+'A') == 'B' || (char)(i+'A') == 'J' || (char)(i+'A') == 'O' || (char)(i+'A') == 'U' || (char)(i+'A') == 'Z') {
					continue;
				}
				for(int j=0; j<aaTable[num][i].length; j++) {
					if((char)(j+'A') == 'B' || (char)(j+'A') == 'J' || (char)(j+'A') == 'O' || (char)(j+'A') == 'U' || (char)(j+'A') == 'Z') {
						continue;
					}
					System.out.print("\t"+((aaTable[num][i][j]/sum)*100));
				}
			}
		}
		
		
	}
	
	public static double[][][] pairCompositionTableFrompFasta (File file, int run, int num) throws IOException {
		
		double[][][] aaTable = new double[run][26][26];
		BufferedReader BR = new BufferedReader(new FileReader(file));
		ArrayList<String> sequences = new ArrayList<String>();
		String line = null;
		
		StringBuilder str = new StringBuilder();
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(str.length() != 0) {
					sequences.add(str.toString());
				}
				str.setLength(0);
			} else {
				str.append(line);
			}
		}
		sequences.add(str.toString());
		
		BR.close();
		
		Random random = new Random();
		for(int i=0; i<run; i++) {
			
			Hashtable<Integer, Boolean> nums = new Hashtable<Integer, Boolean>();
			while(nums.size() < num) {
				int idx = random.nextInt(sequences.size());
				nums.put(idx, true);
			}
			
			final int runIdx = i;
			nums.forEach((idx, nil) -> {
				String peptide = sequences.get(idx);
				for(int j=1; j<peptide.length(); j++) {
					char prevAA = peptide.charAt(j-1);
					char nextAA = peptide.charAt(j);
					aaTable[runIdx][prevAA-'A'][nextAA-'A']++;
				}
			});
		}
		
		return aaTable;
	}
}
