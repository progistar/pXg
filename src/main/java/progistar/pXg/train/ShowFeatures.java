package progistar.pXg.train;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.stream.Stream;

public class ShowFeatures {

	public static void main(String[] args) throws IOException {
		File file = new File("/Users/gistar/resources/Sequences/Human_UP000005640_202204.fasta");
		long[][] aaTable = pairCompositionTableFrompFasta(file);
		
//		File file = new File("/Users/gistar/projects/pXg/Laumont_NatCommun2016/S1.RAW.PEAKS.DecoyOut.pxg");
//		long[][] aaTable = pairCompositionTableFrompXg(file);

//		File file = new File("/Users/gistar/projects/pXg/IEDB/IEDB.host");
//		long[][] aaTable = pairCompositionTableFromIEDB(file);
		
		for(int i=0; i<26; i++) {
			if((char)(i+'A') == 'B' || (char)(i+'A') == 'J' || (char)(i+'A') == 'O' || (char)(i+'A') == 'U' || (char)(i+'A') == 'Z') {
				continue;
			}
			System.out.print("\t"+(char)(i+'A'));
		}
		System.out.println();
		for(int i=0; i<aaTable.length; i++) {
			if((char)(i+'A') == 'B' || (char)(i+'A') == 'J' || (char)(i+'A') == 'O' || (char)(i+'A') == 'U' || (char)(i+'A') == 'Z') {
				continue;
			}
			System.out.print((char)(i+'A'));
			for(int j=0; j<aaTable[i].length; j++) {
				if((char)(j+'A') == 'B' || (char)(j+'A') == 'J' || (char)(j+'A') == 'O' || (char)(j+'A') == 'U' || (char)(j+'A') == 'Z') {
					continue;
				}
				System.out.print("\t"+aaTable[i][j]);
			}
			System.out.println();
		}
		
	}
	
	public static long[][] pairCompositionTableFrompFasta (File file) throws IOException {
		
		long[][] aaTable = new long[26][26];
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
		for(String peptide : sequences) {
			for(int i=1; i<peptide.length(); i++) {
				char prevAA = peptide.charAt(i-1);
				char nextAA = peptide.charAt(i);
				aaTable[prevAA-'A'][nextAA-'A']++;
			}
			
		}
		return aaTable;
	}
	
	public static long[][] pairCompositionTableFrompXg (File file) throws IOException {
		long[][] aaTable = new long[26][26];
		
		Stream<String> text = Files.lines(file.toPath()).parallel();
		
		Hashtable<String, String> checker = new Hashtable<String, String>();
		
		text.parallel().forEach(line -> {
			String[] fields = line.split("\t");
			String label = fields[0];
			String isCanonical = fields[37];
//			String binder = fields[43];
			
			if(label.equalsIgnoreCase("-1") && isCanonical.equalsIgnoreCase("true")) {
				String peptide = fields[21];
				if(checker.get(peptide) == null) {
					for(int i=1; i<peptide.length(); i++) {
						char prevAA = peptide.charAt(i-1);
						char nextAA = peptide.charAt(i);
						aaTable[prevAA-'A'][nextAA-'A']++;
					}
					checker.put(peptide,"");
				}
			}
		});
		
		text.close();
		
		return aaTable;
	}
	
	public static long[][] pairCompositionTableFromIEDB (File file) throws IOException {
		long[][] aaTable = new long[26][26];
		
		Stream<String> text = Files.lines(file.toPath()).parallel();
		
		text.parallel().forEach(line -> {
			String[] fields = line.split("\t");
			if(!fields[0].startsWith("Description")) {
				String peptide = fields[0].split("\\s")[0];
				for(int i=1; i<peptide.length(); i++) {
					char prevAA = peptide.charAt(i-1);
					char nextAA = peptide.charAt(i);

					try {
						aaTable[prevAA-'A'][nextAA-'A']++;
					}catch(Exception e) {
						System.out.println(line);
					}
				}
			}
		});
		
		text.close();
		
		return aaTable;
	}
}
