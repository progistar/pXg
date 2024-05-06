package progistar.pXg.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

import progistar.pXg.data.GenomicSequence;

public class GenTranslationDB {

	public static final String INPUT = "-i";
	public static final String MODE = "-m";
	public static File inputFile = null;
	public static String fileFormat = null;
	public static String mode = "";
	public static int len = 15;

	public static void main(String[] args) throws IOException {
		if(!parseArg(args)) {
			System.exit(1);
		}

		long startTime = System.currentTimeMillis();
		BufferedReader BR = new BufferedReader(new FileReader(inputFile));
		BufferedWriter BW = new BufferedWriter(new FileWriter(inputFile.getAbsolutePath().replace(".sam", ".fasta")));
		String line = null;

		Hashtable<String, Integer> map = new Hashtable<>();
		while((line = BR.readLine()) != null) {
			if(line.startsWith("@")) {
				continue;
			}
			String[] fields = line.split("\t");
			String nucleotide = fields[9];
			// forward
			for(int frame=0; frame<3; frame++) {
				String[] peptides = GenomicSequence.translation(nucleotide, frame).split("X");
				for (String peptide : peptides) {
					int len = peptide.length();

					if(len < 8) {
						continue;
					}

					int start = 0;
					int end = len < 15 ? len : 15;

					while(end <= len) {
						String subStr = peptide.substring(start, end);
						Integer cnt = map.get(subStr);
						if(cnt == null) {
							cnt = 0;
						}
						cnt++;
						map.put(subStr, cnt);
						end++;
						start++;
					}
				}
			}
			// reverse
			for(int frame=0; frame<3; frame++) {
				String[] peptides = GenomicSequence.reverseComplementTranslation(nucleotide, frame).split("X");
				for (String peptide : peptides) {
					int len = peptide.length();

					if(len < 8) {
						continue;
					}

					int start = 0;
					int end = len < 15 ? len : 15;

					while(end <= len) {
						String subStr = peptide.substring(start, end);
						Integer cnt = map.get(subStr);
						if(cnt == null) {
							cnt = 0;
						}
						cnt++;
						map.put(subStr, cnt);
						end++;
						start++;
					}
				}
			}
		}

		// 0 for peptide cnt
		// 1 for aa cnt
		long[] summary = new long[2];

		long maxMemory = Runtime.getRuntime().maxMemory();
		long allocatedMemory = Runtime.getRuntime().totalMemory();
		long freeMemory = Runtime.getRuntime().freeMemory();
		long totalMemory = Runtime.getRuntime().totalMemory();

		System.out.println("maxMem:" +maxMemory/1024);
		System.out.println("allocated: "+allocatedMemory/1024);
		System.out.println("freeMem: "+freeMemory/1024);
		System.out.println("totalMem: "+totalMemory/1024);
		System.out.println("usedMem: "+(totalMemory-freeMemory)/1024);

		map.forEach((peptide, cnt)->{
			try {
				summary[0]++;
				summary[1] += peptide.length();

				BW.append(">"+summary[0]+" "+cnt);
				BW.newLine();
				BW.append(peptide);
				BW.newLine();
			}catch(IOException ioe) {

			}
		});
		long endTime = System.currentTimeMillis();

		System.out.println((endTime-startTime)/1000 +"sec");

		System.out.println("number of peptides: "+summary[0]);
		System.out.println("AA size: "+summary[1]);

		BW.close();
		BR.close();
	}



	public static boolean parseArg(String[] args) {
		for(int i=0; i<args.length; i++) {
			if(args[i].equalsIgnoreCase(INPUT)) {
				inputFile = new File(args[++i]);
				if(inputFile.exists()) {
					if(inputFile.getName().toLowerCase().endsWith(".sam")) {
						fileFormat = "sam";
					}
				}
			} else if(args[i].equalsIgnoreCase(MODE)) {
			 	mode = args[++i];
			}
		}

		// check
		if(inputFile == null || !inputFile.exists() || ! (mode.equalsIgnoreCase("s") || mode.equalsIgnoreCase("t")) || fileFormat == null) {
			System.out.println("Usage: java -Xmx2G -jar -i [input file name] -m [translation mode]");
			System.out.println("-i: sam is possible");
			System.out.println("-m: s for six-frame, t for three-frame");
			return false;
		}
		return true;
	}
}
