package progistar.pXg.constants;

import java.util.Hashtable;

public class RunInfo {

	public static long totalProcessedReads = 0;
	public static long totalProcessedPeptides = 0;
	public static Hashtable<String, Integer> processedChromosomes = new Hashtable<String, Integer>();
	
	
	public static void printProcessedChromosomes () {
		StringBuilder chrList = new StringBuilder();
		processedChromosomes.forEach((chr, value)->{
			chrList.append("|").append(chr).append("=").append(value);
		});
		
		System.out.println(chrList.substring(1));
	}
}
