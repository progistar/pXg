package progistar.pXg.constants;

import java.util.Hashtable;

public class RunInfo {

	public static long totalProcessedReads = 0;
	public static long totalProcessedPeptides = 0;
	public static int[] cutoffReads = null;
	public static double cPSMScoreTreshold = 0;
	public static double ncPSMScoreTreshold = 0;
	public static Hashtable<String, Integer> processedChromosomes = new Hashtable<String, Integer>();
	
	
	public static void printProcessedChromosomes () {
		StringBuilder chrList = new StringBuilder();
		processedChromosomes.forEach((chr, value)->{
			chrList.append("|").append(chr).append("=").append(value);
		});
		
		System.out.println(chrList.substring(1));
	}
	
	public static void printPSMCutoff () {
		System.out.println("PSM FDR cutoff for canonical: "+cPSMScoreTreshold);
		System.out.println("PSM FDR cutoff for noncanonical: "+ncPSMScoreTreshold);
	}
}
