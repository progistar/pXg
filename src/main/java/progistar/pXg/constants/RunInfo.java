package progistar.pXg.constants;

import java.util.Hashtable;

import progistar.pXg.utils.Logger;

public class RunInfo {
	
	// Count PSMs
	public static long[] totalRankPSMs = new long[101];
	
	public static long[] workerProcessedReads = null;
	public static long totalProcessedReads = 0;
	public static long totalProcessedPeptides = 0;
	public static int[] cutoffReads = null;
	public static double cPSMScoreTreshold = 0;
	public static double ncPSMScoreTreshold = 0;
	public static Hashtable<String, Integer> processedChromosomes = new Hashtable<String, Integer>();
	
	// original scans & peptides
	public static int initialScanNum				= 0;
	public static int initialPeptideNum				= 0;
	
	// step1: rank filter
	public static int rankFilterScanNum1			= 0;
	public static int rankFilterPeptideNum1			= 0;
	
	// step2: length filter
	public static int lengthFilterScanNum2			= 0;
	public static int lengthFilterPeptideNum2		= 0;
	
	// step3: xBlock mapping filter (exp read mapping)
	public static int mappingFilterScanNum3			= 0;
	public static int mappingFilterPeptideNum3		= 0;
	
	// step4: p-value
	public static int pvalueFilterScanNum4			= 0;
	public static int pvalueFilterPeptideNum4		= 0;
	
	public static void printProcessedChromosomes () {
		StringBuilder chrList = new StringBuilder();
		processedChromosomes.forEach((chr, value)->{
			chrList.append("|").append(chr).append("=").append(value);
		});
		
		System.out.println(chrList.substring(1));
		
		// append to logger
		Logger.append(chrList.substring(1));
		Logger.newLine();
	}
	
	/**
	 * @deprecated
	 * Do not use after v2.0.0
	 * 
	 */
	public static void printPSMCutoff () {
		System.out.println("Minimum PSM score threshold to accept as canonical PSMs: "+cPSMScoreTreshold);
		System.out.println("Minimum PSM score threshold to accept as noncanonical PSMs: "+ncPSMScoreTreshold);
		
		// append to logger
		Logger.append("Minimum PSM score threshold to accept as canonical PSMs: "+cPSMScoreTreshold);
		Logger.newLine();
		Logger.append("Minimum PSM score threshold to accept as noncanonical PSMs: "+ncPSMScoreTreshold);
		Logger.newLine();
	}
	
	public static void printFilterStat () {
		System.out.println("Step\tScans\tPeptides");
		System.out.println("Initial\t"+initialScanNum+"\t"+initialPeptideNum);
		System.out.println("Rank Filter\t"+rankFilterScanNum1+"\t"+rankFilterPeptideNum1);
		System.out.println("Length Filter\t"+lengthFilterScanNum2+"\t"+lengthFilterPeptideNum2);
		
		System.out.println("Mapped Candidates\t"+mappingFilterScanNum3+"\t"+mappingFilterPeptideNum3);
		if(Parameters.pvalue < 1.00) {
			System.out.println("Significantly Mapped\t"+pvalueFilterScanNum4+"\t"+pvalueFilterPeptideNum4);
		}
		
		
		// append to logger
		Logger.append("Step\tScans\tPeptides");
		Logger.newLine();
		Logger.append("Initial\t"+initialScanNum+"\t"+initialPeptideNum);
		Logger.newLine();
		Logger.append("Rank Filter\t"+rankFilterScanNum1+"\t"+rankFilterPeptideNum1);
		Logger.newLine();
		Logger.append("Length Filter\t"+lengthFilterScanNum2+"\t"+lengthFilterPeptideNum2);
		Logger.newLine();
		Logger.append("Mapped Candidates\t"+mappingFilterScanNum3+"\t"+mappingFilterPeptideNum3);
		Logger.newLine();
		if(Parameters.pvalue < 1.00) {
			Logger.append("Significantly Mapped\t"+pvalueFilterScanNum4+"\t"+pvalueFilterPeptideNum4);
			Logger.newLine();
		}
	}
	
}
