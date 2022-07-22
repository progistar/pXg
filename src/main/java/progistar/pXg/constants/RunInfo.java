package progistar.pXg.constants;

import java.util.Hashtable;

import progistar.pXg.utils.Logger;

public class RunInfo {
	
	// Count PSMs
	// 0 for canonical
	// 1 for noncanonical
	public static long[][] targetRankPSMs = new long[2][101];
	public static long[][] decoyRankPSMs = new long[2][101];
	
	// Count PSMs
	public static long[] totalRankPSMs = new long[101];

	// Count peptides
	public static Hashtable[] targetPeptideHash		= new Hashtable[Parameters.maxPeptLen - Parameters.minPeptLen + 1];
	public static Hashtable[] decoyPeptideHash		= new Hashtable[Parameters.maxPeptLen - Parameters.minPeptLen + 1];
	public static Hashtable[] overlappedPeptideHash	= new Hashtable[Parameters.maxPeptLen - Parameters.minPeptLen + 1];
	
	// target decoy sequence statistics
	public static long[] targetAAFreqs = new long [26];
	public static long[] decoyAAFreqs = new long [26];
	
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
	
	// step5: region filter
	public static int regionFilterScanNum5			= 0;
	public static int regionFilterPeptideNum5		= 0;
	
	// step6: top-score filter
	public static int topscoreFilterScanNum6		= 0;
	public static int topscoreFilterPeptideNum6		= 0;
	
	// step7: FDR
	public static int fdrFilterScanNum7				= 0;
	public static int fdrFilterPeptideNum7			= 0;
	
	public static void printRankPSM () {
		Logger.append("Rank\tTargetCandidate\tDecoyCandidate\tIsCanonical");
		Logger.newLine();
		for(int i=0; i<RunInfo.targetRankPSMs[0].length; i++) {
			if(RunInfo.targetRankPSMs[0][i] == 0 && RunInfo.decoyRankPSMs[0][i] == 0) {
				continue;
			} else {
				Logger.append(i+"\t"+RunInfo.targetRankPSMs[0][i]+"\t"+RunInfo.decoyRankPSMs[0][i]+"\tTRUE");
				Logger.newLine();
			}
		}
		for(int i=0; i<RunInfo.targetRankPSMs[1].length; i++) {
			if(RunInfo.targetRankPSMs[1][i] == 0 && RunInfo.decoyRankPSMs[1][i] == 0) {
				continue;
			} else {
				Logger.append(i+"\t"+RunInfo.targetRankPSMs[1][i]+"\t"+RunInfo.decoyRankPSMs[1][i]+"\tFALSE");
				Logger.newLine();
			}
		}
	}
	
	public static void countTDPeptide (String target, String decoy) {
		int targetLen = target.length();
		
		// count target and decoy peptide
		for(int len=Parameters.minPeptLen; len<=Parameters.maxPeptLen; len++) {
			int lenIdx = len - Parameters.minPeptLen;
			if(RunInfo.targetPeptideHash[lenIdx] == null) {
				RunInfo.targetPeptideHash[lenIdx] = new Hashtable<String, String>();
			}
			
			if(RunInfo.decoyPeptideHash[lenIdx] == null) {
				RunInfo.decoyPeptideHash[lenIdx] = new Hashtable<String, String>();
			}
			
			if(RunInfo.overlappedPeptideHash[lenIdx] == null) {
				RunInfo.overlappedPeptideHash[lenIdx] = new Hashtable<String, String>();
			}
			
			for(int j=len; j<targetLen; j++) {
				String subTarget = target.substring(j-len, j);
				String subDecoy = decoy.substring(j-len, j);
				
				if(!subTarget.contains("X")) {
					RunInfo.targetPeptideHash[lenIdx].put(subTarget, "");
					if(RunInfo.decoyPeptideHash[lenIdx].get(subTarget) != null) {
						RunInfo.overlappedPeptideHash[lenIdx].put(subTarget, "");
					}
				}
				if(!subDecoy.contains("X")) {
					RunInfo.decoyPeptideHash[lenIdx].put(subDecoy, "");
					if(RunInfo.targetPeptideHash[lenIdx].get(subDecoy) != null) {
						RunInfo.overlappedPeptideHash[lenIdx].put(subDecoy, "");
					}
				}
				
			}
		}
		
	}
	
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
	
	public static void printPSMCutoff () {
		System.out.println("Minimum PSM score threshold to accept as canonical PSMs: "+cPSMScoreTreshold);
		System.out.println("Minimum PSM score threshold to accept as noncanonical PSMs: "+ncPSMScoreTreshold);
		
		// append to logger
		Logger.append("Minimum PSM score threshold to accept as canonical PSMs: "+cPSMScoreTreshold);
		Logger.newLine();
		Logger.append("Minimum PSM score threshold to accept as noncanonical PSMs: "+ncPSMScoreTreshold);
		Logger.newLine();
	}
	
	public static void printAAStat () {
		System.out.println("AA\tTargetAA\tDecoyAA");
		Logger.append("AA\tTargetAA\tDecoyAA");
		Logger.newLine();
		for(int i=0; i<26; i++) {
			char aa = Character.valueOf((char) ('A'+i));
			System.out.println(aa+"\t"+RunInfo.targetAAFreqs[i]+"\t"+RunInfo.decoyAAFreqs[i]);
			Logger.append(aa+"\t"+RunInfo.targetAAFreqs[i]+"\t"+RunInfo.decoyAAFreqs[i]);
			Logger.newLine();
		}
	}
	
	public static void printFilterStat () {
		System.out.println("Step\tScans\tPeptides");
		System.out.println("Initial\t"+initialScanNum+"\t"+initialPeptideNum);
		System.out.println("Rank Filter\t"+rankFilterScanNum1+"\t"+rankFilterPeptideNum1);
		System.out.println("Length Filter\t"+lengthFilterScanNum2+"\t"+lengthFilterPeptideNum2);
		
		System.out.println("Mapped Candidates\t"+mappingFilterScanNum3+"\t"+mappingFilterPeptideNum3);
		System.out.println("Significantly Mapped\t"+pvalueFilterScanNum4+"\t"+pvalueFilterPeptideNum4);
		System.out.println("Region Decision\t"+regionFilterScanNum5+"\t"+regionFilterPeptideNum5);
		System.out.println("PSM Decision\t"+topscoreFilterScanNum6+"\t"+topscoreFilterPeptideNum6);
		System.out.println("FDR Estimation\t"+fdrFilterScanNum7+"\t"+fdrFilterPeptideNum7);
		
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
		Logger.append("Significantly Mapped\t"+pvalueFilterScanNum4+"\t"+pvalueFilterPeptideNum4);
		Logger.newLine();
		Logger.append("Region Decision\t"+regionFilterScanNum5+"\t"+regionFilterPeptideNum5);
		Logger.newLine();
		Logger.append("PSM Decision\t"+topscoreFilterScanNum6+"\t"+topscoreFilterPeptideNum6);
		Logger.newLine();
		Logger.append("FDR Estimation\t"+fdrFilterScanNum7+"\t"+fdrFilterPeptideNum7);
		Logger.newLine();
	}
	
}
