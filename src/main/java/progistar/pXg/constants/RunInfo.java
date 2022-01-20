package progistar.pXg.constants;

import java.util.Hashtable;

public class RunInfo {

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
	
	public static void printProcessedChromosomes () {
		StringBuilder chrList = new StringBuilder();
		processedChromosomes.forEach((chr, value)->{
			chrList.append("|").append(chr).append("=").append(value);
		});
		
		System.out.println(chrList.substring(1));
	}
	
	public static void printPSMCutoff () {
		System.out.println("Minimum PSM score threshold to accept as canonical PSMs: "+cPSMScoreTreshold);
		System.out.println("Minimum PSM score threshold to accept as noncanonical PSMs: "+ncPSMScoreTreshold);
	}
	
	public static void printFilterStat () {
		System.out.println("Step\tScans\tPeptides");
		System.out.println("Initial\t"+initialScanNum+"\t"+initialPeptideNum);
		System.out.println("Rank\t"+rankFilterScanNum1+"\t"+rankFilterPeptideNum1);
		System.out.println("Length\t"+lengthFilterScanNum2+"\t"+lengthFilterPeptideNum2);
		
		System.out.println("ExpReads\t"+mappingFilterScanNum3+"\t"+mappingFilterPeptideNum3);
		System.out.println("p-value\t"+pvalueFilterScanNum4+"\t"+pvalueFilterPeptideNum4);
		System.out.println("RegionPanelty\t"+regionFilterScanNum5+"\t"+regionFilterPeptideNum5);
		System.out.println("top-score\t"+topscoreFilterScanNum6+"\t"+topscoreFilterPeptideNum6);
		System.out.println("FDR\t"+fdrFilterScanNum7+"\t"+fdrFilterPeptideNum7);
	}
}
