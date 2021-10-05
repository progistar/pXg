package progistar.pXg.utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Priority {

	// Region Score
	// XXX smaller is better
	// the left most digit: RegionStrVar (CDS = 0, UTR/Noncoding/Frameshift = 1, StrVar = 2)
	// the right most digit: Number of mutations
	public static final Pattern REGION_REG = Pattern.compile("\\([0-9A-Za-z;\\-]*\\)");
	
	/**
	 * 
	 * Smaller is better. <br>
	 * 
	 * @param region
	 * @param mutations
	 * @return
	 */
	public static double getRegionScore (String region, String mutations) {
		double score = Double.MAX_VALUE;
		
		// region
		Matcher matcher = REGION_REG.matcher(region);
		
		if(matcher.find()) {
			score = 0;
			
			region = matcher.group();
			region = region.substring(1, region.length()-1); // remove parenthesis
			
			String[] regions = region.split("\\;");
			
			double mScore = 0;
			double rScore = 0;
			
			// mutation score
			if(!mutations.equalsIgnoreCase("-")) {
				mScore = mutations.split("\\|").length;
			}
			
			// region score
			// F = 5'-UTR
			// T = 3'-UTR
			// N = non-coding such as lncRNA/pseudogenes
			// I = in-frame (!= I means frameshift or non-coding)
			if(regions[0].contains("F") || regions[0].contains("T") || regions[0].contains("N") || regions[2].charAt(0) != 'I') {
				rScore = 10;
			}
			
			// strvar
			// I = intron
			// - = intergenic
			// anti-sense = anti-sense
			
			if(regions[0].contains("I") || regions[0].contains("-") || regions[1].equalsIgnoreCase("anti-sense")) {
				rScore = 20;
			}
			
			score = rScore + mScore;
		}
		
		
		return score;
	}
}
