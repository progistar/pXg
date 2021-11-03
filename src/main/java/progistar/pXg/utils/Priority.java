package progistar.pXg.utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import progistar.pXg.constants.Constants;

public class Priority {

	// Region Score
	// XXX smaller is better
	// the left most digit: RegionStrVar (CDS = 0, UTR/Noncoding/Frameshift = 1, StrVar = 2)
	// the right most digit: Number of mutations
	public static final Pattern REGION_REG = Pattern.compile("\\([0-9A-Za-z;\\-]*\\)");
	
	/**
	 * 
	 * @param region
	 * @param mutations
	 * @return
	 */
	public static double getRegionPenalty (String region, String mutations) {
		double penalty = Double.MAX_VALUE;
		
		// region
		Matcher matcher = REGION_REG.matcher(region);
		
		if(matcher.find()) {
			penalty = 0;
			
			region = matcher.group();
			region = region.substring(1, region.length()-1); // remove parenthesis
			
			String[] regions = region.split("\\;");
			
			double mPenalty = 0;
			double rPenalty = 0;
			
			// mutation score
			if(!mutations.equalsIgnoreCase("-")) {
				mPenalty = mutations.split("\\|").length;
			}
			
			// region score
			// F = 5'-UTR
			// T = 3'-UTR
			// N = non-coding such as lncRNA/pseudogenes
			// I = in-frame (!= I means frameshift or non-coding)
			if(regions[0].contains(Constants.MARK_5UTR+"") || 
					regions[0].contains(Constants.MARK_3UTR+"") || regions[0].contains(Constants.MARK_NCDS+"") 
					|| regions[2].charAt(0) != Constants.IN_FRAME) {
				rPenalty += 10;
			}
			
			// strvar
			// I = intron
			// - = intergenic
			if(regions[0].contains(Constants.MARK_INTRON+"") || regions[0].contains(Constants.MARK_INTERGENIC+"")) {
				rPenalty += 20;
			}
			
			// antisense = antisense
			if(regions[1].equalsIgnoreCase("antisense")) {
				rPenalty += 20;
			}
			
			penalty = rPenalty + mPenalty;
		}
		
		return penalty;
	}
	
	/**
	 * Input region example: ENST0000...(3F5C;sense;N).<br>
	 * 
	 * @param region
	 * @return
	 */
	public static String getRegionEvent (String region) {
		String event = Constants.EVENT_PROTEINCODING;
		
		// region
		Matcher matcher = REGION_REG.matcher(region);
		
		if(matcher.find()) {
			region = matcher.group();
			region = region.substring(1, region.length()-1); // remove parenthesis
			
			String[] regions = region.split("\\;");
			
			// region score
			// F = 5'-UTR
			// T = 3'-UTR
			// N = non-coding such as lncRNA/pseudogenes
			// I = in-frame (!= I means frameshift or non-coding)
			if(regions[0].contains(Constants.MARK_5UTR+"")) {
				event = Constants.EVENT_5UTR;
			} else if(regions[0].contains(Constants.MARK_3UTR+"")) {
				event = Constants.EVENT_3UTR;
			} else if(regions[0].contains(Constants.MARK_NCDS+"")) {
				event = Constants.EVENT_NONCODING;
			} else if(regions[0].contains(Constants.MARK_INTERGENIC+"")) {
				event = Constants.EVENT_INTERGENIC;
			} else if(regions[0].contains(Constants.MARK_INTRON+"")) {
				event = Constants.EVENT_INTRON;
			} else if(regions[2].contains(Constants.OUT_OF_FRAME+"")) {
				event = Constants.EVENT_FRAMESHIFT;
			}
			
			if(regions[1].equalsIgnoreCase("antisense")) {
				event += ";"+Constants.EVENT_ANTISENSE;
			} else {
				event += ";"+Constants.EVENT_SENSE;
			}
		}
		
		
		return event;
	}
}
