package progistar.pXg.data.parser;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import progistar.pXg.data.Cigar;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.utils.IndexConvertor;

/**
 * ready, parseSam, and finish methods must be run consecutively.<br>
 * ready: open the SAM file. <br>
 * parseSam: parsing SAM file one by one. <br>
 * finish: close the SAM file. <br>
 *
 * @author progi
 *
 */
public class SamParser {

	private static final Pattern EACH_CIGAR_REGEX = Pattern.compile("([0-9]+)([MINDSHPX=])");

	private static int QNAME_IDX 		= 0;
	public static int CHR_IDX 			= 2;
	public static int START_POS_IDX 	= 3;
	private static int CIGAR_IDX 		= 5;
	private static int SEQUENCE_IDX 	= 9;
	private static int QUALITY_IDX 		= 10;

	// prevent to generate constructor
	private SamParser () {}

	/**
	 * Read SAM file upto readNum. <br>

	 * @param readNum
	 * @return
	 */
	public static GenomicSequence parseSam (String samRead) {


		String[] fields = samRead.split("\\s");
		// TODO: unmapped reads cannot have genomic position information.
		// In case of unmapped read,  must consider it!
		String qName = fields[QNAME_IDX];
		String chr = fields[CHR_IDX];
		Integer startPosition = Integer.parseInt(fields[START_POS_IDX]);
		String cigarString = fields[CIGAR_IDX];
		String nucleotides = fields[SEQUENCE_IDX];
		String phred33 = fields[QUALITY_IDX];

		// average of phred33 QScore
		int length = phred33.length();
		double meanQScore = 0;
		for(int i=0; i<length; i++) {
			char qChar = phred33.charAt(i);
			if(qChar != '*') {
				meanQScore += (qChar-33);
			}
		}

		meanQScore /= length;

		// Note that
		// Chr of unmapped reads are marked as *
		// From this, we can recognize unmapped reads

		// Cigar has nucleotides and relative positions to the start position.
		ArrayList<Cigar> cigars = parseCigarString(cigarString, nucleotides);

		// find MD string
		String mdStr = null;
		for(int i=SEQUENCE_IDX; i<fields.length; i++) {
			if(fields[i].startsWith("MD:Z:")) {
				mdStr = fields[i].replace("MD:Z:", "");
				break;
			}
		}

		if(mdStr == null) {
			//System.out.println(line);
		}

		int chrIndex = IndexConvertor.chrToIndex(chr);

		// Genomic Sequence
		return new GenomicSequence(qName, chrIndex, startPosition, cigars, mdStr, meanQScore);

	}

	/**
	 * parsing Cigar string and assign nucleotide sequence to each Cigar operation. <br>
	 *
	 * @param cigarString
	 * @param nucleotides
	 * @return
	 */
	private static ArrayList<Cigar> parseCigarString (String cigarString, String nucleotides) {
		Matcher matcher = EACH_CIGAR_REGEX.matcher(cigarString);
		ArrayList<Cigar> results = new ArrayList<>();
		ArrayList<Cigar> filterResults = new ArrayList<>();

		// unmapped reads have no matcher...!
	    while (matcher.find()) {
	      int markerSize = Integer.parseInt(matcher.group(1));
	      char operation = matcher.group(2).charAt(0);

	      results.add(new Cigar(markerSize, operation));
	    }

	    // Unmapped read checker
	    if(results.size() == 0 && cigarString.equalsIgnoreCase("*")) {
	    	// add unmapped read cigar
	    	results.add(new Cigar(nucleotides.length(), '*'));
	    }

	    int ntIndex = 0;
	    int relPos = 0;
	    int[] relativePositions = null;
	    // drop all cigars without MIND
	    for(int i=0; i<results.size(); i++) {
	    	Cigar cigar = results.get(i);
	    	char op = cigar.operation;

	    	switch (op) {
	    	case 'M': // match or mismatch
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize);
	    		ntIndex += cigar.markerSize;

	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos++;
	    		}

	    		cigar.relativePositions = relativePositions;
	    		filterResults.add(cigar);
	    		break;

	    	case 'S': // soft clip
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize);
	    		ntIndex += cigar.markerSize;

	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos; // relPos is not changed... consistent!
	    		}

	    		cigar.relativePositions = relativePositions;
	    		filterResults.add(cigar);
	    		break;

	    	case 'I': // insertion
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize);
	    		ntIndex += cigar.markerSize;

	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos; // relPos is not changed... consistent!
	    		}

	    		cigar.relativePositions = relativePositions;
	    		filterResults.add(cigar);
	    		break;

	    	case 'D': // deletion
	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos++;
	    		}

	    		cigar.relativePositions = relativePositions;
	    		filterResults.add(cigar);
	    		break;

	    	case 'N': // skip (ex> exon junction)
	    		relPos += cigar.markerSize;
	    		filterResults.add(cigar);
	    		break;

	    	case '*': // unmapped
	    		cigar.nucleotides = nucleotides;
	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos++;
	    		}
	    		cigar.relativePositions = relativePositions;

	    		filterResults.add(cigar);
	    		break;
	    	}
	    }

	    return filterResults;
	}
}
