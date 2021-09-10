package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Parameter;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import progistar.pXg.constants.Parameters;
import progistar.pXg.data.Cigar;
import progistar.pXg.data.GenomicAnnotation;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.processor.Mapper;
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
	private static File file = null;
	private static BufferedReader BR = null;
	private static long totalReadCount = 0;
	private static boolean isEndOfFile = false;
	
	// prevent to generate constructor
	private SamParser () {}
	
	private enum FieldIndex {
		QNAME(0), CHR(2), START_POS(3), CIGAR(5), SEQUENCE(9);

		private int value;

		FieldIndex(int value) {
			this.value = value;
		}
	}
	
	/**
	 * Create BufferedReader for "samFilePath" <br>
	 * 
	 * @param samFilePath
	 */
	public static void ready (String samFilePath) {
		try {
			file = new File(samFilePath);
			BR = new BufferedReader(new FileReader(file));
		}catch(IOException ioe) {
			
		}
	}
	
	public static void finish () {
		try {
			if(BR != null) BR.close();
		}catch(IOException ioe) {
			
		}
	}
	
	public static boolean isEndOfFile () {
		return isEndOfFile;
	}
	
	/**
	 * Read SAM file upto readNum. <br>

	 * @param readNum
	 * @return
	 */
	public static ArrayList<GenomicSequence> parseSam (long readNum) {
		assert BR != null;
		
		long startTime = System.currentTimeMillis();
		
		ArrayList<GenomicSequence> gSeqs = new ArrayList<GenomicSequence>();
		
		System.out.print("reading "+file.getName()+"... ("+(totalReadCount+1)+"-"+(totalReadCount+readNum)+")");
		try {
			
			String line = null;
			
			long readCount = 0;
			while((line = BR.readLine()) != null) {
				if(line.startsWith("@")) continue; // skip meta
				
				String[] fields = line.split("\\s");
				
				// TODO: unmapped reads cannot have genomic position information.
				// In case of unmapped read,  must consider it! 
				String qName = fields[FieldIndex.QNAME.value];
				String chr = fields[FieldIndex.CHR.value];
				Integer startPosition = Integer.parseInt(fields[FieldIndex.START_POS.value]);
				String cigarString = fields[FieldIndex.CIGAR.value];
				String nucleotides = fields[FieldIndex.SEQUENCE.value];
				
				// Cigar has nucleotides and relative positions to the start position.
				ArrayList<Cigar> cigars = parseCigarString(cigarString, nucleotides);
				
				// the index for that chr is automatically assigned by auto-increment key.
				IndexConvertor.putChrIndexer(chr);
				byte chrIndex = IndexConvertor.chrToIndex(chr);
				
				// Genomic Sequence
				GenomicSequence gSeq = new GenomicSequence(qName, chrIndex, startPosition, cigars);
				gSeqs.add(gSeq);
				
				readCount ++;
				if(readCount == readNum) break;
			}
			
			totalReadCount += readCount;
			
			// is end of file
			if(line == null) {
				isEndOfFile = true;
			}
		}catch(IOException ioe) {
			
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");
		
		return gSeqs;
		
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
		ArrayList<Cigar> results = new ArrayList<Cigar>();
		
	    while (matcher.find()) {
	      int markerSize = Integer.parseInt(matcher.group(1));
	      char operation = matcher.group(2).charAt(0);

	      results.add(new Cigar(markerSize, operation));
	    }
	    
	    int ntIndex = 0;
	    int relPos = 0;
	    int[] relativePositions = null;
	    for(int i=0; i<results.size(); i++) {
	    	Cigar cigar = results.get(i);
	    	char op = cigar.operation;
	    	
	    	switch (op) {
	    	case 'S': // soft clip
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize);
	    		ntIndex += cigar.markerSize;
	    		break;
	    		
	    	case 'M': // match or mismatch
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize);
	    		ntIndex += cigar.markerSize;
	    		
	    		relativePositions = new int[cigar.nucleotides.length()];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos++;
	    		}
	    		
	    		cigar.relativePositions = relativePositions;
	    		
	    		break;
	    		
	    	case 'I': // insertion
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize); 
	    		ntIndex += cigar.markerSize;
	    		
	    		relativePositions = new int[cigar.nucleotides.length()];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos; // relPos is not changed... consistent!
	    		}
	    		
	    		cigar.relativePositions = relativePositions;
	    		
	    		break;
	    		
	    	case 'D': // deletion
	    		relativePositions = new int[cigar.nucleotides.length()];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos++;
	    		}
	    		
	    		cigar.relativePositions = relativePositions;
	    		
	    		break;	
	    		
	    	case 'N': // skip (ex> exon junction)
	    		relPos += cigar.markerSize;
	    		
	    		break;	
	    		
	    	case '*': // unmapped
	    		//TODO unmapped reads
	    		
	    		break;
	    	}
	    }

	    return results;
	}
}
