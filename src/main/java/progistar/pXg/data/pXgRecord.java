package progistar.pXg.data;

import progistar.pXg.constants.Parameters;
import progistar.pXg.data.parser.pXgParser;

public class pXgRecord {
	public static final int PE = 2;
	private String[] fields = null;
	
	public pXgRecord (String[] fields) {
		this.fields = fields;
	}
	
	public String getHeader () {
		StringBuilder header = new StringBuilder(">pXg");
		
		String id = getID();
		String isCanonical = getValueByFieldName("isCanonical").equalsIgnoreCase("true") ? "Canonical" : "Noncanonical";
		String gn = getValueByFieldName("GeneNames");
		String ev = getValueByFieldName("Events");
		String pe = PE+"";
		String pep = getValueByFieldName("InferredPeptide");
		String exp = getValueByFieldName("Reads");
		String var = getValueByFieldName("Mutations");
		String alt = getValueByFieldName("MutationStatus");
		String gId = getValueByFieldName("GeneIDs");
		String rna = getNucleotideSequence();
		
		header.append("|").append(id)
		.append("|").append(id+"_"+isCanonical)
		.append(" ").append(gId)
		.append(" ").append("GN="+gn)
		.append(" ").append("EV="+ev)
		.append(" ").append("PEP="+pep)
		.append(" ").append("EXP="+exp)
		.append(" ").append("VAR="+var)
		.append(" ").append("ALT="+alt)
		.append(" ").append("RNA="+rna)
		.append(" ").append("PE="+pe);
		
		return header.toString();
	}
	
	public String getID () {
		String id = null;
		
		String genomicLoci = getValueByFieldName("GenomicLoci").replace("|", ",");
		String centerSeuqnece = getValueByFieldName("ObservedNucleotide");
		String strand = getValueByFieldName("Strand");
		
		id = genomicLoci+":"+strand+":"+centerSeuqnece;
		
		return id;
	}
	
	public String getTranslatedSequence () {
		String strand = getValueByFieldName("Strand");
		String nucleotide = getNucleotideSequence();
		if(Parameters.isIncludedFlankSequence) {
			nucleotide = nucleotide.replaceAll("[-\\|]", "");
		} else {
			nucleotide = nucleotide.split("\\|")[1];
		}
		
		
		if(strand.equalsIgnoreCase("+")) {
			return GenomicSequence.translation(nucleotide, 0);
		} else {
			return GenomicSequence.reverseComplementTranslation(nucleotide, 0);
		}
	}
	
	public boolean isCanonical () {
		return getValueByFieldName("isCanonical").equalsIgnoreCase("true") ? true : false;
	}
	
	/**
	 * return 
	 * left-flank|center|right-flank
	 * 
	 * @return
	 */
	private String getNucleotideSequence () {
		StringBuilder sequence = new StringBuilder();
		String centerSeuqnece = getValueByFieldName("ObservedNucleotide").toUpperCase();
		String leftFlank = getValueByFieldName("ObservedLeftFlankNucleotide").toUpperCase();
		String rightFlank = getValueByFieldName("ObservedRightFlankNucleotide").toUpperCase();
		sequence.append(leftFlank);
		sequence.append("|");
		sequence.append(centerSeuqnece);
		sequence.append("|");
		sequence.append(rightFlank);
		
		return sequence.toString();
	}
	

	private String getValueByFieldName (String fieldName) {
		String[] header = pXgParser.header;
		String value = null;
		for(int i=0; i<header.length; i++) {
			if(header[i].equalsIgnoreCase(fieldName)) {
				if(value != null) {
					System.out.println(fieldName+" is duplciated");
				} else {
					value = fields[i];
				}
			}
		}
		return value;
	}
}