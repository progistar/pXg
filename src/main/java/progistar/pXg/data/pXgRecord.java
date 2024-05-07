package progistar.pXg.data;

import progistar.pXg.constants.Parameters;
import progistar.pXg.data.parser.pXgParser;

public class pXgRecord {
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
		String pe = getValueByFieldName("isCanonical").equalsIgnoreCase("true") ? "1" : "2";
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
		String label = getValueByFieldName("Label");

		id = genomicLoci+":"+strand+":"+centerSeuqnece+":"+label;

		return id;
	}

	public boolean hasFastaID () {
		return getValueByFieldName("FastaIDs").equalsIgnoreCase("-") ? false : true;
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


	public String getValueByFieldName (String fieldName) {
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

	public void setValueByFieldName (String fieldName, String value) {
		String[] header = pXgParser.header;
		for(int i=0; i<header.length; i++) {
			if(header[i].equalsIgnoreCase(fieldName)) {
				fields[i] = value;
			}
		}
	}

	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();

		for(int i=0; i<fields.length; i++) {
			if(i != 0) {
				str.append("\t");
			}
			str.append(fields[i]);
		}

		return str.toString();
	}
}
