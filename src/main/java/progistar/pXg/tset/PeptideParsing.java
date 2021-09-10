package progistar.pXg.tset;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import progistar.pXg.constants.Parameters;

public class PeptideParsing {

	public static void main(String[] args) {
		String peptide = "K.+229.163EHS+26.049PDEFIK+229.163DEQNK+229.163.G";
		
		Pattern peptideRegExr = Pattern.compile(Parameters.peptideParserRegExr);
		
		StringBuilder pSeq = new StringBuilder();
		Matcher matcher = peptideRegExr.matcher(peptide);
		
		while(matcher.find()) {
			pSeq.append(matcher.group());
		}
		
		System.out.println(pSeq.toString());
	}
}
