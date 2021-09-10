package progistar.pXg;

import java.io.IOException;

import progistar.pXg.constants.Parameters;
import progistar.pXg.data.PeptideAnnotation;
import progistar.pXg.data.parser.PeptideParser;

public class pMappingApp {

	public static void main(String[] args) throws IOException {
		Parameters.peptideFilePath = "C:\\Users\\progi\\Desktop\\Projects\\pXg\\peptideAnnotation.txt";
		Parameters.peptideColumnIndex = 4 - 1;
		PeptideAnnotation peptideAnnotation = PeptideParser.parseResult(Parameters.peptideFilePath);
		
		peptideAnnotation.writeFastaQuery();
	}
}
