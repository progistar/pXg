package progistar.pXg.tset;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

import progistar.pXg.constants.Parameters;
import progistar.pXg.data.PxGAnnotation;
import progistar.pXg.data.parser.PeptideParser;
import progistar.pXg.data.parser.ResultParser;

public class PriorityTest {
	
	public static void main(String[] args) throws IOException {
		
		PeptideParser.parseResult("C:\\Users\\progi\\Desktop\\Projects\\pXg\\peptideAnnotationM.txt"); // static..!
		
		// read tmp output files
		ArrayList<File> tmpOutputFiles = new ArrayList<File>();
		tmpOutputFiles.add(new File("C:\\Users\\progi\\eclipse-workspace\\pXg\\test\\SRR1925276_Aligned.filter.sorted.sam.worker1.tmp"));

		
		PxGAnnotation pXgA = ResultParser.parseResult(tmpOutputFiles);
		// filter by pvalue
		// among them, use highest-scored PSM
		pXgA.topScoreFilter();
		// filter regions
		pXgA.regionScoreFilter();
	}
	
}
