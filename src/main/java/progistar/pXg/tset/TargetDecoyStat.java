package progistar.pXg.tset;

import java.io.File;
import java.io.IOException;

import progistar.pXg.data.PeptideAnnotation;

public class TargetDecoyStat {

	public static void main(String[] args) throws IOException {
		
		File[] fileList = new File[1];
		fileList[0] = new File("C:\\Users\\progi\\Desktop\\Projects\\pXg\\chr1toy.sam.worker1.tmp");
		
		PeptideAnnotation.parseTmpResults(fileList);
	}
}
