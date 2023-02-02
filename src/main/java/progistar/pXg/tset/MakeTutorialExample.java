package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class MakeTutorialExample {

	public static void main(String[] args) throws IOException {
		File samFile = new File("/Users/gistar/eclipse-workspace/pXg/tutorial/test.sorted.sam");
		
		BufferedReader BR = new BufferedReader(new FileReader(samFile));
		String line = null;
		
		while((line = BR.readLine()) != null) {
			String asterisk = line.split("\\s")[2];
			if(asterisk.equalsIgnoreCase("*")) {
				String sequence = line.split("\\s")[9];
				if(sequence.contains("GGAGAAGTCATCGGAACTCGATGG") && !sequence.contains("CCATCGAGTTCCGATGACTTCTCC")) {
					System.out.println(line);
				}
			} else {
				System.out.println(line);
			}
		}
		
		BR.close();
	}
}
