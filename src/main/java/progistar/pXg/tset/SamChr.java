package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class SamChr {

	private enum FieldIndex {
		QNAME(0), CHR(2), START_POS(3), CIGAR(5), SEQUENCE(9);

		private int value;

		FieldIndex(int value) {
			this.value = value;
		}
	}
	
	public static void main(String[] args) throws IOException {
		String fileName = "C:\\Users\\progi\\eclipse-workspace\\pXg\\test\\chr1toy.sam";
		
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		String line = null;
		
		int count = 1000;
		while((line = BR.readLine()) != null) {
			if(line.startsWith("@")) {
				count--;
			}
			
			System.out.println(line);
			
			if(count == 0 ) break;
		}
		
		BR.close();
	}
}
