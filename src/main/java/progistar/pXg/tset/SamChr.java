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
		String fileName = args[0];
		
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		String line = null;
		
		Hashtable<String, String> chrs = new Hashtable<String, String>();
		long totalReads = 0;
		while((line = BR.readLine()) != null) {
			if(line.startsWith("@")) continue;
			String[] fields = line.split("\\s");
			String chr = fields[FieldIndex.CHR.value];
			chrs.put(chr, "");
			totalReads++;
		}
		
		BR.close();
		
		System.out.println(totalReads);
		chrs.forEach((chr, value) -> {
			System.out.println(chr);
		});
	}
}
