package progistar.thirdparty.sSim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;

public class PeptideLoader {

	// key: peptide
	// value: pool ID
	public static Hashtable<String, String> matchedPeptides = new Hashtable<>();

	public static void loadPeptideList (String peptideList) throws IOException {
		File file = new File(peptideList);

		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;

		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");

			if(matchedPeptides.get(fields[0]) != null) {
				System.out.println("Unexpected case: duplicated peptide was detected... => "+fields[0]+" in "+fields[1]) ;
			} else {
				matchedPeptides.put(fields[0], fields[1]);
			}
		}


		BR.close();
	}
}
